"""
Inception curves: locate the roots of det Q(E/N, p·d) = 0 and track them
across a p·d sweep, giving the (partial) discharge inception voltage
PDIV = V*(p·d).

For each value of p·d the critical E/N satisfying det Q = 0 is found with
the Brent root-finding method (scipy.optimize.brentq) after a logarithmic
sign-change scan.  The corresponding breakdown voltage is

    V* = E* · d = (E/N)* · (p·d) · 1e-16 / (kB · T)

using N = p · 1e5 / (kB·T) and the unit conversion between Townsend
(1e-21 V m²) and bar (1e5 Pa).
"""

import math

import numpy as np
import scipy.optimize

from incept1d.constants import kB as _kB
from incept1d.solver import inception_det

_ROOT_CHECK_TOL = 0.01  # |det Q(root)| / bracket-scale above this triggers a warning
_NAN_SENTINEL = 1e-290  # |value| below this → came from NaN → -1e-300 mapping


def _accept_root(root, EN_a, EN_b, fa, fb, pd, det_fn, mod, p, T):
    """
    Verify that det Q is genuinely near zero at a candidate root.

    Returns True if the root is accepted, False if it should be discarded.

    Two failure modes are distinguished:

    1. NaN-sentinel artefact: f_brentq mapped NaN → -1e-300 at a bracket
       endpoint, creating an artificial sign change.  det_fn(root) = NaN and
       at least one bracket fa/fb equals the sentinel (``|value| <= 1e-290``).
       → Discard and print a warning.

    2. Genuine singularity: both bracket endpoints are finite with opposite
       signs; brentq converges to the true det Q = 0 locus where Q is exactly
       singular so cond(Q) → ∞ and det_fn returns NaN.
       → Accept silently (NaN here is the expected consequence of Q → 0).

    3. Large finite residual: det_fn(root) is finite but ``|det Q|`` is large
       relative to the bracket scale.
       → Accept but print a warning.

    Parameters
    ----------
    root : float        — candidate E/N root returned by brentq
    EN_a, EN_b : float  — bracket endpoints used by brentq
    fa, fb : float      — f_brentq values at the bracket endpoints
                          (must be the FINAL bracket values, not the original
                           proxy-scan values — caller is responsible)
    pd : float          — pressure × gap length in bar·m
    det_fn : callable   — full-resolution det function (used for residual eval)
    mod, p, T           — mechanism, pressure, temperature passed to det_fn
    """
    det_root = det_fn(root, pd, mod, p, T)

    if not np.isfinite(det_root):
        # Distinguish genuine singularity from NaN-sentinel artefact.
        fa_sentinel = abs(fa) <= _NAN_SENTINEL
        fb_sentinel = abs(fb) <= _NAN_SENTINEL
        if fa_sentinel or fb_sentinel:
            print(
                f"  [det Q check] EN = {root:.4f} Td  pd = {pd*1e3:.4g} bar·mm: "
                f"det Q = NaN at root; bracket endpoint is NaN sentinel "
                f"(fa={fa:.3e}, fb={fb:.3e}) — artefact, discarded"
            )
            return False
        # Both bracket endpoints were genuine finite values; NaN at root means
        # Q is exactly singular (cond → ∞, det → 0).  Root is genuine.
        return True

    bracket_scale = max(
        abs(fa) if np.isfinite(fa) else 0.0, abs(fb) if np.isfinite(fb) else 0.0, 1e-300
    )
    if abs(det_root) > _ROOT_CHECK_TOL * bracket_scale:
        print(
            f"  [det Q check] EN = {root:.4f} Td  pd = {pd*1e3:.4g} bar·mm: "
            f"det Q = {det_root:.3e}  (bracket scale {bracket_scale:.3e}) — suspect root"
        )

    return True


def find_all_breakdown_EN(
    pd,
    mod,
    p,
    T,
    EN_lo=10.0,
    EN_hi=3e5,
    n_scan=200,
    first_only=False,
    det_fn=None,
    fast_det_fn=None,
    med_det_fn=None,
    EN_hints=None,
    hint_factor=2.0,
):
    """
    Find E/N values where det Q(E/N, pd) = 0.

    A logarithmic scan over [EN_lo, EN_hi] locates every sign change in
    det Q(EN); scipy.optimize.brentq refines each bracket independently.
    Using n_scan=200 reduces the chance of missing a sign change when multiple
    roots are close together.

    With first_only=True the value returned is always the *globally lowest*
    root in [EN_lo, EN_hi]; all roots are returned only when first_only=False
    (``--all-branches``).

    If EN_hints is supplied and first_only=True, a warm-start pass is
    attempted first: the lowest hint is bracketed by
    [hint/hint_factor, hint*hint_factor] and refined directly with brentq.
    The resulting root is accepted — and the coarse scan skipped — only if it
    also brackets successfully *and* no sign change is found beneath it
    (sampled at the same points-per-decade as the full scan); otherwise the
    full scan runs and decides.  Warm-start is intentionally disabled when
    first_only=False (--all-branches) because new branches can appear at any
    pd step; the coarse scan is required to detect them.

    Parameters
    ----------
    pd : float
        Product of pressure and gap length in bar·m.
    mod : Mechanism
        Loaded mechanism (returned by load_mechanism).
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.
    EN_lo : float
        Lower bound of the E/N search range in Townsend.  Default 10.0 Td.
    EN_hi : float
        Upper bound of the E/N search range in Townsend.  Default 3e5 Td.
    n_scan : int
        Number of points in the coarse scan.  Default 200.
    first_only : bool
        If True, return as soon as the first (lowest E/N) root is found.
    fast_det_fn : callable or None
        Optional cheaper det function used only for the coarse sign-change
        scan.  Brentq refinement always uses det_fn for full accuracy.
        If None, det_fn is used for both scan and refinement.
    med_det_fn : callable or None
        Optional medium-fidelity det function used for the fallback scan when
        fast_det_fn finds no sign changes.  Cheaper than det_fn for scanning
        but more accurate than fast_det_fn.  If None, det_fn is used as the
        fallback.
    EN_hints : list of float or None
        E/N values from the previous pd point to use as warm-start brackets.
        When all hints bracket successfully the full scan is skipped entirely.
        Pass [] or None to disable (first pd point, or after a failed pd).
    hint_factor : float
        Bracket half-width multiplier: each hint is bracketed by
        [hint/hint_factor, hint*hint_factor].  Default 2.0 (one octave each way).

    Returns
    -------
    list of float
        Critical E/N values in Townsend, sorted ascending.
        Empty list if no sign change is found in [EN_lo, EN_hi].
    """
    if det_fn is None:
        det_fn = inception_det
    scan_fn = fast_det_fn if fast_det_fn is not None else det_fn

    def f_brentq(EN):
        # NaN → small negative sentinel: NaN always occurs above the critical
        # E/N where the true sign is −1, so -1e-300 is physically correct.
        val = det_fn(EN, pd, mod, p, T)
        return val if np.isfinite(val) else -1e-300

    def _no_root_below(EN_ref):
        """
        True when the coarse scan finds no sign change in [EN_lo, EN_ref].

        A warm-start root is only the *lowest* root if nothing lies below it.
        The check samples [EN_lo, EN_ref] at the same points-per-decade as the
        full scan over [EN_lo, EN_hi], so accepting a warm-start root is
        exactly as reliable as running the scan it replaces — just cheaper,
        because the interval is shorter.
        """
        if EN_ref <= EN_lo:
            return True
        per_decade = (n_scan - 1) / math.log10(EN_hi / EN_lo)
        n_guard = max(2, int(math.ceil(per_decade * math.log10(EN_ref / EN_lo))) + 1)
        grid = np.logspace(np.log10(EN_lo), np.log10(EN_ref), n_guard)
        vals = np.array([scan_fn(en, pd, mod, p, T) for en in grid])
        vals = np.where(np.isfinite(vals), vals, 0.0)  # NaN → 0, as in _scan_with
        return not np.any(vals[:-1] * vals[1:] < 0.0)

    # Warm-start: only attempted when first_only=True (single-branch mode).
    # With first_only=False (--all-branches) the coarse scan is mandatory
    # because new branches can appear at any pd point; skipping it would
    # silently miss roots not bracketed by existing hints.
    #
    # Even in single-branch mode the hint only says where the *previous* pd
    # point's root was: a new, lower root may have appeared since (this is
    # exactly what happens near the left-branch asymptote, where the lowest
    # root drops by orders of magnitude between neighbouring pd points).
    # Refining the hint alone would then return a root that is not the lowest,
    # and the branch would stay locked onto the wrong one for the rest of the
    # sweep, so the warm-start root is accepted only after _no_root_below
    # confirms that nothing lies underneath it.
    if EN_hints and first_only:
        warm_roots = []
        all_ok = True
        hints_to_try = EN_hints[:1] if first_only else EN_hints
        for en_hint in hints_to_try:
            EN_a = max(EN_lo, en_hint / hint_factor)
            EN_b = min(EN_hi, en_hint * hint_factor)
            fa = f_brentq(EN_a)
            fb = f_brentq(EN_b)
            if fa * fb < 0.0:
                # When the upper endpoint is a NaN sentinel, brentq would
                # converge to the NaN boundary rather than the true root.
                # Narrow EN_b geometrically (in log-EN space) until we land on
                # a finite negative value, giving brentq a clean bracket.
                if abs(fb) <= _NAN_SENTINEL:
                    lo_fn, hi_nan = EN_a, EN_b
                    found = False
                    for _ in range(20):
                        mid = math.sqrt(lo_fn * hi_nan)
                        v = det_fn(mid, pd, mod, p, T)
                        if not np.isfinite(v):
                            hi_nan = mid
                        elif v < 0.0:
                            EN_b, fb = mid, v
                            found = True
                            break
                        else:
                            lo_fn = mid
                    if not found:
                        all_ok = False
                        break
                root = scipy.optimize.brentq(
                    f_brentq, EN_a, EN_b, xtol=1e-8, rtol=1e-14
                )
                if _accept_root(root, EN_a, EN_b, fa, fb, pd, det_fn, mod, p, T):
                    warm_roots.append(float(root))
                else:
                    all_ok = False
                    break
            else:
                all_ok = False
                break
        if all_ok and warm_roots and _no_root_below(min(warm_roots)):
            return sorted(warm_roots)
        # Warm start incomplete, or a lower root exists — fall through to the
        # full scan, which is authoritative.

    EN_scan = np.logspace(np.log10(EN_lo), np.log10(EN_hi), n_scan)

    def _scan_with(sfn):
        """Coarse sign-change scan using sfn; Brentq refinement always uses det_fn."""

        def f_s(EN):
            val = sfn(EN, pd, mod, p, T)
            return val if np.isfinite(val) else 0.0

        D = np.array([f_s(en) for en in EN_scan])
        local_roots = []
        for i in range(len(EN_scan) - 1):
            if D[i] * D[i + 1] < 0.0:
                EN_a, EN_b = EN_scan[i], EN_scan[i + 1]
                fa, fb = f_brentq(EN_a), f_brentq(EN_b)
                if fa * fb >= 0.0:
                    # Proxy scan bracket doesn't hold for the full det; widen by
                    # one scan step on each side and re-scan with the full det.
                    lo = EN_scan[max(0, i - 1)]
                    hi = EN_scan[min(len(EN_scan) - 1, i + 2)]
                    fine = np.logspace(np.log10(lo), np.log10(hi), 30)
                    fine_v = [f_brentq(en) for en in fine]
                    found = False
                    for k in range(len(fine) - 1):
                        if fine_v[k] * fine_v[k + 1] < 0.0:
                            EN_a, EN_b = fine[k], fine[k + 1]
                            fa, fb = fine_v[k], fine_v[k + 1]
                            found = True
                            break
                    if not found:
                        continue
                root = scipy.optimize.brentq(
                    f_brentq, EN_a, EN_b, xtol=1e-8, rtol=1e-14
                )
                if not _accept_root(root, EN_a, EN_b, fa, fb, pd, det_fn, mod, p, T):
                    continue
                local_roots.append(float(root))
                if first_only:
                    break
            elif D[i] > 0.0 and D[i + 1] == 0.0:
                # sfn returned NaN (mapped to 0) at the upper endpoint.  A sign
                # change may be hidden just before the NaN boundary; fine-scan
                # with the full det_fn to locate it.
                fine = np.logspace(np.log10(EN_scan[i]), np.log10(EN_scan[i + 1]), 40)
                fine_v = [det_fn(en, pd, mod, p, T) for en in fine]
                pos_j = None
                for j, v in enumerate(fine_v):
                    if np.isfinite(v) and v > 0.0:
                        pos_j = j
                    elif np.isfinite(v) and v < 0.0 and pos_j is not None:
                        root = scipy.optimize.brentq(
                            f_brentq, fine[pos_j], fine[j], xtol=1e-8, rtol=1e-14
                        )
                        if not _accept_root(
                            root,
                            fine[pos_j],
                            fine[j],
                            fine_v[pos_j],
                            v,
                            pd,
                            det_fn,
                            mod,
                            p,
                            T,
                        ):
                            break
                        local_roots.append(float(root))
                        break
                if first_only and local_roots:
                    break
        return local_roots

    roots = _scan_with(scan_fn)
    if not roots and fast_det_fn is not None:
        # Fast scan found no sign changes; fall back to the full det_fn scan.
        # med_det_fn is intentionally not used here: it almost never brackets
        # roots that fast_det missed, so using it only adds scan overhead before
        # the full fallback (which is always needed anyway).
        roots = _scan_with(det_fn)
    return roots


def compute_inception_curve(
    pd_arr,
    mod,
    p,
    T,
    all_branches=False,
    det_fn=None,
    fast_det_fn=None,
    med_det_fn=None,
):
    """
    Compute inception-curve branches across a p*d sweep.

    For each pd point roots of det Q(E/N, pd) = 0 are found via
    find_all_breakdown_EN.  Branches are tracked by continuity: each new root
    is matched to the existing branch whose most-recent E/N is nearest in
    log space, using a greedy nearest-neighbour assignment.  Unmatched roots
    start new branches.  This correctly handles saddle-node bifurcations where
    two new low-E/N roots appear without disrupting the pre-existing branch.

    Parameters
    ----------
    pd_arr : numpy.ndarray, shape (M,)
        Array of p*d values in bar·m.
    mod : Mechanism
        Loaded mechanism (returned by load_mechanism).
    p : float or numpy.ndarray, shape (M,)
        Gas pressure in bar.  A scalar is used for all points (fixed-p mode);
        an array of length M allows pressure to vary per point (fixed-d mode).
    T : float
        Gas temperature in Kelvin.
    all_branches : bool
        If False (default), only the lowest-E/N root (branch 0) is kept at
        each pd point.  If True, all roots are collected.

    Returns
    -------
    list of dict
        One dict per branch.  Branch 0 is the one whose first root appeared
        earliest (lowest pd); within that, ordered by first-appearance E/N.
        Each dict has keys:

        - ``'idx'`` — integer indices into pd_arr where this branch has a root
        - ``'pd'``, ``'EN'``, ``'V'``, ``'p'``, ``'d'`` — corresponding arrays

        Returns an empty list if no solutions are found anywhere.
    """
    if det_fn is None:
        det_fn = inception_det
    p_arr = np.full_like(pd_arr, p) if np.ndim(p) == 0 else np.asarray(p, dtype=float)
    branches = []
    branch_last_logEN = []  # last log10(EN) for each branch, for continuity tracking
    prev_roots_EN = []  # roots found at the previous pd point (warm-start hints)

    for i, (pd_i, p_i) in enumerate(zip(pd_arr, p_arr)):
        roots = find_all_breakdown_EN(
            pd_i,
            mod,
            p_i,
            T,
            first_only=not all_branches,
            det_fn=det_fn,
            fast_det_fn=fast_det_fn,
            med_det_fn=med_det_fn,
            EN_hints=prev_roots_EN,
        )
        prev_roots_EN = list(roots)
        if not roots:
            continue
        d_i = pd_i / p_i

        if not branches:
            # First pd with roots: create one branch per root, sorted ascending
            for EN in sorted(roots):
                V = EN * pd_i * 1e-21 / (_kB * T) * 1e5
                branches.append(
                    {
                        "idx": [i],
                        "pd": [pd_i],
                        "EN": [EN],
                        "V": [V],
                        "p": [p_i],
                        "d": [d_i],
                    }
                )
                branch_last_logEN.append(np.log10(EN))
        else:
            # Match each new root to the nearest existing branch in log-EN space.
            # Build cost matrix: rows = existing branches, cols = new roots.
            log_roots = np.log10(np.array(sorted(roots)))
            log_last = np.array(branch_last_logEN)
            cost = np.abs(log_last[:, None] - log_roots[None, :])

            assigned_branches = set()
            assigned_roots = set()
            root_to_branch = {}

            # Greedy: repeatedly pick the (branch, root) pair with smallest distance
            flat_order = np.argsort(cost, axis=None)
            for idx in flat_order:
                b, r = divmod(int(idx), len(log_roots))
                if b in assigned_branches or r in assigned_roots:
                    continue
                root_to_branch[r] = b
                assigned_branches.add(b)
                assigned_roots.add(r)
                if len(assigned_roots) == min(len(branches), len(log_roots)):
                    break

            sorted_roots = sorted(roots)
            for r, EN in enumerate(sorted_roots):
                V = EN * pd_i * 1e-21 / (_kB * T) * 1e5
                if r in root_to_branch:
                    b = root_to_branch[r]
                    branches[b]["idx"].append(i)
                    branches[b]["pd"].append(pd_i)
                    branches[b]["EN"].append(EN)
                    branches[b]["V"].append(V)
                    branches[b]["p"].append(p_i)
                    branches[b]["d"].append(d_i)
                    branch_last_logEN[b] = np.log10(EN)
                else:
                    branches.append(
                        {
                            "idx": [i],
                            "pd": [pd_i],
                            "EN": [EN],
                            "V": [V],
                            "p": [p_i],
                            "d": [d_i],
                        }
                    )
                    branch_last_logEN.append(np.log10(EN))

    for br in branches:
        for k in ("pd", "EN", "V", "p", "d"):
            br[k] = np.array(br[k], dtype=float)
        br["idx"] = np.array(br["idx"], dtype=int)

    return branches
