# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Inception curves: locate the roots of det Q(E/N, p·d) = 0 and track them
across a p·d sweep, giving the (partial) discharge inception voltage
PDIV = V*(p·d).

The root finding, warm starts and branch tracking are described in
``docs/source/numerics/rootfinding.rst``.  The criterion whose root is
sought is :func:`incept1d.solver.riccati_criterion` unless a ``det_fn`` is
given; :func:`incept1d.solver.inception_det` evaluates the same condition
as det Q.
"""

import math
import time

import numpy as np
import scipy.optimize

from incept1d.constants import kB as _kB
from incept1d.parallel import parallel_map
from incept1d.solver import riccati_criterion

_ROOT_CHECK_TOL = 0.01  # |det Q(root)| / bracket-scale above this triggers a warning
_NAN_SENTINEL = 1e-290  # |value| below this → came from NaN → -1e-300 mapping


def _accept_root(root, EN_a, EN_b, fa, fb, pd, det_fn, mod, p, T):
    """
    Verify that det Q is genuinely near zero at a candidate root.

    A NaN at the root is expected — Q is exactly singular there — so the
    question is whether the sign change that produced it was real.  A
    bracket endpoint sitting at the NaN sentinel means it was not, and the
    root is discarded with a warning.

    Parameters
    ----------
    root : float
        Candidate E/N root returned by brentq.
    EN_a, EN_b : float
        Bracket endpoints used by brentq.
    fa, fb : float
        Values at those endpoints.  These must be the *final* bracket
        values, not the ones from the proxy scan; the caller is
        responsible for that.
    pd : float
        Pressure × gap length in bar·m.
    det_fn : callable
        Full-resolution determinant, for evaluating the residual.
    mod, p, T :
        Mechanism, pressure and temperature, passed to det_fn.

    Returns
    -------
    bool
        True if the root is accepted.
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


#: How far beyond a proxy bracket the full determinant is searched before the
#: bracket is declared unconfirmed.  A uniform-field proxy misplaces a root of
#: a strongly non-uniform field by up to ~30 %, so this leaves margin on top.
_CONFIRM_FACTOR = 1.6

#: Lowest E/N (Td) searched when the criterion shows inception below EN_lo.
_EN_FLOOR = 0.1


def _has_physical_sign(det_fn):
    """True when *det_fn* is (a partial of) the Riccati criterion."""
    return getattr(det_fn, "func", det_fn) is riccati_criterion


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

    A logarithmic scan brackets every sign change and Brent's method
    refines each bracket.  With ``first_only`` the value returned is
    always the *globally lowest* root in the search range, whether it came
    from the scan or from a warm start.

    Parameters
    ----------
    pd : float
        Product of pressure and gap length in bar·m.
    mod : Mechanism
        Loaded mechanism.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.
    EN_lo, EN_hi : float
        Bounds of the E/N search range in Townsend.
    n_scan : int
        Number of points in the coarse scan.
    first_only : bool
        Keep only the lowest root.  Warm starts apply in this mode only,
        since in all-branches mode a new branch can appear at any pd and
        the scan is needed to find it.
    fast_det_fn, med_det_fn : callable or None
        Cheaper determinants for the coarse scan and its fallback.  Brent
        refinement always uses det_fn.  Both default to det_fn.
    EN_hints : list of float or None
        Roots from the previous pd point, used as warm-start brackets.
        Pass [] or None at the first pd point, or after a failed one.
    hint_factor : float
        Each hint is bracketed by [hint/hint_factor, hint*hint_factor].

    Returns
    -------
    list of float
        Critical E/N values in Townsend, sorted ascending.
        Empty list if no sign change is found in [EN_lo, EN_hi].
    """
    if det_fn is None:
        det_fn = riccati_criterion
    scan_fn = fast_det_fn if fast_det_fn is not None else det_fn

    def f_brentq(EN):
        # NaN → small negative sentinel: NaN always occurs above the critical
        # E/N where the true sign is −1, so -1e-300 is physically correct.
        val = det_fn(EN, pd, mod, p, T)
        return val if np.isfinite(val) else -1e-300

    def _no_root_below(EN_ref):
        """
        True when no sign change lies below EN_ref.

        A warm-start root is the lowest root only if nothing is beneath
        it, and this samples at the same points per decade as the full
        scan, so accepting one is as reliable as running the scan it
        replaces.
        """
        if EN_ref <= EN_lo:
            return True
        per_decade = (n_scan - 1) / math.log10(EN_hi / EN_lo)
        n_guard = max(3, int(math.ceil(per_decade * math.log10(EN_ref / EN_lo))) + 1)
        grid = np.logspace(np.log10(EN_lo), np.log10(EN_ref), n_guard)
        vals = np.array([scan_fn(en, pd, mod, p, T) for en in grid])
        vals = np.where(np.isfinite(vals), vals, 0.0)  # NaN → 0, as in _scan_with
        # The last grid point IS the warm-start root, where det Q is zero to
        # within the root-finder tolerance and its sign is therefore numerical
        # noise.  Comparing it against its neighbour reports a sign change
        # roughly half the time, and a spurious "there is a root below" sends
        # every pd point through the full scan — silently undoing the warm
        # start it was meant to protect.  Test only the intervals that lie
        # strictly below the root; a genuine lower root produces a sign change
        # there, at the same resolution as the scan this replaces.
        below = vals[:-1]
        return not np.any(below[:-1] * below[1:] < 0.0)

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

    # A criterion whose sign is physical (Riccati: g > 0 below inception for
    # every mechanism) tells at a single point whether the lowest root lies
    # below the scan range.  Then no scan above EN_lo can find it; search
    # down instead, where a strongly non-uniform gap puts its gap-averaged
    # inception field.
    if first_only and _has_physical_sign(det_fn):
        f_lo = f_brentq(EN_lo)
        if f_lo < 0.0:
            EN_b = EN_lo
            while EN_b > _EN_FLOOR:
                EN_a = max(EN_b / 2.0, _EN_FLOOR)
                if f_brentq(EN_a) > 0.0:
                    root = scipy.optimize.brentq(
                        f_brentq, EN_a, EN_b, xtol=1e-8, rtol=1e-14
                    )
                    return [float(root)]
                EN_b = EN_a
            return []

    EN_scan = np.logspace(np.log10(EN_lo), np.log10(EN_hi), n_scan)

    def _scan_with(sfn):
        """Coarse sign-change scan using sfn; Brentq refinement always uses det_fn.

        Returns ``(roots, unconfirmed)``.  ``unconfirmed`` is True when a
        bracket reported by *sfn* could not be confirmed with the full
        determinant, which matters when *sfn* is a proxy: the caller must
        then not trust the roots that were found above it.
        """

        def f_s(EN):
            val = sfn(EN, pd, mod, p, T)
            return val if np.isfinite(val) else 0.0

        # Evaluated lazily, bottom up: in first_only mode the scan stops at the
        # first confirmed root instead of paying for every point above it,
        # which matters when sfn is the full criterion.
        D = [f_s(EN_scan[0])]
        local_roots = []
        unconfirmed = False
        for i in range(len(EN_scan) - 1):
            D.append(f_s(EN_scan[i + 1]))
            if D[i] * D[i + 1] < 0.0:
                EN_a, EN_b = EN_scan[i], EN_scan[i + 1]
                fa, fb = f_brentq(EN_a), f_brentq(EN_b)
                if fa * fb >= 0.0:
                    # The proxy bracket does not hold for the full determinant.
                    # The two disagree on where a root sits by as much as tens
                    # of per cent for a strongly non-uniform field, which is
                    # many scan steps, so the confirmation window is a fixed
                    # factor rather than a step count.
                    lo = max(EN_lo, EN_a / _CONFIRM_FACTOR)
                    hi = min(EN_hi, EN_b * _CONFIRM_FACTOR)
                    fine = np.logspace(np.log10(lo), np.log10(hi), 40)
                    fine_v = [f_brentq(en) for en in fine]
                    found = False
                    for k in range(len(fine) - 1):
                        if fine_v[k] * fine_v[k + 1] < 0.0:
                            EN_a, EN_b = fine[k], fine[k + 1]
                            fa, fb = fine_v[k], fine_v[k + 1]
                            found = True
                            break
                    if not found:
                        # Either the proxy invented a root, or the real one lies
                        # outside the window.  Both are only resolvable with the
                        # full determinant, and skipping on would silently
                        # promote a higher root to "lowest".
                        unconfirmed = True
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
        return local_roots, unconfirmed

    roots, unconfirmed = _scan_with(scan_fn)
    if fast_det_fn is not None and (not roots or unconfirmed):
        # Either the fast scan found no sign change at all, or one of its
        # brackets could not be confirmed.  In both cases only the full
        # determinant can settle it; in the second case the roots above the
        # unconfirmed bracket must be discarded, since the lowest of them is
        # not necessarily the lowest root.
        # med_det_fn is intentionally not used here: it almost never brackets
        # roots that fast_det missed, so using it only adds scan overhead before
        # the full fallback (which is always needed anyway).
        roots, _ = _scan_with(det_fn)
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
    progress=None,
    jobs=1,
    verify=None,
):
    """
    Compute inception-curve branches across a p·d sweep.

    Roots are found at each pd point and assigned to branches by
    continuity in log(E/N); unmatched roots start new branches.

    Parameters
    ----------
    pd_arr : numpy.ndarray, shape (M,)
        Array of p*d values in bar·m.
    mod : Mechanism
        Loaded mechanism.
    p : float or numpy.ndarray, shape (M,)
        Gas pressure in bar: a scalar for fixed-pressure mode, or one value
        per pd point for fixed-distance mode.
    T : float
        Gas temperature in Kelvin.
    all_branches : bool
        If False (default), only the lowest-E/N root (branch 0) is kept at
        each pd point.  If True, all roots are collected.
    progress : callable or None
        Called after each pd point as
        ``progress(i, n, pd, p, roots, seconds, check)``, with the E/N roots
        found there (ascending, possibly empty), the time it took and the
        value *verify* returned (None without it), so a caller can report
        the sweep as it runs.  With ``jobs > 1`` the calls come
        in completion order, not pd order.
    jobs : int
        Worker processes for the root searches (see
        :func:`incept1d.parallel.parallel_map`).  The pd points are split
        into contiguous blocks solved concurrently, each in order with
        warm starts, and assigned to branches afterwards in pd order.
    verify : callable or None
        ``verify(pd, p, roots) -> check``, evaluated next to the root search
        (in the same worker) and handed to *progress*; e.g. an independent
        check of the lowest root.

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
        det_fn = riccati_criterion
    p_arr = np.full_like(pd_arr, p) if np.ndim(p) == 0 else np.asarray(p, dtype=float)
    branches = []
    branch_last_logEN = []  # last log10(EN) for each branch, for continuity tracking

    def solve(i, previous=None):
        # previous: this block's result at the pd point before, whose roots
        # warm-start this one (see parallel_map).
        hints = list(previous[0]) if previous is not None else None
        t_start = time.perf_counter()
        roots = find_all_breakdown_EN(
            pd_arr[i],
            mod,
            p_arr[i],
            T,
            first_only=not all_branches,
            det_fn=det_fn,
            fast_det_fn=fast_det_fn,
            med_det_fn=med_det_fn,
            EN_hints=hints,
        )
        roots = sorted(float(r) for r in roots)
        seconds = time.perf_counter() - t_start
        check = verify(pd_arr[i], p_arr[i], roots) if verify is not None else None
        return roots, seconds, check

    def report(i, result):
        if progress is not None:
            progress(i, len(pd_arr), pd_arr[i], p_arr[i], *result)

    # Contiguous blocks of pd points, each solved in order with warm starts;
    # jobs = 1 is a single block, i.e. the plain sequential sweep.
    solved = parallel_map(solve, len(pd_arr), jobs, on_result=report)

    for i, (pd_i, p_i) in enumerate(zip(pd_arr, p_arr)):
        roots = solved[i][0]
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
