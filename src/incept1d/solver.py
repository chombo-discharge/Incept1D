# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Core solver: the augmented ODE, its propagators and the boundary
determinant det Q(λ).

The augmented state is θ = (y, Ψ⁺, Ψ⁻), where y = V n are the species
fluxes and Ψ± the forward and backward photon fluxes.  Integrating
∂_x θ = A_aug θ across the gap gives the propagator M(d), from which the
boundary-condition matrix Q is assembled; inception is det Q(λ) = 0 at
λ = 0.

The model and its boundary conditions are derived in the Theory chapter
of the documentation (``docs/source/theory/``, eq. ``eq_augmented_ode``
and ``eq_det_criterion``); how the propagator and determinant are
actually evaluated is in ``docs/source/numerics/``.  The root of det Q in
E/N at fixed p·d is found in :mod:`incept1d.inception`.
"""

import functools
import itertools
import math

import numpy as np
import scipy.integrate
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

from incept1d.constants import kB as _kB, c_light as _C_LIGHT  # noqa: F401
from incept1d.fields import FieldDistribution  # noqa: F401


def _build_A_aug(EN, d, mod, p, T, lam=0.0):
    """Build the augmented ODE matrix for reduced field EN and gap length d.

    Every photon group is propagated explicitly, however optically thick.
    Folding a thick group into A as the local source 2 b_j c_jᵀ / κ_j keeps
    the photoionization it produces but puts each photoelectron at the
    point of emission, and so removes the upstream seeding that makes
    photoionization a feedback loop: in an inhomogeneous gap the ionization
    zone can be as thin as 1/κ even when κ d ≫ 1.  The stiffness of the
    e^{±κx} photon modes is handled by the determinant instead (see
    :func:`_det_Q`).

    Parameters
    ----------
    EN, d, p, T : float
        Reduced field (Td), gap length (m), pressure (bar), temperature (K).
    mod : Mechanism
        Polarity-resolved mechanism.
    lam : float, optional
        Temporal growth rate λ in s⁻¹; 0 is the inception threshold.

    Returns
    -------
    tuple
        ``(A_aug, N_gamma_eff, aug_mask, n_aug)``.  ``aug_mask`` selects the
        explicitly propagated photon groups — all of them — and is None when
        N_γ = 0.
    """
    n = len(mod.SPECIES)

    R = mod.get_R(EN, p, T)
    V = mod.get_V(EN, p, T)
    # V is always diagonal; multiply R and C by 1/v column-wise via broadcasting
    # instead of forming the full n×n inverse matrix.  If a mechanism ever returns
    # non-diagonal V this will give wrong results — add np.linalg.inv fallback then.
    v_inv = 1.0 / np.diag(V)
    A = R * v_inv[np.newaxis, :]
    if lam != 0.0:
        A = A - lam * np.diag(v_inv)  # (R − λI) V⁻¹

    B = mod.get_B(EN, p, T)  # (N_s, N_gamma)
    C = mod.get_C(EN, p, T) * v_inv[np.newaxis, :]  # (N_gamma, N_s)
    kappa = mod.get_kappa(p, T)  # (N_gamma,)
    if lam != 0.0:
        kappa = kappa + lam / _C_LIGHT  # κ + λ/c
    N_gamma = B.shape[1]

    N_gamma_eff = N_gamma
    aug_mask = np.ones(N_gamma, dtype=bool) if N_gamma > 0 else None
    B_aug, C_aug_ph, kappa_aug = B, C, kappa

    if N_gamma_eff == 0:
        A_aug = A
        n_aug = n
    else:
        D_mat = np.diag(kappa_aug)
        z_gg = np.zeros((N_gamma_eff, N_gamma_eff))
        A_aug = np.block(
            [
                [A, B_aug, B_aug],
                [C_aug_ph, -D_mat, z_gg],
                [-C_aug_ph, z_gg, D_mat],
            ]
        )
        n_aug = n + 2 * N_gamma_eff

    return A_aug, N_gamma_eff, aug_mask, n_aug


def _expm_shifted(M):
    """Matrix exponential of M, shifted so that it cannot overflow.

    M is the pre-formed argument to expm: A_aug * h for the midpoint rule,
    or the full Ω for a Magnus expansion.

    The shift factor is deliberately *not* restored.  It is a positive
    scalar common to every row of Q_d, so it changes the magnitude of
    det Q but not its sign, which is all the root finder uses.
    """
    # The shift is the true largest real eigenvalue rather than a cheaper
    # Gershgorin bound: a bound over-shrinks M, and the error compounds
    # across the N_steps products until Q trips the conditioning guard.
    real_eigs = np.real(np.linalg.eigvals(M))
    lam_max = float(np.max(real_eigs))
    lam = lam_max if lam_max > 0.0 else 0.0
    return scipy.linalg.expm(M - lam * np.eye(M.shape[0]))


DX_N_MIN_DEFAULT = 5
# Measured on sphere-plane gaps (R = 0.32 and 5 mm, both polarities): the
# old 200 / 0.03 put the inception field 0.3–2.3 % off the converged value;
# 1000 / 1e-3 keeps it within 0.1 % at 3–4× the cost.
DX_N_MAX_DEFAULT = 1000
DX_TOL_DEFAULT = 1e-3


def parse_dx_spec(tokens, parser=None):
    """
    Parse a --dx token list into (N_min, N_max, tol).

    Tokens are all optional and fall back to the module defaults.
    N_min = N_max disables adaptation and gives a uniform grid.

    Parameters
    ----------
    tokens : list of str or None
        Raw argparse token list, e.g. ['5', '1000', '1e-3'].
    parser : argparse.ArgumentParser or None
        If given, validation errors are routed through parser.error().
    """

    def _err(msg):
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    N_min, N_max, tol = DX_N_MIN_DEFAULT, DX_N_MAX_DEFAULT, DX_TOL_DEFAULT
    if tokens:
        try:
            if len(tokens) >= 1:
                N_min = int(tokens[0])
            if len(tokens) >= 2:
                N_max = int(tokens[1])
            if len(tokens) >= 3:
                tol = float(tokens[2])
        except ValueError as exc:
            _err(f"--dx: {exc}")
    if N_min < 1:
        _err(f"--dx: N_min must be ≥ 1, got {N_min}")
    if N_max < N_min:
        _err(f"--dx: N_max ({N_max}) must be ≥ N_min ({N_min})")
    if not 0.0 < tol < 1.0:
        _err(f"--dx: tol must be in (0, 1), got {tol}")
    return N_min, N_max, tol


def _midpoint_exponent(A_func, x_lo, x_hi):
    """Exponent of the midpoint rule on [x_lo, x_hi]: A(x_mid) · h."""
    return A_func(0.5 * (x_lo + x_hi)) * (x_hi - x_lo)


def midpoint_propagator(A_func, x_lo, x_hi):
    """
    Midpoint-rule propagator: expm(A(x_mid) * h).

    Parameters
    ----------
    A_func : callable
        ``A_func(x) -> A_aug``.  x is a code coordinate in metres; the
        polarity flip and field evaluation are the caller's business.
    x_lo, x_hi : float
        Subinterval bounds in metres.

    Returns
    -------
    ndarray
        Propagator for the subinterval.
    """
    return _expm_shifted(_midpoint_exponent(A_func, x_lo, x_hi))


def _adaptive_midpoint_segment(
    A_func, x_lo, x_hi, tol, max_depth, P_coarse=None, _level=0, _diag=None
):
    """
    Midpoint exponents for [x_lo, x_hi], refined by recursive step halving.

    Parameters
    ----------
    A_func : callable
        Same signature as for midpoint_propagator.
    x_lo, x_hi : float
        Subinterval bounds in code coordinates (metres).
    tol : float
        Relative Frobenius error at which a segment is accepted.
    max_depth : int
        Remaining recursion depth.  At 0 the coarse estimate is returned
        without a halving check.
    P_coarse : tuple or None
        The parent's half-step ``(exponent, propagator)``, which is this
        call's coarse estimate.  Passing it saves one expm; None at the top
        level.
    _level : int
        Recursion depth, for diagnostics only.
    _diag : list or None
        If given, one dict per halving check is appended:
        ``{'x_lo', 'x_hi', 'level', 'err', 'accepted'}``.  Leaf calls at
        max_depth 0 attempt no halving and are not recorded.

    Returns
    -------
    list of ndarray
        Step exponents covering [x_lo, x_hi] from cathode to anode; the
        propagator is the ordered product of their exponentials.  They are
        returned rather than multiplied so that :func:`_det_Q` can also carry
        them across the gap without forming the product, which is what
        loses the subdominant modes (see :func:`_propagate_compound`).
    """
    if P_coarse is None:
        Om = _midpoint_exponent(A_func, x_lo, x_hi)
        P_coarse = (Om, _expm_shifted(Om))
    if max_depth == 0:
        return [P_coarse[0]]

    x_mid = 0.5 * (x_lo + x_hi)
    Om_l = _midpoint_exponent(A_func, x_lo, x_mid)
    Om_r = _midpoint_exponent(A_func, x_mid, x_hi)
    P_l, P_r = _expm_shifted(Om_l), _expm_shifted(Om_r)
    P_fine = P_r @ P_l

    # The accuracy test compares whole propagators.  They are dominated by
    # the leading mode, which is what the step size has to resolve; the
    # subdominant modes are kept separately, by the subspace propagation.
    norm_ref = np.linalg.norm(P_fine, "fro")
    if norm_ref < 1e-300:
        norm_ref = 1.0
    err = np.linalg.norm(P_fine - P_coarse[1], "fro") / norm_ref

    accepted = err <= tol
    if _diag is not None:
        _diag.append(
            {
                "x_lo": x_lo,
                "x_hi": x_hi,
                "level": _level,
                "err": err,
                "accepted": accepted,
            }
        )

    if accepted:
        return [Om_l, Om_r]

    # Error too large — refine each half, passing the already-computed half-step
    # propagators as the coarse estimates for the sub-calls.
    return _adaptive_midpoint_segment(
        A_func, x_lo, x_mid, tol, max_depth - 1, (Om_l, P_l), _level + 1, _diag
    ) + _adaptive_midpoint_segment(
        A_func, x_mid, x_hi, tol, max_depth - 1, (Om_r, P_r), _level + 1, _diag
    )


def _magnus2_exponent(A_func, x_lo, x_hi):
    """
    Second-order Magnus exponent Ω on [x_lo, x_hi].

    Falls back to the midpoint exponent where the Magnus series is not
    guaranteed to converge; see :func:`magnus2_propagator`.
    """
    h = x_hi - x_lo
    mid = 0.5 * (x_lo + x_hi)
    offset = h / (2.0 * math.sqrt(3.0))
    A1 = A_func(mid - offset)
    A2 = A_func(mid + offset)
    Omega1 = 0.5 * h * (A1 + A2)
    # When ‖Ω₁‖_F > π the Magnus series is not guaranteed to converge; the
    # commutator correction would dominate and produce a wildly wrong exponent.
    # Fall back to the midpoint rule, which uses expm(A_mid·h) directly and
    # remains accurate for any step size via scipy's scaling-and-squaring.
    if np.linalg.norm(Omega1, "fro") > math.pi:
        return _midpoint_exponent(A_func, x_lo, x_hi)
    Omega2 = (math.sqrt(3.0) * h**2 / 12.0) * (A2 @ A1 - A1 @ A2)
    return Omega1 + Omega2


def magnus2_propagator(A_func, x_lo, x_hi):
    """
    Second-order Magnus propagator, with two-point Gauss-Legendre quadrature.

    Where the Magnus series is not guaranteed to converge (∫‖A‖ dx > π) the
    interval falls back to the midpoint rule, so the caller does not have
    to pick a grid fine enough for large p·d.

    Parameters
    ----------
    A_func : callable
        A_func(x) → A_aug (ndarray).  x is a code-coordinate in metres.
    x_lo, x_hi : float
        Subinterval bounds in code coordinates (metres).

    Returns
    -------
    ndarray
        Propagator matrix P ≈ expm(∫_{x_lo}^{x_hi} A dx).
    """
    return _expm_shifted(_magnus2_exponent(A_func, x_lo, x_hi))


# Largest condition number allowed for one substep of the compound
# propagation.  The k-th compound of a propagator P has condition number up
# to cond(P)^k, so substeps are sized by cond(P) ≤ _MAX_STEP_COND^(1/k).  The
# bound is on the condition number, not on the eigenvalue spread: A_aug is
# strongly non-normal (electron and ion speeds differ by orders of magnitude),
# and a step whose eigenvalues span a few e-folds can have cond(P) ~ 1e18.
_MAX_STEP_COND = 1e8


@functools.lru_cache(maxsize=None)
def _compound_index(n, k):
    """
    Index tables for the k-th additive compound of an n×n matrix.

    Returns ``(subsets, position, rows, cols, a, b, sign)``: the k-subsets of
    range(n) in lexicographic order, their positions, and one entry per
    nonzero of the compound, ``K[rows, cols] += sign * Ω[a, b]``.
    """
    subsets = list(itertools.combinations(range(n), k))
    position = {s: i for i, s in enumerate(subsets)}
    rows, cols, a_, b_, sign = [], [], [], [], []
    for subset in subsets:
        for ia, a in enumerate(subset):
            rest = subset[:ia] + subset[ia + 1 :]
            for b in range(n):
                if b in rest:
                    continue
                J = tuple(sorted(rest + (b,)))
                rows.append(position[subset])
                cols.append(position[J])
                a_.append(a)
                b_.append(b)
                sign.append((-1.0) ** (ia + J.index(b)))
    return (
        subsets,
        position,
        np.array(rows),
        np.array(cols),
        np.array(a_),
        np.array(b_),
        np.array(sign),
    )


def _additive_compound(Omega, k, sparse=False):
    """
    The k-th additive compound Ω^[k], the generator of ∧^k exp(Ω).

    exp(Ω^[k]) acts on the Plücker coordinates (the k×k minors) of an
    n×k matrix exactly as exp(Ω) acts on the matrix itself.  Each row has at
    most 1 + k(n − k) nonzeros; with *sparse* it is returned in CSR form.
    """
    _, _, rows, cols, a, b, sign = _compound_index(Omega.shape[0], k)
    size = math.comb(Omega.shape[0], k)
    if sparse:
        # Duplicate (row, col) pairs — the diagonal — are summed on conversion.
        return scipy.sparse.coo_matrix(
            (sign * Omega[a, b], (rows, cols)), shape=(size, size)
        ).tocsr()
    K = np.zeros((size, size))
    np.add.at(K, (rows, cols), sign * Omega[a, b])
    return K


# Above this many Plücker coordinates the compound propagator is applied to
# the coordinate vector (scipy.sparse.linalg.expm_multiply) instead of being
# formed: C(12, 5) = 792 for three explicit photon groups, where a dense
# 792×792 exponential per step dominated the run time.  Like the dense
# exponential, the action only follows coupling paths of Ω^[k], so the
# structural zeros that keep small coordinates accurate are preserved.
_DENSE_COMPOUND_MAX = 120


def _propagate_compound(exponents, Y):
    """
    Carry the Plücker coordinates of span(Y) across the gap.

    Forming the propagator M = ∏ exp(Ω_i) and then selecting rows of it
    fails once the k modes the anode rows see differ in growth by more than
    about 35 e-folds: the rows S M are then parallel to machine precision.
    With λ > 0 that happens at once, because λ V⁻¹ separates the slow ion
    modes by thousands of e-folds.

    Re-orthonormalising the propagated subspace (Godunov–Conte) is not a
    cure here.  Orthogonalisation mixes all components, so a component that
    is transiently tiny — the electrons, attached in the low-field region —
    is buried under the rounding error of the large ion components, and the
    avalanche downstream amplifies that error to O(1).  The k×k minors of
    the propagated n×k matrix obey a *linear* equation, dp/dx = A^[k] p,
    whose coefficients keep every structural zero of A (positive ions feed
    no other species), so each minor is carried with its own relative
    accuracy and the one det Q needs is read off at the anode without
    cancellation (the compound-matrix method).

    Parameters
    ----------
    exponents : list of ndarray
        Step exponents Ω_i, cathode to anode.
    Y : ndarray
        (n_aug, k) matrix whose span is propagated.

    Returns
    -------
    p : ndarray
        Plücker coordinates of ∏ exp(Ω_i) Y, scaled so the largest ``abs(p)`` is 1, ordered
        as :func:`_compound_index`.
    log_scale : float
        log of the factor removed: the true coordinates are p · e^log_scale.
    """
    n, k = Y.shape
    subsets = _compound_index(n, k)[0]
    p = np.array([np.linalg.det(Y[list(rows), :]) for rows in subsets])
    log_scale = math.log(np.abs(p).max())
    p = p / np.abs(p).max()
    for Omega in exponents:
        # Size substeps so the compound propagator stays well conditioned.
        # Splitting is exact: exp(Ω) = exp(Ω/m)^m.
        mu = np.sort(np.real(np.linalg.eigvals(Omega)))[::-1]
        m = max(1, math.ceil((mu[0] - mu[-1]) / math.log(_MAX_STEP_COND) * k))
        cond_max = _MAX_STEP_COND ** (1.0 / k)
        while np.linalg.cond(_expm_shifted(Omega / m)) > cond_max and m < 2**20:
            m *= 2
        # Shift by the largest eigenvalue of Ω^[k], the sum of the k largest
        # of Ω, so that the compound propagator cannot overflow.
        shift = float(mu[:k].sum()) / m
        if len(subsets) <= _DENSE_COMPOUND_MAX:
            Pk = scipy.linalg.expm(
                _additive_compound(Omega / m, k) - shift * np.eye(len(subsets))
            )
            step = Pk.__matmul__
        else:
            Kk = _additive_compound(Omega / m, k, sparse=True) - shift * (
                scipy.sparse.identity(len(subsets), format="csr")
            )
            step = functools.partial(scipy.sparse.linalg.expm_multiply, Kk)
        done = 0
        while done < m:
            q = step(p)
            scale = np.abs(q).max()
            log_scale += shift + math.log(scale)
            q = q / scale
            done += 1
            # Once p is an eigenvector of Pk, each remaining substep only
            # rescales it by the same factor; skipping them is exact and is
            # what keeps a uniform gap at large λ from costing m products.
            # The test must be relative, coordinate by coordinate: a
            # coordinate that is transiently tiny (electrons attached in a
            # low-field region) can still be changing by orders of magnitude
            # while its absolute change is far below any tolerance, and the
            # avalanche downstream would amplify it back to O(1).
            nz = p != 0.0
            if (
                done < m
                and not np.any(q[~nz])
                and np.max(np.abs(q[nz] / p[nz] - 1.0)) < 1e-12
            ):
                log_scale += (m - done) * (shift + math.log(scale))
                p = q
                break
            p = q
    return p, log_scale


def _log_row_norm(row, exponents):
    """
    log ‖row · ∏ exp(Ω_i)‖, carried anode to cathode with a running scale.

    Each exponent is split into substeps spanning at most 30 e-folds, so no
    entry of a substep propagator underflows; the shifts that keep them from
    overflowing are added back, so the result is the true norm.
    """
    r, log_r = np.asarray(row, dtype=float).copy(), 0.0
    n = r.shape[0]
    for Omega in reversed(exponents):
        mu = np.real(np.linalg.eigvals(Omega))
        m = max(1, math.ceil((mu.max() - mu.min()) / 30.0))
        shift = float(mu.max()) / m
        P = scipy.linalg.expm(Omega / m - shift * np.eye(n))
        done = 0
        while done < m:
            q = r @ P
            scale = np.abs(q).max()
            log_r += shift + math.log(scale)
            q = q / scale
            done += 1
            nz = r != 0.0
            if (
                done < m
                and not np.any(q[~nz])
                and np.max(np.abs(q[nz] / r[nz] - 1.0)) < 1e-12
            ):
                log_r += (m - done) * (shift + math.log(scale))
                r = q
                break
            r = q
    return log_r + math.log(np.linalg.norm(r))


def _boundary_rows(EN_cathode, mod, p, T, N_gamma_eff, aug_mask, n_aug):
    """
    The two halves of the boundary-condition matrix Q.

    Returns ``(Q0, S)``: Q0 holds the cathode conditions on θ(0) (electron
    emission, no negative ions, no forward photons), and S selects the
    components that vanish at the anode (positive ions, backward photons),
    so that Q = [Q0; S M].  *mod* must already be polarity-resolved.
    """
    n = len(mod.SPECIES)

    Pi_e = mod.get_Pi_e()
    Pi_plus = mod.get_Pi_plus()
    Pi_minus = mod.get_Pi_minus()
    gplus = mod.get_gamma_plus(EN_cathode, p, T)

    if N_gamma_eff == 0:
        Q0_e = Pi_e + gplus[np.newaxis, :] @ Pi_plus
        return np.vstack([Q0_e, Pi_minus]), Pi_plus

    z_gamma = np.zeros((Pi_e.shape[0], 2 * N_gamma_eff))
    z_plus = np.zeros((Pi_plus.shape[0], 2 * N_gamma_eff))
    z_minus = np.zeros((Pi_minus.shape[0], 2 * N_gamma_eff))
    Pi_e_aug = np.hstack([Pi_e, z_gamma])
    Pi_plus_aug = np.hstack([Pi_plus, z_plus])
    Pi_minus_aug = np.hstack([Pi_minus, z_minus])

    Pi_fwd = np.zeros((N_gamma_eff, n_aug))
    Pi_bck = np.zeros((N_gamma_eff, n_aug))
    for j in range(N_gamma_eff):
        Pi_fwd[j, n + j] = 1.0
        Pi_bck[j, n + N_gamma_eff + j] = 1.0

    g_Psi = mod.get_gamma_Psi(EN_cathode, p, T)[aug_mask]

    Q0_e = Pi_e_aug + gplus[np.newaxis, :] @ Pi_plus_aug - g_Psi[np.newaxis, :] @ Pi_bck
    Q0 = np.vstack([Q0_e, Pi_minus_aug, Pi_fwd])
    return Q0, np.vstack([Pi_plus_aug, Pi_bck])


def _det_Q_norm(Q):
    """
    det of Q with every row scaled to unit norm; NaN if ill-conditioned.

    Row scaling removes the arbitrary scale of M relative to Q0 and leaves
    the sign unchanged.  Must be called inside a numpy errstate block.
    """
    if np.any(~np.isfinite(Q)):
        return np.nan
    row_norms = np.linalg.norm(Q, axis=1, keepdims=True)
    row_norms = np.where(row_norms < 1e-300, 1.0, row_norms)
    Q_norm = Q / row_norms
    sign, logabsdet = np.linalg.slogdet(Q_norm)
    if np.linalg.cond(Q_norm) > 1e14:
        return np.nan
    if not np.isfinite(logabsdet):
        return np.nan
    if sign == 0:
        return 0.0
    return float(sign) * min(np.exp(logabsdet), 1e300)


def _det_Q_compound(exponents, Q0, S):
    """
    Row-normalised det [Q0; S M] by the compound-matrix method.

    Returns the same value as forming M and taking the row-normalised
    determinant (:func:`_det_Q_norm`), but without forming M, so it stays
    exact where the rows S M are parallel to machine precision.  See
    :func:`_det_Q` for the identity it rests on.
    """
    if S.shape[0] + Q0.shape[0] != S.shape[1]:
        raise ValueError(
            f"Q is {Q0.shape[0] + S.shape[0]}x{S.shape[1]}: every species must be "
            "selected by exactly one of Pi_e, Pi_plus, Pi_minus"
        )
    G = Q0 @ Q0.T
    # N0 is the anode-selected coordinates projected onto null(Q0): it
    # varies continuously with E/N, so det T cannot change sign between
    # evaluations.
    N0 = S.T - Q0.T @ np.linalg.solve(G, Q0 @ S.T)
    s_G, ld_G = np.linalg.slogdet(G)
    s_T, ld_T = np.linalg.slogdet(np.hstack([Q0.T, N0]))

    pl, log_scale = _propagate_compound(exponents, N0)
    subsets, position = _compound_index(S.shape[1], N0.shape[1])[:2]
    selected = [int(np.argmax(row)) for row in S]
    # det(S Y) is the minor on rows `selected`, in S's row order.
    order = np.argsort(selected)
    perm_sign = np.linalg.det(np.eye(len(selected))[order])
    coord = perm_sign * pl[position[tuple(sorted(selected))]]
    if coord == 0.0 or not np.isfinite(log_scale):
        return 0.0 if coord == 0.0 else np.nan

    # Row norms of Q.  The rows of S M can be hundreds of e-folds below
    # M's largest entry, so each is carried separately with its own scale.
    log_rows = float(np.sum(np.log(np.linalg.norm(Q0, axis=1))))
    for row in S:
        log_rows += _log_row_norm(row, exponents)
    log_det = ld_G - ld_T + math.log(abs(coord)) + log_scale - log_rows
    sign = s_G * s_T * np.sign(coord)
    # Floor the magnitude: an underflow to 0.0 would read as an exact root,
    # and anything below 1e-290 as the root finders' NaN sentinel.
    return float(sign) * math.exp(min(max(log_det, -644.0), 690.0))


def _det_Q(
    exponents, EN_cathode, mod, p, T, N_gamma_eff, aug_mask, n_aug, resolve=False
):
    """
    Row-normalised det Q for the propagator M = ∏ exp(Ω_i).

    Q = [Q0; S M] (eq. ``eq_Q_system``).  The direct evaluation — form M,
    assemble Q, take the determinant — is exact whenever Q is well
    conditioned, and is cheap, so it is tried first.  Where it is not (the
    anode rows S M parallel to machine precision, typically above inception
    or at λ > 0) the same number is computed by the compound-matrix method
    instead: with N0 spanning null(Q0) and T = [Q0ᵀ | N0],

        det Q · det T = det(Q0 Q0ᵀ) · det(S M N0),

    and det(S M N0) is one Plücker coordinate of M N0, which
    :func:`_propagate_compound` carries across the gap without forming M.
    Everything is kept in logarithms, then divided by the row norms of Q, so
    both routes return the same value.

    The compound route costs 10²–10⁴ times the direct one, so it runs only
    when *resolve* is True.  Otherwise a singular Q returns NaN, which the
    inception scans read as "above threshold" (true at λ = 0, where the
    anode rows become parallel only above inception).

    Must be called inside a numpy errstate(over='ignore', invalid='ignore')
    block.  Returns NaN if the requested routes cannot evaluate it.
    """
    Q0, S = _boundary_rows(EN_cathode, mod, p, T, N_gamma_eff, aug_mask, n_aug)

    # Direct route.  M is carried with a running scale so it cannot overflow;
    # a positive scalar on M does not change the row-normalised determinant.
    M = np.eye(n_aug)
    log_M = 0.0
    for Omega in exponents:
        M = _expm_shifted(Omega) @ M
        scale = np.abs(M).max()
        if not np.isfinite(scale) or scale == 0.0:
            return np.nan
        M /= scale
        log_M += math.log(scale)
    direct = _det_Q_norm(np.vstack([Q0, S @ M]))
    if np.isfinite(direct) or not resolve:
        return direct

    return _det_Q_compound(exponents, Q0, S)


def polarities_equivalent(mod, field_dist):
    """
    True when positive and negative polarity must give the same det Q.

    Swapping the electrodes is a reflection of the gap, which leaves the
    answer unchanged only if the field *and* the electrode surfaces are
    symmetric.  A configuration with per-polarity overrides describes two
    different cathodes, so it is asymmetric even in a uniform field.

    Parameters
    ----------
    mod : Mechanism
        Loaded mechanism, before polarity resolution.
    field_dist : FieldDistribution
        Gap geometry.

    Returns
    -------
    bool
        True if one polarity can stand in for the other.
    """
    return field_dist.is_symmetric and not mod.has_polarity_overrides


def inception_det(
    EN_ref,
    pd,
    mod,
    p,
    T,
    field_dist,
    N_min=DX_N_MIN_DEFAULT,
    N_max=DX_N_MAX_DEFAULT,
    tol=DX_TOL_DEFAULT,
    lam=0.0,
    positive_polarity=True,
    propagator=midpoint_propagator,
    diag_list=None,
    resolve=True,
):
    """
    Evaluate det Q(λ) for the inception boundary-value problem.

    Integrates the augmented ODE across the gap, assembles the
    boundary-condition matrix Q and returns its determinant, which vanishes
    at inception.

    A uniform field is a single exact matrix exponential, so the grid and
    propagator options are then ignored and every setting returns the same
    det Q.  For a non-uniform field see ``docs/source/numerics/``.

    Parameters
    ----------
    EN_ref : float
        Uniform-equivalent reduced field in Townsend.
    pd : float
        Product of pressure and gap length in bar·m.
    mod : Mechanism
        Loaded mechanism, with its configuration already baked in.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.
    field_dist : FieldDistribution
        Gap geometry.
    N_min, N_max : int
        Initial number of segments, and the budget of fine steps the
        adaptive refinement may spend.  Equal values disable adaptation.
    tol : float
        Relative Frobenius error at which a segment is accepted.
    lam : float
        Temporal growth rate λ in s⁻¹; 0 is the inception threshold.
    resolve : bool
        Evaluate det Q even where Q is numerically singular, by the
        compound-matrix method (see :func:`_det_Q`), instead of returning
        NaN.  Needed above threshold, at λ > 0, and wherever an optically
        thick photon group makes the propagator stiff; slower, but only
        used where the direct evaluation fails.  False restores NaN there.
    positive_polarity : bool
        True when the ξ = 0 electrode is the anode.  Ignored for the
        symmetric geometries.
    propagator : callable
        :func:`midpoint_propagator` (adaptive halving) or
        :func:`magnus2_propagator` (uniform N_min grid).  It selects the
        quadrature; how the step exponents are combined is up to
        :func:`_det_Q`.

    Returns
    -------
    float
        det Q(λ).  Inception threshold: det Q = 0.
    """
    with np.errstate(over="ignore", invalid="ignore"):
        exponents, EN_cathode, eff, N_gamma_eff, aug_mask, n_aug = _discretise(
            EN_ref,
            pd,
            mod,
            p,
            T,
            field_dist,
            N_min,
            N_max,
            tol,
            lam,
            positive_polarity,
            propagator,
            diag_list,
        )
        return _det_Q(
            exponents,
            EN_cathode,
            eff,
            p,
            T,
            N_gamma_eff,
            aug_mask,
            n_aug,
            resolve,
        )


def _discretise(
    EN_ref,
    pd,
    mod,
    p,
    T,
    field_dist,
    N_min,
    N_max,
    tol,
    lam,
    positive_polarity,
    propagator,
    diag_list,
):
    """
    Step exponents across the gap and the boundary data both criteria need.

    Returns ``(exponents, EN_cathode, eff, N_gamma_eff, aug_mask, n_aug)``:
    the step exponents Ω_i from cathode to anode (one exact exponent for a
    uniform field; the midpoint or Magnus grid otherwise), the reduced field
    at the cathode, the polarity-resolved mechanism, and the photon
    structure of the augmented system.  :func:`inception_det` and
    :func:`riccati_criterion` share it, so they solve the same discrete
    problem on the same grid.
    """
    eff = mod.resolve(positive_polarity)
    d = pd / p

    # Uniform field: A_aug does not depend on x, so the path-ordered product
    # collapses to a single matrix exponential,
    #
    #     M(d) = exp(A_aug · d),
    #
    # which is the *exact* propagator rather than a quadrature of it — there is
    # nothing for the midpoint rule or the adaptive grid to improve on, and the
    # N_min−1 matrix products are pure roundoff.  Magnus2 also reduces to this
    # (its commutator term vanishes for constant A), so the shortcut applies
    # whichever propagator was requested.  The compound route still splits
    # the single exponent A_aug·d into substeps, but that split is exact
    # (exp(Ω) = exp(Ω/m)^m), not a quadrature.
    if field_dist.field_type == "uniform":
        A_aug, N_gamma_eff, aug_mask, n_aug = _build_A_aug(
            EN_ref, d, eff, p, T, lam=lam
        )
        with np.errstate(over="ignore", invalid="ignore"):
            if diag_list is not None:
                diag_list.append(
                    {
                        "segment": 0,
                        "x_lo": 0.0,
                        "x_hi": d,
                        "halvings": [
                            {
                                "x_lo": 0.0,
                                "x_hi": d,
                                "level": 0,
                                "err": 0.0,
                                "accepted": True,
                                "exact": True,
                            }
                        ],
                    }
                )
            return [A_aug * d], EN_ref, eff, N_gamma_eff, aug_mask, n_aug

    f = field_dist.build(d)
    d_step = d / N_min
    max_depth = max(0, int(math.log2(max(1, N_max // N_min))))

    # A_func encapsulates the coordinate flip and all physics context.
    # Code coordinates run 0 (cathode) → d (anode); for positive polarity the
    # formula coordinate is reversed so xi=0 (sphere/anode) maps to x_code=d.
    def A_func(x_code, _f=f, _d=d, _EN=EN_ref, _pos=positive_polarity):
        x_formula = (_d - x_code) if _pos else x_code
        A_aug, *_ = _build_A_aug(_EN * _f(x_formula / _d), _d, eff, p, T, lam=lam)
        return A_aug

    # Extract augmented-system metadata once from the first step midpoint.
    # Photon structure (aug_mask, N_gamma_eff) is assumed constant across the gap.
    xi_mid0 = 0.5 * d_step / d
    xi_first = (1.0 - xi_mid0) if positive_polarity else xi_mid0
    _, N_gamma_eff, aug_mask, n_aug = _build_A_aug(
        EN_ref * f(xi_first), d, eff, p, T, lam=lam
    )
    exponents = []

    with np.errstate(over="ignore", invalid="ignore"):
        if propagator is midpoint_propagator:
            # Adaptive midpoint: each of the N_min segments is refined by step
            # halving until the relative Frobenius error < tol or max_depth is
            # exhausted.  max_depth = 0 (when N_min = N_max) gives a constant
            # uniform grid with no halving attempted.
            for i in range(N_min):
                seg_diag = [] if diag_list is not None else None
                exponents += _adaptive_midpoint_segment(
                    A_func,
                    i * d_step,
                    (i + 1) * d_step,
                    tol,
                    max_depth,
                    _diag=seg_diag,
                )
                if diag_list is not None:
                    diag_list.append(
                        {
                            "segment": i,
                            "x_lo": i * d_step,
                            "x_hi": (i + 1) * d_step,
                            "halvings": seg_diag,
                        }
                    )
        elif propagator is magnus2_propagator:
            for i in range(N_min):
                exponents.append(
                    _magnus2_exponent(A_func, i * d_step, (i + 1) * d_step)
                )
        else:
            raise ValueError(
                "propagator must be midpoint_propagator or magnus2_propagator"
            )

        xi_cathode = 1.0 if positive_polarity else 0.0
        EN_cathode = EN_ref * f(xi_cathode)
        return exponents, EN_cathode, eff, N_gamma_eff, aug_mask, n_aug


# Largest condition number of one Riccati substep's propagator.  The linear
# fractional map tolerates more than the compound route: on the verification
# cases 1e12 leaves every root unchanged to 1e-9, 1e16 moves one by 2.5e-7
# and 1e24 by 0.4 %; 1e12 needs 1.5x fewer substeps than 1e8.
_RICCATI_STEP_COND = 1e12


def _riccati_g(exponents, Q0, S):
    """
    The inception criterion by the Riccati (reflection) method.

    Split θ into the components that travel toward the anode (electrons,
    negative ions, forward photons: f) and toward the cathode (positive ions,
    backward photons: b, the ones S selects).  The anode condition b(d) = 0
    is carried back to the cathode as b = P f, where the reflection operator
    P obeys

        P′ = A_bf + A_bb P − P A_ff − P A_fb P,   P(d) = 0.

    Integrated from the anode, every component travels in its own direction,
    so the modes that make the propagator M useless — e^{+κx} for backward
    photons, e^{+λx/|v₊|} for ions — decay instead of growing.  Each step
    with constant A is described by its reflection and transmission matrices
    (:func:`_slab_scattering`), built by adding–doubling from thin, well
    conditioned substeps, and P is carried across it exactly by

        P ← R_l + T_b P (I − R_r P)⁻¹ T_f,

    so the result is the same discrete problem det Q solves on the same
    grid.  At the cathode the conditions Q0 θ = 0 with θ_b = P θ_f give

        g = det(Q0_f + Q0_b P(0)),

    which is zero exactly where det Q is.  With one electron row it is
    1 − (loop gain): positive below inception, negative above, and of order
    one; P(0) e_e are the ion and backward-photon fluxes returning to the
    cathode per emitted electron.

    P stays finite below inception: it can only blow up where a sub-slab
    [x, d] is self-sustaining without the cathode, which makes the whole gap
    supercritical.  That shows as the loop gain R_r P (or R_r^A R_l^B inside
    a step) reaching spectral radius 1 (see :func:`_spectral_radius`), and
    returns −1.
    """
    n = Q0.shape[1]
    b = [int(np.argmax(row)) for row in S]
    fw = [i for i in range(n) if i not in b]
    if Q0.shape[0] != len(fw):
        raise ValueError(
            f"Q0 has {Q0.shape[0]} rows for {len(fw)} forward components: every "
            "species must be selected by exactly one of Pi_e, Pi_plus, Pi_minus"
        )
    nf = len(fw)
    P = np.zeros((len(b), nf))
    for Omega in reversed(exponents):
        slab = _slab_scattering(Omega, b, fw)
        if slab is None:
            return -1.0
        R_l, R_r, T_f, T_b = slab
        # Riccati update across the step: with b = P f at its anode side,
        # P ← R_l + T_b P (I − R_r P)⁻¹ T_f.  R_r P is the loop gain between
        # this step and everything downstream of it; once its spectral
        # radius reaches 1 that combination sustains itself: above inception.
        loop = R_r @ P
        if _spectral_radius(loop) >= 1.0:
            return -1.0
        P = R_l + T_b @ P @ np.linalg.solve(np.eye(nf) - loop, T_f)
        if not np.all(np.isfinite(P)):
            return -1.0
    return float(np.linalg.det(Q0[:, fw] + Q0[:, b] @ P))


def _slab_scattering(Omega, b, fw):
    """
    Reflection and transmission of one step with exponent Ω.

    For a slab whose ends are related by θ_left = E θ_right, E = exp(−Ω),
    the fluxes leaving it are given by the fluxes entering it:

        f_right = T_f f_left + R_r b_right,
        b_left  = R_l f_left + T_b b_right,

    with T_f = E_ff⁻¹, R_r = −E_ff⁻¹ E_fb, R_l = E_bf E_ff⁻¹ and
    T_b = E_bb − E_bf E_ff⁻¹ E_fb.  A thin substep, small enough for its
    propagator to be well conditioned, is formed directly; the step is then
    built from m such substeps by repeated squaring with the Redheffer star
    product (the adding–doubling method), in log₂ m combinations instead of
    m.  For a system whose couplings are all non-negative sources these
    matrices are non-negative too, so composing them involves no
    cancellation.

    Returns ``(R_l, R_r, T_f, T_b)``, or None if part of the step is already
    self-sustaining on its own (above inception).
    """
    mu = np.real(np.linalg.eigvals(Omega))
    m = max(1, math.ceil((mu.max() - mu.min()) / math.log(_RICCATI_STEP_COND)))
    E = _expm_shifted(-Omega / m)
    while np.linalg.cond(E) > _RICCATI_STEP_COND and m < 2**40:
        m *= 2
        E = _expm_shifted(-Omega / m)
    # _expm_shifted dropped the factor e^shift from E; restore it in the
    # transmissions (the reflections are ratios and do not see it).
    shift = max(0.0, float(np.max(np.real(np.linalg.eigvals(-Omega / m)))))
    E_ff_inv = np.linalg.inv(E[np.ix_(fw, fw)])
    E_fb, E_bf, E_bb = E[np.ix_(fw, b)], E[np.ix_(b, fw)], E[np.ix_(b, b)]
    R_l = E_bf @ E_ff_inv
    R_r = -E_ff_inv @ E_fb
    T_f = E_ff_inv * math.exp(-shift)
    T_b = (E_bb - R_l @ E_fb) * math.exp(shift)
    unit = (R_l, R_r, T_f, T_b)
    if not all(np.all(np.isfinite(x)) for x in unit):
        return None  # overflowed: see _spectral_radius

    result = None
    while m:
        if m & 1:
            result = unit if result is None else _star(result, unit)
            if result is None:
                return None
        m >>= 1
        if m:
            unit = _star(unit, unit)
            if unit is None:
                return None
    return result


def _spectral_radius(M):
    """
    Largest |eigenvalue| of a loop-gain matrix.

    The test for "self-sustaining" is ρ(loop) ≥ 1, not a sign change of
    det(I − loop): for a system whose couplings are all non-negative
    sources, I − K has a non-negative inverse exactly when ρ(K) < 1, and ρ
    stays above 1 once it has crossed, whereas the determinant changes sign
    back when a second mode crosses too.  Far above inception that happens
    within one step, and a determinant test then reported an isolated,
    spurious "below inception" (g ≈ +1e6 at 153 Td in a coaxial gap whose
    inception field is 44 Td).
    """
    if M.size == 0:
        return 0.0
    if not np.all(np.isfinite(M)):
        # Overflowed: the gain is astronomically large, i.e. above threshold.
        # (Only a system with no feedback at all, whose loop gain is exactly
        # zero, could overflow its transmissions without that being so, and
        # such a system has no inception to find.)
        return np.inf
    return float(np.max(np.abs(np.linalg.eigvals(M))))


def _star(A, B):
    """
    Redheffer star product: slab A (cathode side) followed by slab B.

    Returns None if the pair is self-sustaining: the loop gain R_r^A R_l^B
    between them has spectral radius 1 or more.
    """
    Rl_a, Rr_a, Tf_a, Tb_a = A
    Rl_b, Rr_b, Tf_b, Tb_b = B
    loop = Rr_a @ Rl_b
    if _spectral_radius(loop) >= 1.0:
        return None
    K = np.eye(Tf_a.shape[0]) - loop
    K_Tf = np.linalg.solve(K, Tf_a)
    K_Rr = np.linalg.solve(K, Rr_a)
    out = (
        Rl_a + Tb_a @ Rl_b @ K_Tf,
        Rr_b + Tf_b @ K_Rr @ Tb_b,
        Tf_b @ K_Tf,
        Tb_a @ (Tb_b + Rl_b @ K_Rr @ Tb_b),
    )
    if not all(np.all(np.isfinite(x)) for x in out):
        return None  # overflowed: see _spectral_radius
    return out


def riccati_criterion(
    EN_ref,
    pd,
    mod,
    p,
    T,
    field_dist,
    N_min=DX_N_MIN_DEFAULT,
    N_max=DX_N_MAX_DEFAULT,
    tol=DX_TOL_DEFAULT,
    lam=0.0,
    positive_polarity=True,
    propagator=midpoint_propagator,
    diag_list=None,
):
    """
    Evaluate the inception criterion by the Riccati method.

    A drop-in alternative to :func:`inception_det`, with the same arguments
    and the same step grid, returning g = det(Q0_f + Q0_b P(0)) instead of
    det Q (see :func:`_riccati_g`).  g vanishes exactly where det Q does;
    unlike det Q it is of order one, positive below inception and negative
    above for every mechanism, and cheap to evaluate however optically thick
    the photon groups or however large λ.

    Returns
    -------
    float
        g; −1 where the reflection operator has a pole inside the gap
        (above inception).
    """
    with np.errstate(over="ignore", invalid="ignore"):
        exponents, EN_cathode, eff, N_gamma_eff, aug_mask, n_aug = _discretise(
            EN_ref,
            pd,
            mod,
            p,
            T,
            field_dist,
            N_min,
            N_max,
            tol,
            lam,
            positive_polarity,
            propagator,
            diag_list,
        )
        Q0, S = _boundary_rows(EN_cathode, eff, p, T, N_gamma_eff, aug_mask, n_aug)
        return _riccati_g(exponents, Q0, S)


#: Inception criteria selectable with ``--criterion``; each has the signature
#: of :func:`inception_det` (without ``resolve``) and vanishes at inception.
CRITERIA = {"riccati": riccati_criterion, "detq": inception_det}

#: The criterion used unless ``--criterion`` says otherwise.
CRITERION_DEFAULT = "riccati"


def add_criterion_argument(parser):
    """Register the shared ``--criterion`` option on an argparse *parser*."""
    parser.add_argument(
        "--criterion",
        choices=sorted(CRITERIA),
        default=CRITERION_DEFAULT,
        help=(
            "How the inception condition is evaluated.  'riccati' (default): "
            "the reflection operator carried from the anode, of order one "
            "and cheap however optically thick the photon groups.  'detq': "
            "the boundary determinant det Q, the formulation of the Theory "
            "chapter; slow where the propagator is stiff.  Both solve the "
            "same discrete problem on the same grid."
        ),
    )
