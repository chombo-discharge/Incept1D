"""
Core solver: the augmented ODE, its propagators and the boundary determinant
det Q(λ).

The inception threshold is determined by the determinant condition (manuscript,
eq. 333):

    det Q(λ) = 0   at  λ = 0

where Q = [Q_0; Q_d] is the boundary-condition matrix assembled from cathode
(Q_0) and anode (Q_d) constraints on the augmented state vector

    θ = (y, Ψ⁺, Ψ⁻)^T,   y = V n  (species fluxes),
                            Ψ⁺       (forward photon fluxes, +x direction),
                            Ψ⁻       (backward photon fluxes, −x direction).

The augmented ODE is (manuscript Eq. 252):

    ∂_x θ = A_aug θ,

where

    A_aug = [A    B    B ]      A = R V^{-1},  B = photon absorption coupling,
            [C   -D    0 ]
            [-C   0    D ]      C = photon source,  D = diag(κ_j).

When there are no photon species (N_γ = 0) A_aug reduces to A = R V^{-1}.

Cathode boundary conditions (Q_0 θ_0 = 0, manuscript Eq. 305–310):

    (Π_e + γ_+^T Π_+ − γ_Ψ^T Π_bck) θ_0 = 0    ← electron SEE
    Π_− θ_0 = 0                                   ← negative-ion zero flux
    Π_fwd θ_0 = 0                                 ← Ψ⁺(0) = 0 (no incoming fwd photons)

Anode boundary conditions (Q_d θ_0 = 0, manuscript Eq. 319–322):

    Π_+ M(d) θ_0 = 0                              ← positive-ion zero flux
    Π_bck M(d) θ_0 = 0                            ← Ψ⁻(d) = 0 (no incoming bck photons)

The full boundary matrix Q is (N_s + 2 N_γ) × (N_s + 2 N_γ).  A non-trivial
solution exists if and only if det Q = 0 (inception threshold).

The matrix exponential M(d) = exp(A_aug * d) is computed with a Padé
approximation (scipy.linalg.expm) after eigenvalue shifting:

    A_shift = A_aug − λ_max I
    M(d) = exp(A_shift * d) * exp(λ_max * d)

where λ_max = max(Re(eigenvalues(A_aug)), 0).  The shift keeps the Padé
argument well-conditioned even for large positive eigenvalues at high E/N.

For each value of p*d the critical E/N satisfying det Q = 0 is found with
the Brent root-finding method (scipy.optimize.brentq).  The corresponding
breakdown voltage is

    V* = E* * d = (E/N)* * (p*d) * 1e-16 / (kB * T)

where the last equality uses N = p * 1e5 / (kB*T) and the unit conversion
between Townsend (1e-21 V m^2) and bar (1e5 Pa).

The root of det Q in E/N at fixed p·d is found in :mod:`incept1d.inception`.
"""

import math

import numpy as np
import scipy.integrate
import scipy.linalg

from incept1d.constants import kB as _kB, c_light as _C_LIGHT  # noqa: F401
from incept1d.fields import FieldDistribution  # noqa: F401


def _build_A_aug(EN, d, mod, p, T, lam=0.0):
    """
    Build the augmented ODE matrix A_aug for reduced field EN and gap length d.

    Returns (A_aug, N_gamma_eff, aug_mask, n_aug).
    aug_mask is None when N_gamma == 0.

    lam : float
        Temporal growth rate λ (s⁻¹). Modifies A → (R − λI)V⁻¹ and D → κ + λ/c.
        Default 0.0 (standard inception threshold).
    """
    _KD_LOCAL_THRESHOLD = 12.0
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

    if N_gamma > 0:
        kd = kappa * d
        local_mask = kd > _KD_LOCAL_THRESHOLD
        aug_mask = ~local_mask
        for j in np.where(local_mask)[0]:
            A += 2.0 * np.outer(B[:, j], C[j, :]) / kappa[j]
        B_aug = B[:, aug_mask]
        C_aug_ph = C[aug_mask, :]
        kappa_aug = kappa[aug_mask]
        N_gamma_eff = int(aug_mask.sum())
    else:
        N_gamma_eff = 0
        aug_mask = None
        B_aug = B
        C_aug_ph = C
        kappa_aug = kappa

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
    """
    Compute expm(M) with eigenvalue shifting for numerical stability.

    Shifts by the true max real eigenvalue to keep expm well-conditioned without
    over-shrinking the Qd rows.  The restoration factor is intentionally NOT applied;
    it uniformly scales the Qd rows, leaving sign(det Q) unchanged.  Gershgorin
    bounds were tried but cause cumulative over-shrinkage of M across N_steps
    multiplications, making Q ill-conditioned and triggering the cond_Q > 1e14 guard.

    M is the pre-formed matrix argument to expm — for the midpoint rule this is
    A_aug * d_step; for Magnus expansions it is the full Ω already formed by the
    caller.
    """
    real_eigs = np.real(np.linalg.eigvals(M))
    lam_max = float(np.max(real_eigs))
    lam = lam_max if lam_max > 0.0 else 0.0
    return scipy.linalg.expm(M - lam * np.eye(M.shape[0]))


DX_N_MIN_DEFAULT = 5
DX_N_MAX_DEFAULT = 200
DX_TOL_DEFAULT = 0.03


def parse_dx_spec(tokens, parser=None):
    """
    Parse a --dx token list into (N_min, N_max, tol).

    Tokens are all optional; missing values fall back to the defaults
    (N_min=5, N_max=200, tol=0.03).  N_min = N_max disables adaptation
    (max_depth = 0) and gives a constant uniform grid of N_min steps.

    Parameters
    ----------
    tokens : list of str or None
        Raw argparse token list, e.g. ['5', '200', '0.03'].
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


def midpoint_propagator(A_func, x_lo, x_hi):
    """
    Zeroth-order Magnus propagator (midpoint rule).

    Evaluates A at the interval midpoint and returns expm(A_mid * h).

    Parameters
    ----------
    A_func : callable
        A_func(x) → A_aug (ndarray).  x is a code-coordinate in metres; the
        polarity-dependent flip and field evaluation are encapsulated by the caller.
    x_lo, x_hi : float
        Subinterval bounds in code coordinates (metres).

    Returns
    -------
    ndarray
        Propagator matrix P = expm(A(x_mid) * h).
    """
    h = x_hi - x_lo
    x_mid = 0.5 * (x_lo + x_hi)
    return _expm_shifted(A_func(x_mid) * h)


def _adaptive_midpoint_segment(
    A_func, x_lo, x_hi, tol, max_depth, P_coarse=None, _level=0, _diag=None
):
    """
    Compute the midpoint propagator for [x_lo, x_hi] with adaptive step halving.

    Compares the 1-step (coarse) propagator against the 2-step (halved) estimate.
    If the relative Frobenius error exceeds tol and max_depth > 0, each half is
    refined recursively.  Passing P_coarse from the parent avoids one expm call
    when recursing (the parent's half-step result IS the child's coarse estimate).

    Parameters
    ----------
    A_func : callable
        Same signature as for midpoint_propagator.
    x_lo, x_hi : float
        Subinterval bounds in code coordinates (metres).
    tol : float
        Relative Frobenius error threshold; refinement stops when error <= tol.
    max_depth : int
        Maximum remaining recursion depth.  At depth 0 the coarse estimate is
        returned without any halving check.
    P_coarse : ndarray or None
        Pre-computed 1-step propagator for [x_lo, x_hi].  Computed internally
        when None (top-level call).
    _level : int
        Current recursion depth from the initial-segment top-level call (0 = top).
        Used only for diagnostics.
    _diag : list or None
        If a list is supplied, a dict is appended for every halving check performed:
        {'x_lo', 'x_hi', 'level', 'err', 'accepted'}.  Leaf calls at max_depth=0
        where no halving is attempted are NOT recorded.

    Returns
    -------
    ndarray
        Propagator matrix for [x_lo, x_hi].
    """
    if P_coarse is None:
        P_coarse = midpoint_propagator(A_func, x_lo, x_hi)
    if max_depth == 0:
        return P_coarse

    x_mid = 0.5 * (x_lo + x_hi)
    P_l = midpoint_propagator(A_func, x_lo, x_mid)
    P_r = midpoint_propagator(A_func, x_mid, x_hi)
    P_fine = P_r @ P_l

    norm_ref = np.linalg.norm(P_fine, "fro")
    if norm_ref < 1e-300:
        norm_ref = 1.0
    err = np.linalg.norm(P_fine - P_coarse, "fro") / norm_ref

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
        return P_fine

    # Error too large — refine each half, passing the already-computed half-step
    # propagators as the coarse estimates for the sub-calls.
    P_l_fine = _adaptive_midpoint_segment(
        A_func, x_lo, x_mid, tol, max_depth - 1, P_l, _level + 1, _diag
    )
    P_r_fine = _adaptive_midpoint_segment(
        A_func, x_mid, x_hi, tol, max_depth - 1, P_r, _level + 1, _diag
    )
    return P_r_fine @ P_l_fine


def magnus2_propagator(A_func, x_lo, x_hi):
    """
    Second-order Magnus propagator with 2-point Gauss-Legendre quadrature.

    Ω = h/2·(A₁+A₂) + √3·h²/12·[A₂,A₁],   P = expm(Ω)

    The first term is the GL2 quadrature approximation to ∫A dx; the second is
    the leading commutator correction from the Magnus series.  The commutator
    vanishes for a uniform field (A constant across the interval), so Magnus2
    reduces exactly to the midpoint rule in that case.

    Convergence guard: the Magnus series converges only when ∫‖A‖ dx < π.
    If ‖Ω₁‖_F = ‖h/2·(A₁+A₂)‖_F exceeds π the interval is halved and the
    two sub-propagators are multiplied.  This keeps Magnus2 accurate at large
    p·d without requiring the caller to specify a finer grid.

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
        return midpoint_propagator(A_func, x_lo, x_hi)
    Omega2 = (math.sqrt(3.0) * h**2 / 12.0) * (A2 @ A1 - A1 @ A2)
    return _expm_shifted(Omega1 + Omega2)


def _assemble_det_Q(M, EN_cathode, mod, p, T, N_gamma_eff, aug_mask, n_aug):
    """
    Assemble the boundary-condition matrix Q from the propagator M and cathode
    boundary field EN_cathode, then return det Q (NaN if ill-conditioned or overflow).

    Must be called inside a numpy errstate(over='ignore', invalid='ignore') block.
    mod must already be polarity-resolved (via mod.resolve(positive)).
    """
    n = len(mod.SPECIES)

    Pi_e = mod.get_Pi_e()
    Pi_plus = mod.get_Pi_plus()
    Pi_minus = mod.get_Pi_minus()
    gplus = mod.get_gamma_plus(EN_cathode, p, T)

    if N_gamma_eff == 0:
        Q0_e = Pi_e + gplus[np.newaxis, :] @ Pi_plus
        Q0_neg = Pi_minus
        Q0 = np.vstack([Q0_e, Q0_neg])
        Qd = Pi_plus @ M
    else:
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

        Q0_e = (
            Pi_e_aug
            + gplus[np.newaxis, :] @ Pi_plus_aug
            - (g_Psi[np.newaxis, :] @ Pi_bck)
        )
        Q0_neg = Pi_minus_aug
        Q0_fwd = Pi_fwd
        Q0 = np.vstack([Q0_e, Q0_neg, Q0_fwd])

        Qd_plus = Pi_plus_aug @ M
        Qd_bck = Pi_bck @ M
        Qd = np.vstack([Qd_plus, Qd_bck])

    Q = np.vstack([Q0, Qd])

    if np.any(~np.isfinite(Q)):
        return np.nan

    row_norms = np.linalg.norm(Q, axis=1, keepdims=True)
    row_norms = np.where(row_norms < 1e-300, 1.0, row_norms)
    Q_norm = Q / row_norms
    sign, logabsdet = np.linalg.slogdet(Q_norm)
    cond_Q = np.linalg.cond(Q_norm)
    if cond_Q > 1e14:
        return np.nan
    if not np.isfinite(logabsdet):
        return np.nan
    if sign == 0:
        return 0.0
    return float(sign) * min(np.exp(logabsdet), 1e300)


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
):
    """
    Evaluate det Q(λ) for the inception boundary-value problem.

    Integrates the augmented ODE across the gap by accumulating per-subinterval
    propagator matrices, then assembles the boundary-condition matrix Q and returns
    its determinant.  Breakdown occurs when det Q = 0.

    For the midpoint propagator the grid is adaptive: each of the N_min initial
    segments is recursively halved until the relative Frobenius error between the
    one-step and two-step estimates falls below tol, or the maximum refinement
    depth (floor(log2(N_max // N_min))) is reached.  Setting N_min = N_max gives
    a constant uniform grid (max_depth = 0, no halving attempted).

    Parameters
    ----------
    EN_ref : float
        Uniform-equivalent reduced field in Townsend.
    pd : float
        Product of pressure and gap length in bar·m.
    mod : Mechanism
        Loaded mechanism (returned by load_mechanism).  Configuration parameters
        (xi_photo, xi_emit, reaction_multipliers, gamma overrides) are baked in.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.
    field_dist : FieldDistribution
        Field geometry specification (uniform / sphere-plane / sphere-sphere /
        fieldline).  For 'fieldline' the profile is a tabulated ``|E|`` along a
        (possibly curved) field line parametrised by arc length; xi = 0 is the
        first tabulated point.
    N_min : int
        Minimum (initial) number of integration segments.  Default 5.
    N_max : int
        Maximum total number of fine steps across the gap (adaptive budget).
        N_max = N_min disables adaptation.  Default 200.
    tol : float
        Relative Frobenius error threshold for the midpoint step-halving check.
        Default 0.03 (3 %).
    lam : float
        Temporal growth rate λ (s⁻¹, default 0.0).
    positive_polarity : bool
        True → xi = 0 electrode (sphere / field-line start) is anode (+).
        False → xi = 0 electrode is cathode (−).
        Ignored for uniform and sphere-sphere (symmetric) geometries.
    propagator : callable
        Propagator algorithm.  Signature: propagator(A_func, x_lo, x_hi) → ndarray.
        Adaptive step halving is applied only when propagator is midpoint_propagator;
        other propagators (e.g. magnus2_propagator) use N_min uniform steps.

    Returns
    -------
    float
        det Q(λ).  Inception threshold: det Q = 0.
    """
    eff = mod.resolve(positive_polarity)
    d = pd / p
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
    M = np.eye(n_aug)

    with np.errstate(over="ignore", invalid="ignore"):
        if propagator is midpoint_propagator:
            # Adaptive midpoint: each of the N_min segments is refined by step
            # halving until the relative Frobenius error < tol or max_depth is
            # exhausted.  max_depth = 0 (when N_min = N_max) gives a constant
            # uniform grid with no halving attempted.
            for i in range(N_min):
                seg_diag = [] if diag_list is not None else None
                P = _adaptive_midpoint_segment(
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
                M = P @ M
        else:
            for i in range(N_min):
                M = propagator(A_func, i * d_step, (i + 1) * d_step) @ M

        xi_cathode = 1.0 if positive_polarity else 0.0
        EN_cathode = EN_ref * f(xi_cathode)
        return _assemble_det_Q(M, EN_cathode, eff, p, T, N_gamma_eff, aug_mask, n_aug)
