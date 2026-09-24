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

import math

import numpy as np
import scipy.integrate
import scipy.linalg

from incept1d.constants import kB as _kB, c_light as _C_LIGHT  # noqa: F401
from incept1d.fields import FieldDistribution  # noqa: F401


def _build_A_aug(EN, d, mod, p, T, lam=0.0):
    """Build the augmented ODE matrix for reduced field EN and gap length d.

    Photon groups that are re-absorbed within a fraction of a step are
    folded into A rather than propagated explicitly, so the augmented
    system can be smaller than N_s + 2 N_γ.

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
        explicitly propagated photon groups, and is None when N_γ = 0.
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
DX_N_MAX_DEFAULT = 200
DX_TOL_DEFAULT = 0.03


def parse_dx_spec(tokens, parser=None):
    """
    Parse a --dx token list into (N_min, N_max, tol).

    Tokens are all optional and fall back to the module defaults.
    N_min = N_max disables adaptation and gives a uniform grid.

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
    h = x_hi - x_lo
    x_mid = 0.5 * (x_lo + x_hi)
    return _expm_shifted(A_func(x_mid) * h)


def _adaptive_midpoint_segment(
    A_func, x_lo, x_hi, tol, max_depth, P_coarse=None, _level=0, _diag=None
):
    """
    Midpoint propagator for [x_lo, x_hi], refined by recursive step halving.

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
    P_coarse : ndarray or None
        The parent's half-step propagator, which is this call's coarse
        estimate.  Passing it saves one expm; None at the top level.
    _level : int
        Recursion depth, for diagnostics only.
    _diag : list or None
        If given, one dict per halving check is appended:
        ``{'x_lo', 'x_hi', 'level', 'err', 'accepted'}``.  Leaf calls at
        max_depth 0 attempt no halving and are not recorded.

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
    Second-order Magnus propagator, with two-point Gauss-Legendre quadrature.

    The Magnus series converges only for ∫‖A‖ dx < π, so an interval whose
    quadrature term exceeds that is halved and the two sub-propagators
    multiplied.  The caller therefore does not have to pick a grid fine
    enough for large p·d.

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
    positive_polarity : bool
        True when the ξ = 0 electrode is the anode.  Ignored for the
        symmetric geometries.
    propagator : callable
        ``propagator(A_func, x_lo, x_hi) -> ndarray``.  Adaptive halving is
        applied only to :func:`midpoint_propagator`; anything else runs on
        the uniform N_min grid.

    Returns
    -------
    float
        det Q(λ).  Inception threshold: det Q = 0.
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
    # whichever propagator was requested.  _expm_shifted drops the same scalar
    # factor exp(λ_max·d) that the stepped product drops as N factors of
    # exp(λ_max·h), so det Q is identical, not merely equivalent.
    if field_dist.field_type == "uniform":
        A_aug, N_gamma_eff, aug_mask, n_aug = _build_A_aug(
            EN_ref, d, eff, p, T, lam=lam
        )
        with np.errstate(over="ignore", invalid="ignore"):
            M = _expm_shifted(A_aug * d)
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
            return _assemble_det_Q(M, EN_ref, eff, p, T, N_gamma_eff, aug_mask, n_aug)

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
