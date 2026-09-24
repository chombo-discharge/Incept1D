# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Temporal growth rate λ as a function of voltage above inception.

At the inception voltage the discharge is marginally stable: det Q(0) = 0.
Above it, det Q(λ) = 0 has a root λ > 0, which is the rate at which the
discharge grows in time.  This module finds the inception point first and
then follows λ up the voltage range.

The root finding, and the NaN convention it relies on, are described in
``docs/source/numerics/rootfinding.rst``.

The command-line front end is :mod:`incept1d.cli.growth`.
"""

import numpy as np
import scipy.optimize

from incept1d.constants import kB as _kB
from incept1d.eigenvalues import max_real_eigenvalue as _max_real_eigenvalue
from incept1d.solver import inception_det

# ── Core solver ───────────────────────────────────────────────────────────────


def find_lambda_for_voltage(
    EN_ref,
    pd,
    mod,
    p,
    T,
    field_dist,
    N_min,
    N_max,
    tol,
    propagator,
    lam_scale=None,
    positive_polarity=True,
):
    """
    Find λ* > 0 such that det Q(λ*, EN_ref) = 0.

    Above inception det Q(0) is NaN, because Q is numerically singular; a
    large enough λ damps the solution and restores a finite, positive
    determinant.  The bracket between the two is refined with Brent's
    method.

    Parameters
    ----------
    EN_ref : float
        Reduced field in Townsend.
    pd : float
        Pressure × gap in bar·m.
    mod : Mechanism
        Loaded mechanism.
    p, T : float
        Pressure in bar and temperature in K.
    field_dist : FieldDistribution
        Gap geometry.
    N_min, N_max, tol : int, int, float
        Integration grid, as for :func:`incept1d.solver.inception_det`.
    propagator : callable
        Propagator algorithm.
    lam_scale : float or None
        Approximate λ scale in s⁻¹ for the bracket search.  Estimated from
        the dominant ionization eigenvalue and the electron speed if None.
    positive_polarity : bool
        True when the ξ = 0 electrode is the anode, as for
        :func:`incept1d.solver.inception_det`.

    Returns
    -------
    lam_star : float
        Growth rate in s⁻¹, NaN if no root could be resolved.
    status : str
        'ok', 'suspect', or the reason for failure.
    """

    def _det_raw(lam_val):
        return inception_det(
            EN_ref,
            pd,
            mod,
            p,
            T,
            field_dist=field_dist,
            N_min=N_min,
            N_max=N_max,
            tol=tol,
            lam=lam_val,
            positive_polarity=positive_polarity,
            propagator=propagator,
        )

    # NaN convention: NaN → -1e-300 (= "above inception", same as Inception.py).
    # This is correct for λ < λ* (where Q is near-singular due to overvoltage),
    # AND is used consistently throughout brentq.
    def _f_brentq(lam_val):
        val = _det_raw(lam_val)
        return val if np.isfinite(val) else -1e-300

    # det Q(0) > 0 means V ≤ V*: already sub-threshold at λ=0.
    f0 = _det_raw(0.0)
    if np.isfinite(f0) and f0 > 0.0:
        return 0.0, "below_inception"

    # Estimate the scale for the upper-bracket search.
    if lam_scale is None:
        lam_eig = max(1.0, _max_real_eigenvalue(mod, EN_ref, p, T))
        Vmat = mod.get_V(EN_ref, p, T)
        v_e = abs(float(np.diag(Vmat)[mod.ELECTRON_INDEX]))
        lam_scale = max(lam_eig * v_e, 1e4)

    # Probe around lam_scale (and broader) for a λ_hi where det Q is finite & > 0.
    # This point is in the "sub-threshold" region (λ > λ*) where Q is well-conditioned.
    lam_hi = None
    probes = [1.0, 0.3, 3.0, 0.1, 10.0, 30.0, 0.03, 100.0, 300.0, 1000.0]
    for factor in probes:
        test_lam = lam_scale * factor
        v = _det_raw(test_lam)
        if np.isfinite(v) and v > 0.0:
            lam_hi = test_lam
            break
    if lam_hi is None:
        for log_lam in range(3, 14):
            v = _det_raw(10.0**log_lam)
            if np.isfinite(v) and v > 0.0:
                lam_hi = 10.0**log_lam
                break
    if lam_hi is None:
        return float("nan"), "no_bracket_found"

    # f_brentq(0) ≤ 0  (NaN → -1e-300, or genuinely negative)
    # f_brentq(lam_hi) > 0
    # → valid bracket for brentq.
    lam_star = scipy.optimize.brentq(
        _f_brentq,
        0.0,
        lam_hi,
        xtol=1.0,
        rtol=1e-6,
        maxiter=60,
    )

    # Very small result means V ≈ V* (Q was singular at λ=0 and root is near zero).
    if lam_star < 1.0:
        return 0.0, "ok"

    f_check = _det_raw(lam_star)

    # A root must have det Q ≈ 0.  A *NaN* there means brentq did not find a
    # zero at all: _f_brentq maps NaN to a negative sentinel, so a run of NaN
    # below a positive value looks exactly like a sign change, and brentq
    # converges on the edge of the NaN region instead.
    #
    # That is not hypothetical.  Above roughly 1.2·V* the spectral spread of
    # A_aug·d exceeds the double-precision underflow limit (~709), the
    # subdominant modes of M underflow to exactly zero, Q becomes rank
    # deficient and _assemble_det_Q returns NaN for *every* λ.  The apparent
    # root is then the λ at which a photon group crosses κd = 12 and leaves
    # the augmented block, changing Q's size — which is independent of E/N,
    # so the same λ was reported for every overvoltage.
    #
    # There is no way to recover a growth rate from a determinant that cannot
    # be evaluated, so report the failure rather than the artefact.  This is
    # the same check _accept_root performs in incept1d.inception.
    if not np.isfinite(f_check):
        return float("nan"), "det_Q_unresolved"

    f_ref = _det_raw(lam_hi)
    if np.isfinite(f_ref) and abs(f_check) > 1e-2 * abs(f_ref):
        return lam_star, "suspect"
    return lam_star, "ok"


def compute_lambda_curve(
    EN_star,
    pd,
    mod,
    p,
    T,
    field_dist,
    N_min,
    N_max,
    tol,
    propagator,
    n_voltages,
    v_max_factor,
    positive_polarity=True,
):
    """
    Sweep voltages from V* to v_max_factor·V* and find λ at each point.

    Parameters
    ----------
    EN_star : float
        Inception E/N in Townsend, the lower end of the sweep.
    pd : float
        Pressure × gap in bar·m.
    mod : Mechanism
        Loaded mechanism.
    p, T : float
        Pressure in bar and temperature in K.
    field_dist : FieldDistribution
        Gap geometry.
    N_min, N_max, tol : int, int, float
        Integration grid, as for :func:`incept1d.solver.inception_det`.
    propagator : callable
        Propagator algorithm.
    n_voltages : int
        Number of voltages in the sweep.
    v_max_factor : float
        Upper end of the sweep, as a multiple of V*.
    positive_polarity : bool
        True when the ξ = 0 electrode is the anode.  *EN_star* must be the
        inception field for the same polarity.

    Returns
    -------
    list of tuple
        One ``(V_kV, V_ratio, EN_ref, lam_star, tau_ns, nu_ion, status)``
        per voltage: the voltage in kV and as a multiple of V*, the reduced
        field, the growth rate in s⁻¹, the e-folding time in ns, the
        ionization frequency scale in s⁻¹, and the per-point status.
    """
    V_star = EN_star * pd * 1e-16 / (_kB * T)  # inception voltage in V
    voltages_V = np.geomspace(V_star, v_max_factor * V_star, n_voltages)

    results = []
    for i, V in enumerate(voltages_V):
        EN_ref = V * _kB * T / (pd * 1e-16)  # Townsend
        V_kV = V * 1e-3
        V_ratio = V / V_star

        # Ionisation frequency scale: max real eigenvalue of A = R V⁻¹, times v_e.
        lam_eig = _max_real_eigenvalue(mod, EN_ref, p, T)  # m⁻¹
        Vmat = mod.get_V(EN_ref, p, T)
        v_e = abs(float(np.diag(Vmat)[mod.ELECTRON_INDEX]))  # m/s
        nu_ion = max(lam_eig, 0.0) * v_e  # s⁻¹
        lam_scale = max(nu_ion, 1e4)

        if i == 0:
            # V = V* by construction; λ* = 0 by definition.
            lam_star, status = 0.0, "ok"
        else:
            lam_star, status = find_lambda_for_voltage(
                EN_ref,
                pd,
                mod,
                p,
                T,
                field_dist,
                N_min,
                N_max,
                tol,
                propagator,
                lam_scale=lam_scale,
                positive_polarity=positive_polarity,
            )

        tau_ns = (
            1e9 / lam_star
            if (np.isfinite(lam_star) and lam_star > 0.0)
            else float("nan")
        )
        results.append((V_kV, V_ratio, EN_ref, lam_star, tau_ns, nu_ion, status))

        tag = f"  [{status}]" if status != "ok" else ""
        if i == 0:
            print(f"  V = {V_kV:8.4f} kV  (1.000×V*)  λ = 0  (inception){tag}")
        elif np.isfinite(lam_star):
            print(
                f"  V = {V_kV:8.4f} kV  ({V_ratio:.3f}×V*)  "
                f"λ = {lam_star:.4e} s⁻¹  τ = {tau_ns:.3f} ns{tag}"
            )
        else:
            print(f"  V = {V_kV:8.4f} kV  ({V_ratio:.3f}×V*)  λ = NaN{tag}")

    return results


# ── File output ───────────────────────────────────────────────────────────────
