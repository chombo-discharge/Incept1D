# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Temporal growth rate λ as a function of voltage above inception.

At the inception voltage the discharge is marginally stable: det Q(0) = 0.
Above it, det Q(λ) = 0 has a root λ > 0, which is the rate at which the
discharge grows in time.  This module finds the inception point first and
then follows λ up the voltage range.

Why the growth rate is the *largest* real root of det Q(λ), and how it is
bracketed, is described in ``docs/source/numerics/rootfinding.rst``.

The command-line front end is :mod:`incept1d.cli.growth`.
"""

import numpy as np
import scipy.optimize

from incept1d.constants import kB as _kB
from incept1d.eigenvalues import max_real_eigenvalue as _max_real_eigenvalue
from incept1d.solver import inception_det

# ── Core solver ───────────────────────────────────────────────────────────────


def peak_ionization_frequency(mod, EN_ref, pd, p, T, field_dist, n_samples=201):
    """
    The fastest local net ionization rate anywhere in the gap, in s⁻¹.

    The largest real eigenvalue of A = R V⁻¹ times the electron speed,
    evaluated at the peak of the field profile.  In a non-uniform gap the
    gap-averaged field can be below the ionization threshold while the
    field at the electrode is far above it, so the average would give 0.

    Parameters
    ----------
    mod : Mechanism
        Loaded mechanism.
    EN_ref : float
        Uniform-equivalent reduced field in Townsend.
    pd : float
        Pressure × gap in bar·m.
    p, T : float
        Pressure in bar and temperature in K.
    field_dist : FieldDistribution
        Gap geometry.
    n_samples : int
        Points at which the profile f(ξ) is sampled, endpoints included.

    Returns
    -------
    float
        ν_ion ≥ 0 in s⁻¹.
    """
    f = field_dist.build(pd / p)
    EN_peak = EN_ref * max(f(xi) for xi in np.linspace(0.0, 1.0, n_samples))
    lam_eig = _max_real_eigenvalue(mod, EN_peak, p, T)  # m⁻¹
    v_e = abs(float(np.diag(mod.get_V(EN_peak, p, T))[mod.ELECTRON_INDEX]))
    return max(lam_eig, 0.0) * v_e


# Ratio between successive λ in the downward scan for the largest root.  Two
# roots closer than this factor can be missed; roots from different feedback
# loops (ion transit, photons) are usually orders of magnitude apart.
_SCAN_FACTOR = 3.0


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
    Find the growth rate λ* > 0 of the dominant mode at E/N = EN_ref.

    det Q(λ) = 0 has many roots — one per mode of the gap, including
    complex ones.  Every coupling of the linear model is a non-negative
    source, so the gain of each feedback loop decreases with λ; the dominant
    mode is then the unique real λ* at which the strongest loop gain is 1,
    and no root, real or complex, lies above it.  λ* is therefore the
    *largest* real root, and it is bracketed from above: det Q is evaluated
    at ten times the fastest local ionization rate — well above any rate at
    which the gap as a whole can grow — and then at successively smaller λ
    until its sign changes, and that bracket is refined with Brent's method.
    Searching upward from λ = 0 instead can stop at a slower mode.

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
        The fastest local ionization rate in s⁻¹, from which the scan
        starts; :func:`peak_ionization_frequency` if None.
    positive_polarity : bool
        True when the ξ = 0 electrode is the anode, as for
        :func:`incept1d.solver.inception_det`.

    Returns
    -------
    lam_star : float
        Growth rate in s⁻¹, NaN if no root could be resolved.
    status : str
        'ok'; 'below_inception'; 'suspect' if det Q does not change sign
        across λ* (a jump rather than a zero); or the reason for failure.
    """

    def _det(lam_val):
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
            resolve=True,
        )

    # det Q(0) > 0 means V ≤ V*: already sub-threshold at λ = 0.
    f0 = _det(0.0)
    if not np.isfinite(f0):
        return float("nan"), "det_Q_unresolved"
    if f0 > 0.0:
        return 0.0, "below_inception"
    if f0 == 0.0:
        return 0.0, "ok"

    if lam_scale is None:
        lam_scale = peak_ionization_frequency(mod, EN_ref, pd, p, T, field_dist)
    lam_scale = max(lam_scale, 1e4)

    # Start above the fastest local ionization rate — no mode can outgrow
    # the fastest local multiplication — and confirm det Q > 0 there.
    lam_hi, f_hi = 10.0 * lam_scale, None
    for _ in range(20):
        f_hi = _det(lam_hi)
        if not np.isfinite(f_hi):
            return float("nan"), "det_Q_unresolved"
        if f_hi > 0.0:
            break
        lam_hi *= 10.0
    else:
        return float("nan"), "no_bracket_found"

    # Scan down to the first sign change: the largest root.
    lam_lo = lam_hi
    while True:
        lam_lo = lam_hi / _SCAN_FACTOR
        if lam_lo < 1.0:
            lam_lo, f_lo = 0.0, f0
            break
        f_lo = _det(lam_lo)
        if not np.isfinite(f_lo):
            return float("nan"), "det_Q_unresolved"
        if f_lo <= 0.0:
            break
        lam_hi, f_hi = lam_lo, f_lo

    lam_star = scipy.optimize.brentq(
        _det, lam_lo, lam_hi, xtol=1e-3, rtol=1e-9, maxiter=100
    )

    # Very small result means V ≈ V*.
    if lam_star < 1.0:
        return 0.0, "ok"

    # A root must be a crossing, not a jump: det Q has to change sign across
    # a small neighbourhood of λ*.  (Brent's method also converges on a
    # discontinuity, e.g. where a photon group crosses κd = 12 and leaves the
    # augmented system, changing the size of Q.)
    below = _det(lam_star * (1.0 - 1e-4))
    above = _det(lam_star * (1.0 + 1e-4))
    if not (below < 0.0 < above):
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
        peak-field ionization frequency (:func:`peak_ionization_frequency`)
        in s⁻¹, and the per-point status.
    """
    V_star = EN_star * pd * 1e-16 / (_kB * T)  # inception voltage in V
    voltages_V = np.geomspace(V_star, v_max_factor * V_star, n_voltages)

    results = []
    for i, V in enumerate(voltages_V):
        EN_ref = V * _kB * T / (pd * 1e-16)  # Townsend
        V_kV = V * 1e-3
        V_ratio = V / V_star

        nu_ion = peak_ionization_frequency(mod, EN_ref, pd, p, T, field_dist)

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
                lam_scale=nu_ion,
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
