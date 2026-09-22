# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Temporal growth rate λ vs. voltage above inception.

For a fixed geometry (pd, p, T, field), finds the inception voltage V* by solving
det Q(λ=0, EN*) = 0, then for each voltage V ∈ [V*, F·V*] finds λ* > 0 such that
det Q(λ*, EN) = 0.  The result is the discharge growth rate as a function of applied
voltage, with V* as the natural normalisation point.

The root det Q(λ) = 0 is found for fixed EN by:

  1. Evaluating det Q at λ=0 to establish the reference sign.
  2. Expanding an upper bracket geometrically (×10 per step) from 1 s⁻¹ until a
     sign change is detected.
  3. Refining with scipy.optimize.brentq.

NaN returns from inception_det (Q nearly singular near the root) are mapped to the
opposite-sign sentinel so brentq converges through the singularity.

The command-line front end is :mod:`incept1d.cli.growth`.
"""

import numpy as np
import scipy.optimize

from incept1d.constants import kB as _kB
from incept1d.eigenvalues import max_real_eigenvalue as _max_real_eigenvalue
from incept1d.solver import inception_det

# ── Core solver ───────────────────────────────────────────────────────────────


def find_lambda_for_voltage(
    EN_ref, pd, mod, p, T, field_dist, N_min, N_max, tol, propagator, lam_scale=None
):
    """
    Find λ* > 0 such that det Q(λ*, EN_ref) = 0.

    For EN_ref > EN* (overvoltage), det Q(λ=0) is NaN because Q is near-singular
    (same convention as Inception.py: NaN = "above inception").  The root λ* is
    found by locating a λ_hi where det Q is finite and positive (sub-threshold),
    then bracketing [0, λ_hi] with brentq (NaN → -1e-300 sentinel throughout).

    Parameters
    ----------
    EN_ref     : float   — reduced field in Townsend
    pd         : float   — pressure × gap in bar·m
    mod        : Mechanism
    p          : float   — pressure in bar
    T          : float   — temperature in K
    field_dist : FieldDistribution
    N_min, N_max, tol : int, int, float
    propagator : callable
    lam_scale  : float or None
        Approximate λ scale in s⁻¹ for the upper-bracket search.  If None,
        estimated from the dominant ionisation eigenvalue × electron speed.

    Returns
    -------
    lam_star : float   — growth rate in s⁻¹ (NaN on failure)
    status   : str     — 'ok', 'suspect', or reason for failure
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
            positive_polarity=True,
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
    f_ref = _det_raw(lam_hi)
    if np.isfinite(f_check) and np.isfinite(f_ref) and abs(f_check) > 1e-2 * abs(f_ref):
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
):
    """
    Sweep voltages from V* to v_max_factor·V* and find λ at each point.

    Parameters
    ----------
    EN_star     : float   — inception E/N in Townsend
    pd          : float   — pressure × gap in bar·m
    ...

    Returns
    -------
    list of (V_kV, V_ratio, EN_ref, lam_star, tau_ns, nu_ion, status)
        V_kV    — voltage in kV
        V_ratio — V / V*
        EN_ref  — reduced field in Townsend
        lam_star — growth rate in s⁻¹
        tau_ns  — e-folding time in ns
        nu_ion  — ionisation frequency scale in s⁻¹ (max eigenvalue × v_e)
        status  — 'ok' or diagnostic string
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
