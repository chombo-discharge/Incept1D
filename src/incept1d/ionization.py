# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Ionization integrals across the gap.

Two integrals are available for a given voltage and geometry: the classical
∫ max(α−η, 0) dx, and ∫ max(Re λ_max(R V⁻¹), 0) dx, which replaces the
swarm coefficients with the leading eigenvalue of the transport matrix and
so accounts for the ion chemistry the first one ignores.  Both are
approximations to the full criterion of :mod:`incept1d.solver`, and exist
to be compared against it.

The command-line front end is :mod:`incept1d.cli.ionization`.
"""

import numpy as np

from incept1d.eigenvalues import max_real_eigenvalue as _max_real_eigenvalue
from incept1d.fields import FieldDistribution


def aed_integral(EN_ref, p_val, d_val, mod, T, field_dist: FieldDistribution, N: int):
    """Return ∫ max(α−η, 0) dx in m⁻¹, or NaN on overflow.

    Parameters
    ----------
    EN_ref : float
        E/N in Td at the reference point.
    p_val, d_val, T : float
        Pressure in bar, gap distance in metres, temperature in K.
    mod : Mechanism
        Loaded mechanism; must expose ``alpha`` and ``eta``.
    field_dist : FieldDistribution
        Gap geometry.
    N : int
        Number of midpoint-rule integration steps.
    """
    try:
        if field_dist.field_type == "uniform":
            val = (mod.alpha(EN_ref, p_val, T) - mod.eta(EN_ref, p_val, T)) * d_val
            result = max(0.0, float(val))
            return result if np.isfinite(result) else float("nan")

        f = field_dist.build(d_val)

        def _net(xi):
            en = EN_ref * f(xi)
            return mod.alpha(en, p_val, T) - mod.eta(en, p_val, T)

        if not field_dist.is_monotone_decreasing:
            # General profile (sphere-sphere, tabulated field line): the active
            # region need not start at xi = 0 nor be a single interval, so use
            # plain midpoint quadrature of max(α−η, 0) over the whole gap.
            xis = (np.arange(N) + 0.5) / N
            diff = np.array([_net(xi) for xi in xis])
            result = float(np.sum(np.maximum(0.0, diff)) * d_val / N)
            return result if np.isfinite(result) else float("nan")

        # If there is no net ionisation even at the sphere surface, return 0 immediately.
        if _net(0.0) <= 0.0:
            return 0.0

        # Locate xi_cross ∈ (0, 1] where α(E) = η(E).  The field is monotonically
        # decreasing from sphere (xi=0) to plane/far-sphere (xi=1), so _net is also
        # decreasing and the active region is exactly [0, xi_cross].
        # Restricting integration to this interval avoids the partial-cell error that
        # arises when a coarse quadrature cell straddles the threshold.
        if _net(1.0) > 0.0:
            xi_cross = 1.0  # active all the way to the far electrode
        else:
            lo, hi = 0.0, 1.0
            for _ in range(52):  # 52 bisections → ~1e-15 precision in xi
                mid = 0.5 * (lo + hi)
                if _net(mid) > 0.0:
                    lo = mid
                else:
                    hi = mid
            xi_cross = 0.5 * (lo + hi)

        # Midpoint-rule quadrature over [0, xi_cross·d].  Every cell is within the
        # active region so all contributions are positive — no max(0,·) needed.
        xis = (np.arange(N) + 0.5) / N * xi_cross
        EN_arr = EN_ref * np.array([f(xi) for xi in xis])
        ds = xi_cross * d_val / N
        diff = np.array(
            [mod.alpha(en, p_val, T) - mod.eta(en, p_val, T) for en in EN_arr]
        )
        result = float(np.sum(diff) * ds)

    except (OverflowError, ValueError):
        return float("nan")

    return result if np.isfinite(result) else float("nan")


def eig_integral(EN_ref, p_val, d_val, mod, T, field_dist: FieldDistribution, N: int):
    """Return ∫ max(Re λ_max(R V⁻¹), 0) dx in m⁻¹, or NaN on overflow.

    Parameters
    ----------
    EN_ref : float
        E/N in Td at the reference point.
    p_val, d_val, T : float
        Pressure in bar, gap distance in metres, temperature in K.
    mod : Mechanism
        Loaded mechanism; must expose ``get_R`` and ``get_V``.
    field_dist : FieldDistribution
        Gap geometry.
    N : int
        Number of midpoint-rule integration steps.
    """
    try:
        if field_dist.field_type == "uniform":
            val = _max_real_eigenvalue(mod, EN_ref, p_val, T) * d_val
            result = max(0.0, float(val))
            return result if np.isfinite(result) else float("nan")

        f = field_dist.build(d_val)

        def _lmax(xi):
            return _max_real_eigenvalue(mod, EN_ref * f(xi), p_val, T)

        if not field_dist.is_monotone_decreasing:
            # General profile: see aed_integral.
            xis = (np.arange(N) + 0.5) / N
            eigs = np.array([_lmax(xi) for xi in xis])
            result = float(np.sum(np.maximum(0.0, eigs)) * d_val / N)
            return result if np.isfinite(result) else float("nan")

        if _lmax(0.0) <= 0.0:
            return 0.0

        if _lmax(1.0) > 0.0:
            xi_cross = 1.0
        else:
            lo, hi = 0.0, 1.0
            for _ in range(52):
                mid = 0.5 * (lo + hi)
                if _lmax(mid) > 0.0:
                    lo = mid
                else:
                    hi = mid
            xi_cross = 0.5 * (lo + hi)

        xis = (np.arange(N) + 0.5) / N * xi_cross
        EN_arr = EN_ref * np.array([f(xi) for xi in xis])
        ds = xi_cross * d_val / N
        eigs = np.array([_max_real_eigenvalue(mod, en, p_val, T) for en in EN_arr])
        result = float(np.sum(eigs) * ds)

    except (OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan")

    return result if np.isfinite(result) else float("nan")


# ------------------------------------------------------------------
# Data-file helpers
# ------------------------------------------------------------------
