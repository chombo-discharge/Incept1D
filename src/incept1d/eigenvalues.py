# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Eigenvalues of the local transport matrix A = R V⁻¹.

A is formed at each E/N and all its eigenvalues are computed.  A positive
real eigenvalue means a spatially growing mode: net ionization exceeds
attachment, and the electron flux amplifies across the gap.  This is a
local diagnostic — no gap integration, no boundary conditions — so it says
where growth is possible, not whether the gap breaks down.

The command-line front end is :mod:`incept1d.cli.eigenvalues`.
"""

import numpy as np
import scipy.optimize

from incept1d.constants import kB as _kB


def max_real_eigenvalue(mod, EN, p, T):
    """Return Re(λ_max) of A = R V^{-1} at a single E/N value [m^{-1}]."""
    R = mod.get_R(EN, p, T)
    V = mod.get_V(EN, p, T)
    return float(np.max(np.real(np.linalg.eigvals(R @ np.linalg.inv(V)))))


def _track_step(prev, curr):
    """
    Reorder *curr* eigenvalues to best continue the tracks in *prev*.

    Uses the Hungarian algorithm to find the permutation of *curr* that
    minimises total squared distance in the complex plane from *prev*.
    """
    cost = np.abs(prev[:, np.newaxis] - curr[np.newaxis, :]) ** 2
    _, col_ind = scipy.optimize.linear_sum_assignment(cost)
    return curr[col_ind]


def compute_eigenvalues(mod, EN_range, p, T):
    """
    Compute tracked eigenvalues of A = R V⁻¹ for each value of E/N.

    Eigenvalues are matched between adjacent E/N so that a track follows
    one mode rather than jumping when two eigenvalues cross.  Track 0 is
    the most-growing mode at the first E/N.

    Parameters
    ----------
    mod : Mechanism
        Loaded mechanism.
    EN_range : array-like, shape (M,)
        E/N values in Townsend.
    p, T : float
        Gas pressure in bar and temperature in K.

    Returns
    -------
    numpy.ndarray, shape (M, N), dtype complex
        Tracked eigenvalues.
    """
    n = len(mod.SPECIES)
    n_EN = len(EN_range)
    eigvals = np.zeros((n_EN, n), dtype=complex)

    # Scan high→low so the ionisation mode (unambiguously the largest real
    # eigenvalue at high E/N) seeds track 0, avoiding misassignment to the two
    # structural zero eigenvalues that exist because columns 1 and 2 of R are
    # identically zero (N2+, O2+ have no off-diagonal source reactions).
    for i, EN in enumerate(EN_range[::-1]):
        R = mod.get_R(EN, p, T)
        V = mod.get_V(EN, p, T)
        A = R @ np.linalg.inv(V)
        eigs = np.linalg.eigvals(A)

        if i == 0:
            eigs = eigs[np.argsort(-np.real(eigs))]
        else:
            eigs = _track_step(eigvals[i - 1], eigs)

        eigvals[i] = eigs

    return eigvals[::-1]


def compute_pressure_scan(mod, EN_range, pressures, T, eig_index):
    """
    Compute Re(λ_j / N) vs E/N for each pressure, for one track j.

    Parameters
    ----------
    mod : Mechanism
    EN_range : array-like, shape (M,)
        E/N values in Townsend.
    pressures : array-like, shape (P,)
        Pressures in bar.
    T : float
        Gas temperature in K.
    eig_index : int
        Which eigenvalue track to extract.

    Returns
    -------
    numpy.ndarray, shape (M, P)
        Re(λ_{eig_index} / N) for each (E/N, pressure) pair, in m².
    """
    n_EN = len(EN_range)
    n_p = len(pressures)
    result = np.zeros((n_EN, n_p))

    for pi, p in enumerate(pressures):
        N_density = p * 1e5 / (_kB * T)
        ev = compute_eigenvalues(mod, EN_range, p, T)
        result[:, pi] = np.real(ev[:, eig_index]) / N_density

    return result
