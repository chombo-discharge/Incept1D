# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Eigenvalues of the local transport matrix A = R V^{-1}.

For a given reaction mechanism the matrix A = R V^{-1} is formed at each value
of E/N and all eigenvalues are computed.  A positive real eigenvalue at a given
E/N means the discharge has a spatially growing mode (net ionisation exceeds
attachment) and the electron flux amplifies as it travels from cathode to anode.

Two operating modes are available:

Default mode
------------
Plot all N eigenvalue tracks (Re λ_j / N) vs E/N for a fixed pressure.
One subplot per configuration; subplots share the x-axis.

Pressure-scan mode  (--pressure-scan)
--------------------------------------
Plot the leading eigenvalue (or a user-selected track) vs E/N for several
log-spaced pressures.  One subplot per configuration, one line per pressure.

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
    Compute tracked eigenvalues of A = R V^{-1} for each value of E/N.

    At each E/N the transport matrix A = R(EN) @ inv(V(EN)) is formed and all
    N eigenvalues are computed.  Eigenvalues are continuously tracked across
    E/N using the Hungarian algorithm to prevent spurious reordering jumps.
    The initial ordering is by descending real part (most-growing mode in
    slot 0).

    Parameters
    ----------
    mod : Mechanism
        Loaded mechanism (from Inception.load_mechanism).
    EN_range : array-like, shape (M,)
        E/N values in Townsend.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.

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
    Compute Re(λ_{eig_index} / N) vs E/N for each pressure.

    Parameters
    ----------
    mod : Mechanism
    EN_range : array-like, shape (M,)
    pressures : array-like, shape (P,)
        Pressures in bar.
    T : float
        Gas temperature in Kelvin.
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
