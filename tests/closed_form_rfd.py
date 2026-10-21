# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Self-sustained threshold of the relativistic feedback model by normal modes.

Transcribed from the steady state of the model equations in
``docs/source/configuration/examples/rfd.rst`` (eq. ``eq_rfd_model``) and the
normal-mode expansion documented there (eq. ``eq_rfd_normal_modes``), in the
paper's own variables: densities n = (n_r, n_pos, n_up, n_down) in z, with
dn/dz = M n.  Like :mod:`closed_form` this shares no code with
:mod:`incept1d`, and takes the coefficients as numbers.

The steady state is a sum of modes v_j e^{kappa_j z} (M v_j = kappa_j v_j).
With nothing entering, n_r(0) = n_up(0) = n_pos(L) = n_down(L) = 0, and a
non-trivial combination of the modes exists where the 4x4 matrix of those
four boundary values has a zero determinant.  That determinant depends on how
each v_j is scaled and ordered; divided by det [v_1 ... v_4] it does not, and
is real.
"""

import numpy as np
import scipy.optimize


def mode_matrix(
    inv_lam,
    P_bb,
    P_cp_up,
    P_cp_down,
    P_pp,
    inv_lam_pos,
    P_br_up,
    P_br_down,
    inv_lam_up,
    inv_lam_down,
):
    """
    M of dn/dz = M n for n = (n_r, n_pos, n_up, n_down), per metre.

    Rows are the steady states of the four model equations: the RRE and the
    forward photon equations as written, the positron and backward photon
    equations (which carry -dn/dz) multiplied through by -1.
    """
    return np.array(
        [
            [inv_lam, P_bb, P_cp_up, P_cp_down],
            [0.0, inv_lam_pos, -P_pp, 0.0],
            [P_br_up, 0.0, -inv_lam_up, 0.0],
            [-P_br_down, 0.0, 0.0, inv_lam_down],
        ]
    )


def boundary_determinant(M, L):
    """
    det A / det [v_1 ... v_4], A the boundary values of the normal modes.

    Zero where a self-sustained steady state exists.  Scaled by
    exp(-sum max(Re kappa_j, 0) L), a positive factor, to stay of order one.
    """
    kappa, W = np.linalg.eig(M)
    # Column j of A scaled by s_j = exp(-max(Re kappa_j, 0) L), applied
    # before the exponential so that nothing overflows; det A gains the
    # factor prod s_j, the scale in the docstring.
    shift = np.maximum(kappa.real, 0.0) * L
    s = np.exp(-shift)
    at_L = np.exp(kappa * L - shift)
    A = np.vstack(
        [
            W[0] * s,  # n_r(0)
            W[2] * s,  # n_up(0)
            W[1] * at_L,  # n_pos(L)
            W[3] * at_L,  # n_down(L)
        ]
    )
    return float(np.real(np.linalg.det(A) / np.linalg.det(W)))


def threshold(M_of_E, L, E_lo, E_hi, n_scan=400):
    """
    Lowest E in [E_lo, E_hi] where :func:`boundary_determinant` changes sign.

    Parameters
    ----------
    M_of_E : callable
        E -> M, e.g. built from :func:`mode_matrix`.
    L : float
        Length of the high-field region [m].
    """
    f = lambda E: boundary_determinant(M_of_E(E), L)  # noqa: E731
    grid = np.linspace(E_lo, E_hi, n_scan)
    values = [f(E) for E in grid]
    for a, b, fa, fb in zip(grid, grid[1:], values, values[1:]):
        if fa * fb < 0.0:
            return scipy.optimize.brentq(f, a, b, xtol=1e-10, rtol=1e-14)
    raise RuntimeError(f"no self-sustained threshold in [{E_lo}, {E_hi}] at L = {L}")
