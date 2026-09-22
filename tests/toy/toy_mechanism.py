# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Three-species toy mechanism with a closed-form inception condition.

This is the reference against which the solver is verified.  It implements the
reduced model of ``docs/source/theory/inceptioncriterion.rst``: electrons, one
positive ion and one negative ion, coupled by ionization, attachment and
detachment, with ion-induced secondary emission at the cathode and no photon
feedback.  With

    alpha(E/N)  ionization coefficient        [m^-1]
    eta(E/N)    attachment coefficient        [m^-1]
    delta       detachment coefficient        [m^-1]
    gamma       cathode secondary emission    [-]

the transport matrix reduces to

    A = R V^-1 = [[alpha - eta, 0,      delta],
                  [alpha,       0,      0    ],
                  [eta,         0,     -delta]]

whose electron/negative-ion block has eigenvalues lambda_pm, and the inception
condition det Q(0) = 0 has the closed form :eq:`eq_generalized_paschen`.  That
expression is implemented independently in ``tests/closed_form.py`` -- keep the
two derivations separate so a mistake in one does not hide a mistake in the
other.

The coefficients are deliberately simple analytic functions of E/N (Townsend's
form for alpha, constant eta and delta) rather than table lookups, so the
mechanism loads instantly and the closed form stays exact.

Species order: ``["e", "M+", "M-"]``; the coordinate convention is the one
documented in ``mechanisms/air/air_pancheshnyi.py`` (cathode at x = 0, anode at
x = d, electrons and negative ions drift toward +x, positive ions toward -x).
"""

import numpy as np

from incept1d.constants import kB

SPECIES = ["e", "M+", "M-"]
ELECTRON_INDEX = 0

# Defaults; a companion config.py may override them through pre_exec_vars().
_ALPHA_A = globals().get("_ALPHA_A", 1.0e4)  # m^-1 bar^-1
_ALPHA_B = globals().get("_ALPHA_B", 300.0)  # Td
_ETA_0 = globals().get("_ETA_0", 20.0)  # m^-1 bar^-1
_DELTA_0 = globals().get("_DELTA_0", 5.0)  # m^-1 bar^-1
_GAMMA0 = globals().get("_GAMMA0", 1.0e-2)  # -
_MU_E_N = globals().get("_MU_E_N", 4.0e23)  # (m V s)^-1
_MU_ION_N = globals().get("_MU_ION_N", 2.0e21)  # (m V s)^-1


def _N(p, T):
    """Neutral number density [m^-3] at pressure p [bar], temperature T [K]."""
    return p * 1e5 / (kB * T)


def alpha(EN, p=1.0, T=293.0):
    """Townsend ionization coefficient [m^-1]: alpha = A p exp(-B / (E/N))."""
    if EN <= 0.0:
        return 0.0
    return _ALPHA_A * p * np.exp(-_ALPHA_B / EN)


def eta(EN, p=1.0, T=293.0):
    """Attachment coefficient [m^-1]; constant, so alpha = eta has one root."""
    return _ETA_0 * p


def delta(EN, p=1.0, T=293.0):
    """Detachment coefficient [m^-1] of the negative ion."""
    return _DELTA_0 * p


def drift_speed(EN, p=1.0, T=293.0):
    """Electron drift speed [m/s] (magnitude)."""
    return _MU_E_N * EN * 1e-21


def ion_drift_speed(EN, p=1.0, T=293.0):
    """Ion drift speed [m/s] (magnitude)."""
    return _MU_ION_N * EN * 1e-21


def get_V(EN, p=1.0, T=293.0):
    """
    Diagonal drift-velocity matrix V [m/s].

    Signs follow the coordinate convention: electrons (+x), positive ions
    (-x), negative ions (+x).
    """
    ve = drift_speed(EN, p, T)
    vi = ion_drift_speed(EN, p, T)
    return np.diag([ve, -vi, vi])


def get_R(EN, p=1.0, T=293.0, multipliers=None):
    """
    Reaction matrix R [s^-1] chosen so that A = R V^-1 is the reduced model.

    ``multipliers`` accepts the same reaction-string keys as a real mechanism;
    the three supported keys are ``"e -> 2e + M+"``, ``"e -> M-"`` and
    ``"M- -> e"``.
    """
    m = {k.replace(" ", ""): v for k, v in (multipliers or {}).items()}
    a = alpha(EN, p, T) * m.get("e->2e+M+", 1.0)
    e_att = eta(EN, p, T) * m.get("e->M-", 1.0)
    d_det = delta(EN, p, T) * m.get("M->e", m.get("M--->e", 1.0))

    ve = drift_speed(EN, p, T)
    vi = ion_drift_speed(EN, p, T)

    nu_ion, nu_att, nu_det = a * ve, e_att * ve, d_det * vi

    R = np.zeros((3, 3))
    R[0, 0] = nu_ion - nu_att  # e from e
    R[1, 0] = nu_ion  # M+ from e
    R[2, 0] = nu_att  # M- from e
    R[0, 2] = nu_det  # e from M-
    R[2, 2] = -nu_det  # M- loss
    return R


def get_Pi_e():
    """Row selector for the electron flux, shape (1, 3)."""
    P = np.zeros((1, 3))
    P[0, 0] = 1.0
    return P


def get_Pi_plus():
    """Row selector for the positive-ion flux, shape (1, 3)."""
    P = np.zeros((1, 3))
    P[0, 1] = 1.0
    return P


def get_Pi_minus():
    """Row selector for the negative-ion flux, shape (1, 3)."""
    P = np.zeros((1, 3))
    P[0, 2] = 1.0
    return P


def get_gamma_plus(EN, p=1.0, T=293.0):
    """Ion-induced secondary emission yield, shape (1,).  Field independent."""
    return np.array([_GAMMA0])


def get_gamma_plus_with(EN, p, T, gamma0=None, gamma1=None, eref=None, beta=None):
    """As get_gamma_plus, with the per-polarity override used by Mechanism."""
    g0 = _GAMMA0 if gamma0 is None else gamma0
    return np.array([g0])


# ---------------------------------------------------------------------------
# No photon feedback: the photon blocks are zero-width.
# ---------------------------------------------------------------------------


def get_B(EN, p=1.0, T=293.0):
    """Photon absorption coupling, shape (3, 0) -- no photon groups."""
    return np.zeros((3, 0))


def get_C(EN, p=1.0, T=293.0):
    """Photon source, shape (0, 3) -- no photon groups."""
    return np.zeros((0, 3))


def get_kappa(p=1.0, T=293.0):
    """Photon absorption coefficients, shape (0,) -- no photon groups."""
    return np.zeros(0)


def get_gamma_Psi(EN, p=1.0, T=293.0):
    """Photon-induced emission yields, shape (0,) -- no photon groups."""
    return np.zeros(0)
