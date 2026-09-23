# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Classical Townsend mechanism for helium, argon and air.

This is the textbook model behind Paschen's law, expressed in the interface
Incept1D expects.  There are no tabulated cross sections and no swarm data:
the ionization coefficient is the two-parameter Townsend form

    alpha / p = A exp(-B p / E),

equivalently, in the reduced variables this code uses,

    alpha = A p exp(-B / (E/N) * (N/p)^-1)  ->  see alpha() below,

and the cathode releases ``gamma`` electrons per arriving positive ion.  With
no attachment, no detachment and no photon feedback, the inception condition
reduces exactly to

    alpha d = ln(1 + 1/gamma)                        (eq_standard_paschen)

whose solution is the classical Paschen curve

    U = B (pd) / [ ln(A pd) - ln(ln(1 + 1/gamma)) ],

with minimum

    U_min  = e B / A * ln(1 + 1/gamma)
    (pd)_min = e / A * ln(1 + 1/gamma).

The point of shipping it is that these closed forms are known, so the example
doubles as a validation of the solver against something that can be checked by
hand.  It is *not* a quantitative model of a real discharge: the Townsend
coefficients are fits valid over a limited E/p range, and air in particular is
electronegative, so the real left branch and the real high-pd behaviour need
``mechanisms/air/``.

Coefficients
------------
A [cm^-1 Torr^-1] and B [V cm^-1 Torr^-1] are the standard values tabulated by
Raizer, *Gas Discharge Physics* (Springer, 1991), Table 4.1, valid over the
stated E/p range:

    gas      A      B      E/p validity [V cm^-1 Torr^-1]
    helium   3      34      20 -  150
    argon   12     180     100 -  600
    air     15     365     100 -  800

They are converted to SI internally.  ``gamma`` is not a gas constant -- it
depends on the cathode material and its surface state -- so it is a
configuration parameter, defaulting to 0.01, a representative value for a
metal cathode.

Species and coordinates
-----------------------
Two species, ``["e", "M+"]``, with the convention of
``mechanisms/air/air_pancheshnyi.py``: cathode at x = 0, anode at x = d,
electrons drift toward +x and positive ions toward -x.
"""

import os
import sys

import numpy as np

from incept1d.constants import kB

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

SPECIES = ["e", "M+"]
ELECTRON_INDEX = 0

# --- gas data ---------------------------------------------------------------

#: Townsend coefficients, Raizer Table 4.1.  (A [cm^-1 Torr^-1],
#: B [V cm^-1 Torr^-1], E/p validity range [V cm^-1 Torr^-1]).
GASES = {
    "helium": (3.0, 34.0, (20.0, 150.0)),
    "argon": (12.0, 180.0, (100.0, 600.0)),
    "air": (15.0, 365.0, (100.0, 800.0)),
}

_TORR_PER_BAR = 750.061682704  # 1 bar in Torr
_CM_PER_M = 100.0

# Injected by config.py before this module executes.
GAS = globals().get("GAS", "air")
_GAMMA0 = globals().get("_GAMMA0", 1.0e-2)
#: Ion mobility * N [(m V s)^-1]; only the ratio to the electron mobility
#: matters for the inception condition, so one representative value is used
#: for every gas.
_MU_ION_N = globals().get("_MU_ION_N", 2.0e21)
_MU_E_N = globals().get("_MU_E_N", 4.0e23)

if GAS not in GASES:
    raise ValueError(f"Unknown gas {GAS!r}; choose one of {sorted(GASES)}")

_A_cm_torr, _B_cm_torr, _EP_RANGE = GASES[GAS]

#: A in [m^-1 bar^-1] and B in [V m^-1 bar^-1].
A_SI = _A_cm_torr * _CM_PER_M * _TORR_PER_BAR
B_SI = _B_cm_torr * _CM_PER_M * _TORR_PER_BAR


def _N(p, T):
    """Neutral number density [m^-3]."""
    return p * 1e5 / (kB * T)


def E_field(EN, p=1.0, T=293.0):
    """Electric field [V/m] from the reduced field EN [Td]."""
    return EN * 1e-21 * _N(p, T)


def alpha(EN, p=1.0, T=293.0):
    """
    Townsend ionization coefficient [m^-1].

        alpha = A p exp(-B p / E)

    Returns 0 at zero field rather than underflowing.
    """
    E = E_field(EN, p, T)
    if E <= 0.0:
        return 0.0
    return A_SI * p * np.exp(-B_SI * p / E)


def eta(EN, p=1.0, T=293.0):
    """Attachment coefficient [m^-1].  Zero: this is the non-attaching model."""
    return 0.0


def reduced_field_validity(p=1.0, T=293.0):
    """
    The E/N range [Td] over which the tabulated fit is valid.

    The Townsend two-parameter form is a fit, not a law; outside this window
    the curve it produces is extrapolation.  Converted from the E/p range in
    :data:`GASES`.
    """
    lo, hi = _EP_RANGE  # V cm^-1 Torr^-1
    scale = _CM_PER_M * _TORR_PER_BAR * 1e21 / (1e5 / (kB * T))
    return lo * scale, hi * scale


def drift_speed(EN, p=1.0, T=293.0):
    """Electron drift speed [m/s]."""
    return _MU_E_N * EN * 1e-21


def ion_drift_speed(EN, p=1.0, T=293.0):
    """Positive-ion drift speed [m/s]."""
    return _MU_ION_N * EN * 1e-21


def get_V(EN, p=1.0, T=293.0):
    """Diagonal drift-velocity matrix [m/s]; positive ions move toward -x."""
    return np.diag([drift_speed(EN, p, T), -ion_drift_speed(EN, p, T)])


def get_R(EN, p=1.0, T=293.0, multipliers=None):
    """
    Reaction matrix [s^-1] for the single reaction ``e -> 2e + M+``.

    ``multipliers`` accepts the key ``"e -> 2e + M+"``.
    """
    m = {k.replace(" ", ""): v for k, v in (multipliers or {}).items()}
    a = alpha(EN, p, T) * m.get("e->2e+M+", 1.0)
    nu = a * drift_speed(EN, p, T)
    R = np.zeros((2, 2))
    R[0, 0] = nu  # net electron production
    R[1, 0] = nu  # one positive ion per ionization
    return R


def get_Pi_e():
    """Row selector for the electron flux, shape (1, 2)."""
    P = np.zeros((1, 2))
    P[0, 0] = 1.0
    return P


def get_Pi_plus():
    """Row selector for the positive-ion flux, shape (1, 2)."""
    P = np.zeros((1, 2))
    P[0, 1] = 1.0
    return P


def get_Pi_minus():
    """No negative ions in this model: shape (0, 2)."""
    return np.zeros((0, 2))


def get_gamma_plus(EN, p=1.0, T=293.0):
    """Ion-induced secondary emission yield, shape (1,).  Field independent."""
    return np.array([_GAMMA0])


def get_gamma_plus_with(EN, p, T, gamma0=None, gamma1=None, eref=None, beta=None):
    """As :func:`get_gamma_plus`, with the per-polarity override."""
    return np.array([_GAMMA0 if gamma0 is None else gamma0])


# --- no photon feedback ------------------------------------------------------


def get_B(EN, p=1.0, T=293.0):
    """Photon absorption coupling, shape (2, 0)."""
    return np.zeros((2, 0))


def get_C(EN, p=1.0, T=293.0):
    """Photon source, shape (0, 2)."""
    return np.zeros((0, 2))


def get_kappa(p=1.0, T=293.0):
    """Photon absorption coefficients, shape (0,)."""
    return np.zeros(0)


def get_gamma_Psi(EN, p=1.0, T=293.0):
    """Photon-induced emission yields, shape (0,)."""
    return np.zeros(0)


# --- closed forms, for the documentation and the tests -----------------------


def paschen_voltage(pd_bar_m, gamma=None):
    """
    Closed-form breakdown voltage [V] at ``pd`` [bar m].

        U = B (pd) / [ln(A pd) - ln(ln(1 + 1/gamma))]

    Returns NaN below the left-branch asymptote, where the denominator is not
    positive and no solution exists.
    """
    g = _GAMMA0 if gamma is None else gamma
    denom = np.log(A_SI * pd_bar_m) - np.log(np.log1p(1.0 / g))
    if denom <= 0.0:
        return float("nan")
    return B_SI * pd_bar_m / denom


def paschen_minimum(gamma=None):
    """
    Closed-form minimum of the Paschen curve.

    Returns ``(pd_min [bar m], U_min [V])`` with

        (pd)_min = e / A  ln(1 + 1/gamma),
        U_min    = e B / A ln(1 + 1/gamma).
    """
    g = _GAMMA0 if gamma is None else gamma
    L = np.log1p(1.0 / g)
    pd_min = np.e * L / A_SI
    return pd_min, np.e * B_SI * L / A_SI
