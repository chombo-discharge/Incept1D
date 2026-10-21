# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Relativistic feedback discharge (RFD) in air: the fluid model of Færder et al.

Relativistic runaway electrons (RREs) avalanche along the field and emit
bremsstrahlung photons forward and backward.  Forward photons make new RREs
(Compton scattering, photoabsorption) and positrons (pair production); the
positrons run back against the field and make RREs by Bhabha scattering, and
so do the backward photons.  The discharge is self-sustained when this
feedback closes inside the high-field region on its own, with nothing coming
in at either end.  The model and every coefficient below are from
[Faerder2026b]_ (main text, eqs. 3-6, and its supplementary material,
eqs. S1-S9, where they were fitted to GEANT4 simulations).

Species and photon groups
-------------------------
Two species, ``["rre", "e+"]``, and two photon groups, ``["forward",
"backward"]`` (``PHOTON_DIRECTIONS``).  Coordinates follow
``mechanisms/air/pancheshnyi/air_pancheshnyi.py``: x = 0 is where the RREs
enter the high-field region (the paper's z = 0) and x = d its far end.  RREs
drift toward +x at c~ = 0.89 c and positrons toward -x at c~, so the RRE is
the "electron" of the interface and the positron its "positive ion": no
positron enters at x = d, and nothing is emitted at x = 0
(``gamma_+ = gamma_Psi = 0``).  The photon fluxes are Psi = c~ n_gamma.

Field and density
-----------------
The paper gives every coefficient at sea-level density, where all its
results are, and scales only the avalanche thresholds with the air density
relative to sea level, n~ (E_t = 284 kV/m n~).  Off sea level this module
does what the authors' code does (doi:10.5281/zenodo.22983026): the
avalanche length is U_R / (E - n~ E_d), and every other coefficient is its
sea-level function of the field E itself, times n~.  The sea-level number
density N0 that sets n~ = N / N0 is the authors' 2.688e25 m^-3
(``config.py``); at it, n~ = 1 and every coefficient is the paper's.

Coefficients [1/m, at n~ = 1; E in V/m]
-----------------------------------------
    1/lambda      = (E - E_d) / U_R,  U_R = 7300 kV,
                    E_d = 276 kV/m for E >= E_t = 284 kV/m, else 320 kV/m
    P_bb          = fit(E; 4.325e-3, 7.084e-10, 305.4e3, 427.6e3, 1.1),
                    0 for E <= 305.4 kV/m                            (S1)
    P_cp,fwd      = fit(E; 2.468e-3, 1.028e-8, -623.9e3, 1452.9e3, 5.5)  (S3)
    P_cp,bck      = fit(E; 1.019e-3, 1.062e-8, -10.7e3, 874.9e3, 4.1)    (S4)
    P_pp,fwd      = 8.580e-5                                         (S5)
    P_br,fwd      = 1.442e-2,  P_br,bck = 9.666e-4                   (S6)
    1/lambda_a,pos = 127 V/m^2 / (E - 259.1 kV/m)                    (S9)
    lambda_a,fwd  = 307.64 m,  lambda_a,bck = 92.22 m (R = 5 km)    (S7)

with fit(E; P, s, E_d, E_s, r) = (P + s E)(1 - exp(-((E - E_d)/E_s)^r)).
The photon attenuation lengths depend on the horizontal radius R of the
high-field region; ``radius.json`` holds the paper's two cases (S7, S8).

Below 305.4 kV/m there is no Bhabha production, so positrons feed nothing
back, and 1/lambda_a,pos is held at its 305.4 kV/m value there instead of
diverging at 259.1 kV/m.  This changes nothing the positrons can affect.
"""

import numpy as np

from incept1d.constants import c_light, kB

SPECIES = ["rre", "e+"]
ELECTRON_INDEX = 0
PHOTON_DIRECTIONS = ["forward", "backward"]

# Injected by config.py before this module executes.
#: Sea-level number density [m^-3] that n~ = N / N0 is relative to.
_N0 = globals().get("_N0", 2.688e25)
#: Photon attenuation lengths at n~ = 1 [m], forward and backward (S7).
_LAMBDA_A_FWD = globals().get("_LAMBDA_A_FWD", 307.64)
_LAMBDA_A_BCK = globals().get("_LAMBDA_A_BCK", 92.22)

#: Fluid velocity of RREs and positrons [m/s].
C_TILDE = 0.89 * c_light

# --- coefficients at n~ = 1, as functions of E [V/m] ---------------------------

U_R = 7.3e6  # V
E_T = 284e3  # V/m, RREA threshold
_E_D_ABOVE, _E_D_BELOW = 276e3, 320e3  # V/m, avalanche length offsets

#: (P_inf [1/m], s_inf [1/V], E_d [V/m], E_s [V/m], r) of S1, S3, S4.
_FIT_BB = (4.325e-3, 7.084e-10, 305.4e3, 427.6e3, 1.1)
_FIT_CP_FWD = (2.468e-3, 1.028e-8, -623.9e3, 1452.9e3, 5.5)
_FIT_CP_BCK = (1.019e-3, 1.062e-8, -10.7e3, 874.9e3, 4.1)

P_PP_FWD = 8.580e-5  # 1/m (S5)
P_BR_FWD = 1.442e-2  # 1/m (S6)
P_BR_BCK = 9.666e-4  # 1/m (S6)

_A_POS = 127.0  # V/m^2 (S9)
_E_D_POS = 259.1e3  # V/m (S9)


def _fit(E, P_inf, s_inf, E_d, E_s, r):
    """(P + s E)(1 - exp(-((E - E_d)/E_s)^r)) for E > E_d, else 0."""
    if E <= E_d:
        return 0.0
    return (P_inf + s_inf * E) * -np.expm1(-(((E - E_d) / E_s) ** r))


def inv_avalanche_length(E):
    """1/lambda [1/m] of the RRE avalanche at n~ = 1."""
    return (E - (_E_D_ABOVE if E >= E_T else _E_D_BELOW)) / U_R


def P_bb(E):
    """RREs per metre per positron, Bhabha scattering (S1)."""
    return _fit(E, *_FIT_BB)


def P_cp_fwd(E):
    """RREs per metre per forward photon, Compton + photoabsorption (S3)."""
    return _fit(E, *_FIT_CP_FWD)


def P_cp_bck(E):
    """RREs per metre per backward photon, Compton + photoabsorption (S4)."""
    return _fit(E, *_FIT_CP_BCK)


def inv_positron_attenuation(E):
    """1/lambda_a,pos [1/m] (S9), held below the Bhabha threshold."""
    return _A_POS / (max(E, _FIT_BB[2]) - _E_D_POS)


# --- field and density -------------------------------------------------------


def _N(p, T):
    """Neutral number density [m^-3] at p [bar], T [K]."""
    return p * 1e5 / (kB * T)


def relative_density(p=1.0, T=293.0):
    """n~ = N / N0."""
    return _N(p, T) / _N0


def E_field(EN, p=1.0, T=293.0):
    """Electric field [V/m] from the reduced field EN [Td]."""
    return EN * 1e-21 * _N(p, T)


def _avalanche(E, n):
    """1/lambda [1/m] at field E, density n: (E - n E_d) / U_R, branch at n E_t."""
    return n * inv_avalanche_length(E / n)


# --- interface ----------------------------------------------------------------


def get_V(EN, p=1.0, T=293.0):
    """Diagonal velocity matrix [m/s]: RREs toward +x, positrons toward -x."""
    return np.diag([C_TILDE, -C_TILDE])


def get_R(EN, p=1.0, T=293.0, multipliers=None):
    """
    Reaction matrix [1/s] acting on the densities (n_rre, n_e+).

    Per unit length times c~: RREs avalanche at 1/lambda and are made by
    positrons at P_bb; positrons are lost at 1/lambda_a,pos.
    """
    n = relative_density(p, T)
    E = E_field(EN, p, T)
    return C_TILDE * np.array(
        [
            [_avalanche(E, n), n * P_bb(E)],
            [0.0, -n * inv_positron_attenuation(E)],
        ]
    )


def get_B(EN, p=1.0, T=293.0):
    """
    Production by the photon fluxes, shape (2, 2): RREs from both groups,
    positrons from the forward one by pair production.
    """
    n = relative_density(p, T)
    E = E_field(EN, p, T)
    return n * np.array([[P_cp_fwd(E), P_cp_bck(E)], [P_PP_FWD, 0.0]])


def get_C(EN, p=1.0, T=293.0):
    """Bremsstrahlung source per RRE density, shape (2, 2) [1/s]."""
    n = relative_density(p, T)
    return C_TILDE * n * np.array([[P_BR_FWD, 0.0], [P_BR_BCK, 0.0]])


def get_kappa(p=1.0, T=293.0):
    """Photon attenuation coefficients [1/m], forward and backward group."""
    return relative_density(p, T) / np.array([_LAMBDA_A_FWD, _LAMBDA_A_BCK])


def get_Pi_e():
    """Row selector for the RRE flux, shape (1, 2)."""
    return np.array([[1.0, 0.0]])


def get_Pi_plus():
    """Row selector for the positron flux, shape (1, 2)."""
    return np.array([[0.0, 1.0]])


def get_Pi_minus():
    """Nothing drifts toward +x but RREs: shape (0, 2)."""
    return np.zeros((0, 2))


def get_gamma_plus(EN, p=1.0, T=293.0):
    """No emission at x = 0 by arriving positrons."""
    return np.zeros(1)


def get_gamma_plus_with(EN, p, T, gamma0=None, gamma1=None, eref=None, beta=None):
    """As :func:`get_gamma_plus`; there is no yield to override."""
    return np.zeros(1)


def get_gamma_Psi(EN, p=1.0, T=293.0):
    """No emission at x = 0 by arriving photons."""
    return np.zeros(2)
