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
The paper gives every coefficient at the field magnitude E and the air
density relative to sea level, n~, with every threshold field scaling as n~
(E_t = 284 kV/m n~).  Here each coefficient is the paper's function of the
sea-level-equivalent field E / n~, times n~, so that it depends only on E/N
and the inception curve only on p d, as for every other mechanism.  At
n~ = 1 this is the paper's expression exactly.  The sea-level number density
N0 that sets n~ = N / N0, and with it the conversion between E/N and E,
defaults to 1 atm and 293.15 K (``config.py``).

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
    lambda_a,fwd  = 307.64 m,  lambda_a,bck = 92.22 m (R >> L)       (S7)

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
_N0 = globals().get("_N0", 101325.0 / (kB * 293.15))
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


def sea_level_field(EN):
    """E / n~ [V/m], the field at sea-level density with the same E/N."""
    return EN * 1e-21 * _N0


def E_field(EN, p=1.0, T=293.0):
    """Electric field [V/m] from the reduced field EN [Td]."""
    return EN * 1e-21 * _N(p, T)


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
    E = sea_level_field(EN)
    return (
        C_TILDE
        * n
        * np.array(
            [
                [inv_avalanche_length(E), P_bb(E)],
                [0.0, -inv_positron_attenuation(E)],
            ]
        )
    )


def get_B(EN, p=1.0, T=293.0):
    """
    Production by the photon fluxes, shape (2, 2): RREs from both groups,
    positrons from the forward one by pair production.
    """
    n = relative_density(p, T)
    E = sea_level_field(EN)
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


# --- normal modes, for the documentation and the example ---------------------


def mode_matrix(E, n=1.0):
    """
    M of dn/dz = M n for the steady state, n = (n_r, n_pos, n_up, n_down).

    The model equations with the time derivatives dropped, in the paper's
    density variables, at field E [V/m] and relative density n.
    """
    Er = E / n
    return n * np.array(
        [
            [inv_avalanche_length(Er), P_bb(Er), P_cp_fwd(Er), P_cp_bck(Er)],
            [0.0, inv_positron_attenuation(Er), -P_PP_FWD, 0.0],
            [P_BR_FWD, 0.0, -1.0 / _LAMBDA_A_FWD, 0.0],
            [-P_BR_BCK, 0.0, 0.0, 1.0 / _LAMBDA_A_BCK],
        ]
    )


def normal_mode_threshold(L, n=1.0, E_lo=285e3, E_hi=6e6, n_scan=2000):
    """
    Self-sustained field E_s [V/m] of a region of length L [m], by normal modes.

    The steady state is a sum of modes v_j e^{kappa_j z}.  With nothing
    entering, n_r(0) = n_up(0) = n_pos(L) = n_down(L) = 0, which a
    non-trivial combination satisfies where the determinant of the modes'
    boundary values vanishes.  Divided by det [v_1 ... v_4] that determinant
    is real and does not depend on how the modes are scaled or ordered.
    Returns the lowest sign change in [E_lo, E_hi] (n~ E_lo ... n~ E_hi), or
    NaN if there is none.
    """
    import scipy.optimize

    def g(E):
        kappa, W = np.linalg.eig(mode_matrix(E, n))
        shift = np.maximum(kappa.real, 0.0) * L
        s, at_L = np.exp(-shift), np.exp(kappa * L - shift)
        A = np.vstack([W[0] * s, W[2] * s, W[1] * at_L, W[3] * at_L])
        return float(np.real(np.linalg.det(A) / np.linalg.det(W)))

    grid = n * np.linspace(E_lo, E_hi, n_scan)
    values = [g(E) for E in grid]
    for a, b, ga, gb in zip(grid, grid[1:], values, values[1:]):
        if ga * gb < 0.0:
            return scipy.optimize.brentq(g, a, b, xtol=1e-10, rtol=1e-14)
    return float("nan")
