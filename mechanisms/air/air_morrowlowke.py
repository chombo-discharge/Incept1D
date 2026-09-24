# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Morrow–Lowke transport model for air, with Zheleznyak photoionization.

The swarm data are the analytic fits of Morrow and Lowke, J. Phys. D: Appl.
Phys. 30, 614 (1997), Appendix: one lumped positive ion, one lumped negative
ion, impact ionization, two- and three-body attachment, and no detachment.
Photoionization and cathode secondary emission are the same as in
``air_pancheshnyi.py`` (the two-stream fit to the Zheleznyak absorption
curve from ``zheleznyak.py``), except that the photon source is driven by
the Morrow–Lowke ionization frequency α|W_e|.

Omitted from the published model, because the inception criterion is the
*linearised* drift-reaction problem:

* Electron-ion and ion-ion recombination — quadratic in the densities, so
  they vanish on linearisation about the neutral gas;
* Electron diffusion — the solver is drift-reaction only
  (see docs/source/theory/overview.rst, eq_drift_reaction).

Units
-----
The published fits use E/N in V cm², N in cm⁻³ and velocities in cm/s.  They
are evaluated in those units inside the ``_ml_*`` helpers and converted to SI
at the public boundary: E/N in Td (1 Td = 1e-17 V cm²), α and η in m⁻¹, drift
velocities in m/s, rates in s⁻¹.

Coordinate convention
---------------------
As in ``air_pancheshnyi.py``: x increases from cathode (x = 0) to anode
(x = d), E points in the −x direction, EN = |E|/N > 0, and

    V[i,i] = -sgn(Z_i) * |W_i|

so electrons and negative ions drift towards +x and positive ions towards −x.

Species index
-------------
0  e      electron
1  M+     positive ion (lumped)
2  M-     negative ion (lumped)
"""

import os
import sys
import argparse
import math

import numpy as np
import matplotlib.pyplot as plt

from incept1d.constants import kB
from incept1d.reactions import (
    compile_reactions as _compile_reactions,
    build_R_from_compiled as _build_R_fast,
)

# Directory of this file: the Zheleznyak helper module is resolved relative to it.
_HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Photoionization constants — identical to air_pancheshnyi.py
# ---------------------------------------------------------------------------

_XI_EXC = globals().get("_XI_EXC", 0.6)  # Excitation+emission efficiency ν_exc
_XI_IONI = globals().get("_XI_IONI", 0.1)  # Photoionization efficiency ξ
_XI_EMIT = globals().get("_XI_EMIT", 0.1)  # Photoemission efficiency ξ
_PQ_BAR = globals().get("_PQ_BAR", 30.0 / 750.064)  # Quenching pressure: 30 Torr → bar

# Ion secondary-emission coefficients — identical to air_pancheshnyi.py
_GAMMA0 = globals().get("_GAMMA0", 1e-3)  # base SEE yield
_GAMMA1 = globals().get("_GAMMA1", 0.0)  # exponential prefactor
_EREF = globals().get("_EREF", 170e7)  # reference field [V/m]
_BETA = globals().get("_BETA", 1.0)  # field scaling exponent

# Mutable module state — set by init_photoionization() at import time
_N_GAMMA = None  # number of photon groups
_kappa_SI = None  # kappa/pO2 [m⁻¹ Pa⁻¹], shape (N_γ,)
_g_groups = None  # photon group fractions, shape (N_γ,), sum = 1
_CONE_FACTOR = None  # ΔΩ/(4π) = (1 − cos θ_cone) / 2


def init_photoionization(ngroups=3, cone_angle_deg=45.0):
    """
    Initialise the two-stream photoionization model.

    Same model as ``air_pancheshnyi.init_photoionization``: absorption
    coefficients κ_j and group fractions g_j from
    ``zheleznyak.fit_twostream``, and the cone factor
    ΔΩ/(4π) = (1 − cos θ) / 2.

    Parameters
    ----------
    ngroups : int
        Number of photon groups (default 3).
    cone_angle_deg : float
        Half-opening angle of the emission cone in degrees (default 45.0).
    """
    global _N_GAMMA, _kappa_SI, _g_groups, _CONE_FACTOR
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    from zheleznyak import fit_twostream

    kappa, g, _ = fit_twostream(ngroups)
    _N_GAMMA = ngroups
    _kappa_SI = kappa
    _g_groups = g
    _CONE_FACTOR = (1.0 - math.cos(cone_angle_deg * math.pi / 180.0)) / 2.0


init_photoionization()

# ---------------------------------------------------------------------------
# Gas composition and species
# ---------------------------------------------------------------------------

xO2 = 0.2  # O2 mole fraction: sets the photon absorption, as in air_pancheshnyi
SPECIES = ["e", "M+", "M-"]
ELECTRON_INDEX = 0
_N_SPECIES = len(SPECIES)

# Reference gas density of the published ion mobilities [cm^-3].  Morrow and
# Lowke give the ion drift velocities as W = mu * E at atmospheric density;
# mu * N is the pressure-independent quantity.
_N0_CM3 = 2.5e19


def _make_N(p=1.0, T=293.0):
    """
    Neutral number density N = p / (kB T) in m⁻³, with p in bar.

    Parameters
    ----------
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in K.
    """
    return p * 1e5 / (kB * T)


# ---------------------------------------------------------------------------
# Morrow–Lowke fits in the published units (E/N in V cm², N in cm⁻³)
# ---------------------------------------------------------------------------


def _ml_We(en):
    """Electron drift speed |W_e| in cm/s; *en* = E/N in V cm²."""
    if en > 2e-15:
        return 7.4e21 * en + 7.1e6
    if en > 1e-16:
        return 1.03e22 * en + 1.3e6
    if en > 2.6e-17:
        return 7.2973e21 * en + 1.63e6
    return 6.87e22 * en + 3.38e4


def _ml_alpha_N(en):
    """Reduced ionization coefficient α/N in cm²; *en* = E/N in V cm²."""
    if en > 1.5e-15:
        return 2e-16 * math.exp(-7.248e-15 / en)
    return 6.619e-17 * math.exp(-5.593e-15 / en)


def _ml_eta2_N(en):
    """
    Reduced two-body attachment coefficient η₂/N in cm²; *en* = E/N in V cm².

    The low-field branch of the published fit turns negative below
    E/N ≈ 4.75e-16 V cm² (47.5 Td); it is clipped at zero there.
    """
    if en > 1.05e-15:
        return 8.889e-5 * en + 2.567e-19
    return max(0.0, 6.089e-4 * en - 2.893e-19)


def _ml_eta3_N2(en):
    """Three-body attachment η₃/N² in cm⁵; *en* = E/N in V cm²."""
    return 4.7778e-59 * en ** (-1.2749)


def _ml_mu_minus(en):
    """Negative-ion mobility in cm²/(V s) at N = _N0_CM3; *en* = E/N in V cm²."""
    return 2.7 if en > 5e-16 else 1.86


_ML_MU_PLUS = 2.34  # Positive-ion mobility in cm²/(V s) at N = _N0_CM3

# ---------------------------------------------------------------------------
# Transport and Townsend coefficients (SI, E/N in Td)
# ---------------------------------------------------------------------------


def ElectronDriftVelocity(EN):
    """
    Electron drift speed |W_e|.

    Parameters
    ----------
    EN : float
        Reduced electric field in Td.

    Returns
    -------
    float
        |W_e| in m/s.
    """
    return 1e-2 * _ml_We(EN * 1e-17)


def alpha(EN, p=1.0, T=293.0):
    """
    Townsend ionization coefficient α.

    Parameters
    ----------
    EN : float
        Reduced electric field in Td.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in K.

    Returns
    -------
    float
        α in m⁻¹.
    """
    return _ml_alpha_N(EN * 1e-17) * 1e-4 * _make_N(p, T)


def eta2(EN, p=1.0, T=293.0):
    """Two-body attachment coefficient η₂ in m⁻¹ (E/N in Td, p in bar, T in K)."""
    return _ml_eta2_N(EN * 1e-17) * 1e-4 * _make_N(p, T)


def eta3(EN, p=1.0, T=293.0):
    """Three-body attachment coefficient η₃ in m⁻¹ (E/N in Td, p in bar, T in K)."""
    N_cm3 = _make_N(p, T) * 1e-6
    return _ml_eta3_N2(EN * 1e-17) * N_cm3**2 * 1e2


def eta(EN, p=1.0, T=293.0):
    """
    Total attachment coefficient η = η₂ + η₃.

    Parameters
    ----------
    EN : float
        Reduced electric field in Td.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in K.

    Returns
    -------
    float
        η in m⁻¹.
    """
    return eta2(EN, p, T) + eta3(EN, p, T)


def _ion_muN(mu_cm2):
    """Reduced mobility μN in m⁻¹ V⁻¹ s⁻¹ from a mobility in cm²/(V s) at N0."""
    return mu_cm2 * 1e-4 * _N0_CM3 * 1e6


# ---------------------------------------------------------------------------
# Reaction list and matrix interface
# ---------------------------------------------------------------------------
# Each rate is a frequency, coefficient × |W_e|, in s⁻¹.  Recombination
# (e + M+, M- + M+) is quadratic in the charged densities and drops out of
# the linearised inception problem, so it is not listed.

REACTIONS = [
    ("e + M -> 2e + M+", lambda EN, p, T: alpha(EN, p, T) * ElectronDriftVelocity(EN)),
    ("e + M -> M-", lambda EN, p, T: eta2(EN, p, T) * ElectronDriftVelocity(EN)),
    ("e + 2M -> M- + M", lambda EN, p, T: eta3(EN, p, T) * ElectronDriftVelocity(EN)),
]

_COMPILED_REACTIONS = _compile_reactions(REACTIONS, SPECIES)


def get_R(EN, p=1.0, T=293.0, multipliers=None):
    """
    Build the 3x3 reaction matrix R in s⁻¹.

    Parameters
    ----------
    EN : float
        Reduced electric field in Td.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in K.
    multipliers : dict, optional
        Per-reaction multipliers keyed by reaction string.

    Returns
    -------
    numpy.ndarray, shape (3, 3)
    """
    return _build_R_fast(_COMPILED_REACTIONS, _N_SPECIES, EN, p, T, multipliers)


def get_V(EN, p=1.0, T=293.0):
    """
    Build the diagonal drift-velocity matrix V in m/s.

    V[i,i] = -sgn(Z_i) |W_i|: electrons |W_e| from the Morrow–Lowke fit,
    ions μN · E/N with the published mobilities scaled to μN at N0.

    Parameters
    ----------
    EN : float
        Reduced electric field in Td.
    p, T : float
        Unused: every velocity here is a function of E/N alone.

    Returns
    -------
    numpy.ndarray, shape (3, 3)
    """
    en_Vm2 = EN * 1e-21
    return np.diag(
        [
            ElectronDriftVelocity(EN),
            -_ion_muN(_ML_MU_PLUS) * en_Vm2,
            _ion_muN(_ml_mu_minus(EN * 1e-17)) * en_Vm2,
        ]
    )


def get_Pi_e():
    """Electron row-selection matrix Π_e, shape (1, 3)."""
    M = np.zeros((1, _N_SPECIES))
    M[0, ELECTRON_INDEX] = 1.0
    return M


def get_Pi_plus():
    """Positive-ion row-selection matrix Π_+, shape (1, 3)."""
    M = np.zeros((1, _N_SPECIES))
    M[0, 1] = 1.0
    return M


def get_Pi_minus():
    """Negative-ion row-selection matrix Π_-, shape (1, 3)."""
    M = np.zeros((1, _N_SPECIES))
    M[0, 2] = 1.0
    return M


def get_gamma_plus(EN, p, T):
    """
    Ion-induced secondary emission yield γ_+, shape (1,).

    γ(E) = γ₀ + γ₁ exp(−E_ref / (β E)), as in air_pancheshnyi.py.
    """
    return get_gamma_plus_with(EN, p, T)


def get_gamma_plus_with(EN, p, T, gamma0=None, gamma1=None, eref=None, beta=None):
    """Like get_gamma_plus, with optional per-call overrides (None → default)."""
    g0 = gamma0 if gamma0 is not None else _GAMMA0
    g1 = gamma1 if gamma1 is not None else _GAMMA1
    er = eref if eref is not None else _EREF
    bt = beta if beta is not None else _BETA
    E = EN * _make_N(p, T) * 1e-21
    g = np.zeros(1)
    if E > 0.0:
        g[0] = g0 + g1 * np.exp(-er / (bt * E))
    return g


# ---------------------------------------------------------------------------
# Photoionization interface — two-stream Zheleznyak model
# ---------------------------------------------------------------------------


def get_kappa(p=1.0, T=293.0):
    """Absorption coefficients κ_j = (κ_j/p_O2) · x_O2 · p in m⁻¹, shape (N_γ,)."""
    return _kappa_SI * xO2 * p * 1e5


def get_B(EN, p=1.0, T=293.0):
    """
    Photon-to-species coupling B, shape (3, N_γ), in m⁻¹.

    Photoionization produces one electron and one positive ion:
    B[e, j] = B[M+, j] = ξ κ_j.
    """
    kappa = get_kappa(p, T)
    B = np.zeros((_N_SPECIES, _N_GAMMA))
    B[0, :] = _XI_IONI * kappa
    B[1, :] = _XI_IONI * kappa
    return B


def get_C(EN, p=1.0, T=293.0):
    """
    Photon-source coupling C_raw, shape (N_γ, 3), in s⁻¹.

    C_raw[j, e] = (ΔΩ/4π) ν_exc p_q/(p + p_q) ν_i g_j, with the ionization
    frequency ν_i = α |W_e| of the Morrow–Lowke fit.  The solver
    right-multiplies by V⁻¹.
    """
    nu_i = alpha(EN, p, T) * ElectronDriftVelocity(EN)
    rho = _XI_EXC * _PQ_BAR / (p + _PQ_BAR) * nu_i
    C = np.zeros((_N_GAMMA, _N_SPECIES))
    C[:, ELECTRON_INDEX] = _CONE_FACTOR * rho * _g_groups
    return C


def get_gamma_Psi(EN, p=1.0, T=293.0):
    """Photoemission yield per backward photon, shape (N_γ,)."""
    return _XI_EMIT * np.ones(_N_GAMMA)


# ---------------------------------------------------------------------------
# __main__: plot the Townsend coefficients and the electron drift velocity
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot the Morrow–Lowke α, η and |W_e| against E/N."
    )
    parser.add_argument("--p", type=float, default=1.0, help="Pressure in bar.")
    parser.add_argument("--T", type=float, default=293.0, help="Temperature in K.")
    args = parser.parse_args()

    EN_arr = np.logspace(1, 3, 500)
    a = np.array([alpha(e, args.p, args.T) for e in EN_arr])
    e2 = np.array([eta2(e, args.p, args.T) for e in EN_arr])
    e3 = np.array([eta3(e, args.p, args.T) for e in EN_arr])
    We = np.array([ElectronDriftVelocity(e) for e in EN_arr])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9))
    ax1.loglog(EN_arr, a, label=r"$\alpha$")
    ax1.loglog(EN_arr, e2, label=r"$\eta_2$")
    ax1.loglog(EN_arr, e3, label=r"$\eta_3$")
    ax1.set_ylabel(r"Townsend coefficient (m$^{-1}$)")
    ax1.set_title(f"Morrow–Lowke air,  p = {args.p} bar,  T = {args.T} K")
    ax1.legend()
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax2.loglog(EN_arr, We)
    ax2.set_xlabel("E/N (Td)")
    ax2.set_ylabel(r"$|W_e|$ (m/s)")
    ax2.grid(True, which="both", ls="--", alpha=0.4)
    plt.tight_layout()
    plt.show()
