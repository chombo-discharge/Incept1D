# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Reaction mechanism for dry air (80% N2, 20% O2).

This module defines all rate coefficients, transport properties, and matrix
interface functions needed for 1-D drift-diffusion discharge modelling of air.
The public interface (get_R, get_V, get_Pi_e, get_Pi_plus, get_Pi_minus,
get_gamma_plus, get_B, get_C, get_kappa, get_gamma_Psi,
SPECIES, ELECTRON_INDEX) is consumed by Eigenvalues.py and Inception.py.

Coordinate convention
---------------------
The spatial variable x increases from cathode (x = 0) to anode (x = d).
The electric field E points in the -x direction (from anode toward cathode,
consistent with E = -dV/dx for a cathode at lower potential).  The signed
field is therefore E_signed = -|E|, but EN = |E|/N > 0 throughout.  The
drift-velocity matrix V is constructed so that

    V[i,i] = -sgn(Z_i) * mu_i * |E|

giving electrons (Z = -1) a positive velocity toward the anode and positive
ions (Z = +1) a negative velocity toward the cathode.  With A = R V^{-1},
positive eigenvalues of A correspond to spatially growing solutions (net
ionisation exceeds attachment).

Species index
-------------
0  e          electron
1  N2+        molecular nitrogen cation
2  O2+        molecular oxygen cation
3  O-         atomic oxygen anion
4  O2-(exc)   vibrationally excited molecular oxygen anion
5  O2-        molecular oxygen anion
6  O3-*       excited ozone anion
7  O3-        ozone anion
"""

import os
import argparse

import numpy as np
import math
import matplotlib.pyplot as plt

from incept1d.mechanism import load_helper
from incept1d.constants import kB, Q
from incept1d.reactions import (
    compile_reactions as _compile_reactions,
    build_R_from_compiled as _build_R_fast,
)

# Directory of this file: data tables (BOLSIG+ output, mobilities) and the
# Zheleznyak helper module are resolved relative to it.
_HERE = os.path.dirname(os.path.abspath(__file__))
# Shared by the air mechanisms: the photoionization fit, and the LXCat data.
_AIR = os.path.dirname(_HERE)
_LXCAT = os.path.join(_AIR, "lxcat")

# ---------------------------------------------------------------------------
# Photoionization constants — two-stream model (Zheleznyak absorption curve)
# ---------------------------------------------------------------------------

_XI_EXC = globals().get("_XI_EXC", 0.6)  # Excitation+emission efficiency ν_exc
_XI_IONI = globals().get("_XI_IONI", 0.1)  # Photoionization efficiency ξ
_XI_EMIT = globals().get("_XI_EMIT", 0.1)  # Photoemission efficiency ξ
_PQ_BAR = globals().get("_PQ_BAR", 30.0 / 750.064)  # Quenching pressure: 30 Torr → bar

# Ion secondary-emission coefficients — configurable via globals().get() injection
_GAMMA0 = globals().get("_GAMMA0", 1e-6)  # base SEE yield (N2+, O2+)
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

    Calls Zheleznyak.fit_twostream to obtain absorption coefficients κ_j and
    group fractions g_j, then computes the spatial cone factor ΔΩ/(4π) from
    the half-opening angle.  Called automatically at module import with
    defaults; call again to change the model parameters at run time.

    Parameters
    ----------
    ngroups : int
        Number of photon groups (default 3).
    cone_angle_deg : float
        Half-opening angle of the emission cone in degrees (default 45.0).
        Translated to ΔΩ/(4π) = (1 − cos θ) / 2.
    """
    global _N_GAMMA, _kappa_SI, _g_groups, _CONE_FACTOR
    fit_twostream = load_helper(os.path.join(_AIR, "zheleznyak.py")).fit_twostream

    kappa, g, _ = fit_twostream(ngroups)
    _N_GAMMA = ngroups
    _kappa_SI = kappa  # kappa/pO2 [m^-1 Pa^-1] from Zheleznyak fit
    _g_groups = g
    cone_rad = cone_angle_deg * math.pi / 180.0
    _CONE_FACTOR = (1.0 - math.cos(cone_rad)) / 2.0


init_photoionization()

# ---------------------------------------------------------------------------
# Gas composition
# ---------------------------------------------------------------------------

xO2 = 0.21  # Mole fraction of O2
xN2 = 0.79  # Mole fraction of N2

# ---------------------------------------------------------------------------
# Ion transport
# ---------------------------------------------------------------------------

_N_1bar = 1e5 / (kB * 300.0)  # Neutral density at 1 bar, 300 K [m^-3]
ion_muN = 2e-4 * _N_1bar  # N2+, O2+ reduced mobility mu*N [m^-1 V^-1 s^-1]
_muN_Om = 1.2e22  # O-           reduced mobility mu*N [m^-1 V^-1 s^-1]
m_O3 = 48 * 1.66053906660e-27  # Ozone (O3) mass [kg]

# ---------------------------------------------------------------------------
# Species list and electron index
# ---------------------------------------------------------------------------

SPECIES = ["e", "N2+", "O2+", "O-", "O2-(exc)", "O2-", "O3-*", "O3-"]
ELECTRON_INDEX = 0
_N_SPECIES = len(SPECIES)

# ---------------------------------------------------------------------------
# BOLSIG+ output file parser
# ---------------------------------------------------------------------------


def _load_bolsig(path):
    """
    Parse a BOLSIG+ output file and return six two-column tables.

    The file is expected to have two sections separated by CR line endings:

    Section 1 — transport coefficients (20 fixed columns):
        col 1 = E/N (Td), col 2 = A1 mean energy (eV),
        col 3 = A2 mobility*N (1/m/V/s), col 4 = A6 diffusion*N (1/m/s)

    Section 2 — per-process rate coefficients (variable columns):
        A legend maps each Cx label to a species and process type.
        All N2 Ionization columns are summed to give k1, all O2 Ionization
        columns to give k2, and all O2 Attachment columns to give k3.

    Returns
    -------
    tuple of six ndarray, each shaped (n_pts, 2) with columns [E/N, value]:
        energy_table, mobility_table, diffusion_table,
        n2ionize_table, o2ionize_table, o2dissociate_table
    """
    with open(path, "rb") as f:
        raw = f.read()
    lines = [ln.decode("ascii", errors="replace").strip() for ln in raw.split(b"\r")]

    # --- Section 1 ---
    sec1_rows, in_sec1 = [], False
    for ln in lines:
        if "R#" in ln and "A1" in ln and "A2" in ln:
            in_sec1 = True
            continue
        if in_sec1:
            if not ln or "Rate coefficients" in ln:
                break
            vals = ln.split()
            if vals and vals[0].lstrip("-").isdigit():
                sec1_rows.append([float(v) for v in vals])
    sec1 = np.array(sec1_rows)
    EN1 = sec1[:, 1]
    energy = sec1[:, 2]
    muN = sec1[:, 3]
    DN = sec1[:, 4]

    # --- Section 2 legend ---
    legend = {}  # cx_num (int) -> (species, proc_type)
    in_legend = False
    sec2_hdr_idx = None
    for i, ln in enumerate(lines):
        if "Rate coefficients (m3/s)" in ln and "Inverse" not in ln:
            in_legend = True
            continue
        if in_legend:
            if "R#" in ln and "E/N" in ln:
                sec2_hdr_idx = i
                in_legend = False
                break
            parts = ln.split()
            if parts and parts[0].startswith("C") and parts[0][1:].isdigit():
                cx = int(parts[0][1:])
                sp = parts[1] if len(parts) > 1 else ""
                pt = parts[2] if len(parts) > 2 else ""
                legend[cx] = (sp, pt)

    # --- Section 2 data ---
    sec2_rows = []
    for ln in lines[sec2_hdr_idx + 1 :]:
        if not ln or "Inverse" in ln:
            break
        vals = ln.split()
        if not vals or not vals[0].lstrip("-").isdigit():
            break
        sec2_rows.append([float(v) for v in vals])
    sec2 = np.array(sec2_rows)
    EN2 = sec2[:, 1]

    # col index for process Cx in sec2: R#=0, E/N=1, Energy=2, C1=3 → Cx at x+2
    k1_t = np.zeros(len(sec2))
    k2_t = np.zeros(len(sec2))
    k3_t = np.zeros(len(sec2))
    for cx, (sp, pt) in legend.items():
        col = cx + 2
        if col >= sec2.shape[1]:
            continue
        if sp == "N2" and pt == "Ionization":
            k1_t += sec2[:, col]
        elif sp == "O2" and pt == "Ionization":
            k2_t += sec2[:, col]
        elif sp == "O2" and pt == "Attachment":
            k3_t += sec2[:, col]

    return (
        np.column_stack([EN1, energy]),
        np.column_stack([EN1, muN]),
        np.column_stack([EN1, DN]),
        np.column_stack([EN2, k1_t]),
        np.column_stack([EN2, k2_t]),
        np.column_stack([EN2, k3_t]),
    )


def _load_lxcat_mobility(path):
    """
    Parse an LXCat-format ion mobility file and return a table of
    [EN_Td, muN] where muN = mu*N in m^-1 V^-1 s^-1.

    Ko (cm^2/(V*s)) from the file is converted via muN = Ko * 1e-4 * N0,
    where N0 = 2.6868e25 m^-3 is the Loschmidt constant (0 degC, 1 atm).
    Data between the first pair of '----' separator lines is read;
    np.interp on the resulting table clamps out-of-range queries to the
    nearest boundary value.
    """
    _N0 = 2.6868e25  # Loschmidt constant [m^-3]
    rows = []
    in_data = False
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("----"):
                if not in_data:
                    in_data = True
                else:
                    break
                continue
            if in_data:
                parts = stripped.split()
                if len(parts) == 2:
                    EN_td = float(parts[0])
                    Ko = float(parts[1])  # cm^2/(V*s)
                    rows.append([EN_td, Ko * 1e-4 * _N0])
    return np.array(rows)


# ---------------------------------------------------------------------------
# Cross-section and transport data tables (loaded once at import time)
# ---------------------------------------------------------------------------

BOLSIG_FILE = globals().get("BOLSIG_FILE", os.path.join(_LXCAT, "phelps.txt"))

(
    energy_table,
    mobility_table,
    diffusion_table,
    n2ionize_table,
    o2ionize_table,
    o2dissociate_table,
) = _load_bolsig(BOLSIG_FILE)

_o2m_mobility_table = _load_lxcat_mobility(os.path.join(_LXCAT, "o2m_mobility.txt"))
_o3m_mobility_table = _load_lxcat_mobility(os.path.join(_LXCAT, "o3m_mobility.txt"))


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _make_N(p=1.0, T=293.0):
    """
    Compute the neutral number density from pressure and temperature.

    Uses the ideal gas law N = p / (kB * T) with p expressed in bar
    (1 bar = 1e5 Pa).

    Parameters
    ----------
    p : float
        Gas pressure in bar.  Default 1.0 bar.
    T : float
        Gas temperature in Kelvin.  Default 293.0 K.

    Returns
    -------
    float
        Neutral number density N in m^-3.
    """
    return p * 1e5 / (kB * T)


# ---------------------------------------------------------------------------
# Transport coefficient functions
# ---------------------------------------------------------------------------


def ElectronMeanEnergy(EN):
    """
    Return the mean electron energy as a function of the reduced electric field.

    Interpolated from the tabulated EEDF data in energy.dat.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend (1 Td = 1e-21 V m^2).

    Returns
    -------
    float or ndarray
        Mean electron energy in eV.
    """
    return np.interp(EN, energy_table[:, 0], energy_table[:, 1])


def ElectronTemperature(EN):
    """
    Return the electron temperature in Kelvin.

    Derived from the mean energy via T_e = (2/3) * <epsilon> * Q / kB,
    which assumes a Maxwellian EEDF.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        Electron temperature in Kelvin.
    """
    return 2.0 * ElectronMeanEnergy(EN) * Q / (3.0 * kB)


def ElectronMobility(EN):
    """
    Return the reduced electron mobility mu_e * N.

    Interpolated from mobility.dat.  The product mu_e * N is nearly independent
    of gas pressure, making it the natural quantity to tabulate.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        mu_e * N in m^-1 V^-1 s^-1.
    """
    return np.interp(EN, mobility_table[:, 0], mobility_table[:, 1])


def O2mMobility(EN):
    """
    Return the reduced O2- (and O2-(exc)) mobility mu * N.

    Interpolated from o2m_mobility.txt (Viehland/LXCat).  Queries outside
    the table range [5, 250] Td are clamped to the nearest boundary value.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        mu * N in m^-1 V^-1 s^-1.
    """
    return np.interp(EN, _o2m_mobility_table[:, 0], _o2m_mobility_table[:, 1])


def O3mMobility(EN):
    """
    Return the reduced O3- (and O3-*) mobility mu * N.

    Interpolated from o3m_mobility.txt (Viehland/LXCat).  Queries outside
    the table range [5, 200] Td are clamped to the nearest boundary value.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        mu * N in m^-1 V^-1 s^-1.
    """
    return np.interp(EN, _o3m_mobility_table[:, 0], _o3m_mobility_table[:, 1])


def ElectronDiffusion(EN):
    """
    Return the reduced electron diffusion coefficient D_e * N.

    Interpolated from diffusion.dat.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        D_e * N in m^-1 s^-1.
    """
    return np.interp(EN, diffusion_table[:, 0], diffusion_table[:, 1])


# ---------------------------------------------------------------------------
# Rate coefficient functions k1 – k12
# ---------------------------------------------------------------------------


def k1(EN):
    """
    Rate coefficient for electron-impact ionisation of N2.

    Reaction: e + N2 -> e + e + N2+

    Summed from all N2 Ionization processes in BOLSIG_FILE.  Multiply by
    the N2 number density [N2] = xN2 * N to obtain the volumetric rate in s^-1.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        Rate coefficient in m^3 s^-1.
    """
    return np.interp(EN, n2ionize_table[:, 0], n2ionize_table[:, 1])


def k2(EN):
    """
    Rate coefficient for electron-impact ionisation of O2.

    Reaction: e + O2 -> e + e + O2+

    Summed from all O2 Ionization processes in BOLSIG_FILE.  Multiply by
    [O2] = xO2 * N.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        Rate coefficient in m^3 s^-1.
    """
    return np.interp(EN, o2ionize_table[:, 0], o2ionize_table[:, 1])


def k3(EN):
    """
    Rate coefficient for dissociative attachment of electrons to O2.

    Reaction: e + O2 -> O- + O  (plus any other O2 attachment channels in the database)

    Summed from all O2 Attachment processes in BOLSIG_FILE.  Multiply by
    [O2] = xO2 * N.

    Parameters
    ----------
    EN : float or array-like
        Reduced electric field in Townsend.

    Returns
    -------
    float or ndarray
        Rate coefficient in m^3 s^-1.
    """
    return np.interp(EN, o2dissociate_table[:, 0], o2dissociate_table[:, 1])


def k45(EN, T=293.0):
    """
    Three-body rate coefficient for O2- formation.

    Reaction: e + O2 + O2 -> O2- + O2

    Formula from Kossyi et al.  The three-body nature means the returned
    coefficient must be multiplied by [O2]^2 = (xO2 * N)^2 to obtain
    the volumetric rate in s^-1.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    T : float
        Gas temperature in Kelvin.  Default 293.0 K.

    Returns
    -------
    float
        Three-body rate coefficient in m^6 s^-1.
    """
    Te = ElectronTemperature(EN)
    return (
        1.4e-29
        * (300.0 / Te)
        * math.exp(-600.0 / T)
        * math.exp(700.0 * (Te - T) / (Te * T))
        * 1e-12
    )


def k4(EN, T=293.0):
    """
    Rate coefficient for electron attachment forming vibrationally excited O2-.

    Reaction: e + O2 -> O2-(exc)

    Reverse-engineered from the Kossyi three-body rate k45 using the
    quasi-steady-state balance between O2-(exc) autodetachment (k5) and
    collisional de-excitation (k6):  k4 = k45 * k5 / k6.  Multiply by
    [O2] = xO2 * N.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    T : float
        Gas temperature in Kelvin.  Default 293.0 K.

    Returns
    -------
    float
        Rate coefficient in m^3 s^-1.
    """
    return k45(EN, T) * k5(EN) / k6(EN)


def k5(EN):
    """
    Autodetachment rate from O2-(exc).

    Reaction: O2-(exc) -> e + O2

    This is a unimolecular rate (not a rate coefficient).  It is constant
    and independent of E/N.  The EN parameter is accepted for interface
    consistency only.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend (unused).

    Returns
    -------
    float
        Autodetachment rate in s^-1.
    """
    return 1e10


def k6(EN):
    """
    Rate coefficient for collisional de-excitation of O2-(exc) by O2.

    Reaction: O2-(exc) + O2 -> O2- + O2

    Constant and independent of E/N.  Multiply by [O2] = xO2 * N.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend (unused).

    Returns
    -------
    float
        Rate coefficient in m^3 s^-1.
    """
    return 1e-15


def k7(EN):
    """
    Rate coefficient for collisional detachment from O2-.

    Reaction: O2- + M -> O2 + O2 + e

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.

    Returns
    -------
    float
        Rate coefficient in m^3 s^-1.
    """
    return 1.24e-17 * math.exp(-pow(179.0 / (8.8 + EN), 2))


def k8(EN):
    """
    Rate coefficient for ion conversion from O- to O2-.

    Reaction: O- + O2 -> O + O2-

    Exponential-Gaussian formula from Kossyi et al.  Multiply by [O2] = xO2 * N.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.

    Returns
    -------
    float
        Rate coefficient in m^3 s^-1.
    """
    return 6.96e-17 * math.exp(-pow(198.0 / (5.6 + EN), 2))


def k9(EN):
    """
    Rate coefficient for excited ozone-ion formation from O-.

    Reaction: O- + O2 -> O3-*

    The three-body association constant 1.1e-42 m^6 s^-1 is pre-multiplied
    by k10/k11 to account for the branching between autodetachment and
    collisional stabilisation, folding the three-body dependence into an
    effective two-body rate.  Multiply by [O2] = xO2 * N.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.

    Returns
    -------
    float
        Effective two-body rate coefficient in m^3 s^-1.
    """
    return 1.1e-42 * math.exp(-pow(EN / 65.0, 2)) * k10(EN) / k11(EN)


def k10(EN):
    """
    Autodetachment rate from O3-*.

    Reaction: O3-* -> O- + O2

    Unimolecular rate, constant and independent of E/N.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend (unused).

    Returns
    -------
    float
        Autodetachment rate in s^-1.
    """
    return 1e10


def k11(EN):
    """
    Rate coefficient for collisional stabilisation of O3-*.

    Reaction: O3-* + M -> O3- + M*

    Constant and independent of E/N.  Multiply by the total neutral density N.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend (unused).

    Returns
    -------
    float
        Rate coefficient in m^3 s^-1.
    """
    return 1e-15


def k12(EN, T=293.0):
    """
    Rate coefficient for ozone-ion detachment by O2.

    Reaction: O3- + O2 -> O- + O2

    The rate depends on temperature through the mean ion kinetic energy,
    which includes both thermal (3/2 kB T) and directed drift contributions
    (pi/4 * m_O3 * u_drift^2).  Multiply by [O2] = xO2 * N.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    T : float
        Gas temperature in Kelvin.  Default 293.0 K.

    Returns
    -------
    float
        Rate coefficient in m^3 s^-1.
    """
    u = O3mMobility(EN) * EN * 1e-21  # O3- drift speed [m/s]
    # Mean ion energy [eV].  Kept for the disabled rate expression below.
    eps = (1.5 * kB * T + 0.25 * math.pi * m_O3 * u**2) / Q  # noqa: F841

    #    return 1E-18 * math.exp(-1.5 / eps)
    return 0.0


def alpha(EN, p=1.0, T=293.0):
    """
    Return the Townsend first ionisation coefficient alpha.

    Computed as alpha = (k1 [N2] + k2 [O2]) / (mu_e E), where
    E = EN * N * 1e-21 and N = _make_N(p, T).

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    p : float
        Gas pressure in bar.  Default 1.0 bar.
    T : float
        Gas temperature in Kelvin.  Default 293.0 K.

    Returns
    -------
    float
        Townsend ionisation coefficient in m^-1.
    """
    N = _make_N(p, T)
    K1 = k1(EN) * xN2 * N
    K2 = k2(EN) * xO2 * N
    return (K1 + K2) / (ElectronMobility(EN) * EN * 1e-21)


def eta(EN, p=1.0, T=293.0):
    """
    Return the effective Townsend attachment coefficient eta.

    Treats O2-(exc) as a permanent electron sink (autodetachment ignored),
    so the full K4 rate contributes to electron loss alongside K3.  This
    gives the net electron loss rate per electron divided by the drift
    velocity, consistent with disabling the O2-(exc) autodetachment channel (k5).

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    p : float
        Gas pressure in bar.  Default 1.0 bar.
    T : float
        Gas temperature in Kelvin.  Default 293.0 K.

    Returns
    -------
    float
        Effective Townsend attachment coefficient in m^-1.
    """
    N = _make_N(p, T)
    K3 = k3(EN) * xO2 * N
    K4 = k4(EN, T) * xO2 * N
    return (K3 + K4) / (ElectronMobility(EN) * EN * 1e-21)


# ---------------------------------------------------------------------------
# Reaction list and matrix interface — required by Eigenvalues.py and Inception.py
# ---------------------------------------------------------------------------

REACTIONS = [
    ("e + N2 -> 2e + N2+", lambda EN, p, T: k1(EN) * xN2 * _make_N(p, T)),
    ("e + O2 -> 2e + O2+", lambda EN, p, T: k2(EN) * xO2 * _make_N(p, T)),
    ("e + O2 -> O- + O", lambda EN, p, T: k3(EN) * xO2 * _make_N(p, T)),
    ("e + O2 -> O2-(exc)", lambda EN, p, T: k4(EN, T) * xO2 * _make_N(p, T)),
    ("O2-(exc) -> e + O2", lambda EN, p, T: k5(EN)),
    ("O2-(exc) + O2 -> O2- + O2", lambda EN, p, T: k6(EN) * xO2 * _make_N(p, T)),
    ("O2- + O2 -> e + 2O2", lambda EN, p, T: k7(EN) * _make_N(p, T)),
    ("O- + O2 -> O + O2-", lambda EN, p, T: k8(EN) * xO2 * _make_N(p, T)),
    ("O- + O2 -> O3-*", lambda EN, p, T: k9(EN) * xO2 * _make_N(p, T)),
    ("O3-* -> O- + O2", lambda EN, p, T: k10(EN)),
    ("O3-* + M -> O3- + M", lambda EN, p, T: k11(EN) * _make_N(p, T)),
    ("O3- + O2 -> O- + 2O2", lambda EN, p, T: k12(EN, T) * xO2 * _make_N(p, T)),
]

_COMPILED_REACTIONS = _compile_reactions(REACTIONS, SPECIES)
_N_SPECIES_INT = len(SPECIES)


def get_R(EN, p=1.0, T=293.0, multipliers=None):
    """
    Build and return the 8x8 reaction matrix R.

    R[i, j] is the rate [s^-1] at which one carrier of species j produces
    (positive) or destroys (negative, diagonal) one particle of species i.
    The matrix is assembled automatically from REACTIONS by Reactions.build_R.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    p : float
        Gas pressure in bar.  Default 1.0 bar.
    T : float
        Gas temperature in Kelvin.  Default 293.0 K.
    multipliers : dict, optional
        Per-reaction scalar multipliers keyed by reaction string (spaces
        optional, e.g. `"e+N2->2e+N2+"`).  Missing keys default to 1.0.
        Use 0.0 to disable a reaction, >1.0 to amplify.  An unknown key
        triggers a printed warning and is ignored.

    Returns
    -------
    numpy.ndarray, shape (8, 8)
        Reaction matrix R with entries in s^-1.
    """
    return _build_R_fast(_COMPILED_REACTIONS, _N_SPECIES_INT, EN, p, T, multipliers)


def get_V(EN, p=1.0, T=293.0):
    """
    Build and return the diagonal drift-velocity matrix V.

    The drift velocity of species i is V[i,i] = -sgn(Z_i) * mu_i * |E|,
    where the minus sign encodes the coordinate convention (x from cathode
    to anode, E in the -x direction).  With EN = |E|/N in Townsend and
    using mu_i * N = muiN:

        V[i,i] = -sgn(Z_i) * muiN * EN * 1e-21   [m/s]

    Electrons (Z = -1) obtain positive velocity toward the anode; positive
    ions (Z = +1) obtain negative velocity toward the cathode; negative ions
    (Z = -1) obtain positive velocity toward the anode.

    The electron reduced mobility is field-dependent (BOLSIG+ table).  Ion
    reduced mobilities are species-specific:
        N2+, O2+      → ion_muN (constant)
        O-            → _muN_Om (constant)
        O2-(exc), O2- → O2mMobility(EN) (interpolated from Viehland/LXCat)
        O3-*, O3-     → O3mMobility(EN) (interpolated from Viehland/LXCat)
    p and T are accepted for interface consistency but do not affect V
    in this implementation (mu*N is pressure-independent).

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    p : float
        Gas pressure in bar.  Default 1.0 bar (unused in this mechanism).
    T : float
        Gas temperature in Kelvin.  Default 293.0 K (unused in this mechanism).

    Returns
    -------
    numpy.ndarray, shape (8, 8)
        Diagonal drift-velocity matrix V in m/s.
    """
    Z_diag = np.array([-1.0, +1.0, +1.0, -1.0, -1.0, -1.0, -1.0, -1.0])
    muN = np.array(
        [
            ElectronMobility(EN),  # 0: e
            ion_muN,  # 1: N2+
            ion_muN,  # 2: O2+
            _muN_Om,  # 3: O-
            O2mMobility(EN),  # 4: O2-(exc)
            O2mMobility(EN),  # 5: O2-
            O3mMobility(EN),  # 6: O3-*
            O3mMobility(EN),  # 7: O3-
        ]
    )
    return np.diag(-Z_diag * muN * EN * 1e-21)


# ---------------------------------------------------------------------------
# Row-selection matrices — required by the determinant inception criterion
# ---------------------------------------------------------------------------


def get_Pi_e():
    """
    Return the electron row-selection matrix Π_e.

    Π_e[0, ELECTRON_INDEX] = 1, all other entries zero.  Used in the cathode
    boundary condition: the electron flux row of the boundary matrix Q_0.

    Returns
    -------
    numpy.ndarray, shape (1, 8)
        One-hot row matrix selecting the electron species.
    """
    M = np.zeros((1, _N_SPECIES))
    M[0, ELECTRON_INDEX] = 1.0
    return M


def get_Pi_plus():
    """
    Return the positive-ion row-selection matrix Π_+.

    Selects N2+ (index 1) and O2+ (index 2).  Row 0 → N2+, row 1 → O2+.
    Used to impose the zero positive-ion flux condition at the anode.

    Returns
    -------
    numpy.ndarray, shape (2, 8)
        Two one-hot rows, one per positive-ion species.
    """
    M = np.zeros((2, _N_SPECIES))
    M[0, 1] = 1.0  # N2+
    M[1, 2] = 1.0  # O2+
    return M


def get_Pi_minus():
    """
    Return the negative-ion row-selection matrix Π_-.

    Selects O- (3), O2-(exc) (4), O2- (5), O3-* (6), O3- (7).
    Used to impose the zero negative-ion flux condition at the cathode.

    Returns
    -------
    numpy.ndarray, shape (5, 8)
        Five one-hot rows, one per negative-ion species.
    """
    M = np.zeros((5, _N_SPECIES))
    for k, idx in enumerate([3, 4, 5, 6, 7]):
        M[k, idx] = 1.0
    return M


def get_gamma_plus(EN, p, T):
    """
    Return the secondary emission coefficient vector γ_+ for positive ions.

    γ_+[i] is the mean number of secondary electrons emitted per incident
    particle of positive-ion species i, in the same row order as get_Pi_plus()
    (N2+ first, O2+ second).  The field-dependent formula is identical to the
    legacy get_gamma, restricted to the two cation species:

        γ_i(E) = γ₀ + γ₁ * exp(−E_ref / (β * E))

    where E = EN * N * 1e-21 V/m.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.

    Returns
    -------
    numpy.ndarray, shape (2,)
        Secondary emission yields for [N2+, O2+].
    """
    N = _make_N(p, T)
    E = EN * N * 1e-21
    g = np.zeros(2)
    if E > 0.0:
        factor = np.exp(-_EREF / (_BETA * E))
        g[0] = _GAMMA0 + _GAMMA1 * factor  # N2+
        g[1] = _GAMMA0 + _GAMMA1 * factor  # O2+
    return g


def get_gamma_plus_with(EN, p, T, gamma0=None, gamma1=None, eref=None, beta=None):
    """Like get_gamma_plus but with optional per-call overrides of the four
    gamma parameters.  None → use the module-level defaults (_GAMMA0 etc.)."""
    g0 = gamma0 if gamma0 is not None else _GAMMA0
    g1 = gamma1 if gamma1 is not None else _GAMMA1
    er = eref if eref is not None else _EREF
    bt = beta if beta is not None else _BETA
    N = _make_N(p, T)
    E = EN * N * 1e-21
    g = np.zeros(2)
    if E > 0.0:
        factor = np.exp(-er / (bt * E))
        g[0] = g[1] = g0 + g1 * factor
    return g


# ---------------------------------------------------------------------------
# Photoionization interface — three-group Eddington (SP₁) model
# Bourdon et al. (2007) PSST 16 656
# ---------------------------------------------------------------------------


def get_B(EN, p=1.0, T=293.0):
    """
    Return the photon-to-species coupling matrix B, shape (N_s, N_γ).

    In the two-stream model the augmented ODE is ∂_x y = A y + B(Ψ⁺ + Ψ⁻),
    where Ψ⁺ and Ψ⁻ are photon fluxes [m⁻²s⁻¹].  The photoionisation rate
    per unit volume is ξ · κ_j · Ψ_j, giving:

        B[i,j] = β_j[i] · ξ · κ_j   [m⁻¹]

    UV photons from the N₂ second-positive system ionise O₂, producing equal
    numbers of electrons (row 0) and O₂⁺ (row 2).  All other rows are zero.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend (unused; retained for interface
        consistency).
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin (unused; retained for interface consistency).

    Returns
    -------
    numpy.ndarray, shape (N_s, N_γ)
    """
    kappa = get_kappa(p, T)
    B = np.zeros((_N_SPECIES, _N_GAMMA))
    B[0, :] = _XI_IONI * kappa  # electrons
    B[2, :] = _XI_IONI * kappa  # O₂⁺
    return B


def get_C(EN, p=1.0, T=293.0):
    """
    Return the photon-source coupling matrix C_raw, shape (N_γ, N_s).

    C_raw[j, i] couples species-i density n_i to the source term in the
    two-stream equation for photon group j.  Inception.py right-multiplies
    by V⁻¹ to convert to flux space before inserting into A_aug.

    From the manuscript (Eq. 253):

        C = (ΔΩ/4π) η g u_eᵀ V⁻¹

    so C_raw (before the V⁻¹ factor) has a single non-zero column at the
    electron index:

        C_raw[j, 0] = (ΔΩ/4π) · η · g_j   > 0

    where η = ν_exc · (p_q/(p+p_q)) · K_ion is the total photon emission
    rate per electron and ΔΩ/(4π) = (1 − cos θ_cone) / 2.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.

    Returns
    -------
    numpy.ndarray, shape (N_γ, N_s)
        C_raw acting on species densities n.  The caller must right-multiply
        by V⁻¹ to obtain the flux-space matrix needed by A_aug.
    """
    N = _make_N(p, T)
    xi_eff = _XI_EXC * _PQ_BAR / (p + _PQ_BAR)
    K_ion = k1(EN) * xN2 * N + k2(EN) * xO2 * N
    eta = xi_eff * K_ion
    C = np.zeros((_N_GAMMA, _N_SPECIES))
    C[:, ELECTRON_INDEX] = _CONE_FACTOR * eta * _g_groups
    return C


def get_kappa(p=1.0, T=293.0):
    """
    Return the two-stream absorption coefficient vector κ, shape (N_γ,).

    κ_j = _kappa_SI_j · p_O₂  [m⁻¹]

    where _kappa_SI [m⁻¹ Pa⁻¹] are the Zheleznyak-fitted pressure-reduced
    absorption coefficients (set by init_photoionization) and p_O₂ is the
    O₂ partial pressure in Pa.  κ_j
    appears in the two-stream ODE:

        ∂_x Ψ_j⁺ = −κ_j Ψ_j⁺ + source
        ∂_x Ψ_j⁻ = +κ_j Ψ_j⁻ − source

    Parameters
    ----------
    p : float
        Gas pressure in bar.  κ scales linearly with p_O₂ = x_O₂ · p · 1e5.
    T : float
        Gas temperature in Kelvin (unused; retained for interface consistency).

    Returns
    -------
    numpy.ndarray, shape (N_γ,)
        Absorption coefficients [m⁻¹].
    """
    pO2_Pa = xO2 * p * 1e5
    return _kappa_SI * pO2_Pa


def get_gamma_Psi(EN, p=1.0, T=293.0):
    """
    Return the photon secondary emission coefficient vector γ_Ψ, shape (N_γ,).

    γ_Ψ[j] is the mean number of secondary electrons emitted per backward
    photon of group j incident on the cathode.

    Parameters
    ----------
    EN : float
        Reduced electric field in Townsend (unused).
    p : float
        Gas pressure in bar (unused).
    T : float
        Gas temperature in Kelvin (unused).

    Returns
    -------
    numpy.ndarray, shape (N_γ,)
    """
    return _XI_EMIT * np.ones(_N_GAMMA)


# ---------------------------------------------------------------------------
# __main__: plot R matrix entries and write transport data file
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Plot all non-zero off-diagonal entries of the reaction matrix R "
            "as a function of E/N on a log-log scale, and write a transport "
            "data file written only when --write_data is supplied."
        )
    )
    parser.add_argument(
        "--p", type=float, default=1.0, help="Gas pressure in bar (default: 1.0)."
    )
    parser.add_argument(
        "--T",
        type=float,
        default=293.0,
        help="Gas temperature in Kelvin (default: 293.0).",
    )
    parser.add_argument(
        "--write_data",
        type=str,
        default=None,
        metavar="FILE",
        help=(
            "Write tabulated transport data to FILE.  "
            "If omitted, no data file is written."
        ),
    )
    parser.add_argument(
        "--ngroups",
        type=int,
        default=3,
        help="Number of photon groups for two-stream photoionization (default: 3).",
    )
    parser.add_argument(
        "--cone-angle",
        dest="cone_angle",
        type=float,
        default=2.0,
        help=(
            "Half-opening angle of the photon emission cone in degrees (default: 45.0).  "
            "Translated to the spatial factor ΔΩ/(4π) = (1 − cos θ) / 2."
        ),
    )
    args = parser.parse_args()

    init_photoionization(args.ngroups, args.cone_angle)

    p_val = args.p
    T_val = args.T
    N_val = _make_N(p_val, T_val)

    # Human-readable labels for each non-zero off-diagonal R entry
    _RATE_LABELS = {
        (1, 0): r"$K_1$: e + N$_2$ $\to$ 2e + N$_2^+$",
        (2, 0): r"$K_2$: e + O$_2$ $\to$ 2e + O$_2^+$",
        (3, 0): r"$K_3$: e + O$_2$ $\to$ O$^-$ + O",
        (4, 0): r"$K_4$: e + O$_2$ $\to$ O$_2^-$(exc)",
        (0, 4): r"$K_5$: O$_2^-$(exc) $\to$ e + O$_2$",
        (5, 4): r"$K_6$: O$_2^-$(exc) + O$_2$ $\to$ O$_2^-$ + O$_2$",
        (0, 5): r"$K_7$: O$_2^-$ + O$_2$ $\to$ e + 2O$_2$",
        (5, 3): r"$K_8$: O$^-$ + O$_2$ $\to$ O + O$_2^-$",
        (6, 3): r"$K_9$: O$^-$ + O$_2$ $\to$ O$_3^{-*}$",
        (3, 6): r"$K_{10}$: O$_3^{-*}$ $\to$ O$^-$ + O$_2$",
        (7, 6): r"$K_{11}$: O$_3^{-*}$ + M $\to$ O$_3^-$ + M",
        (3, 7): r"$K_{12}$: O$_3^-$ + O$_2$ $\to$ O$^-$ + O$_2$",
    }

    numPts = 500
    EN_arr = np.logspace(1, 3, numPts)  # 10 – 1000 Td

    # Accumulate off-diagonal R entries across the E/N range
    offdiag_data = {}  # (i, j) -> ndarray
    for en_idx, en in enumerate(EN_arr):
        R = get_R(en, p_val, T_val)
        for i in range(_N_SPECIES):
            for j in range(_N_SPECIES):
                if i == j:
                    continue
                val = R[i, j]
                if val != 0.0:
                    key = (i, j)
                    if key not in offdiag_data:
                        offdiag_data[key] = np.zeros(numPts)
                    offdiag_data[key][en_idx] = val

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10))

    for (i, j), label in _RATE_LABELS.items():
        vals = offdiag_data.get((i, j))
        if vals is not None:
            nonzero = vals > 0
            if np.any(nonzero):
                ax1.loglog(EN_arr[nonzero], vals[nonzero], label=label)

    ax1.set_xlabel("E/N (Td)")
    ax1.set_ylabel("Rate (s$^{-1}$)")
    ax1.set_title(f"Reaction rates — air_2body.py,  p = {p_val} bar,  T = {T_val} K")
    ax1.legend(loc="best", fontsize=8, ncol=2)
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.set_ylim(bottom=1e-6)

    # --- Second subplot: α, η, and α−η ---
    alpha_arr = np.array([alpha(en, p_val, T_val) for en in EN_arr])
    eta_arr = np.array([eta(en, p_val, T_val) for en in EN_arr])
    net_arr = alpha_arr - eta_arr

    ax2.loglog(EN_arr, alpha_arr, label=r"$\alpha$")
    ax2.loglog(EN_arr, eta_arr, label=r"$\eta$")

    pos_mask = net_arr > 0
    neg_mask = net_arr < 0
    if np.any(pos_mask):
        ax2.loglog(EN_arr[pos_mask], net_arr[pos_mask], "k-", label=r"$\alpha - \eta$")
    if np.any(neg_mask):
        ax2.loglog(
            EN_arr[neg_mask],
            -net_arr[neg_mask],
            "k--",
            label=r"$\eta - \alpha$  ($\eta > \alpha$ region)",
        )

    ax2.set_xlabel("E/N (Td)")
    ax2.set_ylabel(r"Townsend coefficient (m$^{-1}$)")
    ax2.legend(loc="best", fontsize=8)
    ax2.grid(True, which="both", ls="--", alpha=0.4)

    plt.tight_layout()
    plt.show()

    if args.write_data is not None:
        # Evaluate rates on a denser grid (10-100 Td) for the data file
        numPts_file = 1500
        EN_file = np.logspace(1, 2, numPts_file)

        K1_f = np.array([k1(en) * xN2 * N_val for en in EN_file])
        K2_f = np.array([k2(en) * xO2 * N_val for en in EN_file])
        K3_f = np.array([k3(en) * xO2 * N_val for en in EN_file])
        K4_f = np.array([k4(en, T_val) * xO2 * N_val for en in EN_file])
        K5_f = np.array([k5(en) for en in EN_file])
        K6_f = np.array([k6(en) * xO2 * N_val for en in EN_file])
        K7_f = np.array([k7(en) * xO2 * N_val for en in EN_file])
        K8_f = np.array([k8(en) * xO2 * N_val for en in EN_file])
        K9_f = np.array([k9(en) * xO2 * N_val for en in EN_file])
        K10_f = np.array([k10(en) for en in EN_file])
        K11_f = np.array([k11(en) * N_val for en in EN_file])
        K12_f = np.array([k12(en, T_val) * xO2 * N_val for en in EN_file])
        header = (
            f"Transport data for N2/O2 air (80% N2, 20% O2) "
            f"at {p_val} bar, {T_val} K\n"
            "-------------------------------------------------------------------\n"
            "Column  1: E/N (Td)\n"
            "Column  2: K1  = k1*[N2]   (s^-1)  e + N2 -> 2e + N2+\n"
            "Column  3: K2  = k2*[O2]   (s^-1)  e + O2 -> 2e + O2+\n"
            "Column  4: K3  = k3*[O2]   (s^-1)  e + O2 -> O- + O\n"
            "Column  5: K4  = k4*[O2]   (s^-1)  e + O2 -> O2-(exc)\n"
            "Column  6: K5  = k5        (s^-1)  O2-(exc) -> e + O2\n"
            "Column  7: K6  = k6*[O2]   (s^-1)  O2-(exc) + O2 -> O2- + O2\n"
            "Column  8: K7  = k7*[O2]   (s^-1)  O2- + O2 -> e + 2 O2\n"
            "Column  9: K8  = k8*[O2]   (s^-1)  O- + O2 -> O + O2-\n"
            "Column 10: K9  = k9*[O2]   (s^-1)  O- + O2 -> O3-*\n"
            "Column 11: K10 = k10       (s^-1)  O3-* -> O- + O2\n"
            "Column 12: K11 = k11*N     (s^-1)  O3-* + M -> O3- + M*\n"
            "Column 13: K12 = k12*[O2]  (s^-1)  O3- + O2 -> O- + O2\n"
            "-------------------------------------------------------------------\n"
        )
        np.savetxt(
            args.write_data,
            np.column_stack(
                (
                    EN_file,
                    K1_f / p_val,
                    K2_f / p_val,
                    K3_f / p_val,
                    K4_f / p_val,
                    K5_f,
                    K6_f / p_val,
                    K7_f / p_val,
                    K8_f / p_val,
                    K9_f / p_val,
                    K10_f,
                    K11_f / p_val,
                    K12_f / p_val,
                )
            ),
            header=header,
        )
        print(f"Wrote {args.write_data}")
