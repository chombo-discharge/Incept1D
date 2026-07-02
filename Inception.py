"""
Inception.py — compute and plot the inception breakdown curve det Q(λ=0) = 0.

The inception threshold is determined by the determinant condition (manuscript,
eq. 333):

    det Q(λ) = 0   at  λ = 0

where Q = [Q_0; Q_d] is the boundary-condition matrix assembled from cathode
(Q_0) and anode (Q_d) constraints on the augmented state vector

    θ = (y, Ψ⁺, Ψ⁻)^T,   y = V n  (species fluxes),
                            Ψ⁺       (forward photon fluxes, +x direction),
                            Ψ⁻       (backward photon fluxes, −x direction).

The augmented ODE is (manuscript Eq. 252):

    ∂_x θ = A_aug θ,

where

    A_aug = [A    B    B ]      A = R V^{-1},  B = photon absorption coupling,
            [C   -D    0 ]
            [-C   0    D ]      C = photon source,  D = diag(κ_j).

When there are no photon species (N_γ = 0) A_aug reduces to A = R V^{-1}.

Cathode boundary conditions (Q_0 θ_0 = 0, manuscript Eq. 305–310):

    (Π_e + γ_+^T Π_+ − γ_Ψ^T Π_bck) θ_0 = 0    ← electron SEE
    Π_− θ_0 = 0                                   ← negative-ion zero flux
    Π_fwd θ_0 = 0                                 ← Ψ⁺(0) = 0 (no incoming fwd photons)

Anode boundary conditions (Q_d θ_0 = 0, manuscript Eq. 319–322):

    Π_+ M(d) θ_0 = 0                              ← positive-ion zero flux
    Π_bck M(d) θ_0 = 0                            ← Ψ⁻(d) = 0 (no incoming bck photons)

The full boundary matrix Q is (N_s + 2 N_γ) × (N_s + 2 N_γ).  A non-trivial
solution exists if and only if det Q = 0 (inception threshold).

The matrix exponential M(d) = exp(A_aug * d) is computed with a Padé
approximation (scipy.linalg.expm) after eigenvalue shifting:

    A_shift = A_aug − λ_max I
    M(d) = exp(A_shift * d) * exp(λ_max * d)

where λ_max = max(Re(eigenvalues(A_aug)), 0).  The shift keeps the Padé
argument well-conditioned even for large positive eigenvalues at high E/N.

For each value of p*d the critical E/N satisfying det Q = 0 is found with
the Brent root-finding method (scipy.optimize.brentq).  The corresponding
breakdown voltage is

    V* = E* * d = (E/N)* * (p*d) * 1e-16 / (kB * T)

where the last equality uses N = p * 1e5 / (kB*T) and the unit conversion
between Townsend (1e-21 V m^2) and bar (1e5 Pa).

Usage
-----
    python Inception.py <mechanism> [--p P [P ...]] [--d D [D ...]]
                           [--T T] [--xi-photo tag=val]
                           [--xi-emit tag=val]

Arguments
---------
mechanism
    Path to mechanism Python file.  Example: Data/Air/Air.py.
--p
    One or more gas pressures in bar (default: 1.0).  Fixed-p mode: for each
    pressure, gap length d = pd/p varies across the sweep.  A curve is
    produced for every (pressure, modifier) combination.
--d
    One or more gap lengths in mm.  Fixed-d mode: for each distance, pressure
    p = pd/d varies across the sweep.  A curve is produced for every
    (distance, modifier) combination.  May be combined with --p.
--T
    Gas temperature in Kelvin (default: 293.0).
CONFIG.json
    One or more JSON configuration files for the mechanism. Each file may
    contain a single configuration object or a list of objects under a
    "configurations" key.  All configurations across all files are run in
    sequence.  If omitted, a single baseline configuration with default
    parameters is used.
"""

import os
import json
import argparse
import functools
import importlib.util
import math
from typing import Optional

import numpy as np
import scipy.integrate
import scipy.linalg
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as _mticker

from FieldDistributions import (
    FieldDistribution,
    add_field_argument,
    parse_field_spec,
)

from Constants import kB as _kB, c_light as _C_LIGHT


class Mechanism:
    """Loaded mechanism module with configuration baked in.

    Wraps a bare mechanism module and applies all configuration parameters
    (xi_photo, xi_emit, reaction_multipliers, gamma overrides) directly inside
    the interface methods.  Callers never need a separate modifier object.

    Call resolve(positive) to obtain a polarity-resolved copy with the
    appropriate cathode SEE overrides applied for that polarity.
    """

    def __init__(
        self,
        _mod,
        label="",
        reaction_multipliers=None,
        xi_photo=1.0,
        xi_emit=1.0,
        gamma0=None,
        gamma1=None,
        eref=None,
        beta=None,
        pos_override=None,
        neg_override=None,
    ):
        self._mod = _mod
        self.label = label
        self._reaction_multipliers = reaction_multipliers or {}
        self._xi_photo = float(xi_photo)
        self._xi_emit = float(xi_emit)
        self._gamma0 = gamma0
        self._gamma1 = gamma1
        self._eref = eref
        self._beta = beta
        self._pos_override = pos_override  # raw dict or None
        self._neg_override = neg_override  # raw dict or None
        self.SPECIES = _mod.SPECIES
        self.ELECTRON_INDEX = _mod.ELECTRON_INDEX
        if hasattr(_mod, "alpha"):
            self.alpha = lambda EN, p, T: _mod.alpha(EN, p, T)
        if hasattr(_mod, "eta"):
            self.eta = lambda EN, p, T: _mod.eta(EN, p, T)

    def resolve(self, positive: bool) -> "Mechanism":
        """Return a polarity-resolved copy with cathode SEE overrides applied."""
        ovr = self._pos_override if positive else self._neg_override
        if ovr is None:
            return self
        return Mechanism(
            self._mod,
            self.label,
            reaction_multipliers=self._reaction_multipliers,
            xi_photo=ovr.get("xi_photo", self._xi_photo),
            xi_emit=ovr.get("xi_emit", self._xi_emit),
            gamma0=ovr.get("gamma0", self._gamma0),
            gamma1=ovr.get("gamma1", self._gamma1),
            eref=ovr.get("eref", self._eref),
            beta=ovr.get("beta", self._beta),
        )

    def get_R(self, EN, p, T):
        return self._mod.get_R(EN, p, T, multipliers=self._reaction_multipliers)

    def get_V(self, EN, p, T):
        return self._mod.get_V(EN, p, T)

    def get_B(self, EN, p, T):
        return self._mod.get_B(EN, p, T) * self._xi_photo

    def get_C(self, EN, p, T):
        return self._mod.get_C(EN, p, T)

    def get_kappa(self, p, T):
        return self._mod.get_kappa(p, T)

    def get_gamma_plus(self, EN, p, T):
        return self._mod.get_gamma_plus_with(
            EN,
            p,
            T,
            gamma0=self._gamma0,
            gamma1=self._gamma1,
            eref=self._eref,
            beta=self._beta,
        )

    def get_gamma_Psi(self, EN, p, T):
        return self._mod.get_gamma_Psi(EN, p, T) * self._xi_emit

    def get_Pi_e(self):
        return self._mod.get_Pi_e()

    def get_Pi_plus(self):
        return self._mod.get_Pi_plus()

    def get_Pi_minus(self):
        return self._mod.get_Pi_minus()

    def init_photoionization(self, ngroups=3, cone_angle_deg=2.0):
        if hasattr(self._mod, "init_photoionization"):
            self._mod.init_photoionization(ngroups, cone_angle_deg)


# Required attributes that every mechanism module must expose.
_REQUIRED_ATTRS = [
    "SPECIES",
    "ELECTRON_INDEX",
    "get_R",
    "get_V",
    "get_Pi_e",
    "get_Pi_plus",
    "get_Pi_minus",
    "get_gamma_plus",
    "get_gamma_plus_with",
    "get_B",
    "get_C",
    "get_kappa",
    "get_gamma_Psi",
]


def _read_json_configs(json_paths):
    """Read one or more JSON config files and return a flat list of raw dicts.

    Each file may contain: a dict with a "configurations" key (list), a bare
    list of dicts, or a single dict.  Each returned dict has '_cfg_dir' set to
    the directory of the JSON file so that Config.from_dict can resolve
    relative paths (e.g. cross_sections) correctly.
    """
    results = []
    for path in json_paths:
        abs_path = os.path.abspath(path)
        cfg_dir = os.path.dirname(abs_path)
        with open(abs_path) as fh:
            data = json.load(fh)
        if isinstance(data, dict) and "configurations" in data:
            items = data["configurations"]
        elif isinstance(data, list):
            items = data
        else:
            items = [data]
        for item in items:
            d = dict(item)
            d.setdefault("_cfg_dir", cfg_dir)
            results.append(d)
    return results


def load_mechanism(path, config_dict=None):
    """
    Load a reaction mechanism module from a Python source file.

    The module must expose the standard interface defined in _REQUIRED_ATTRS.

    Parameters
    ----------
    path : str
        Absolute or relative path to the mechanism .py file.
    config_dict : dict or None
        Raw configuration dict (e.g. parsed from JSON).  If given, a Config
        object is constructed from it by the companion Config.py that lives
        alongside the mechanism file.  Pass None for a baseline configuration.

    Returns
    -------
    Mechanism
    """
    abs_path = os.path.abspath(path)
    mech_dir = os.path.dirname(abs_path)
    if not os.path.isfile(abs_path):
        raise FileNotFoundError(f"Mechanism file not found: {abs_path}")

    # Discover Config.py in the same directory as the mechanism.
    cfg_py = os.path.join(mech_dir, "Config.py")
    if os.path.isfile(cfg_py):
        _cspec = importlib.util.spec_from_file_location("_mechconfig", cfg_py)
        _cm = importlib.util.module_from_spec(_cspec)
        _cspec.loader.exec_module(_cm)
        config = _cm.Config.from_dict(config_dict or {}, mech_dir)
    elif config_dict:
        raise ImportError(
            f"No Config.py found in {mech_dir} but a config_dict was supplied"
        )
    else:
        config = None

    spec = importlib.util.spec_from_file_location("mechanism", abs_path)
    mod = importlib.util.module_from_spec(spec)
    if config is not None:
        for attr, val in config.pre_exec_vars().items():
            mod.__dict__[attr] = val
    spec.loader.exec_module(mod)
    missing = [a for a in _REQUIRED_ATTRS if not hasattr(mod, a)]
    if missing:
        raise AttributeError(
            f"Mechanism '{path}' is missing required attributes: {missing}"
        )
    if config is not None:
        config.post_exec_init(mod)
    label = config.label if config is not None else ""
    params = config.mechanism_params() if config is not None else {}
    return Mechanism(mod, label, **params)


def _build_A_aug(EN, d, mod, p, T, lam=0.0):
    """
    Build the augmented ODE matrix A_aug for reduced field EN and gap length d.

    Returns (A_aug, N_gamma_eff, aug_mask, n_aug).
    aug_mask is None when N_gamma == 0.

    lam : float
        Temporal growth rate λ (s⁻¹). Modifies A → (R − λI)V⁻¹ and D → κ + λ/c.
        Default 0.0 (standard inception threshold).
    """
    _KD_LOCAL_THRESHOLD = 12.0
    n = len(mod.SPECIES)

    R = mod.get_R(EN, p, T)
    V = mod.get_V(EN, p, T)
    # V is always diagonal; multiply R and C by 1/v column-wise via broadcasting
    # instead of forming the full n×n inverse matrix.  If a mechanism ever returns
    # non-diagonal V this will give wrong results — add np.linalg.inv fallback then.
    v_inv = 1.0 / np.diag(V)
    A = R * v_inv[np.newaxis, :]
    if lam != 0.0:
        A = A - lam * np.diag(v_inv)  # (R − λI) V⁻¹

    B = mod.get_B(EN, p, T)  # (N_s, N_gamma)
    C = mod.get_C(EN, p, T) * v_inv[np.newaxis, :]  # (N_gamma, N_s)
    kappa = mod.get_kappa(p, T)  # (N_gamma,)
    if lam != 0.0:
        kappa = kappa + lam / _C_LIGHT  # κ + λ/c
    N_gamma = B.shape[1]

    if N_gamma > 0:
        kd = kappa * d
        local_mask = kd > _KD_LOCAL_THRESHOLD
        aug_mask = ~local_mask
        for j in np.where(local_mask)[0]:
            A += 2.0 * np.outer(B[:, j], C[j, :]) / kappa[j]
        B_aug = B[:, aug_mask]
        C_aug_ph = C[aug_mask, :]
        kappa_aug = kappa[aug_mask]
        N_gamma_eff = int(aug_mask.sum())
    else:
        N_gamma_eff = 0
        aug_mask = None
        B_aug = B
        C_aug_ph = C
        kappa_aug = kappa

    if N_gamma_eff == 0:
        A_aug = A
        n_aug = n
    else:
        D_mat = np.diag(kappa_aug)
        z_gg = np.zeros((N_gamma_eff, N_gamma_eff))
        A_aug = np.block(
            [
                [A, B_aug, B_aug],
                [C_aug_ph, -D_mat, z_gg],
                [-C_aug_ph, z_gg, D_mat],
            ]
        )
        n_aug = n + 2 * N_gamma_eff

    return A_aug, N_gamma_eff, aug_mask, n_aug


def _expm_shifted(M):
    """
    Compute expm(M) with eigenvalue shifting for numerical stability.

    Shifts by the true max real eigenvalue to keep expm well-conditioned without
    over-shrinking the Qd rows.  The restoration factor is intentionally NOT applied;
    it uniformly scales the Qd rows, leaving sign(det Q) unchanged.  Gershgorin
    bounds were tried but cause cumulative over-shrinkage of M across N_steps
    multiplications, making Q ill-conditioned and triggering the cond_Q > 1e14 guard.

    M is the pre-formed matrix argument to expm — for the midpoint rule this is
    A_aug * d_step; for Magnus expansions it is the full Ω already formed by the
    caller.
    """
    real_eigs = np.real(np.linalg.eigvals(M))
    lam_max = float(np.max(real_eigs))
    lam = lam_max if lam_max > 0.0 else 0.0
    return scipy.linalg.expm(M - lam * np.eye(M.shape[0]))


_DX_N_MIN_DEFAULT = 5
_DX_N_MAX_DEFAULT = 200
_DX_TOL_DEFAULT = 0.03


def parse_dx_spec(tokens, parser=None):
    """
    Parse a --dx token list into (N_min, N_max, tol).

    Tokens are all optional; missing values fall back to the defaults
    (N_min=5, N_max=200, tol=0.03).  N_min = N_max disables adaptation
    (max_depth = 0) and gives a constant uniform grid of N_min steps.

    Parameters
    ----------
    tokens : list of str or None
        Raw argparse token list, e.g. ['5', '200', '0.03'].
    parser : argparse.ArgumentParser or None
        If given, validation errors are routed through parser.error().
    """

    def _err(msg):
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    N_min, N_max, tol = _DX_N_MIN_DEFAULT, _DX_N_MAX_DEFAULT, _DX_TOL_DEFAULT
    if tokens:
        try:
            if len(tokens) >= 1:
                N_min = int(tokens[0])
            if len(tokens) >= 2:
                N_max = int(tokens[1])
            if len(tokens) >= 3:
                tol = float(tokens[2])
        except ValueError as exc:
            _err(f"--dx: {exc}")
    if N_min < 1:
        _err(f"--dx: N_min must be ≥ 1, got {N_min}")
    if N_max < N_min:
        _err(f"--dx: N_max ({N_max}) must be ≥ N_min ({N_min})")
    if not 0.0 < tol < 1.0:
        _err(f"--dx: tol must be in (0, 1), got {tol}")
    return N_min, N_max, tol


def midpoint_propagator(A_func, x_lo, x_hi):
    """
    Zeroth-order Magnus propagator (midpoint rule).

    Evaluates A at the interval midpoint and returns expm(A_mid * h).

    Parameters
    ----------
    A_func : callable
        A_func(x) → A_aug (ndarray).  x is a code-coordinate in metres; the
        polarity-dependent flip and field evaluation are encapsulated by the caller.
    x_lo, x_hi : float
        Subinterval bounds in code coordinates (metres).

    Returns
    -------
    ndarray
        Propagator matrix P = expm(A(x_mid) * h).
    """
    h = x_hi - x_lo
    x_mid = 0.5 * (x_lo + x_hi)
    return _expm_shifted(A_func(x_mid) * h)


def _adaptive_midpoint_segment(
    A_func, x_lo, x_hi, tol, max_depth, P_coarse=None, _level=0, _diag=None
):
    """
    Compute the midpoint propagator for [x_lo, x_hi] with adaptive step halving.

    Compares the 1-step (coarse) propagator against the 2-step (halved) estimate.
    If the relative Frobenius error exceeds tol and max_depth > 0, each half is
    refined recursively.  Passing P_coarse from the parent avoids one expm call
    when recursing (the parent's half-step result IS the child's coarse estimate).

    Parameters
    ----------
    A_func : callable
        Same signature as for midpoint_propagator.
    x_lo, x_hi : float
        Subinterval bounds in code coordinates (metres).
    tol : float
        Relative Frobenius error threshold; refinement stops when error <= tol.
    max_depth : int
        Maximum remaining recursion depth.  At depth 0 the coarse estimate is
        returned without any halving check.
    P_coarse : ndarray or None
        Pre-computed 1-step propagator for [x_lo, x_hi].  Computed internally
        when None (top-level call).
    _level : int
        Current recursion depth from the initial-segment top-level call (0 = top).
        Used only for diagnostics.
    _diag : list or None
        If a list is supplied, a dict is appended for every halving check performed:
        {'x_lo', 'x_hi', 'level', 'err', 'accepted'}.  Leaf calls at max_depth=0
        where no halving is attempted are NOT recorded.

    Returns
    -------
    ndarray
        Propagator matrix for [x_lo, x_hi].
    """
    if P_coarse is None:
        P_coarse = midpoint_propagator(A_func, x_lo, x_hi)
    if max_depth == 0:
        return P_coarse

    x_mid = 0.5 * (x_lo + x_hi)
    P_l = midpoint_propagator(A_func, x_lo, x_mid)
    P_r = midpoint_propagator(A_func, x_mid, x_hi)
    P_fine = P_r @ P_l

    norm_ref = np.linalg.norm(P_fine, "fro")
    if norm_ref < 1e-300:
        norm_ref = 1.0
    err = np.linalg.norm(P_fine - P_coarse, "fro") / norm_ref

    accepted = err <= tol
    if _diag is not None:
        _diag.append(
            {
                "x_lo": x_lo,
                "x_hi": x_hi,
                "level": _level,
                "err": err,
                "accepted": accepted,
            }
        )

    if accepted:
        return P_fine

    # Error too large — refine each half, passing the already-computed half-step
    # propagators as the coarse estimates for the sub-calls.
    P_l_fine = _adaptive_midpoint_segment(
        A_func, x_lo, x_mid, tol, max_depth - 1, P_l, _level + 1, _diag
    )
    P_r_fine = _adaptive_midpoint_segment(
        A_func, x_mid, x_hi, tol, max_depth - 1, P_r, _level + 1, _diag
    )
    return P_r_fine @ P_l_fine


def magnus2_propagator(A_func, x_lo, x_hi):
    """
    Second-order Magnus propagator with 2-point Gauss-Legendre quadrature.

    Ω = h/2·(A₁+A₂) + √3·h²/12·[A₂,A₁],   P = expm(Ω)

    The first term is the GL2 quadrature approximation to ∫A dx; the second is
    the leading commutator correction from the Magnus series.  The commutator
    vanishes for a uniform field (A constant across the interval), so Magnus2
    reduces exactly to the midpoint rule in that case.

    Convergence guard: the Magnus series converges only when ∫‖A‖ dx < π.
    If ‖Ω₁‖_F = ‖h/2·(A₁+A₂)‖_F exceeds π the interval is halved and the
    two sub-propagators are multiplied.  This keeps Magnus2 accurate at large
    p·d without requiring the caller to specify a finer grid.

    Parameters
    ----------
    A_func : callable
        A_func(x) → A_aug (ndarray).  x is a code-coordinate in metres.
    x_lo, x_hi : float
        Subinterval bounds in code coordinates (metres).

    Returns
    -------
    ndarray
        Propagator matrix P ≈ expm(∫_{x_lo}^{x_hi} A dx).
    """
    h = x_hi - x_lo
    mid = 0.5 * (x_lo + x_hi)
    offset = h / (2.0 * math.sqrt(3.0))
    A1 = A_func(mid - offset)
    A2 = A_func(mid + offset)
    Omega1 = 0.5 * h * (A1 + A2)
    # When ‖Ω₁‖_F > π the Magnus series is not guaranteed to converge; the
    # commutator correction would dominate and produce a wildly wrong exponent.
    # Fall back to the midpoint rule, which uses expm(A_mid·h) directly and
    # remains accurate for any step size via scipy's scaling-and-squaring.
    if np.linalg.norm(Omega1, "fro") > math.pi:
        return midpoint_propagator(A_func, x_lo, x_hi)
    Omega2 = (math.sqrt(3.0) * h**2 / 12.0) * (A2 @ A1 - A1 @ A2)
    return _expm_shifted(Omega1 + Omega2)


def _assemble_det_Q(M, EN_cathode, mod, p, T, N_gamma_eff, aug_mask, n_aug):
    """
    Assemble the boundary-condition matrix Q from the propagator M and cathode
    boundary field EN_cathode, then return det Q (NaN if ill-conditioned or overflow).

    Must be called inside a numpy errstate(over='ignore', invalid='ignore') block.
    mod must already be polarity-resolved (via mod.resolve(positive)).
    """
    n = len(mod.SPECIES)

    Pi_e = mod.get_Pi_e()
    Pi_plus = mod.get_Pi_plus()
    Pi_minus = mod.get_Pi_minus()
    gplus = mod.get_gamma_plus(EN_cathode, p, T)

    if N_gamma_eff == 0:
        Q0_e = Pi_e + gplus[np.newaxis, :] @ Pi_plus
        Q0_neg = Pi_minus
        Q0 = np.vstack([Q0_e, Q0_neg])
        Qd = Pi_plus @ M
    else:
        z_gamma = np.zeros((Pi_e.shape[0], 2 * N_gamma_eff))
        z_plus = np.zeros((Pi_plus.shape[0], 2 * N_gamma_eff))
        z_minus = np.zeros((Pi_minus.shape[0], 2 * N_gamma_eff))
        Pi_e_aug = np.hstack([Pi_e, z_gamma])
        Pi_plus_aug = np.hstack([Pi_plus, z_plus])
        Pi_minus_aug = np.hstack([Pi_minus, z_minus])

        Pi_fwd = np.zeros((N_gamma_eff, n_aug))
        Pi_bck = np.zeros((N_gamma_eff, n_aug))
        for j in range(N_gamma_eff):
            Pi_fwd[j, n + j] = 1.0
            Pi_bck[j, n + N_gamma_eff + j] = 1.0

        g_Psi = mod.get_gamma_Psi(EN_cathode, p, T)[aug_mask]

        Q0_e = (
            Pi_e_aug
            + gplus[np.newaxis, :] @ Pi_plus_aug
            - (g_Psi[np.newaxis, :] @ Pi_bck)
        )
        Q0_neg = Pi_minus_aug
        Q0_fwd = Pi_fwd
        Q0 = np.vstack([Q0_e, Q0_neg, Q0_fwd])

        Qd_plus = Pi_plus_aug @ M
        Qd_bck = Pi_bck @ M
        Qd = np.vstack([Qd_plus, Qd_bck])

    Q = np.vstack([Q0, Qd])

    if np.any(~np.isfinite(Q)):
        return np.nan

    row_norms = np.linalg.norm(Q, axis=1, keepdims=True)
    row_norms = np.where(row_norms < 1e-300, 1.0, row_norms)
    Q_norm = Q / row_norms
    sign, logabsdet = np.linalg.slogdet(Q_norm)
    cond_Q = np.linalg.cond(Q_norm)
    if cond_Q > 1e14:
        return np.nan
    if not np.isfinite(logabsdet):
        return np.nan
    if sign == 0:
        return 0.0
    return float(sign) * min(np.exp(logabsdet), 1e300)


def inception_det(
    EN_ref,
    pd,
    mod,
    p,
    T,
    field_dist,
    N_min=_DX_N_MIN_DEFAULT,
    N_max=_DX_N_MAX_DEFAULT,
    tol=_DX_TOL_DEFAULT,
    lam=0.0,
    positive_polarity=True,
    propagator=midpoint_propagator,
    diag_list=None,
):
    """
    Evaluate det Q(λ) for the inception boundary-value problem.

    Integrates the augmented ODE across the gap by accumulating per-subinterval
    propagator matrices, then assembles the boundary-condition matrix Q and returns
    its determinant.  Breakdown occurs when det Q = 0.

    For the midpoint propagator the grid is adaptive: each of the N_min initial
    segments is recursively halved until the relative Frobenius error between the
    one-step and two-step estimates falls below tol, or the maximum refinement
    depth (floor(log2(N_max // N_min))) is reached.  Setting N_min = N_max gives
    a constant uniform grid (max_depth = 0, no halving attempted).

    Parameters
    ----------
    EN_ref : float
        Uniform-equivalent reduced field in Townsend.
    pd : float
        Product of pressure and gap length in bar·m.
    mod : Mechanism
        Loaded mechanism (returned by load_mechanism).  Configuration parameters
        (xi_photo, xi_emit, reaction_multipliers, gamma overrides) are baked in.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.
    field_dist : FieldDistribution
        Field geometry specification (uniform / sphere-plane / sphere-sphere).
    N_min : int
        Minimum (initial) number of integration segments.  Default 5.
    N_max : int
        Maximum total number of fine steps across the gap (adaptive budget).
        N_max = N_min disables adaptation.  Default 200.
    tol : float
        Relative Frobenius error threshold for the midpoint step-halving check.
        Default 0.03 (3 %).
    lam : float
        Temporal growth rate λ (s⁻¹, default 0.0).
    positive_polarity : bool
        True → high-field electrode is anode (+).
        False → high-field electrode is cathode (−).
        Ignored for uniform and sphere-sphere (symmetric) geometries.
    propagator : callable
        Propagator algorithm.  Signature: propagator(A_func, x_lo, x_hi) → ndarray.
        Adaptive step halving is applied only when propagator is midpoint_propagator;
        other propagators (e.g. magnus2_propagator) use N_min uniform steps.

    Returns
    -------
    float
        det Q(λ).  Inception threshold: det Q = 0.
    """
    eff = mod.resolve(positive_polarity)
    d = pd / p
    f = field_dist.build(d)
    d_step = d / N_min
    max_depth = max(0, int(math.log2(max(1, N_max // N_min))))

    # A_func encapsulates the coordinate flip and all physics context.
    # Code coordinates run 0 (cathode) → d (anode); for positive polarity the
    # formula coordinate is reversed so xi=0 (sphere/anode) maps to x_code=d.
    def A_func(x_code, _f=f, _d=d, _EN=EN_ref, _pos=positive_polarity):
        x_formula = (_d - x_code) if _pos else x_code
        A_aug, *_ = _build_A_aug(_EN * _f(x_formula / _d), _d, eff, p, T, lam=lam)
        return A_aug

    # Extract augmented-system metadata once from the first step midpoint.
    # Photon structure (aug_mask, N_gamma_eff) is assumed constant across the gap.
    xi_mid0 = 0.5 * d_step / d
    xi_first = (1.0 - xi_mid0) if positive_polarity else xi_mid0
    _, N_gamma_eff, aug_mask, n_aug = _build_A_aug(
        EN_ref * f(xi_first), d, eff, p, T, lam=lam
    )
    M = np.eye(n_aug)

    with np.errstate(over="ignore", invalid="ignore"):
        if propagator is midpoint_propagator:
            # Adaptive midpoint: each of the N_min segments is refined by step
            # halving until the relative Frobenius error < tol or max_depth is
            # exhausted.  max_depth = 0 (when N_min = N_max) gives a constant
            # uniform grid with no halving attempted.
            for i in range(N_min):
                seg_diag = [] if diag_list is not None else None
                P = _adaptive_midpoint_segment(
                    A_func,
                    i * d_step,
                    (i + 1) * d_step,
                    tol,
                    max_depth,
                    _diag=seg_diag,
                )
                if diag_list is not None:
                    diag_list.append(
                        {
                            "segment": i,
                            "x_lo": i * d_step,
                            "x_hi": (i + 1) * d_step,
                            "halvings": seg_diag,
                        }
                    )
                M = P @ M
        else:
            for i in range(N_min):
                M = propagator(A_func, i * d_step, (i + 1) * d_step) @ M

        xi_cathode = 1.0 if positive_polarity else 0.0
        EN_cathode = EN_ref * f(xi_cathode)
        return _assemble_det_Q(M, EN_cathode, eff, p, T, N_gamma_eff, aug_mask, n_aug)


_ROOT_CHECK_TOL = 0.01  # |det Q(root)| / bracket-scale above this triggers a warning
_NAN_SENTINEL = 1e-290  # |value| below this → came from NaN → -1e-300 mapping


def _accept_root(root, EN_a, EN_b, fa, fb, pd, det_fn, mod, p, T):
    """
    Verify that det Q is genuinely near zero at a candidate root.

    Returns True if the root is accepted, False if it should be discarded.

    Two failure modes are distinguished:

    1. NaN-sentinel artefact: f_brentq mapped NaN → -1e-300 at a bracket
       endpoint, creating an artificial sign change.  det_fn(root) = NaN and
       at least one bracket fa/fb equals the sentinel (|value| ≤ 1e-290).
       → Discard and print a warning.

    2. Genuine singularity: both bracket endpoints are finite with opposite
       signs; brentq converges to the true det Q = 0 locus where Q is exactly
       singular so cond(Q) → ∞ and det_fn returns NaN.
       → Accept silently (NaN here is the expected consequence of Q → 0).

    3. Large finite residual: det_fn(root) is finite but |det Q| is large
       relative to the bracket scale.
       → Accept but print a warning.

    Parameters
    ----------
    root : float        — candidate E/N root returned by brentq
    EN_a, EN_b : float  — bracket endpoints used by brentq
    fa, fb : float      — f_brentq values at the bracket endpoints
                          (must be the FINAL bracket values, not the original
                           proxy-scan values — caller is responsible)
    pd : float          — pressure × gap length in bar·m
    det_fn : callable   — full-resolution det function (used for residual eval)
    mod, p, T           — mechanism, pressure, temperature passed to det_fn
    """
    det_root = det_fn(root, pd, mod, p, T)

    if not np.isfinite(det_root):
        # Distinguish genuine singularity from NaN-sentinel artefact.
        fa_sentinel = abs(fa) <= _NAN_SENTINEL
        fb_sentinel = abs(fb) <= _NAN_SENTINEL
        if fa_sentinel or fb_sentinel:
            print(
                f"  [det Q check] EN = {root:.4f} Td  pd = {pd*1e3:.4g} bar·mm: "
                f"det Q = NaN at root; bracket endpoint is NaN sentinel "
                f"(fa={fa:.3e}, fb={fb:.3e}) — artefact, discarded"
            )
            return False
        # Both bracket endpoints were genuine finite values; NaN at root means
        # Q is exactly singular (cond → ∞, det → 0).  Root is genuine.
        return True

    bracket_scale = max(
        abs(fa) if np.isfinite(fa) else 0.0, abs(fb) if np.isfinite(fb) else 0.0, 1e-300
    )
    if abs(det_root) > _ROOT_CHECK_TOL * bracket_scale:
        print(
            f"  [det Q check] EN = {root:.4f} Td  pd = {pd*1e3:.4g} bar·mm: "
            f"det Q = {det_root:.3e}  (bracket scale {bracket_scale:.3e}) — suspect root"
        )

    return True


def find_all_breakdown_EN(
    pd,
    mod,
    p,
    T,
    EN_lo=10.0,
    EN_hi=3e5,
    n_scan=200,
    first_only=False,
    det_fn=None,
    fast_det_fn=None,
    med_det_fn=None,
    EN_hints=None,
    hint_factor=2.0,
):
    """
    Find E/N values where det Q(E/N, pd) = 0.

    A logarithmic scan over [EN_lo, EN_hi] locates every sign change in
    det Q(EN); scipy.optimize.brentq refines each bracket independently.
    Using n_scan=200 reduces the chance of missing a sign change when multiple
    roots are close together.

    If EN_hints is supplied and first_only=True, a warm-start pass is
    attempted first: the lowest hint is bracketed by
    [hint/hint_factor, hint*hint_factor] and refined directly with brentq,
    bypassing the coarse scan entirely.  Falls back to the full scan if the
    hint fails to bracket.  Warm-start is intentionally disabled when
    first_only=False (--all-branches) because new branches can appear at any
    pd step; the coarse scan is required to detect them.

    Parameters
    ----------
    pd : float
        Product of pressure and gap length in bar·m.
    mod : Mechanism
        Loaded mechanism (returned by load_mechanism).
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.
    EN_lo : float
        Lower bound of the E/N search range in Townsend.  Default 10.0 Td.
    EN_hi : float
        Upper bound of the E/N search range in Townsend.  Default 3e5 Td.
    n_scan : int
        Number of points in the coarse scan.  Default 200.
    first_only : bool
        If True, return as soon as the first (lowest E/N) root is found.
    fast_det_fn : callable or None
        Optional cheaper det function used only for the coarse sign-change
        scan.  Brentq refinement always uses det_fn for full accuracy.
        If None, det_fn is used for both scan and refinement.
    med_det_fn : callable or None
        Optional medium-fidelity det function used for the fallback scan when
        fast_det_fn finds no sign changes.  Cheaper than det_fn for scanning
        but more accurate than fast_det_fn.  If None, det_fn is used as the
        fallback.
    EN_hints : list of float or None
        E/N values from the previous pd point to use as warm-start brackets.
        When all hints bracket successfully the full scan is skipped entirely.
        Pass [] or None to disable (first pd point, or after a failed pd).
    hint_factor : float
        Bracket half-width multiplier: each hint is bracketed by
        [hint/hint_factor, hint*hint_factor].  Default 2.0 (one octave each way).

    Returns
    -------
    list of float
        Critical E/N values in Townsend, sorted ascending.
        Empty list if no sign change is found in [EN_lo, EN_hi].
    """
    if det_fn is None:
        det_fn = inception_det
    scan_fn = fast_det_fn if fast_det_fn is not None else det_fn

    def f_brentq(EN):
        # NaN → small negative sentinel: NaN always occurs above the critical
        # E/N where the true sign is −1, so -1e-300 is physically correct.
        val = det_fn(EN, pd, mod, p, T)
        return val if np.isfinite(val) else -1e-300

    # Warm-start: only safe when first_only=True (single-branch mode).
    # With first_only=False (--all-branches) the coarse scan is mandatory
    # because new branches can appear at any pd point; skipping it would
    # silently miss roots not bracketed by existing hints.
    if EN_hints and first_only:
        warm_roots = []
        all_ok = True
        hints_to_try = EN_hints[:1] if first_only else EN_hints
        for en_hint in hints_to_try:
            EN_a = max(EN_lo, en_hint / hint_factor)
            EN_b = min(EN_hi, en_hint * hint_factor)
            fa = f_brentq(EN_a)
            fb = f_brentq(EN_b)
            if fa * fb < 0.0:
                # When the upper endpoint is a NaN sentinel, brentq would
                # converge to the NaN boundary rather than the true root.
                # Narrow EN_b geometrically (in log-EN space) until we land on
                # a finite negative value, giving brentq a clean bracket.
                if abs(fb) <= _NAN_SENTINEL:
                    lo_fn, hi_nan = EN_a, EN_b
                    found = False
                    for _ in range(20):
                        mid = math.sqrt(lo_fn * hi_nan)
                        v = det_fn(mid, pd, mod, p, T)
                        if not np.isfinite(v):
                            hi_nan = mid
                        elif v < 0.0:
                            EN_b, fb = mid, v
                            found = True
                            break
                        else:
                            lo_fn = mid
                    if not found:
                        all_ok = False
                        break
                root = scipy.optimize.brentq(
                    f_brentq, EN_a, EN_b, xtol=1e-8, rtol=1e-14
                )
                if _accept_root(root, EN_a, EN_b, fa, fb, pd, det_fn, mod, p, T):
                    warm_roots.append(float(root))
                else:
                    all_ok = False
                    break
            else:
                all_ok = False
                break
        if all_ok and warm_roots:
            return sorted(warm_roots)
        # Warm start incomplete — fall through to full scan.

    EN_scan = np.logspace(np.log10(EN_lo), np.log10(EN_hi), n_scan)

    def _scan_with(sfn):
        """Coarse sign-change scan using sfn; Brentq refinement always uses det_fn."""

        def f_s(EN):
            val = sfn(EN, pd, mod, p, T)
            return val if np.isfinite(val) else 0.0

        D = np.array([f_s(en) for en in EN_scan])
        local_roots = []
        for i in range(len(EN_scan) - 1):
            if D[i] * D[i + 1] < 0.0:
                EN_a, EN_b = EN_scan[i], EN_scan[i + 1]
                fa, fb = f_brentq(EN_a), f_brentq(EN_b)
                if fa * fb >= 0.0:
                    # Proxy scan bracket doesn't hold for the full det; widen by
                    # one scan step on each side and re-scan with the full det.
                    lo = EN_scan[max(0, i - 1)]
                    hi = EN_scan[min(len(EN_scan) - 1, i + 2)]
                    fine = np.logspace(np.log10(lo), np.log10(hi), 30)
                    fine_v = [f_brentq(en) for en in fine]
                    found = False
                    for k in range(len(fine) - 1):
                        if fine_v[k] * fine_v[k + 1] < 0.0:
                            EN_a, EN_b = fine[k], fine[k + 1]
                            fa, fb = fine_v[k], fine_v[k + 1]
                            found = True
                            break
                    if not found:
                        continue
                root = scipy.optimize.brentq(
                    f_brentq, EN_a, EN_b, xtol=1e-8, rtol=1e-14
                )
                if not _accept_root(root, EN_a, EN_b, fa, fb, pd, det_fn, mod, p, T):
                    continue
                local_roots.append(float(root))
                if first_only:
                    break
            elif D[i] > 0.0 and D[i + 1] == 0.0:
                # sfn returned NaN (mapped to 0) at the upper endpoint.  A sign
                # change may be hidden just before the NaN boundary; fine-scan
                # with the full det_fn to locate it.
                fine = np.logspace(np.log10(EN_scan[i]), np.log10(EN_scan[i + 1]), 40)
                fine_v = [det_fn(en, pd, mod, p, T) for en in fine]
                pos_j = None
                for j, v in enumerate(fine_v):
                    if np.isfinite(v) and v > 0.0:
                        pos_j = j
                    elif np.isfinite(v) and v < 0.0 and pos_j is not None:
                        root = scipy.optimize.brentq(
                            f_brentq, fine[pos_j], fine[j], xtol=1e-8, rtol=1e-14
                        )
                        if not _accept_root(
                            root,
                            fine[pos_j],
                            fine[j],
                            fine_v[pos_j],
                            v,
                            pd,
                            det_fn,
                            mod,
                            p,
                            T,
                        ):
                            break
                        local_roots.append(float(root))
                        break
                if first_only and local_roots:
                    break
        return local_roots

    roots = _scan_with(scan_fn)
    if not roots and fast_det_fn is not None:
        # Fast scan found no sign changes; fall back to the full det_fn scan.
        # med_det_fn is intentionally not used here: it almost never brackets
        # roots that fast_det missed, so using it only adds scan overhead before
        # the full fallback (which is always needed anyway).
        roots = _scan_with(det_fn)
    return roots


def compute_paschen_curve(
    pd_arr,
    mod,
    p,
    T,
    all_branches=False,
    det_fn=None,
    fast_det_fn=None,
    med_det_fn=None,
):
    """
    Compute Paschen branches across a p*d sweep.

    For each pd point roots of det Q(E/N, pd) = 0 are found via
    find_all_breakdown_EN.  Branches are tracked by continuity: each new root
    is matched to the existing branch whose most-recent E/N is nearest in
    log space, using a greedy nearest-neighbour assignment.  Unmatched roots
    start new branches.  This correctly handles saddle-node bifurcations where
    two new low-E/N roots appear without disrupting the pre-existing branch.

    Parameters
    ----------
    pd_arr : numpy.ndarray, shape (M,)
        Array of p*d values in bar·m.
    mod : Mechanism
        Loaded mechanism (returned by load_mechanism).
    p : float or numpy.ndarray, shape (M,)
        Gas pressure in bar.  A scalar is used for all points (fixed-p mode);
        an array of length M allows pressure to vary per point (fixed-d mode).
    T : float
        Gas temperature in Kelvin.
    all_branches : bool
        If False (default), only the lowest-E/N root (branch 0) is kept at
        each pd point.  If True, all roots are collected.

    Returns
    -------
    list of dict
        One dict per branch.  Branch 0 is the one whose first root appeared
        earliest (lowest pd); within that, ordered by first-appearance E/N.
        Each dict has keys:

        - ``'idx'`` — integer indices into pd_arr where this branch has a root
        - ``'pd'``, ``'EN'``, ``'V'``, ``'p'``, ``'d'`` — corresponding arrays

        Returns an empty list if no solutions are found anywhere.
    """
    if det_fn is None:
        det_fn = inception_det
    p_arr = np.full_like(pd_arr, p) if np.ndim(p) == 0 else np.asarray(p, dtype=float)
    branches = []
    branch_last_logEN = []  # last log10(EN) for each branch, for continuity tracking
    prev_roots_EN = []  # roots found at the previous pd point (warm-start hints)

    for i, (pd_i, p_i) in enumerate(zip(pd_arr, p_arr)):
        roots = find_all_breakdown_EN(
            pd_i,
            mod,
            p_i,
            T,
            first_only=not all_branches,
            det_fn=det_fn,
            fast_det_fn=fast_det_fn,
            med_det_fn=med_det_fn,
            EN_hints=prev_roots_EN,
        )
        prev_roots_EN = list(roots)
        if not roots:
            continue
        d_i = pd_i / p_i

        if not branches:
            # First pd with roots: create one branch per root, sorted ascending
            for EN in sorted(roots):
                V = EN * pd_i * 1e-21 / (_kB * T) * 1e5
                branches.append(
                    {
                        "idx": [i],
                        "pd": [pd_i],
                        "EN": [EN],
                        "V": [V],
                        "p": [p_i],
                        "d": [d_i],
                    }
                )
                branch_last_logEN.append(np.log10(EN))
        else:
            # Match each new root to the nearest existing branch in log-EN space.
            # Build cost matrix: rows = existing branches, cols = new roots.
            log_roots = np.log10(np.array(sorted(roots)))
            log_last = np.array(branch_last_logEN)
            cost = np.abs(log_last[:, None] - log_roots[None, :])

            assigned_branches = set()
            assigned_roots = set()
            root_to_branch = {}

            # Greedy: repeatedly pick the (branch, root) pair with smallest distance
            flat_order = np.argsort(cost, axis=None)
            for idx in flat_order:
                b, r = divmod(int(idx), len(log_roots))
                if b in assigned_branches or r in assigned_roots:
                    continue
                root_to_branch[r] = b
                assigned_branches.add(b)
                assigned_roots.add(r)
                if len(assigned_roots) == min(len(branches), len(log_roots)):
                    break

            sorted_roots = sorted(roots)
            for r, EN in enumerate(sorted_roots):
                V = EN * pd_i * 1e-21 / (_kB * T) * 1e5
                if r in root_to_branch:
                    b = root_to_branch[r]
                    branches[b]["idx"].append(i)
                    branches[b]["pd"].append(pd_i)
                    branches[b]["EN"].append(EN)
                    branches[b]["V"].append(V)
                    branches[b]["p"].append(p_i)
                    branches[b]["d"].append(d_i)
                    branch_last_logEN[b] = np.log10(EN)
                else:
                    branches.append(
                        {
                            "idx": [i],
                            "pd": [pd_i],
                            "EN": [EN],
                            "V": [V],
                            "p": [p_i],
                            "d": [d_i],
                        }
                    )
                    branch_last_logEN.append(np.log10(EN))

    for br in branches:
        for k in ("pd", "EN", "V", "p", "d"):
            br[k] = np.array(br[k], dtype=float)
        br["idx"] = np.array(br["idx"], dtype=int)

    return branches


def main():
    """
    Parse command-line arguments, solve for breakdown E/N across a p*d sweep,
    and display a two-panel Paschen curve figure.

    The p*d sweep range and resolution are controlled by --pd-min, --pd-max,
    and --pd-num (defaults: 1e-2 to 1e3 bar·mm, 50 points).

    Panel 1 — V* vs p*d (log-log):
        Shows the breakdown voltage.  The classical Paschen minimum appears
        as the lowest point on each curve.

    Panel 2 — E/N* vs p*d (log-log):
        Shows the critical reduced electric field at breakdown for each curve.

    One curve is produced for every (pressure, configuration) combination.

    A summary table is printed to stdout for each (p, tag, pd, E/N*, V*) triple
    where a solution was found.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Compute and plot the Paschen breakdown curve det Q(E/N, pd) = 0.  "
            "Sweeps pd from 1e-4 to 1e4 bar·mm in two modes: fixed pressure "
            "(--p) or fixed gap distance (--d).  Both may be combined."
        )
    )
    parser.add_argument(
        "mechanism",
        help="Path to mechanism Python file (e.g. Data/Air/Air.py).",
    )
    parser.add_argument(
        "configs",
        nargs="*",
        metavar="CONFIG.json",
        help=(
            "One or more JSON configuration files for the mechanism.  "
            "Each file may contain a single configuration object or a list "
            'of objects under a "configurations" key.  All configurations '
            "across all files are run in sequence.  If omitted, a single "
            "baseline configuration with default parameters is used."
        ),
    )
    parser.add_argument(
        "--p",
        type=float,
        nargs="+",
        default=[],
        metavar="P",
        help="Fixed-p mode: pressure(s) in bar (default: 1.0 when --d is not given).",
    )
    parser.add_argument(
        "--d",
        type=float,
        nargs="+",
        default=[],
        metavar="D",
        help="Fixed-d mode: gap distance(s) in mm.",
    )
    parser.add_argument(
        "--pd-min",
        type=float,
        default=1e-2,
        metavar="PD_MIN",
        help="Minimum p*d in bar·mm (default: 1e-2).",
    )
    parser.add_argument(
        "--pd-max",
        type=float,
        default=1e3,
        metavar="PD_MAX",
        help="Maximum p*d in bar·mm (default: 1e3).",
    )
    parser.add_argument(
        "--pd-num",
        type=int,
        default=50,
        metavar="PD_NUM",
        help="Number of logarithmically spaced p*d grid points (default: 50).",
    )
    parser.add_argument(
        "--T",
        type=float,
        default=293.0,
        help="Gas temperature in Kelvin (default: 293.0).",
    )
    parser.add_argument(
        "--write-to-file",
        type=str,
        default=None,
        metavar="FILE",
        help="Write all curves to a tab-separated file.",
    )
    parser.add_argument(
        "--save-subplots",
        action="store_true",
        default=False,
        help="Save each subplot as a separate PDF (no titles).",
    )
    parser.add_argument(
        "--all-branches",
        action="store_true",
        default=False,
        help=(
            "Compute and plot all branches of the breakdown curve.  "
            "By default only branch 1 (lowest E/N root) is computed."
        ),
    )
    add_field_argument(parser)
    parser.add_argument(
        "--dx",
        nargs="*",
        default=None,
        metavar="SPEC",
        help=(
            "Adaptive integration stepping: N_min [N_max [tol]].  "
            "N_min (default 5): minimum number of integration segments.  "
            "N_max (default 200): maximum total fine steps; N_max = N_min gives "
            "a constant uniform grid with no adaptive refinement.  "
            "tol (default 0.03): relative Frobenius error threshold for the "
            "midpoint step-halving check (fraction, not percent).  "
            "Example: --dx 10 400 0.01  sets N_min=10, N_max=400, tol=1%%."
        ),
    )
    parser.add_argument(
        "--lam",
        type=float,
        default=0.0,
        metavar="LAM",
        help=(
            "Temporal growth rate λ in s⁻¹ for the generalised inception criterion "
            "det Q(λ) = 0 (default: 0.0 = standard inception threshold).  "
            "λ > 0 → growing discharge (lower breakdown voltage); "
            "λ < 0 → decaying discharge (higher breakdown voltage)."
        ),
    )
    parser.add_argument(
        "--plot-separate-branches",
        action="store_true",
        default=False,
        help=(
            "Give each branch its own line style (solid / dashed / dotted / "
            "dash-dot) and legend entry ('branch 1', 'branch 2', …).  "
            "By default all branches share the same style and a single legend entry."
        ),
    )
    parser.add_argument(
        "--plot-ionization-integral",
        action="store_true",
        default=False,
        help=(
            "Overlay the ionization integral ∫max(α−η,0)dx on a second y-axis "
            "in the voltage subplot.  Only segments where α > η contribute.  "
            "Requires the mechanism to expose both alpha and eta."
        ),
    )
    parser.add_argument(
        "--streamer-criterion",
        dest="streamer_criterion",
        type=float,
        default=None,
        metavar="C",
        help=(
            "Solve for the streamer criterion ∫max(α−η,0)dx = C.  "
            "C must be > 0.  Plots the streamer curve on both panels and "
            "adds it to --write-to-file output.  "
            "Requires the mechanism to expose alpha and eta."
        ),
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        default=False,
        help="Skip the matplotlib figure entirely (useful for batch/scripted runs).",
    )
    parser.add_argument(
        "--method",
        choices=["midpoint", "magnus2"],
        default="midpoint",
        help=(
            "Propagator algorithm for the path-ordered matrix exponential.  "
            "'midpoint' (default): zeroth-order Magnus / midpoint rule.  "
            "'magnus2': second-order Magnus with 2-point Gauss-Legendre quadrature; "
            "reduces to midpoint for uniform fields."
        ),
    )
    args = parser.parse_args()

    if args.streamer_criterion is not None and args.streamer_criterion <= 0:
        parser.error("--streamer-criterion value must be > 0")

    if not args.p and not args.d:
        args.p = [1.0]

    _field_dist = parse_field_spec(args.field, parser)
    _N_min, _N_max, _dx_tol = parse_dx_spec(args.dx, parser)

    _propagator = (
        midpoint_propagator if args.method == "midpoint" else magnus2_propagator
    )

    mech_dir = os.path.dirname(os.path.abspath(args.mechanism))
    mech_name = os.path.basename(args.mechanism)

    # Read JSON config files as raw dicts; no gas-specific class needed here.
    raw_dicts = _read_json_configs(args.configs) if args.configs else [{}]

    # Determine n (number of species) from the first config.
    _mod0 = load_mechanism(args.mechanism, raw_dicts[0])
    n = len(_mod0.SPECIES)
    del _mod0

    if args.streamer_criterion is not None:
        _mod_check = load_mechanism(args.mechanism, raw_dicts[0])
        if not (hasattr(_mod_check, "alpha") and hasattr(_mod_check, "eta")):
            parser.error(
                "--streamer-criterion requires the mechanism to expose alpha and eta"
            )
        del _mod_check

    # Fixed pd sweep
    pd_arr = np.logspace(
        np.log10(args.pd_min * 1e-3), np.log10(args.pd_max * 1e-3), args.pd_num
    )

    # Build curve specs: (label, p_arr, d_arr, mod)
    # Each config dict gets its own freshly loaded Mechanism instance.
    curve_specs = []
    for cfg_dict in raw_dicts:
        mod = load_mechanism(args.mechanism, cfg_dict)
        cfg_label = mod.label
        if len(mod.SPECIES) != n:
            raise ValueError(
                f"Config '{cfg_label}' has {len(mod.SPECIES)} species; "
                f"expected {n} (from first config)."
            )
        for p in args.p:
            curve_specs.append(
                (
                    f"{cfg_label}, p={p} bar",
                    np.full_like(pd_arr, p),
                    pd_arr / p,
                    mod,
                )
            )
        for d_mm in args.d:
            d_m = d_mm * 1e-3
            curve_specs.append(
                (
                    f"{cfg_label}, d={d_mm} mm",
                    pd_arr / d_m,
                    np.full_like(pd_arr, d_m),
                    mod,
                )
            )

    # Compute alpha=eta crossover E/N from the first config's module (for annotation).
    _en_cross = {}
    if curve_specs:
        _ref_mod = curve_specs[0][3]
        if hasattr(_ref_mod, "alpha") and hasattr(_ref_mod, "eta"):
            EN_scan = np.logspace(np.log10(10.0), np.log10(1e5), 200)
            for p in args.p:
                f_cross = lambda EN, _p=p: _ref_mod.alpha(
                    EN, _p, args.T
                ) - _ref_mod.eta(EN, _p, args.T)
                f_vals = np.array([f_cross(en) for en in EN_scan])
                sign_changes = np.where(f_vals[:-1] * f_vals[1:] < 0)[0]
                if len(sign_changes):
                    i = sign_changes[0]
                    _en_cross[p] = scipy.optimize.brentq(
                        f_cross, EN_scan[i], EN_scan[i + 1], xtol=1e-3, rtol=1e-5
                    )
                else:
                    _en_cross[p] = None

    _first_mod = curve_specs[0][3] if curve_specs else None
    _plot_aed = (
        args.plot_ionization_integral
        and _first_mod is not None
        and hasattr(_first_mod, "alpha")
        and hasattr(_first_mod, "eta")
    )

    def _aed_integral(EN_ref, p_val, d_val, _mod):
        """Compute ∫max(α−η,0)dx respecting the actual field profile."""
        if _field_dist.field_type == "uniform":
            return (
                max(
                    0.0,
                    _mod.alpha(EN_ref, p_val, args.T) - _mod.eta(EN_ref, p_val, args.T),
                )
                * d_val
            )
        f = _field_dist.build(d_val)
        xis = (np.arange(_N_min) + 0.5) / _N_min
        EN_arr = EN_ref * np.array([f(xi) for xi in xis])
        ds = d_val / _N_min
        diff = np.array(
            [
                _mod.alpha(en, p_val, args.T) - _mod.eta(en, p_val, args.T)
                for en in EN_arr
            ]
        )
        return float(np.sum(np.maximum(0.0, diff)) * ds)

    def _polarity_desc(polarity):
        if _field_dist.field_type == "uniform":
            return polarity  # "positive" / "negative"
        return f"sphere={polarity}"  # "sphere=positive" / "sphere=negative"

    def _branch0_grids(branches):
        V_g = np.full_like(pd_arr, np.nan)
        EN_g = np.full_like(pd_arr, np.nan)
        p_g = np.full_like(pd_arr, np.nan)
        d_g = np.full_like(pd_arr, np.nan)
        if branches:
            br0 = branches[0]
            V_g[br0["idx"]] = br0["V"]
            EN_g[br0["idx"]] = br0["EN"]
            p_g[br0["idx"]] = br0["p"]
            d_g[br0["idx"]] = br0["d"]
        return V_g, EN_g, p_g, d_g

    # Collects (label, p_arr, d_arr, EN_star, V_star) for every solved curve.
    _file_records = []
    _streamer_records = []  # (label, p_arr, d_arr, EN_arr, V_arr) per streamer setting
    _all_V_pos = {}  # label -> full 200-point V array, positive polarity
    _all_V_neg = {}  # label -> full 200-point V array, negative polarity
    _curve_color = (
        {}
    )  # label -> matplotlib line colour (for consistent colouring on ax3)

    if not args.no_plot:
        plt.rcParams["font.size"] += 2

        # Set up figure — panels: V*, E/N*, optionally modifier ratio, optionally field-error %
        _show_ratio = len(raw_dicts) > 1
        _n_panels = 2 + int(_show_ratio)
        fig, _axes = plt.subplots(1, _n_panels, figsize=(7 * _n_panels, 6))
        ax1, ax2 = _axes[0], _axes[1]
        _next_ax = 2
        if _show_ratio:
            ax3 = _axes[_next_ax]
            _next_ax += 1
        else:
            ax3 = None
        _field_str = _field_dist.label
        _method_str = "" if args.method == "midpoint" else f",  {args.method}"
        _suptitle = fig.suptitle(
            f"Paschen Curve  —  {mech_name},  T = {args.T} K,  {_field_str}{_method_str}",
            fontsize=13,
        )

        _AED_ALPHA = 0.4
        ax1b = ax1.twinx() if _plot_aed else None
        if ax1b is not None:
            ax1b.set_yscale("symlog", linthresh=1)
            ax1b.yaxis.set_major_locator(
                _mticker.SymmetricalLogLocator(
                    linthresh=1, base=10, subs=[1.0, 2.0, 5.0]
                )
            )
            ax1b.yaxis.set_major_formatter(
                _mticker.LogFormatter(minor_thresholds=(np.inf, np.inf))
            )
            ax1b.set_ylabel(
                r"$\int_0^d \max(\alpha-\eta,\,0)\,\mathrm{d}x$", alpha=_AED_ALPHA
            )
            ax1b.spines["right"].set_alpha(_AED_ALPHA)
            ax1b.tick_params(axis="y", colors=(0, 0, 0, _AED_ALPHA))

            def _aed_format_coord(x, y, _a1=ax1, _a1b=ax1b):
                _, y1 = _a1.transData.inverted().transform(
                    _a1b.transData.transform((x, y))
                )
                return f"pd = {x:.3g} bar·mm    U = {y1:.4g} kV    (α−η)d = {y:.4g}"

            ax1b.format_coord = _aed_format_coord

        _markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
        _BRANCH_LS = ["-", "--", ":", "-."]
    else:
        ax1 = ax2 = ax3 = ax1b = None
        _markers = _BRANCH_LS = None

    for curve_idx, (label, p_arr, _d_arr, mod) in enumerate(curve_specs):

        if _field_dist.field_type == "sphere-sphere":
            max_dR = np.max(_d_arr) / _field_dist.sphere_R
            if max_dR > 4.0:
                print(
                    f"  Warning [{label}]: max d/R = {max_dR:.3g} > 4 — "
                    f"sphere-sphere field approximation may be inaccurate."
                )

        # Determine det functions for both polarities.
        # fast_det_fn  — uniform field, always N=1 (coarse sign-change scan)
        # med_det_fn   — N=min(5,N), modest resolution (fallback scan)
        # det_fn       — full field_dist accuracy (Brentq)
        if _field_dist.field_type == "uniform":
            det_pos = det_neg = functools.partial(
                inception_det,
                field_dist=_field_dist,
                N_min=_N_min,
                N_max=_N_max,
                tol=_dx_tol,
                lam=args.lam,
                positive_polarity=True,
                propagator=_propagator,
            )
            fast_det_pos = fast_det_neg = None
            med_det_pos = med_det_neg = None
        else:
            _fast_fd = FieldDistribution("uniform")
            _med_fd = FieldDistribution(_field_dist.field_type, _field_dist.sphere_R)
            _N_med = min(
                5, _N_min
            )  # cheap constant scan; N_med = N_med disables adaptation
            det_pos = functools.partial(
                inception_det,
                field_dist=_field_dist,
                N_min=_N_min,
                N_max=_N_max,
                tol=_dx_tol,
                lam=args.lam,
                positive_polarity=True,
                propagator=_propagator,
            )
            det_neg = (
                det_pos
                if _field_dist.is_symmetric
                else functools.partial(
                    inception_det,
                    field_dist=_field_dist,
                    N_min=_N_min,
                    N_max=_N_max,
                    tol=_dx_tol,
                    lam=args.lam,
                    positive_polarity=False,
                    propagator=_propagator,
                )
            )
            fast_det_pos = fast_det_neg = functools.partial(
                inception_det,
                field_dist=_fast_fd,
                N_min=1,
                N_max=1,
                tol=_dx_tol,
                lam=args.lam,
                propagator=_propagator,
            )
            med_det_pos = functools.partial(
                inception_det,
                field_dist=_med_fd,
                N_min=_N_med,
                N_max=_N_med,
                tol=_dx_tol,
                lam=args.lam,
                positive_polarity=True,
                propagator=_propagator,
            )
            med_det_neg = (
                med_det_pos
                if _field_dist.is_symmetric
                else functools.partial(
                    inception_det,
                    field_dist=_med_fd,
                    N_min=_N_med,
                    N_max=_N_med,
                    tol=_dx_tol,
                    lam=args.lam,
                    positive_polarity=False,
                    propagator=_propagator,
                )
            )

        print(f"\nSolving Paschen curve: {label} ({_polarity_desc('positive')})")
        branches_pos = compute_paschen_curve(
            pd_arr,
            mod,
            p_arr,
            args.T,
            all_branches=args.all_branches,
            det_fn=det_pos,
            fast_det_fn=fast_det_pos,
            med_det_fn=med_det_pos,
        )
        if _field_dist.is_symmetric:
            branches_neg = branches_pos  # symmetric; reuse same object
        else:
            print(f"\nSolving Paschen curve: {label} ({_polarity_desc('negative')})")
            branches_neg = compute_paschen_curve(
                pd_arr,
                mod,
                p_arr,
                args.T,
                all_branches=args.all_branches,
                det_fn=det_neg,
                fast_det_fn=fast_det_neg,
                med_det_fn=med_det_neg,
            )

        if not branches_pos and not branches_neg:
            print(f"  [{label}] No breakdown found for any pd value.")
            _all_V_pos[label] = np.full_like(pd_arr, np.nan)
            _all_V_neg[label] = np.full_like(pd_arr, np.nan)
            continue

        polarity_pairs = [
            ("positive", branches_pos, "-"),
            ("negative", branches_neg, "--"),
        ]
        for polarity, branches, ls_pol in polarity_pairs:
            if not branches:
                continue
            pol_label = f"{label} ({_polarity_desc(polarity)})"
            for b_idx, br in enumerate(branches):
                order = np.argsort(br["pd"])
                br_pd = br["pd"][order]
                br_EN = br["EN"][order]
                br_V = br["V"][order]
                br_p = br["p"][order]
                br_d = br["d"][order]

                if not args.no_plot:
                    markevery = [0, len(br_pd) - 1]
                    if args.plot_separate_branches:
                        ls = _BRANCH_LS[b_idx % len(_BRANCH_LS)]
                        b_label = f"{pol_label} (branch {b_idx + 1})"
                    else:
                        ls = ls_pol
                        b_label = pol_label if b_idx == 0 else "_nolegend_"
                    is_first = polarity == "positive" and b_idx == 0
                    color_kw = {} if is_first else {"color": _curve_color[label]}
                    plot_kw = dict(
                        label=b_label,
                        marker=_markers[curve_idx % len(_markers)],
                        markevery=markevery,
                        markersize=6,
                        ls=ls,
                        **color_kw,
                    )
                    ax1.loglog(br_pd * 1e3, br_V / 1000, **plot_kw)
                    if is_first:
                        _curve_color[label] = ax1.get_lines()[-1].get_color()
                    col = _curve_color[label]
                    if ax1b is not None:
                        aed = np.array(
                            [
                                (
                                    _aed_integral(en, p, d, mod)
                                    if np.isfinite(en)
                                    else np.nan
                                )
                                for en, p, d in zip(br_EN, br_p, br_d)
                            ]
                        )
                        ax1b.plot(
                            br_pd * 1e3,
                            aed,
                            color=col,
                            alpha=_AED_ALPHA,
                            ls=ls_pol,
                            lw=1.5,
                            zorder=1,
                        )
                    ax2.loglog(
                        br_pd * 1e3,
                        br_EN,
                        **{**plot_kw, "color": col, "label": b_label},
                    )

                # Print summary table for this branch
                print(
                    f"\n  {'pd (bar·mm)':>14}  {'p (bar)':>10}  {'d (mm)':>8}  "
                    f"{'E/N (Td)':>12}  {'U (kV)':>12}  {'E (V/m)':>14}  "
                    f"[{pol_label} (branch {b_idx + 1})]"
                )
                print("  " + "-" * 90)
                step = max(1, len(br_pd) // 20)
                for j in range(0, len(br_pd), step):
                    print(
                        f"  {br_pd[j]*1e3:>14.4e}  "
                        f"{br_p[j]:>10.4g}  "
                        f"{br_d[j]*1e3:>8.4g}  "
                        f"{br_EN[j]:>12.4f}  "
                        f"{br_V[j]/1000:>12.4f}  "
                        f"{br_V[j]/br_d[j]:>14.4e}"
                    )

        # Map branch 0 onto the full pd grid for file output and ratio plot
        V_pos, EN_pos, p_grid, d_grid = _branch0_grids(branches_pos)
        V_neg, EN_neg, _, _ = _branch0_grids(branches_neg)

        _all_V_pos[label] = V_pos
        _all_V_neg[label] = V_neg

        _file_records.append(
            (
                f"{label} ({_polarity_desc('positive')})",
                p_grid,
                d_grid,
                EN_pos,
                V_pos,
                mod,
            )
        )
        _file_records.append(
            (
                f"{label} ({_polarity_desc('negative')})",
                p_grid,
                d_grid,
                EN_neg,
                V_neg,
                mod,
            )
        )

    if args.streamer_criterion is not None:
        _C = args.streamer_criterion
        _EN_sc = np.logspace(1.0, np.log10(3e5), 200)

        for _cs_label, _p_arr_cs, _d_arr_cs, _mod_cs in curve_specs:
            _s_label = f"Streamer (C={_C}), {_cs_label}"
            _EN_s = np.full_like(pd_arr, np.nan)
            _V_s = np.full_like(pd_arr, np.nan)

            print(f"\nSolving streamer criterion: {_s_label}")
            for _i, (_pd_i, _p_i, _d_i) in enumerate(zip(pd_arr, _p_arr_cs, _d_arr_cs)):
                _fvals = np.array(
                    [_aed_integral(_en, _p_i, _d_i, _mod_cs) - _C for _en in _EN_sc]
                )
                _idx = np.where(_fvals[:-1] * _fvals[1:] < 0)[0]
                if _idx.size == 0:
                    continue
                _k = _idx[0]
                try:
                    _root = scipy.optimize.brentq(
                        lambda _en, __p=_p_i, __d=_d_i, __m=_mod_cs: _aed_integral(
                            _en, __p, __d, __m
                        )
                        - _C,
                        _EN_sc[_k],
                        _EN_sc[_k + 1],
                        xtol=1e-6,
                        rtol=1e-10,
                    )
                    _EN_s[_i] = _root
                    _V_s[_i] = _root * _pd_i * 1e-21 / (_kB * args.T) * 1e5
                except ValueError:
                    pass

            _streamer_records.append((_s_label, _p_arr_cs, _d_arr_cs, _EN_s, _V_s))

            _smask = np.isfinite(_EN_s)
            if not np.any(_smask):
                print(f"  [{_s_label}] No solution found for any pd value.")
                continue

            if not args.no_plot:
                _sm_kw = dict(
                    color="k",
                    ls="-.",
                    lw=1.5,
                    marker="x",
                    markersize=5,
                    markevery=max(1, int(_smask.sum()) // 10),
                    label=_s_label,
                )
                ax1.loglog(pd_arr[_smask] * 1e3, _V_s[_smask] / 1000, **_sm_kw)
                ax2.loglog(pd_arr[_smask] * 1e3, _EN_s[_smask], **_sm_kw)

            print(
                f"\n  {'pd (bar·mm)':>14}  {'p (bar)':>10}  {'d (mm)':>8}  "
                f"{'E/N (Td)':>12}  {'U (kV)':>12}  {'E (V/m)':>14}  [{_s_label}]"
            )
            print("  " + "-" * 90)
            _sstep = max(1, int(_smask.sum()) // 20)
            _sindices = np.where(_smask)[0][::_sstep]
            for _j in _sindices:
                print(
                    f"  {pd_arr[_j]*1e3:>14.4e}  "
                    f"{_p_arr_cs[_j]:>10.4g}  "
                    f"{_d_arr_cs[_j]*1e3:>8.4g}  "
                    f"{_EN_s[_j]:>12.4f}  "
                    f"{_V_s[_j]/1000:>12.4f}  "
                    f"{_V_s[_j]/_d_arr_cs[_j]:>14.4e}"
                )

    if args.write_to_file and _file_records:
        import datetime

        SEP = "\t"

        # Build flat column list: (name, value_fn) where value_fn(j) -> float
        columns = [("pd_bar_mm", lambda j: pd_arr[j] * 1e3)]
        for lbl, p_arr_c, d_arr_c, EN_c, V_c, _rec_mod in _file_records:
            columns += [
                (f"p_bar[{lbl}]", lambda j, a=p_arr_c: a[j]),
                (f"d_mm[{lbl}]", lambda j, a=d_arr_c: a[j] * 1e3),
                (f"EN_Td[{lbl}]", lambda j, a=EN_c: a[j]),
                (f"U_kV[{lbl}]", lambda j, a=V_c: a[j] / 1e3),
                (
                    f"E_Vm[{lbl}]",
                    lambda j, a=V_c, b=d_arr_c: (
                        a[j] / b[j] if np.isfinite(a[j]) else np.nan
                    ),
                ),
            ]

        if _plot_aed:
            for lbl, p_arr_c, d_arr_c, EN_c, V_c, _rec_mod in _file_records:
                columns.append(
                    (
                        f"ionization_integral[{lbl}]",
                        lambda j, _EN=EN_c, _p=p_arr_c, _d=d_arr_c, _m=_rec_mod: (
                            _aed_integral(_EN[j], _p[j], _d[j], _m)
                            if np.isfinite(_EN[j])
                            else np.nan
                        ),
                    )
                )

        # Alpha=eta crossover columns (one group per pressure where a root exists)
        for p_val, en_c in _en_cross.items():
            if en_c is None:
                continue
            lbl = f"alpha=eta, p={p_val} bar"
            columns += [
                (f"p_bar[{lbl}]", lambda j, _p=p_val: _p),
                (f"d_mm[{lbl}]", lambda j, _p=p_val: pd_arr[j] / _p * 1e3),
                (f"EN_Td[{lbl}]", lambda j, _en=en_c: _en),
                (
                    f"U_kV[{lbl}]",
                    lambda j, _en=en_c, _T=args.T: _en
                    * pd_arr[j]
                    * 1e-16
                    / (_kB * _T)
                    / 1e3,
                ),
                (
                    f"E_Vm[{lbl}]",
                    lambda j, _en=en_c, _p=p_val, _T=args.T: _en
                    * _p
                    * 1e5
                    / (_kB * _T)
                    * 1e-21,
                ),
            ]

        # Streamer criterion columns
        for _slbl, _sp_arr, _sd_arr, _sEN, _sV in _streamer_records:
            columns += [
                (f"EN_Td[{_slbl}]", lambda j, a=_sEN: a[j]),
                (f"U_kV[{_slbl}]", lambda j, a=_sV: a[j] / 1e3),
                (
                    f"E_Vm[{_slbl}]",
                    lambda j, a=_sV, b=_sd_arr: (
                        a[j] / b[j] if np.isfinite(a[j]) else np.nan
                    ),
                ),
            ]

        # Column width: wide enough for every header name plus a 2-char margin,
        # and at least 14 to fit a 12-char scientific-notation value.
        W = max(14, max(len(name) for name, _ in columns) + 2)

        import subprocess
        import sys as _sys

        try:
            _repo_dir = os.path.dirname(os.path.abspath(__file__))
            _git_hash = (
                subprocess.check_output(
                    ["git", "rev-parse", "--short", "HEAD"],
                    cwd=_repo_dir,
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .strip()
            )
            _dirty = (
                subprocess.check_output(
                    [
                        "git",
                        "status",
                        "--porcelain",
                        os.path.abspath(__file__),
                        os.path.abspath(args.mechanism),
                    ],
                    cwd=_repo_dir,
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .strip()
            )
            git_str = _git_hash + (" (dirty)" if _dirty else "")
        except Exception:
            git_str = "unavailable"

        with open(args.write_to_file, "w") as fh:
            fh.write("# --- METADATA ---\n")
            fh.write(
                f"# Date:    {datetime.datetime.now().isoformat(timespec='seconds')}\n"
            )
            fh.write(f"# Git:     {git_str}\n")
            fh.write(f"# Command: {' '.join(_sys.argv)}\n")
            fh.write("# ---\n")
            fh.write(f"# Mechanism:   {mech_name}\n")
            fh.write(f"# Temperature: {args.T} K\n")
            p_str = ", ".join(f"{p} bar" for p in args.p)
            fh.write(f"# Pressures:   {p_str}\n")
            cfg_labels = ", ".join(d.get("label", "Baseline") for d in raw_dicts)
            fh.write(f"# Configs:     {cfg_labels}\n")
            fh.write(f"# Field type:  {_field_dist.field_type}\n")
            fh.write(f"# Lambda:      {args.lam} s⁻¹\n")
            fh.write(f"# Method:      {args.method}\n")
            fh.write(
                f"# Stepping:    N_min={_N_min}, N_max={_N_max}, tol={_dx_tol:.3g}\n"
            )
            if _field_dist.field_type != "uniform":
                fh.write(f"# Sphere R:    {_field_dist.sphere_R*1e3:.4g} mm\n")
            if _field_dist.field_type == "sphere-plane":
                fh.write(
                    f"# Polarity:    sphere=positive → sphere is anode (+),  "
                    f"sphere=negative → sphere is cathode (−)\n"
                )
            if _streamer_records:
                fh.write(f"# Streamer C:  {args.streamer_criterion}\n")
            fh.write("#\n")

            # Per-column descriptions
            for col_idx, (name, _) in enumerate(columns, start=1):
                fh.write(f"# Column {col_idx}: {name}\n")
            fh.write("#\n")

            # Data rows — one per pd point
            for j in range(len(pd_arr)):
                vals = [fn(j) for _, fn in columns]
                fh.write(SEP.join(f"{v:<{W}.6e}" for v in vals) + "\n")

        print(f"\nResults written to: {args.write_to_file}")

    if not args.no_plot:
        if _en_cross:
            multi_p = len(args.p) > 1
            for p, en_c in _en_cross.items():
                if en_c is None:
                    continue
                lbl = (
                    rf"$\alpha=\eta$, p={p} bar ({en_c:.1f} Td)"
                    if multi_p
                    else rf"$\alpha=\eta$  ({en_c:.1f} Td)"
                )
                ax2.axhline(en_c, color="k", ls=":", lw=1.5, label=lbl)

        if ax3 is not None and _all_V_pos:
            _ref_cfg_label = raw_dicts[0].get("label", "Baseline")
            n_settings = len(args.p) + len(args.d)
            for curve_idx, (label, _p_arr, _d_arr, _cm) in enumerate(curve_specs):
                # Extract config label and p/d setting from the curve label.
                # Label format: "<config.label>, p=X bar" or "<config.label>, d=X mm"
                for _sep in [", p=", ", d="]:
                    if _sep in label:
                        cfg_lbl, setting = label.split(_sep, 1)
                        setting = _sep.lstrip(", ") + setting  # e.g. "p=1.0 bar"
                        break
                else:
                    continue
                if cfg_lbl == _ref_cfg_label:
                    continue
                # Find the matching reference (first config) curve.
                ref_label = f"{_ref_cfg_label}, {setting}"
                if ref_label not in _all_V_pos or label not in _all_V_pos:
                    continue
                V_k = _all_V_pos[label]
                V_b = _all_V_pos[ref_label]
                with np.errstate(invalid="ignore", divide="ignore"):
                    ratio = V_k / V_b
                mask_r = np.isfinite(ratio)
                if not np.any(mask_r):
                    continue
                ratio_label = f"{cfg_lbl} ({setting})" if n_settings > 1 else cfg_lbl
                marker = _markers[curve_idx % len(_markers)]
                markevery = max(1, mask_r.sum() // 15)
                ax3.semilogx(
                    pd_arr[mask_r] * 1e3,
                    ratio[mask_r],
                    label=ratio_label,
                    color=_curve_color.get(label),
                    marker=marker,
                    markevery=markevery,
                    markersize=6,
                )
            ax3.axhline(1.0, color="k", ls="--", lw=1.0)
            ax3.set_xlabel("p·d  (bar·mm)")
            ax3.set_ylabel(rf"$U / U_{{\mathrm{{{_ref_cfg_label}}}}}$")
            ax3.set_title(f"Voltage ratio vs. {_ref_cfg_label!r}")
            ax3.legend(loc="best", framealpha=1.0)
            ax3.grid(True, which="both", ls="--", alpha=0.4)

        ax1.set_xlabel("p·d  (bar·mm)")
        ax1.set_ylabel("U  (kV)")
        ax1.set_title("Breakdown voltage")
        ax1.legend(loc="upper left", framealpha=1.0)
        ax1.grid(True, which="both", ls="--", alpha=0.4)

        ax2.set_xlabel("p·d  (bar·mm)")
        ax2.set_ylabel("E/N  (Td)")
        ax2.set_title("Critical reduced field (average)")
        ax2.legend(loc="best", framealpha=1.0)
        ax2.grid(True, which="both", ls="--", alpha=0.4)

        plt.tight_layout()

        if args.save_subplots:
            import matplotlib.transforms as _mtrans

            mech_stem = os.path.splitext(mech_name)[0]

            # Panels: (primary_ax, twin_ax_or_None, output_filename)
            panels = [
                (ax1, None, f"{mech_stem}_U.pdf"),
                (ax2, None, f"{mech_stem}_EN.pdf"),
            ]
            if ax3 is not None:
                panels.append((ax3, None, f"{mech_stem}_ratio.pdf"))

            # Hide all titles before saving
            _suptitle.set_visible(False)
            _saved_titles = {ax: ax.get_title() for ax, *_ in panels}
            for ax, *_ in panels:
                ax.set_title("")

            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()

            for ax_p, ax_t, fname in panels:
                axes_in_panel = [a for a in (ax_p, ax_t) if a is not None]
                bboxes = [a.get_tightbbox(renderer) for a in axes_in_panel]
                bbox_px = _mtrans.Bbox.union(bboxes)
                bbox_in = bbox_px.transformed(fig.dpi_scale_trans.inverted())
                fig.savefig(fname, bbox_inches=bbox_in)
                print(f"Saved: {fname}")

            # Restore titles for the interactive window
            _suptitle.set_visible(True)
            for ax, *_ in panels:
                ax.set_title(_saved_titles[ax])

        plt.show()


if __name__ == "__main__":
    main()
