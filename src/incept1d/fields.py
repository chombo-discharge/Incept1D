# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Field geometry abstraction for gap-discharge computations.

Provides:
  FieldDistribution   — geometry spec independent of gap length d; .build(d) → f(xi)
  add_field_argument  — add --field SPEC to an argparse.ArgumentParser
  parse_field_spec    — parse a --field token list → (FieldDistribution, N)

Low-level functions (moved here from Inception.py):
  _sphere_sphere_axial_field   — exact bispherical image-charge series
  _compute_n_steps             — minimum N for ≤ field_tol relative change in first step
  _n_steps_sphere_cached       — LRU-cached version of _compute_n_steps for sphere geometry
  load_fieldline               — read a tabulated ``|E|(s)`` along a (curved) field line

The standalone field plotter is ``incept1d field`` (:mod:`incept1d.cli.field`).
"""

import argparse
import dataclasses
import functools
import math
import os
from typing import Callable, Optional

import numpy as np

# ── Bispherical image-charge series ──────────────────────────────────────────


def _sphere_sphere_axial_field(xi, alpha, a_focal, d):
    """
    Exact normalised on-axis field between two equal spheres at ±U/2.

    Bispherical image-charge series:
        a  = R · sinh(α),   z_n = a · coth(n·α),   q_n = ½a / sinh(n·α)
        f(ξ) = d · Σ_{n=1}^∞ q_n · [1/(z_n−z)² + 1/(z_n+z)²],   z = d(½−ξ)

    Returns f(ξ) satisfying ∫₀¹ f(ξ) dξ = 1 by construction.  Vectorised
    over n; truncated at n = 30/α (relative tail < 1e-13) or na = 700.
    """
    if alpha < 1e-14:
        return 1.0
    z = d * (0.5 - xi)
    n_max = min(100000, max(50, int(30.0 / alpha) + 10))
    na = np.arange(1, n_max + 1, dtype=float) * alpha
    na = na[na < 700.0]
    sinh_na = np.sinh(na)
    z_n = a_focal * np.cosh(na) / sinh_na
    q_n = 0.5 * a_focal / sinh_na
    return float(np.sum(q_n * (1.0 / (z_n - z) ** 2 + 1.0 / (z_n + z) ** 2))) * d


# ── Step-count helpers (kept for backward compatibility / plot_sphere_plane.py) ──


def _compute_n_steps(field_func, field_tol, max_steps=2000):
    """Return minimum N so that the relative field change across the first step ≤ field_tol."""
    f0 = field_func(0.0)
    if f0 <= 0.0 or field_tol >= 1.0:
        return 1
    lo, hi = 1, max_steps
    if (f0 - field_func(1.0 / hi)) / f0 > field_tol:
        return hi
    while lo < hi:
        mid = (lo + hi) // 2
        if (f0 - field_func(1.0 / mid)) / f0 <= field_tol:
            hi = mid
        else:
            lo = mid + 1
    return lo


@functools.lru_cache(maxsize=1024)
def _n_steps_sphere_cached(alpha, is_plane, field_tol, max_steps):
    """N_steps for sphere geometry — depends only on shape parameter alpha, not scale."""
    sinh_a = math.sinh(float(alpha))
    cosh_a = math.cosh(float(alpha))
    a_norm = sinh_a / (2.0 * (cosh_a - 1.0))
    if is_plane:

        def _f(xi):
            return _sphere_sphere_axial_field(0.5 * xi, float(alpha), a_norm, 1.0)

    else:

        def _f(xi):
            return _sphere_sphere_axial_field(xi, float(alpha), a_norm, 1.0)

    return _compute_n_steps(_f, field_tol, max_steps)


# ── Tabulated field line ──────────────────────────────────────────────────────

_LENGTH_UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "µm": 1e-6}


def load_fieldline(path, length_unit="m"):
    """
    Read a tabulated electric field along a (possibly curved) field line.

    The file is a plain numeric table (whitespace- or comma-separated; header
    and ``#`` comment lines are skipped).  The column layout is inferred from
    the number of numeric columns:

    ========  ===============================================================
    columns   interpretation
    ========  ===============================================================
    2         ``s, |E|`` — arc length and field magnitude
    4         ``x, y, z, |E|`` — position and field magnitude
    6         ``x, y, z, Ex, Ey, Ez`` — position and field vector
    ========  ===============================================================

    For the 4- and 6-column layouts the arc length is the cumulative
    point-to-point Euclidean distance.  Rows are used in file order: the first
    row defines ``ξ = 0`` and the last row ``ξ = 1``.  Only ``|E|`` enters the
    1-D model, so the sign/direction of the field vector is irrelevant.

    Parameters
    ----------
    path : str
        Path to the data file.
    length_unit : str
        Unit of the length columns: one of ``'m'``, ``'cm'``, ``'mm'``,
        ``'um'``.  Field units are irrelevant (the profile is normalised).

    Returns
    -------
    s : ndarray
        Arc length in metres, starting at 0, strictly increasing.
    E : ndarray
        Field magnitude (arbitrary units, ≥ 0) at each ``s``.
    """
    if length_unit not in _LENGTH_UNITS:
        raise ValueError(
            f"Unknown length unit {length_unit!r}; use one of "
            f"{sorted(_LENGTH_UNITS)}"
        )
    scale = _LENGTH_UNITS[length_unit]

    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            tokens = line.replace(",", " ").replace(";", " ").split()
            try:
                rows.append([float(t) for t in tokens])
            except ValueError:
                continue  # header / non-numeric line
    if not rows:
        raise ValueError(f"No numeric data found in {path}")
    ncol = len(rows[0])
    if any(len(r) != ncol for r in rows):
        raise ValueError(f"Inconsistent number of columns in {path}")
    data = np.asarray(rows, dtype=float)

    if ncol == 2:
        s = data[:, 0] * scale
        E = data[:, 1]
    elif ncol == 4:
        s = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(data[:, :3], axis=0), axis=1)))
        )
        s *= scale
        E = data[:, 3]
    elif ncol == 6:
        s = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(data[:, :3], axis=0), axis=1)))
        )
        s *= scale
        E = np.linalg.norm(data[:, 3:6], axis=1)
    else:
        raise ValueError(
            f"{path}: expected 2 (s,|E|), 4 (x,y,z,|E|) or 6 (x,y,z,Ex,Ey,Ez) "
            f"numeric columns, found {ncol}"
        )

    E = np.abs(E)
    s = s - s[0]
    if len(s) < 2:
        raise ValueError(f"{path}: need at least two points along the field line")
    if np.any(np.diff(s) < 0.0):
        raise ValueError(f"{path}: arc length must be monotonically increasing")
    # Drop duplicate positions (zero-length segments) to keep interpolation sane.
    keep = np.concatenate(([True], np.diff(s) > 0.0))
    s, E = s[keep], E[keep]
    if s[-1] <= 0.0:
        raise ValueError(f"{path}: field line has zero length")
    if not np.all(np.isfinite(E)) or np.all(E == 0.0):
        raise ValueError(f"{path}: field magnitude must be finite and not all zero")
    return s, E


# ── Field distribution class ──────────────────────────────────────────────────


@dataclasses.dataclass(eq=False)
class FieldDistribution:
    """
    Field geometry specification, independent of gap length d and grid resolution.

    Call .build(d) to obtain the normalised field callable f(xi) for a specific
    gap length.  The number of integration steps N is a separate concern and is
    not stored here — it is supplied by the caller (CLI arg, adaptive scheme, etc.).

    Attributes
    ----------
    field_type : str
        One of 'uniform', 'sphere-plane', 'sphere-sphere', 'fieldline'.
    sphere_R : float or None
        Sphere radius in metres.  None unless field_type is a sphere geometry.
    fieldline_xi, fieldline_f : ndarray or None
        For 'fieldline': tabulated normalised arc length ξ = s/L ∈ [0,1] and
        normalised field ``f(ξ) = |E|/⟨|E|⟩`` with ∫₀¹ f dξ = 1 (trapezoidal).
    fieldline_length : float or None
        For 'fieldline': total arc length L of the tabulated line in metres.
    fieldline_path : str or None
        For 'fieldline': source file (for labels / output headers).
    """

    field_type: str
    sphere_R: Optional[float] = None
    fieldline_xi: Optional[np.ndarray] = None
    fieldline_f: Optional[np.ndarray] = None
    fieldline_length: Optional[float] = None
    fieldline_path: Optional[str] = None

    @classmethod
    def from_fieldline(cls, path, length_unit="m"):
        """
        Build a 'fieldline' distribution from a tabulated ``|E|(s)`` file.

        The profile is parametrised by normalised arc length ξ = s/L and
        normalised so that ∫₀¹ f(ξ) dξ = 1 (trapezoidal rule on the data
        grid, consistent with the linear interpolation used in :meth:`build`).
        Hence E_ref = V/L is the uniform-equivalent field along the line, and
        ``V = ∫|E| ds`` is the voltage drop along the field line.

        See :func:`load_fieldline` for the accepted file layouts.
        """
        s, E = load_fieldline(path, length_unit)
        L = float(s[-1])
        xi = s / L
        _trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy ≥2.0 / <2.0
        mean_E = float(_trapz(E, xi))
        return cls(
            field_type="fieldline",
            fieldline_xi=xi,
            fieldline_f=E / mean_E,
            fieldline_length=L,
            fieldline_path=path,
        )

    @property
    def is_symmetric(self) -> bool:
        """True when both polarities give the same inception result."""
        return self.field_type in ("uniform", "sphere-sphere")

    @property
    def is_monotone_decreasing(self) -> bool:
        """
        True when f(ξ) is known to decrease monotonically from ξ = 0 to ξ = 1.

        Sphere-plane (and trivially uniform) profiles have this property;
        sphere-sphere has maxima at both ends, and a tabulated field line may
        have any shape.  Used by quadrature helpers that would otherwise
        assume a single active region starting at ξ = 0.
        """
        return self.field_type in ("uniform", "sphere-plane")

    @property
    def label(self) -> str:
        """Human-readable description for plot titles and file headers."""
        if self.field_type == "uniform":
            return "uniform field"
        if self.field_type == "fieldline":
            return (
                f"field line {os.path.basename(self.fieldline_path)}, "
                f"L = {self.fieldline_length * 1e3:.4g} mm"
            )
        R_mm = self.sphere_R * 1e3
        return f"{self.field_type}, R = {R_mm:.4g} mm"

    def build(self, d: float) -> Callable[[float], float]:
        """
        Return the normalised field callable f(xi) for gap length d.

        Parameters
        ----------
        d : float
            Gap length in metres.

        Returns
        -------
        f : callable
            f(xi) → normalised field at fractional position xi ∈ [0,1].
            xi = 0 is the sphere surface (field maximum); xi = 1 is the plane
            / far sphere.  Satisfies ∫₀¹ f(xi) dxi = 1.
            For uniform field, f(xi) = 1.0 everywhere.
            For 'fieldline', f is linear interpolation of the tabulated profile
            (xi = 0 is the first row of the data file) and is independent of
            d: a gap length d ≠ L corresponds to the same electrode arrangement
            scaled geometrically by d/L, which leaves f(xi) unchanged.
        """
        if self.field_type == "uniform":
            return lambda xi: 1.0

        if self.field_type == "fieldline":
            xi_tab, f_tab = self.fieldline_xi, self.fieldline_f

            def f_fl(xi, _xi=xi_tab, _f=f_tab):
                return float(np.interp(xi, _xi, _f))

            return f_fl

        if self.field_type == "sphere-plane":
            alpha = float(np.arccosh(1.0 + d / self.sphere_R))
            a = self.sphere_R * np.sinh(alpha)

            def f_sp(xi, _alpha=alpha, _a=a, _d=d):
                return _sphere_sphere_axial_field(0.5 * xi, _alpha, _a, 2.0 * _d)

            return f_sp

        # sphere-sphere
        alpha = float(np.arccosh(1.0 + d / (2.0 * self.sphere_R)))
        a = self.sphere_R * np.sinh(alpha)

        def f_ss(xi, _alpha=alpha, _a=a, _d=d):
            return _sphere_sphere_axial_field(xi, _alpha, _a, _d)

        return f_ss


# ── CLI utilities (shared --field argument) ─────────────────────────────────────────────────────────────


def add_field_argument(parser: argparse.ArgumentParser) -> None:
    """Add --field SPEC argument to an argparse.ArgumentParser."""
    parser.add_argument(
        "--field",
        nargs="+",
        default=["uniform"],
        metavar="SPEC",
        help=(
            "Field distribution.  'uniform' (default): spatially uniform field.  "
            "'sphere-plane R_mm': sphere-plane gap with sphere radius R_mm in mm.  "
            "'sphere-sphere R_mm': symmetric sphere-sphere gap.  "
            "'fieldline FILE [UNIT]': tabulated |E| along a (curved) field line "
            "read from FILE — numeric columns 's |E|', 'x y z |E|' or "
            "'x y z Ex Ey Ez' (CSV or whitespace; header lines skipped); "
            "UNIT is the length unit of the file (m, cm, mm, um; default m).  "
            "xi = 0 is the first row.  "
            "Grid resolution is set separately via --dx."
        ),
    )


def parse_field_spec(spec: list, parser=None) -> "FieldDistribution":
    """
    Parse a --field token list into a FieldDistribution.

    Parameters
    ----------
    spec : list of str
        The raw token list from argparse (e.g. ['sphere-plane', '50']).
    parser : argparse.ArgumentParser or None
        If given, error messages are routed through parser.error(); otherwise
        a ValueError is raised.

    Returns
    -------
    FieldDistribution
    """

    def _err(msg):
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    field_type = spec[0]
    if field_type == "uniform":
        return FieldDistribution(field_type="uniform")

    if field_type in ("sphere-plane", "sphere-sphere"):
        if len(spec) < 2:
            _err(f"--field {field_type} requires a sphere radius in mm")
        sphere_R = float(spec[1]) * 1e-3  # mm → m
        return FieldDistribution(field_type=field_type, sphere_R=sphere_R)

    if field_type == "fieldline":
        if len(spec) < 2:
            _err("--field fieldline requires a data file path")
        unit = spec[2] if len(spec) > 2 else "m"
        try:
            return FieldDistribution.from_fieldline(spec[1], unit)
        except (OSError, ValueError) as exc:
            _err(f"--field fieldline: {exc}")

    _err(
        f"Unknown field type: {field_type!r}. "
        "Use 'uniform', 'sphere-plane', 'sphere-sphere', or 'fieldline'."
    )


# ── Standalone field plotter ──────────────────────────────────────────────────
