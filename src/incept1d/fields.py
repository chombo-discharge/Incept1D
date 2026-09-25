# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Gap geometry: the normalised field profile f(ξ) on ξ ∈ [0, 1].

Everything downstream of this module sees a gap only through f, normalised
so that ∫f dξ = 1, so the geometry — uniform, sphere-plane, sphere-sphere,
coaxial cylinders, or a tabulated field line — is invisible to the solvers.
:class:`FieldDistribution` is the geometry itself, independent of gap
length; :meth:`FieldDistribution.build` turns it into f for a specific
gap.  :func:`add_field_argument` and :func:`parse_field_spec` implement
the ``--field`` option shared by every subcommand.

The standalone field plotter is ``incept1d field`` (:mod:`incept1d.cli.field`).
"""

import argparse
import dataclasses
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


# ── Coaxial cylinders ────────────────────────────────────────────────────────


def _coaxial_field(xi, ratio):
    """
    Exact normalised radial field between coaxial cylinders, b/a = ratio.

    E(r) = U / (r ln(b/a)) with r = a + ξ(b − a), so that

        f(ξ) = (ρ − 1) / ((1 + ξ(ρ − 1)) ln ρ),   ρ = b/a,

    and ∫₀¹ f dξ = 1 exactly.  ξ = 0 is the inner conductor.  f depends on
    the radius ratio alone, and tends to 1 as ρ → 1.
    """
    t = ratio - 1.0
    if t < 1e-12:
        return 1.0
    return t / ((1.0 + xi * t) * math.log1p(t))


# ── Step-count helpers (kept for backward compatibility / plot_sphere_plane.py) ──


# ── Tabulated field line ──────────────────────────────────────────────────────

_LENGTH_UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "µm": 1e-6}


def load_fieldline(path, length_unit="m"):
    """
    Read a tabulated electric field along a (possibly curved) field line.

    A plain numeric table of 2, 4 or 6 columns; the layout is inferred from
    the column count and documented in
    ``docs/source/modules/fielddistributions.rst``.  Rows are used in file
    order, so the first row is ξ = 0 and the last is ξ = 1.

    The single position column of a 2-column file may be an arc length or a
    coordinate, in either direction; both are read the same way, since the
    arc length is accumulated from it.  That only describes the line
    correctly when the line is straight, or when the column is an arc length
    already — one coordinate of a *curved* line understates the distance
    travelled, and such a line needs the 4- or 6-column layout.

    Parameters
    ----------
    path : str
        Path to the data file.
    length_unit : str
        Unit of the length columns: ``'m'``, ``'cm'``, ``'mm'`` or
        ``'um'``.  The field units are irrelevant, since only the shape of
        the profile survives normalisation.

    Returns
    -------
    s : ndarray
        Arc length in metres, starting at 0, strictly increasing.
    E : ndarray
        Field magnitude (arbitrary units, ≥ 0) at each ``s``.
    reading : str
        How the position columns were read, for the caller to report.
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
        # The first column is either an arc length already or a position
        # along the line, and the file does not say which.  It does not have
        # to: accumulating |Δ| reads an increasing arc length back unchanged,
        # and a decreasing column can only be a coordinate, so the two
        # readings differ only where the arc-length one is impossible.
        x = data[:, 0]
        dx = np.diff(x)
        if np.any(dx > 0.0) and np.any(dx < 0.0):
            raise ValueError(
                f"{path}: the position column runs both up and down, so it is "
                f"neither an arc length nor a coordinate along a straight "
                f"line.  Sort the rows along the line, or use the 4- or "
                f"6-column layout, which gets the arc length from the "
                f"positions and handles a curved line correctly."
            )
        s = np.concatenate(([0.0], np.cumsum(np.abs(dx))))
        s *= scale
        E = data[:, 1]
        reading = (
            "decreasing coordinate"
            if np.any(dx < 0.0)
            else "arc length or increasing coordinate"
        )
    elif ncol == 4:
        s = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(data[:, :3], axis=0), axis=1)))
        )
        s *= scale
        E = data[:, 3]
        reading = "cumulative distance between positions"
    elif ncol == 6:
        s = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(data[:, :3], axis=0), axis=1)))
        )
        s *= scale
        E = np.linalg.norm(data[:, 3:6], axis=1)
        reading = "cumulative distance between positions"
    else:
        raise ValueError(
            f"{path}: expected 2 (s,|E|), 4 (x,y,z,|E|) or 6 (x,y,z,Ex,Ey,Ez) "
            f"numeric columns, found {ncol}"
        )

    E = np.abs(E)
    s = s - s[0]
    if len(s) < 2:
        raise ValueError(f"{path}: need at least two points along the field line")
    # Every layout now accumulates distance, so s is non-decreasing by
    # construction and the direction the line was traced in does not matter.
    # Drop duplicate positions (zero-length segments) to keep interpolation sane.
    keep = np.concatenate(([True], np.diff(s) > 0.0))
    s, E = s[keep], E[keep]
    if s[-1] <= 0.0:
        raise ValueError(f"{path}: field line has zero length")
    if not np.all(np.isfinite(E)) or np.all(E == 0.0):
        raise ValueError(f"{path}: field magnitude must be finite and not all zero")
    return s, E, reading


# ── Field distribution class ──────────────────────────────────────────────────


@dataclasses.dataclass(eq=False)
class FieldDistribution:
    """
    Gap geometry, independent of gap length and grid resolution.

    Call :meth:`build` to obtain f(ξ) for a specific gap length.  The
    number of integration steps is the caller's business and is
    deliberately not stored here.

    Attributes
    ----------
    field_type : str
        One of 'uniform', 'sphere-plane', 'sphere-sphere', 'coaxial',
        'fieldline'.
    sphere_R : float or None
        Sphere radius in metres.  None unless field_type is a sphere geometry.
    coax_a, coax_b : float or None
        For 'coaxial': inner and outer conductor radii in metres, a < b.
    fieldline_xi, fieldline_f : ndarray or None
        For 'fieldline': tabulated normalised arc length ξ = s/L ∈ [0,1] and
        normalised field ``f(ξ) = |E|/⟨|E|⟩`` with ∫₀¹ f dξ = 1 (trapezoidal).
    fieldline_length : float or None
        For 'fieldline': total arc length L of the tabulated line in metres.
    fieldline_reading : str or None
        For 'fieldline': how the position columns were read, so that the
        commands can say it rather than leave the reader to guess.
    fieldline_integral : float or None
        For 'fieldline': ``∫|E| ds`` along the line, in the file's own field
        units times metres.  It is a voltage only if the file tabulates
        ``|E|`` in V/m, which the file does not say and this module cannot
        know, so it is reported as a units check and never used as one.
    fieldline_applied_voltage : float or None
        For 'fieldline': the excitation the line was computed at, in volts,
        as declared by the caller.  The solve sees only the shape of the
        profile, so this does not affect the result; it is what the
        inception voltage is reported relative to.
    fieldline_path : str or None
        For 'fieldline': source file (for labels / output headers).
    """

    field_type: str
    sphere_R: Optional[float] = None
    coax_a: Optional[float] = None
    coax_b: Optional[float] = None
    fieldline_xi: Optional[np.ndarray] = None
    fieldline_f: Optional[np.ndarray] = None
    fieldline_length: Optional[float] = None
    fieldline_reading: Optional[str] = None
    fieldline_integral: Optional[float] = None
    fieldline_applied_voltage: Optional[float] = None
    fieldline_path: Optional[str] = None

    @classmethod
    def from_fieldline(cls, path, length_unit="m"):
        """
        Build a 'fieldline' distribution from a tabulated ``|E|(s)`` file.

        The profile is parametrised by normalised arc length ξ = s/L and
        normalised so that ∫₀¹ f dξ = 1, which makes the gap length the arc
        length and the reference field the mean field along the line.  The
        absolute scale of the tabulated field therefore drops out: the same
        line exported at 100 kV and at 200 kV gives identical results.

        See :func:`load_fieldline` for the accepted file layouts.
        """
        s, E, reading = load_fieldline(path, length_unit)
        L = float(s[-1])
        xi = s / L
        _trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy ≥2.0 / <2.0
        mean_E = float(_trapz(E, xi))
        return cls(
            field_type="fieldline",
            fieldline_xi=xi,
            fieldline_f=E / mean_E,
            fieldline_length=L,
            fieldline_reading=reading,
            # ∫|E| ds = L ∫|E| dξ = L ⟨|E|⟩, in the field units of the file.
            fieldline_integral=mean_E * L,
            fieldline_path=path,
        )

    @property
    def is_symmetric(self) -> bool:
        """
        True when the geometry is unchanged by swapping the two electrodes.

        This is a statement about the field alone.  Whether the two
        polarities give the same answer also depends on the electrode
        surfaces; see :func:`incept1d.solver.polarities_equivalent`.
        """
        return self.field_type in ("uniform", "sphere-sphere")

    def polarity_label(self, polarity: str) -> str:
        """
        Name *polarity* (``"positive"`` or ``"negative"``) for this geometry.

        Says which electrode is the anode in positive polarity: the sphere,
        the inner conductor, or the first tabulated point.
        """
        if self.field_type == "uniform":
            return polarity
        if self.field_type == "fieldline":
            return f"start={polarity}"
        if self.field_type == "coaxial":
            return f"inner={polarity}"
        return f"sphere={polarity}"

    @property
    def is_monotone_decreasing(self) -> bool:
        """
        True when f(ξ) is known to decrease monotonically from ξ = 0 to ξ = 1.

        Sphere-plane (and trivially uniform) profiles have this property;
        sphere-sphere has maxima at both ends, and a tabulated field line may
        have any shape.  Used by quadrature helpers that would otherwise
        assume a single active region starting at ξ = 0.
        """
        return self.field_type in ("uniform", "sphere-plane", "coaxial")

    @property
    def fixed_gap_length(self) -> Optional[float]:
        """
        Gap length in metres when the geometry itself fixes it, else None.

        A tabulated field line is pinned at its arc length, and a coaxial
        arrangement at ``b − a``.  For these a pd sweep is a pressure sweep,
        since any other d is the arrangement at a different size.
        """
        if self.field_type == "fieldline":
            return self.fieldline_length
        if self.field_type == "coaxial":
            return self.coax_b - self.coax_a
        return None

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
        if self.field_type == "coaxial":
            return (
                f"coaxial, a = {self.coax_a * 1e3:.4g} mm, "
                f"b = {self.coax_b * 1e3:.4g} mm"
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
            ``f(xi)`` at fractional position ξ ∈ [0, 1], with ∫₀¹ f dξ = 1.
            ξ = 0 is the high-field electrode.  A tabulated profile is
            linearly interpolated and, like the coaxial profile, does not
            depend on d, since f is invariant under a geometric rescaling
            of the arrangement.
        """
        if self.field_type == "uniform":
            return lambda xi: 1.0

        if self.field_type == "fieldline":
            xi_tab, f_tab = self.fieldline_xi, self.fieldline_f

            def f_fl(xi, _xi=xi_tab, _f=f_tab):
                return float(np.interp(xi, _xi, _f))

            return f_fl

        if self.field_type == "coaxial":
            ratio = self.coax_b / self.coax_a

            def f_cx(xi, _ratio=ratio):
                return _coaxial_field(xi, _ratio)

            return f_cx

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
            "'coaxial A_mm B_mm': coaxial cylinders with inner radius A_mm and "
            "outer radius B_mm in mm (gap b - a, inner conductor at xi = 0).  "
            "'fieldline FILE [UNIT]': tabulated |E| along a (curved) field line "
            "read from FILE — numeric columns 's |E|', 'x y z |E|' or "
            "'x y z Ex Ey Ez' (CSV or whitespace; header lines skipped); "
            "UNIT is the length unit of the file (m, cm, mm, um; default m).  "
            "xi = 0 is the first row.  "
            "Grid resolution is set separately via --dx."
        ),
    )
    parser.add_argument(
        "--fieldline-voltage",
        type=float,
        default=None,
        metavar="U_KV",
        help=(
            "Excitation in kV at which a tabulated field line was computed.  "
            "Only the shape of the profile enters the solve, so this does not "
            "change the result; it is what the inception voltage is reported "
            "relative to.  Without it no such ratio is reported, since the "
            "field units of the file are unknown."
        ),
    )


def parse_field_spec(
    spec: list, parser=None, applied_voltage_kv=None
) -> "FieldDistribution":
    """
    Parse a --field token list into a FieldDistribution.

    Parameters
    ----------
    spec : list of str
        The raw token list from argparse (e.g. ['sphere-plane', '50']).
    parser : argparse.ArgumentParser or None
        If given, error messages are routed through parser.error(); otherwise
        a ValueError is raised.
    applied_voltage_kv : float or None
        Value of ``--fieldline-voltage``.  Meaningful only for a tabulated
        field line; supplying it for any other geometry is an error, since
        there the voltage is an output rather than a property of the input.

    Returns
    -------
    FieldDistribution
    """

    def _err(msg):
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    field_type = spec[0]
    if applied_voltage_kv is not None and field_type != "fieldline":
        _err(
            "--fieldline-voltage applies only to '--field fieldline'; for "
            "every other geometry the voltage is a result, not an input"
        )
    if field_type == "uniform":
        return FieldDistribution(field_type="uniform")

    if field_type in ("sphere-plane", "sphere-sphere"):
        if len(spec) < 2:
            _err(f"--field {field_type} requires a sphere radius in mm")
        sphere_R = float(spec[1]) * 1e-3  # mm → m
        return FieldDistribution(field_type=field_type, sphere_R=sphere_R)

    if field_type == "coaxial":
        if len(spec) != 3:
            _err("--field coaxial requires an inner and an outer radius in mm")
        try:
            a, b = float(spec[1]) * 1e-3, float(spec[2]) * 1e-3  # mm → m
        except ValueError:
            _err(f"--field coaxial: radii must be numbers, got {spec[1:]}")
        if not 0.0 < a < b:
            _err("--field coaxial: need 0 < inner radius < outer radius")
        return FieldDistribution(field_type="coaxial", coax_a=a, coax_b=b)

    if field_type == "fieldline":
        if len(spec) < 2:
            _err("--field fieldline requires a data file path")
        unit = spec[2] if len(spec) > 2 else "m"
        try:
            fd = FieldDistribution.from_fieldline(spec[1], unit)
        except (OSError, ValueError) as exc:
            _err(f"--field fieldline: {exc}")
        if applied_voltage_kv is not None:
            if applied_voltage_kv <= 0.0:
                _err("--fieldline-voltage must be positive")
            fd.fieldline_applied_voltage = applied_voltage_kv * 1e3
        return fd

    _err(
        f"Unknown field type: {field_type!r}. "
        "Use 'uniform', 'sphere-plane', 'sphere-sphere', 'coaxial', or "
        "'fieldline'."
    )


# ── Standalone field plotter ──────────────────────────────────────────────────
