# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Gap geometry: the normalised field profile f(ξ) on ξ ∈ [0, 1].

Everything downstream of this module sees a gap only through f, normalised
so that ∫f dξ = 1, so the geometry — uniform, sphere-plane, sphere-sphere,
hyperboloid-plane, coaxial cylinders, or a tabulated field line — is
invisible to the solvers.  A protrusion on the ξ = 0 electrode
(``--protrusion``) multiplies the profile by its on-axis enhancement and
normalises it again, so the applied voltage is unchanged.
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
import sys
from typing import Callable, Optional

import numpy as np
from scipy.integrate import quad

from incept1d.protrusions import Protrusion, parse_protrusion_spec

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


# ── Hyperboloid-plane ────────────────────────────────────────────────────────


def _hyperboloid_plane_field(xi, k):
    """
    Exact normalised on-axis field between a hyperboloid tip and a plane.

    In prolate spheroidal coordinates with foci at ±a on the axis, the plane
    and the hyperboloid of revolution are both coordinate surfaces, and the
    potential depends on one coordinate alone.  With the tip at distance d
    from the plane and tip radius of curvature r,

        a = √(d(d + r)),   k = d/a = √(d/(d + r)),

    the on-axis field at height z above the plane is
    E(z) = U a / ((a² − z²) artanh k), and with z = d(1 − ξ)

        f(ξ) = k / ((1 − k²(1 − ξ)²) artanh k),

    so that ∫₀¹ f dξ = 1 exactly.  ξ = 0 is the tip, where
    E = 2U a / (r d ln((1 + k)/(1 − k))) ≈ 2U / (r ln(4d/r)) for r ≪ d
    [Coelho1971]_.  f depends on d/r alone, and tends to 1 as k → 0.
    """
    if k < 1e-7:
        return 1.0
    u = k * (1.0 - xi)
    return k / ((1.0 - u * u) * math.atanh(k))


# ── Protrusions ──────────────────────────────────────────────────────────────

#: A protrusion taller than this fraction of a length scale of the gap it sits
#: in (electrode radius, gap length, ...) is outside the validity of the
#: superposition in :meth:`FieldDistribution.build`; the commands say so.
PROTRUSION_SCALE_LIMIT = 0.1


# ── Tabulated field line ──────────────────────────────────────────────────────

_LENGTH_UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "µm": 1e-6}
_VOLTAGE_UNITS = {"mV": 1e-3, "V": 1.0, "kV": 1e3, "MV": 1e6}


def parse_field_unit(text):
    """
    The value in V/m of one unit of a field column, from ``[MULT*]VOLT/LENGTH``.

    VOLT is one of mV, V, kV, MV and LENGTH one of m, cm, mm, um; MULT is an
    optional positive factor for a column in scaled units.  ``"kV/mm"`` is
    1e6, ``"V/cm"`` 100, and ``"1e3*V/m"`` 1e3.

    Raises
    ------
    ValueError
        If *text* is not of that form.
    """
    mult, _, unit = text.rpartition("*")
    volt, sep, length = unit.partition("/")
    if not sep or volt not in _VOLTAGE_UNITS or length not in _LENGTH_UNITS:
        raise ValueError(
            f"field unit {text!r} is not [MULT*]VOLT/LENGTH with VOLT one of "
            f"{', '.join(_VOLTAGE_UNITS)} and LENGTH one of m, cm, mm, um"
        )
    try:
        factor = float(mult) if mult else 1.0
    except ValueError:
        raise ValueError(f"field unit {text!r}: multiplier {mult!r} is not a number")
    if not (factor > 0.0 and math.isfinite(factor)):
        raise ValueError(f"field unit {text!r}: the multiplier must be positive")
    return factor * _VOLTAGE_UNITS[volt] / _LENGTH_UNITS[length]


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


# ── Field-following panels ─────────────────────────────────────────────────

#: Largest relative change of the field across one initial segment.
FIELD_STEP_REL = 0.2


def field_following_edges(f, n_min, rel=FIELD_STEP_REL, max_edges=4097):
    """
    Initial segment edges in ξ ∈ [0, 1], refined where the field changes.

    Starts from *n_min* equal segments and halves every segment across
    which f changes by more than *rel* relative to its larger end, until
    none does.  The adaptive step halving that follows judges a segment by
    comparing one midpoint step with two half steps, and cannot see a
    feature that falls between those sample points: a thin high-field layer
    at a small electrode (a wire of radius a in a gap of 100 a, say) can
    then be missed altogether, and with it the avalanche.  A field
    following start places sample points in that layer.  For a uniform or
    gently varying field it returns the *n_min* equal segments unchanged.

    Shared by the solver's grid and the ionization integrals of
    :mod:`incept1d.ionization`, which use the edges as quadrature panels.
    """
    edges = list(np.linspace(0.0, 1.0, n_min + 1))
    values = [f(x) for x in edges]
    changed = True
    while changed and len(edges) < max_edges:
        changed = False
        new_edges, new_values = [edges[0]], [values[0]]
        for a, b, fa, fb in zip(edges, edges[1:], values, values[1:]):
            if abs(fa - fb) > rel * max(abs(fa), abs(fb)):
                mid = 0.5 * (a + b)
                new_edges.append(mid)
                new_values.append(f(mid))
                changed = True
            new_edges.append(b)
            new_values.append(fb)
        edges, values = new_edges, new_values
    return np.array(edges)


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
        One of 'uniform', 'sphere-plane', 'sphere-sphere',
        'hyperboloid-plane', 'coaxial', 'fieldline'.
    sphere_R : float or None
        Sphere radius in metres.  None unless field_type is a sphere geometry.
    tip_R : float or None
        For 'hyperboloid-plane': radius of curvature of the tip in metres.
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
        units times metres.  A voltage only once the field unit is known
        (:attr:`fieldline_line_voltage`).
    fieldline_field_unit : str or None
        For 'fieldline': the unit of the field column as declared
        (``"kV/mm"``, ``"1e3*V/m"``, ...), or None when it was not.
    fieldline_applied_voltage : float or None
        For 'fieldline': the excitation the line was computed at, in volts:
        declared with ``--fieldline-voltage``, or else ``∫|E| ds`` when the
        field unit is known.  The solve sees only the shape of the profile,
        so this does not affect the result; it is what the inception voltage
        is reported relative to.
    fieldline_voltage_source : str or None
        For 'fieldline': where :attr:`fieldline_applied_voltage` came from,
        for the reports.
    fieldline_path : str or None
        For 'fieldline': source file (for labels / output headers).
    reversed : bool
        The profile is read from the other electrode, f(ξ) → f(1 − ξ), so
        that ξ = 0 — where a protrusion sits, and the electrode the polarity
        labels name — is the plane, the outer conductor, or the last
        tabulated point.
    protrusion : Protrusion or None
        A protrusion on the ξ = 0 electrode (:mod:`incept1d.protrusions`),
        or None for a smooth electrode.
    """

    field_type: str
    sphere_R: Optional[float] = None
    tip_R: Optional[float] = None
    coax_a: Optional[float] = None
    coax_b: Optional[float] = None
    fieldline_xi: Optional[np.ndarray] = None
    fieldline_f: Optional[np.ndarray] = None
    fieldline_length: Optional[float] = None
    fieldline_reading: Optional[str] = None
    fieldline_integral: Optional[float] = None
    fieldline_field_unit: Optional[str] = None
    fieldline_applied_voltage: Optional[float] = None
    fieldline_voltage_source: Optional[str] = None
    fieldline_path: Optional[str] = None
    reversed: bool = False
    protrusion: Optional[Protrusion] = None

    @classmethod
    def from_fieldline(cls, path, length_unit="m", field_unit=None):
        """
        Build a 'fieldline' distribution from a tabulated ``|E|(s)`` file.

        The profile is parametrised by normalised arc length ξ = s/L and
        normalised so that ∫₀¹ f dξ = 1, which makes the gap length the arc
        length and the reference field the mean field along the line.  The
        absolute scale of the tabulated field therefore drops out: the same
        line exported at 100 kV and at 200 kV gives identical results.

        *field_unit* (``[MULT*]VOLT/LENGTH``, :func:`parse_field_unit`)
        declares the unit of the field column.  It leaves the profile alone
        and makes ``∫|E| ds`` a voltage, which then serves as the excitation.

        See :func:`load_fieldline` for the accepted file layouts.
        """
        s, E, reading = load_fieldline(path, length_unit)
        if field_unit is not None:
            parse_field_unit(field_unit)  # validate before building
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
            fieldline_field_unit=field_unit,
            fieldline_path=path,
        )

    @property
    def fieldline_line_voltage(self) -> Optional[float]:
        """``∫|E| ds`` in volts, when the field unit of the file is known."""
        if self.fieldline_field_unit is None:
            return None
        return self.fieldline_integral * parse_field_unit(self.fieldline_field_unit)

    @property
    def is_uniform(self) -> bool:
        """
        True when f ≡ 1: a uniform gap without a protrusion.

        The solvers take a constant-field shortcut on this, so test it
        rather than ``field_type == 'uniform'``, which a protrusion on a
        plate also has.
        """
        return self.field_type == "uniform" and self.protrusion is None

    @property
    def is_symmetric(self) -> bool:
        """
        True when the geometry is unchanged by swapping the two electrodes.

        This is a statement about the field alone.  Whether the two
        polarities give the same answer also depends on the electrode
        surfaces; see :func:`incept1d.solver.polarities_equivalent`.
        """
        if self.protrusion is not None:
            return False
        return self.field_type in ("uniform", "sphere-sphere")

    def polarity_label(self, polarity: str) -> str:
        """
        Name *polarity* (``"positive"`` or ``"negative"``) for this geometry.

        Says which electrode, the one at ξ = 0, is the anode in positive
        polarity: the sphere, the hyperboloid tip, the inner conductor, the
        first tabulated point, or, reversed, the plane, the outer conductor
        or the last point; on a uniform field, the electrode with the
        protrusion.
        """
        if self.field_type == "uniform" and self.protrusion is not None:
            return f"protrusion={polarity}"
        if self.field_type == "uniform":
            return polarity
        names = {
            "fieldline": ("start", "end"),
            "coaxial": ("inner", "outer"),
            "hyperboloid-plane": ("tip", "plane"),
            "sphere-plane": ("sphere", "plane"),
            "sphere-sphere": ("sphere", "sphere"),
        }[self.field_type]
        return f"{names[self.reversed]}={polarity}"

    @property
    def is_monotone_decreasing(self) -> bool:
        """
        True when f(ξ) is known to decrease monotonically from ξ = 0 to ξ = 1.

        Sphere-plane, hyperboloid-plane, coaxial (and trivially uniform)
        profiles have this property;
        sphere-sphere has maxima at both ends, and a tabulated field line may
        have any shape.  A protrusion keeps it, since its enhancement also
        falls monotonically away from the tip.  A reversed profile rises
        instead.  Used by quadrature helpers that would otherwise assume a
        single active region starting at ξ = 0.
        """
        if self.reversed:
            return False
        return self.field_type in (
            "uniform",
            "sphere-plane",
            "hyperboloid-plane",
            "coaxial",
        )

    @property
    def fixed_gap_length(self) -> Optional[float]:
        """
        Gap length in metres when the geometry itself fixes it, else None.

        A tabulated field line is pinned at its arc length, and a coaxial
        arrangement at ``b − a``.  For these a pd sweep is a pressure sweep,
        since any other d is the arrangement at a different size.  A
        protrusion shortens the gap by its height, since d is measured from
        its tip.
        """
        if self.field_type == "fieldline":
            D = self.fieldline_length
        elif self.field_type == "coaxial":
            D = self.coax_b - self.coax_a
        else:
            return None
        return D - (self.protrusion.height if self.protrusion else 0.0)

    @property
    def label(self) -> str:
        """Human-readable description for plot titles and file headers."""
        label = self._background_label
        if self.reversed:
            label += f" (reversed, xi = 0 at the {self.electrode_name})"
        if self.protrusion is None:
            return label
        return f"{label}, {self.protrusion.label}"

    @property
    def electrode_name(self) -> str:
        """The electrode at ξ = 0, as the polarity labels name it (``"plane"``, ...)."""
        return self.polarity_label("positive").split("=")[0]

    @property
    def _background_label(self) -> str:
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
        if self.field_type == "hyperboloid-plane":
            return f"hyperboloid-plane, r = {self.tip_R * 1e3:.4g} mm"
        R_mm = self.sphere_R * 1e3
        return f"{self.field_type}, R = {R_mm:.4g} mm"

    def build(self, d: float) -> Callable[[float], float]:
        """
        Return the normalised field callable f(xi) for gap length d.

        Parameters
        ----------
        d : float
            Gap length in metres: from the ξ = 0 electrode, or from the tip
            of its protrusion, to the other electrode.

        Returns
        -------
        f : callable
            ``f(xi)`` at fractional position ξ ∈ [0, 1], with ∫₀¹ f dξ = 1.
            ξ = 0 is the high-field electrode.  A tabulated profile is
            linearly interpolated and, like the coaxial profile, does not
            depend on d, since f is invariant under a geometric rescaling
            of the arrangement.

        Notes
        -----
        With a protrusion of height h the electrode surface is at
        D = d + h from the other electrode.  The profile is the background
        built for D, read from the tip onwards, times the protrusion's
        on-axis enhancement :meth:`Protrusion.enhancement`, and normalised
        again over the tip-to-electrode path — which keeps ∫E dx equal to
        the applied voltage.  That is exact on a uniform background and a
        local approximation on a curved one.
        """
        if self.protrusion is None:
            return self._background(d)

        pr = self.protrusion
        h = pr.height
        D = d + h
        f_bg = self._background(D)
        g = pr.enhancement

        def raw(xi, _d=d, _h=h, _D=D):
            s = xi * _d
            return f_bg((_h + s) / _D) * g(s)

        # The enhancement varies on the tip radius of curvature and on h;
        # break the quadrature at a few multiples of each.
        rho = pr.tip_radius
        breaks = sorted(
            {
                x / d
                for L in (rho, h)
                for x in (0.1 * L, L, 10.0 * L, 100.0 * L)
                if 0.0 < x < d
            }
        )
        edges = [0.0, *breaks, 1.0]
        norm = sum(
            quad(raw, a0, a1, epsabs=0.0, epsrel=1e-11, limit=500)[0]
            for a0, a1 in zip(edges[:-1], edges[1:])
        )

        def f_pr(xi, _norm=norm):
            return raw(xi) / _norm

        return f_pr

    def _background(self, d: float) -> Callable[[float], float]:
        """The profile of the electrodes alone, normalised over the gap d."""
        f = self._forward(d)
        if not self.reversed:
            return f

        def f_rev(xi, _f=f):
            return _f(1.0 - xi)

        return f_rev

    def _forward(self, d: float) -> Callable[[float], float]:
        """:meth:`_background` from the electrode the geometry names first."""
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

        if self.field_type == "hyperboloid-plane":
            k = math.sqrt(d / (d + self.tip_R))

            def f_hp(xi, _k=k):
                return _hyperboloid_plane_field(xi, _k)

            return f_hp

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
            "'hyperboloid-plane R_mm': hyperboloid of revolution above a plane, "
            "tip radius of curvature R_mm in mm (exact on-axis field).  "
            "'coaxial A_mm B_mm': coaxial cylinders with inner radius A_mm and "
            "outer radius B_mm in mm (gap b - a, inner conductor at xi = 0).  "
            "'fieldline FILE [LENGTH] [FIELD]': tabulated |E| along a (curved) "
            "field line read from FILE — numeric columns 's |E|', 'x y z |E|' "
            "or 'x y z Ex Ey Ez' (CSV or whitespace; header lines skipped); "
            "LENGTH is the length unit of the file (m, cm, mm, um; default m) "
            "and FIELD the unit of the field column, [MULT*]VOLT/LENGTH such "
            "as kV/mm, V/m or 1e3*V/cm, which makes the line integral the "
            "applied voltage.  "
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
            "relative to.  Taken from the line integral when the field unit "
            "is given in --field fieldline; without either, no such ratio is "
            "reported."
        ),
    )
    parser.add_argument(
        "--reverse-field",
        action="store_true",
        default=False,
        help=(
            "Read the profile from the other electrode, f(xi) -> f(1 - xi), so "
            "that xi = 0 -- the electrode a --protrusion sits on and the "
            "polarity labels name -- is the plane (sphere-plane, "
            "hyperboloid-plane), the outer conductor (coaxial) or the last "
            "tabulated point (fieldline).  No effect on a symmetric gap."
        ),
    )
    parser.add_argument(
        "--protrusion",
        nargs="+",
        default=None,
        metavar="SPEC",
        help=(
            "Protrusion on the xi = 0 electrode.  'spheroid H_mm B_mm': a "
            "half-spheroid of height H_mm and base radius B_mm (needle for "
            "H > B, hemisphere for H = B, flat bump for H < B).  "
            "'cone H_mm R_mm ANGLE_deg': a cone of half-angle ANGLE_deg "
            "with its apex rounded to radius R_mm.  'rod H_mm R_mm': a "
            "cylinder of radius R_mm with a hemispherical cap.  The gap d "
            "is then measured from its tip, and the profile is renormalised "
            "so that the voltage is unchanged.  Exact on a uniform field; "
            "on a curved electrode it needs H small against the radius, "
            "which the command checks."
        ),
    )


def field_from_args(args, parser=None) -> "FieldDistribution":
    """
    The :class:`FieldDistribution` described by the options that
    :func:`add_field_argument` registers, with its notes printed.
    """
    return parse_field_spec(
        args.field,
        parser,
        applied_voltage_kv=args.fieldline_voltage,
        reverse=args.reverse_field,
        protrusion=args.protrusion,
    )


def _fieldline_table(fd):
    """The tabulated (ξ, f), in the order the solve reads them (reversal applied)."""
    xi, f = fd.fieldline_xi, fd.fieldline_f
    if fd.reversed:
        return 1.0 - xi[::-1], f[::-1]
    return xi, f


def _fieldline_scale(fd):
    """
    Length over which a tabulated profile changes near ξ = 0, f/|df/ds|.

    The slope is a least-squares fit through the first few points, so a
    single noisy sample does not decide it.  Infinite for a flat start.
    """
    xi, f = _fieldline_table(fd)
    n = min(5, len(xi))
    slope = np.polyfit(xi[:n] * fd.fieldline_length, f[:n], 1)[0]
    if slope == 0.0:
        return math.inf
    return float(f[0] / abs(slope))


def protrusion_notes(fd, d_min=None) -> list:
    """
    Say where a protrusion is outside the validity of its model.

    The superposition in :meth:`FieldDistribution.build` needs the
    protrusion small against the length scale of the electrode it sits on
    (sphere radius, tip radius, inner or outer radius, or the scale over
    which a tabulated field changes) and against the gap, since the closed form
    leaves the other electrode an equipotential only to O((h/d)³).  A field
    line that is stronger at its far end probably has the protrusion on the
    wrong electrode.

    Parameters
    ----------
    fd : FieldDistribution
    d_min : float or None
        Smallest gap length (m) the run uses; None skips the gap check.

    Returns
    -------
    list of str
        One sentence per concern; empty when there is none or no protrusion.
    """
    if fd.protrusion is None:
        return []
    h = fd.protrusion.height
    notes = []
    scale = None
    # On a reversed plane the field varies over the gap, which the gap
    # check below covers.
    if fd.field_type == "sphere-sphere" or (
        fd.field_type == "sphere-plane" and not fd.reversed
    ):
        scale, what = fd.sphere_R, "the sphere radius R"
    elif fd.field_type == "hyperboloid-plane" and not fd.reversed:
        scale, what = fd.tip_R, "the tip radius r"
    elif fd.field_type == "coaxial" and not fd.reversed:
        scale, what = fd.coax_a, "the inner radius a"
    elif fd.field_type == "coaxial":
        scale, what = fd.coax_b, "the outer radius b"
    elif fd.field_type == "fieldline":
        scale = _fieldline_scale(fd)
        what = "the length f/|df/ds| over which the tabulated field changes"
        _, f = _fieldline_table(fd)
        if f[-1] > f[0]:
            notes.append(
                "the field line is stronger at its far end than at xi = 0, "
                "where the protrusion sits; --reverse-field puts the other "
                "end at xi = 0 if the protrusion belongs there"
            )
    if scale is not None and h > PROTRUSION_SCALE_LIMIT * scale:
        notes.append(
            f"protrusion height h = {h * 1e3:.4g} mm is {h / scale:.3g} of "
            f"{what} = {scale * 1e3:.4g} mm; the enhancement is superposed on "
            f"the background field, which is accurate only for h/scale below "
            f"about {PROTRUSION_SCALE_LIMIT:g}"
        )
    if d_min is not None and h > PROTRUSION_SCALE_LIMIT * d_min:
        notes.append(
            f"protrusion height h = {h * 1e3:.4g} mm is {h / d_min:.3g} of the "
            f"smallest gap d = {d_min * 1e3:.4g} mm; the closed form leaves the "
            f"other electrode an equipotential only to O((h/d)^3), accurate "
            f"for h/d below about {PROTRUSION_SCALE_LIMIT:g}"
        )
    return notes


def print_protrusion_notes(fd, d_min=None, gap_only=False) -> None:
    """
    Print :func:`protrusion_notes` to stderr, one ``note:`` line each.

    :func:`parse_field_spec` prints the notes on the geometry, and on the gap
    where the geometry fixes it; a command that chooses the gap itself calls
    this with its smallest gap and *gap_only* once it knows it.
    """
    if gap_only and (fd.protrusion is None or fd.fixed_gap_length):
        return
    notes = protrusion_notes(fd, d_min)
    if gap_only:  # the gap note comes last, after the geometry ones
        notes = notes[len(protrusion_notes(fd)) :]
    for note in notes:
        print(f"note: {note}.", file=sys.stderr)


def parse_field_spec(
    spec: list,
    parser=None,
    applied_voltage_kv=None,
    reverse=False,
    protrusion=None,
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
    reverse : bool
        Value of ``--reverse-field``: read the profile from the other
        electrode.  A symmetric gap is unchanged by it, which is noted.
    protrusion : list of str or None
        The raw ``--protrusion`` tokens (e.g. ['spheroid', '0.5', '0.1']).
        Where the geometry fixes the gap, the notes of
        :func:`protrusion_notes` are printed here; elsewhere the gap is not
        known yet and the caller prints them with its smallest gap.

    Returns
    -------
    FieldDistribution
    """

    def _err(msg):
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    fd = _parse_background(spec, _err, applied_voltage_kv)
    if reverse and fd.is_symmetric:
        print(
            f"note: --field {fd.field_type} is symmetric, so --reverse-field "
            f"changes nothing.",
            file=sys.stderr,
        )
    elif reverse:
        fd.reversed = True
    if protrusion is None:
        return fd

    try:
        fd.protrusion = parse_protrusion_spec(protrusion)
    except ValueError as exc:
        _err(f"--protrusion: {exc}")
    h = fd.protrusion.height
    fixed = fd.fixed_gap_length
    if fixed is not None and fixed <= 0.0:
        _err(
            f"--protrusion: height {h * 1e3:.4g} mm leaves no gap in the "
            f"{fd.field_type} geometry ({(fixed + h) * 1e3:.4g} mm)"
        )
    print_protrusion_notes(fd, fixed)
    return fd


def _parse_background(spec, _err, applied_voltage_kv):
    """The electrodes alone, from the ``--field`` tokens."""
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

    if field_type == "hyperboloid-plane":
        if len(spec) != 2:
            _err("--field hyperboloid-plane requires a tip radius in mm")
        try:
            tip_R = float(spec[1]) * 1e-3  # mm → m
        except ValueError:
            _err(
                f"--field hyperboloid-plane: tip radius must be a number, got {spec[1]!r}"
            )
        if not tip_R > 0.0:
            _err("--field hyperboloid-plane: tip radius must be positive")
        return FieldDistribution(field_type="hyperboloid-plane", tip_R=tip_R)

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
        if len(spec) > 4:
            _err("--field fieldline takes FILE [LENGTH] [FIELD]")
        # A field unit has a '/', a length unit does not; either may be given.
        fields = [t for t in spec[2:] if "/" in t]
        lengths = [t for t in spec[2:] if "/" not in t]
        if len(fields) > 1 or len(lengths) > 1:
            _err("--field fieldline takes at most one length and one field unit")
        try:
            fd = FieldDistribution.from_fieldline(
                spec[1], lengths[0] if lengths else "m", fields[0] if fields else None
            )
        except (OSError, ValueError) as exc:
            _err(f"--field fieldline: {exc}")
        line_voltage = fd.fieldline_line_voltage
        if applied_voltage_kv is not None:
            if applied_voltage_kv <= 0.0:
                _err("--fieldline-voltage must be positive")
            fd.fieldline_applied_voltage = applied_voltage_kv * 1e3
            fd.fieldline_voltage_source = "--fieldline-voltage"
            if line_voltage is not None and not math.isclose(
                line_voltage, applied_voltage_kv * 1e3, rel_tol=0.02
            ):
                print(
                    f"note: ∫|E| ds = {line_voltage / 1e3:.6g} kV with |E| in "
                    f"{fd.fieldline_field_unit}, but --fieldline-voltage is "
                    f"{applied_voltage_kv:g} kV; the declared excitation is "
                    f"used.  Check the units, the column, or whether the line "
                    f"spans the whole gap.",
                    file=sys.stderr,
                )
        elif line_voltage is not None:
            fd.fieldline_applied_voltage = line_voltage
            fd.fieldline_voltage_source = (
                f"∫|E| ds with |E| in {fd.fieldline_field_unit}"
            )
        return fd

    _err(
        f"Unknown field type: {field_type!r}. "
        "Use 'uniform', 'sphere-plane', 'sphere-sphere', 'hyperboloid-plane', "
        "'coaxial', or 'fieldline'."
    )


# ── Standalone field plotter ──────────────────────────────────────────────────
