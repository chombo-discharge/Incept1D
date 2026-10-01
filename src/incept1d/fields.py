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
import heapq
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

    Parameters
    ----------
    text : str
        The unit as written on the command line.

    Returns
    -------
    float
        One unit of the column, in V/m.

    Raises
    ------
    ValueError
        If *text* is not of that form.
    """
    mult, star, unit = text.rpartition("*")
    volt, sep, length = unit.partition("/")
    if not sep or volt not in _VOLTAGE_UNITS or length not in _LENGTH_UNITS:
        raise ValueError(
            f"field unit {text!r} is not [MULT*]VOLT/LENGTH with VOLT one of "
            f"{', '.join(_VOLTAGE_UNITS)} and LENGTH one of m, cm, mm, um"
        )
    try:
        factor = float(mult) if star else 1.0
    except ValueError:
        raise ValueError(f"field unit {text!r}: multiplier {mult!r} is not a number")
    if not math.isfinite(factor):
        raise ValueError(f"field unit {text!r}: the multiplier must be finite")
    if not factor > 0.0:
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

#: Field (in units of the mean) below which a change is judged absolutely
#: rather than relative to the field itself.
FIELD_STEP_FLOOR = 1e-3


def field_following_edges(
    f, n_min, rel=FIELD_STEP_REL, max_edges=4097, nodes=None, floor=FIELD_STEP_FLOOR
):
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

    Comparing the ends of a segment cannot see a feature *inside* it — a
    narrow peak halfway along a tabulated field line.  Where f is known at
    *nodes* (the tabulated points, :meth:`FieldDistribution.nodes`), the
    values at the nodes inside a segment are compared as well, so such a
    peak is refined like any other change.

    Shared by the solver's grid and the ionization integrals of
    :mod:`incept1d.ionization`, which use the edges as quadrature panels.

    Parameters
    ----------
    f : callable
        The normalised profile f(ξ), of mean 1.
    n_min : int
        Number of equal segments to start from.
    rel : float
        Largest change of f across one segment, relative to the larger of
        its values there.
    max_edges : int
        Largest number of edges returned.  Where more segments would need
        splitting than that leaves room for, those across which f changes
        most are split first, so a noisy table or a field that falls to zero
        cannot run away with the cost.
    nodes : array_like or None
        Positions in (0, 1) where f is tabulated, if it is.
    floor : float
        A change is judged relative to at least this value, so that a field
        falling to zero does not ask for ever finer segments next to it.

    Returns
    -------
    ndarray
        Increasing edges from 0 to 1.
    """
    nodes = np.sort(np.asarray([] if nodes is None else nodes, dtype=float))
    node_values = np.array([f(x) for x in nodes])

    def variation(a, b, fa, fb):
        i0 = np.searchsorted(nodes, a, side="right")
        i1 = np.searchsorted(nodes, b, side="left")
        vals = [fa, fb, *node_values[i0:i1]]
        hi, lo = max(vals), min(vals)
        return (hi - lo) / max(abs(hi), abs(lo), floor)

    # Split the segment across which f changes most, one at a time, so that
    # the cap, if reached, has been spent where the field changes most.
    # Without the cap the result does not depend on the order: a segment is
    # split exactly when it varies by more than rel.
    xs = np.linspace(0.0, 1.0, n_min + 1)
    fs = [f(x) for x in xs]
    edges = set(xs)
    heap = [
        (-variation(a, b, fa, fb), a, b, fa, fb)
        for a, b, fa, fb in zip(xs, xs[1:], fs, fs[1:])
    ]
    heapq.heapify(heap)
    while heap and len(edges) < max_edges:
        neg, a, b, fa, fb = heapq.heappop(heap)
        if -neg <= rel:
            break
        mid = 0.5 * (a + b)
        fm = f(mid)
        edges.add(mid)
        heapq.heappush(heap, (-variation(a, mid, fa, fm), a, mid, fa, fm))
        heapq.heappush(heap, (-variation(mid, b, fm, fb), mid, b, fm, fb))
    return np.array(sorted(edges))


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

    def nodes(self, d: float) -> Optional[np.ndarray]:
        """
        Where a tabulated profile has its data points, in the ξ of :meth:`build`.

        The positions inside (0, 1) of the tabulated points of a field line,
        after reversal and after the shift to the tip of a protrusion; None
        for an analytic profile.  Between them the background is linear, so
        they are what a grid or a quadrature must not step over
        (:func:`field_following_edges`).

        Parameters
        ----------
        d : float
            Gap length in metres, as passed to :meth:`build`.

        Returns
        -------
        ndarray or None
            Increasing positions in (0, 1), or None.
        """
        if self.field_type != "fieldline":
            return None
        xi = self.fieldline_xi
        if self.reversed:
            xi = 1.0 - xi[::-1]
        if self.protrusion is not None:
            h = self.protrusion.height
            xi = (xi * (d + h) - h) / d
        return xi[(xi > 0.0) & (xi < 1.0)]

    @property
    def fieldline_line_voltage(self) -> Optional[float]:
        """``∫|E| ds`` in volts, when the field unit of the file is known."""
        if self.fieldline_field_unit is None:
            return None
        return self.fieldline_integral * parse_field_unit(self.fieldline_field_unit)

    @property
    def is_uniform(self) -> bool:
        """
        True when f ≡ 1, that is, for a uniform gap without a protrusion.

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
        or the last point; on a uniform or sphere-sphere gap, the electrode with the
        protrusion.
        """
        # Where both electrodes would carry the same name, the protrusion
        # tells them apart.
        if self.field_type in ("uniform", "sphere-sphere") and self.protrusion:
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
            ξ = 0 is the electrode the geometry names first (sphere, tip,
            inner conductor, first tabulated point), or the other one when
            :attr:`reversed`.  A tabulated profile is linearly interpolated
            and, like the coaxial profile, does not depend on d without a
            protrusion, since f is invariant under a geometric rescaling of
            the arrangement; a protrusion keeps its size as d changes.

        Notes
        -----
        With a protrusion of height h the electrode surface is at
        D = d + h from the other electrode.  The profile is the background
        built for D, read from the tip onwards, times the protrusion's
        on-axis enhancement :meth:`incept1d.protrusions.Protrusion.enhancement`, and normalised
        again over the tip-to-electrode path — which keeps ∫E dx equal to
        the applied voltage.  That is exact on a uniform background and a
        local approximation on a curved one.
        """
        if self.protrusion is None:
            return self._background(d)
        # The normalisation is a quadrature, and the solvers build the same
        # gap many times over a root search.
        return self._cached(("build", d), lambda: self._build_protruded(d))

    def grid_edges(self, d: float, n_min: int, max_edges: int = 4097) -> np.ndarray:
        """
        :func:`field_following_edges` of :meth:`build` for the gap *d*, cached.

        The edges depend on the profile alone, not on the field strength, so
        a root search over E/N at one gap computes them once.

        Parameters
        ----------
        d : float
            Gap length in metres.
        n_min : int
            Number of equal segments to start from.
        max_edges : int
            Largest number of edges.

        Returns
        -------
        ndarray
            Increasing edges in ξ from 0 to 1.
        """
        return self._cached(
            ("edges", d, n_min, max_edges),
            lambda: field_following_edges(
                self.build(d), n_min, max_edges=max_edges, nodes=self.nodes(d)
            ),
        )

    def _cached(self, key, compute):
        """*compute()*, remembered under *key* and the current geometry."""
        geometry = (
            self.field_type,
            self.sphere_R,
            self.tip_R,
            self.coax_a,
            self.coax_b,
            id(self.fieldline_f),
            self.reversed,
            self.protrusion,
        )
        cache = self.__dict__.setdefault("_cache", {})
        full = (geometry, *key)
        if full not in cache:
            cache[full] = compute()
        return cache[full]

    def _build_protruded(self, d: float) -> Callable[[float], float]:
        """:meth:`build` with a protrusion, uncached."""
        pr = self.protrusion
        h = pr.height
        D = d + h
        f_bg = self._background(D)
        g = pr.enhancement

        def raw(xi, _d=d, _h=h, _D=D):
            s = xi * _d
            return f_bg((_h + s) / _D) * g(s)

        # The enhancement varies on the tip radius of curvature and on h;
        # break the quadrature at a few multiples of each, and at the kinks
        # of a tabulated background.
        rho = pr.tip_radius
        breaks = sorted(
            {
                x / d
                for L in (rho, h)
                for x in (0.1 * L, L, 10.0 * L, 100.0 * L)
                if 0.0 < x < d
            }
            | set(self.nodes(d) if self.field_type == "fieldline" else ())
        )
        edges = np.array([0.0, *breaks, 1.0])
        # A panel much narrower than its distance from the tip sees a smooth
        # enhancement and, on a tabulated background, a linear field:
        # Simpson's rule is exact enough there and costs two evaluations,
        # which keeps a table of thousands of points cheap.  Near the tip,
        # adaptive quadrature.
        widths = np.diff(edges) * d
        smooth = widths < 0.02 * (edges[:-1] * d + rho)
        at_edges = {}

        def raw_at(x):
            if x not in at_edges:
                at_edges[x] = raw(x)
            return at_edges[x]

        norm = 0.0
        for a0, a1, easy in zip(edges[:-1], edges[1:], smooth):
            if easy:
                mid = raw(0.5 * (a0 + a1))
                norm += (a1 - a0) * (raw_at(a0) + 4.0 * mid + raw_at(a1)) / 6.0
            else:
                norm += quad(raw, a0, a1, epsabs=0.0, epsrel=1e-11, limit=500)[0]

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

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments, with ``field``, ``fieldline_voltage``,
        ``reverse_field`` and ``protrusion``.
    parser : argparse.ArgumentParser or None
        If given, errors are routed through ``parser.error()``.

    Returns
    -------
    FieldDistribution
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


def describe_fixed_gap(fd) -> Optional[str]:
    """
    The gap a fixed-size geometry imposes, in words, for the commands.

    For a coaxial arrangement or a tabulated field line, ``"d = L = 20 mm
    (the arc length of line.csv)"``, or with a protrusion ``"d = L - h =
    20 - 1 = 19 mm (... less the protrusion height, measured from its
    tip)"``; None where the gap is free.

    Parameters
    ----------
    fd : FieldDistribution

    Returns
    -------
    str or None
    """
    if fd.fixed_gap_length is None:
        return None
    h_mm = fd.protrusion.height * 1e3 if fd.protrusion else 0.0
    L_mm = fd.fixed_gap_length * 1e3
    what = (
        f"the arc length of {os.path.basename(fd.fieldline_path)}"
        if fd.field_type == "fieldline"
        else "b - a"
    )
    if h_mm:
        return (
            f"d = L - h = {L_mm + h_mm:.4g} - {h_mm:.4g} = {L_mm:.4g} mm "
            f"({what} less the protrusion height, measured from its tip)"
        )
    return f"d = L = {L_mm:.4g} mm ({what})"


def fixed_gap_warning(fd, d_mm) -> Optional[str]:
    """
    Say how gaps *d_mm* (mm) other than the one a fixed geometry imposes are
    solved, or None if they all match it (or the gap is free).

    Such a d solves the arrangement scaled so that its gap is d: the whole
    arrangement by d/L, or, with a protrusion, the electrodes by (d + h)/L
    and the protrusion at its given size.

    Parameters
    ----------
    fd : FieldDistribution
    d_mm : list of float
        Gap lengths in mm.

    Returns
    -------
    str or None
        One sentence for stderr, or None.
    """
    if fd.fixed_gap_length is None:
        return None
    L_mm = fd.fixed_gap_length * 1e3
    off = [d for d in d_mm if abs(d / L_mm - 1.0) > 1e-6]
    if not off:
        return None
    scaled = (
        "the electrodes are solved scaled by (d + h)/L, and the protrusion at "
        "its given size"
        if fd.protrusion
        else "the arrangement is solved scaled by d/L"
    )
    return (
        f"--d {', '.join(f'{d:g}' for d in off)} mm differs from "
        f"{describe_fixed_gap(fd)}; {scaled}.  Use d = {L_mm:.4g} mm for the "
        f"geometry as given."
    )


def field_header_lines(fd) -> list:
    """
    The geometry lines of a ``--write-to-file`` header, without the leading ``#``.

    Protrusion, reversal, electrode dimensions, the polarity convention, and
    for a tabulated line its source, arc length, line integral and the
    applied voltage, so that a results file says what it was computed for.
    Shared by the commands that write results.

    Parameters
    ----------
    fd : FieldDistribution

    Returns
    -------
    list of str
        One header line each, in the ``Key:  value`` layout of the headers.
    """
    lines = []
    if fd.protrusion is not None:
        lines.append(
            f"Protrusion:  {fd.protrusion.label}, on the xi = 0 electrode; "
            f"beta on a uniform field = {fd.protrusion.beta:.6g}; d is measured "
            f"from its tip"
        )
    if fd.reversed:
        lines.append(
            f"Reversed:    xi = 0 is the {fd.electrode_name} (--reverse-field)"
        )
    if fd.sphere_R is not None:
        lines.append(f"Sphere R:    {fd.sphere_R * 1e3:.4g} mm")
    if fd.field_type == "hyperboloid-plane":
        lines.append(f"Tip R:       {fd.tip_R * 1e3:.4g} mm")
    if fd.field_type == "coaxial":
        lines.append(
            f"Radii:       a = {fd.coax_a * 1e3:.6g} mm, b = {fd.coax_b * 1e3:.6g} mm"
        )
    if not fd.is_symmetric:
        el = fd.electrode_name
        what = {
            "inner": "inner conductor",
            "outer": "outer conductor",
            "start": "first tabulated point",
            "end": "last tabulated point",
            "protrusion": "electrode with the protrusion",
        }.get(el, el)
        lines.append(
            f"Polarity:    {el}=positive → {what} is anode (+),  "
            f"{el}=negative → {what} is cathode (−)"
        )
    if fd.field_type == "fieldline":
        lines.append(f"Field line:  {fd.fieldline_path}")
        lines.append(f"Arc length:  {fd.fieldline_length * 1e3:.6g} mm")
        u_line = fd.fieldline_line_voltage
        if u_line is not None:
            lines.append(
                f"∫|E| ds:     {u_line / 1e3:.6g} kV (|E| in {fd.fieldline_field_unit})"
            )
        else:
            lines.append(
                f"∫|E| ds:     {fd.fieldline_integral:.6g} (field units of the file "
                f"× m; a voltage only if the file tabulates |E| in V/m)"
            )
        if fd.fieldline_applied_voltage is not None:
            lines.append(
                f"U_applied:   {fd.fieldline_applied_voltage / 1e3:.6g} kV "
                f"({fd.fieldline_voltage_source}; only the shape of the profile "
                f"enters the solve, so U*/U_applied is the factor the excitation "
                f"must be scaled by to reach inception)"
            )
    return lines


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
    # Both the height and the footprint must be small against the scale:
    # a flat bump as wide as the sphere is no local perturbation either.
    base = fd.protrusion.base_radius
    size, size_what = (h, "height h") if h >= base else (base, "base radius")
    if scale is not None and size > PROTRUSION_SCALE_LIMIT * scale:
        notes.append(
            f"protrusion {size_what} = {size * 1e3:.4g} mm is "
            f"{size / scale:.3g} of {what} = {scale * 1e3:.4g} mm; the "
            f"enhancement is superposed on the background field, which is "
            f"accurate only while the protrusion is below about "
            f"{PROTRUSION_SCALE_LIMIT:g} of it"
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

    Parameters
    ----------
    fd : FieldDistribution
    d_min : float or None
        Smallest gap length in metres, for the gap note.
    gap_only : bool
        Print only the gap note, and nothing where the geometry fixes the
        gap (it was printed when the options were read).
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
        if len(spec) != 1:
            _err(f"--field uniform takes no values, got {spec[1:]}")
        return FieldDistribution(field_type="uniform")

    if field_type in ("sphere-plane", "sphere-sphere"):
        if len(spec) != 2:
            _err(f"--field {field_type} requires a sphere radius in mm")
        try:
            sphere_R = float(spec[1]) * 1e-3  # mm → m
        except ValueError:
            _err(f"--field {field_type}: radius must be a number, got {spec[1]!r}")
        if not (sphere_R > 0.0 and math.isfinite(sphere_R)):
            _err(f"--field {field_type}: radius must be positive")
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
