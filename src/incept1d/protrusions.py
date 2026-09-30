# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Protrusions on an electrode: the on-axis field enhancement E(z)/E₀.

A protrusion is an axisymmetric conductor standing on a grounded plane in a
uniform field E₀.  What the 1-D model needs from it is the enhancement
g(s) = E/E₀ on the axis at distance s above its tip;
:meth:`incept1d.fields.FieldDistribution.build` multiplies the background
profile by it and normalises again.  Three shapes are provided:

* :class:`Spheroid` — a half-spheroid, in closed form;
* :class:`Cone` — a cone with a spherically rounded tip;
* :class:`Rod` — a cylinder with a hemispherical cap.

The cone and the rod have no closed form and are solved once, on
construction, with the charge simulation method (:func:`_csm_solve`): ring
charges inside the surface, mirrored in the plane, whose strengths make the
surface an equipotential at a set of collocation points.
"""

import dataclasses
import math

import numpy as np
from scipy.special import ellipkm1

# ── Spheroid (closed form) ───────────────────────────────────────────────────

_SERIES_TERMS = 60  # |ratio| ≤ 1/4 below, so the tail is < 4⁻⁶⁰


def _spheroid_series(r):
    """Σ_{k≥1} r^{k−1}/(2k+1) and Σ_{k≥1} 2k r^{k−1}/(2k+1), for abs(r) ≤ 1/4."""
    s0 = s1 = 0.0
    rk = 1.0
    for k in range(1, _SERIES_TERMS + 1):
        s0 += rk / (2 * k + 1)
        s1 += 2 * k * rk / (2 * k + 1)
        rk *= r
    return s0, s1


def _spheroid_axial_field(s, h, b):
    """
    Exact on-axis field above a conducting half-spheroid on a plane, E/E₀.

    The half-spheroid has height *h* along the field and base radius *b*,
    and stands on a grounded plane in the uniform field E₀; with its
    mirror image it is a whole spheroid in a uniform field.  In spheroidal
    coordinates with focal distance c = √|h² − b²| the perturbation is the
    exterior harmonic Q₁, which on the axis at height z = h + s gives

        E(z)/E₀ = 1 + N(z)/D,   τ = (h² − b²)/h²,   w = h/z,
        D = Σ_{k≥1} τ^{k−1}/(2k+1),
        N = w³ Σ_{k≥1} 2k/(2k+1) (τw²)^{k−1},

    one expression for prolate (τ > 0), spherical (τ = 0: 1 + 2h³/z³) and
    oblate (τ < 0) shapes.  The series are summed where they converge fast
    (abs(τ) ≤ 1/4, abs(τ)w² ≤ 1/4) and replaced by their closed forms elsewhere,
    with Q₁(x) = x artanh(1/x) − 1 (prolate, x = z/c) or
    q(x) = 1 − x arccot x (oblate, x = z/c):

        prolate:  D = Q₁(h/c)/τ,     N = −Q₁′(z/c)/τ^{3/2},
        oblate:   D = q(h/c)/abs(τ),    N = −q′(z/c)/abs(τ)^{3/2}.

    At the tip (s = 0) this is the field-enhancement factor
    β = 1/((ξ₀² − 1)(ξ₀ artanh(1/ξ₀) − 1)), ξ₀ = h/c, of a prolate
    spheroid, 3 for a hemisphere, and → 1 for a flat disc.

    Parameters
    ----------
    s : float
        Distance above the tip, ≥ 0, in the units of *h* and *b*.
    h, b : float
        Height and base radius, > 0.
    """
    z = h + s
    tau = (h - b) * (h + b) / (h * h)
    w = h / z
    if abs(tau) <= 0.25:
        D, _ = _spheroid_series(tau)
    else:
        x0 = 1.0 / math.sqrt(abs(tau))
        if tau > 0.0:
            D = (x0 * math.atanh(1.0 / x0) - 1.0) / tau
        else:
            D = (1.0 - x0 * math.atan2(1.0, x0)) / -tau
    r = tau * w * w
    if abs(r) <= 0.25:
        _, S1 = _spheroid_series(r)
        N = w**3 * S1
    elif tau > 0.0:
        # x − 1 from s and h − c = b²/(h + c), which stays accurate at the
        # tip of a needle where x → 1.
        c = math.sqrt((h - b) * (h + b))
        xm1 = (s + b * b / (h + c)) / c
        x = 1.0 + xm1
        dQ = 0.5 * math.log1p(2.0 / xm1) - x / (xm1 * (x + 1.0))
        N = -dQ / tau**1.5
    else:
        c = math.sqrt((b - h) * (b + h))
        x = z / c
        dq = -math.atan2(1.0, x) + x / (1.0 + x * x)
        N = -dq / (-tau) ** 1.5
    return 1.0 + N / D


def spheroid_enhancement(h, b):
    """
    Field-enhancement factor β at the tip of a half-spheroid on a plane.

    β = E(tip)/E₀ for height *h* and base radius *b* (any common unit); see
    :func:`_spheroid_axial_field`.  3 for a hemisphere,
    ≈ (h/b)²/(ln(2h/b) − 1) for a needle (h ≫ b), → 1 for a flat disc.
    """
    return _spheroid_axial_field(0.0, h, b)


# ── Charge simulation method ─────────────────────────────────────────────────

#: Collocation spacing as a fraction of the local length scale, and the
#: largest ratio between neighbouring spacings.
_CSM_SPACING = 0.3
_CSM_GROWTH = 1.15
#: Inward offset of each ring charge from its collocation point, in units of
#: the local spacing (the "assignment factor").
_CSM_ASSIGNMENT = 1.0
#: Largest surface potential the solution may leave, relative to the applied
#: potential E₀z it cancels there (E₀r_tip near the base).
_CSM_RESIDUAL_LIMIT = 2e-3
#: Largest number of collocation points (the matrix is dense).
_CSM_MAX_POINTS = 4000


def _ring_potential(rho, z, a, z0):
    """
    Potential of a ring of unit charge, radius *a* at height *z0*, at (ρ, z).

    φ = (2/π) K(m) / √((ρ + a)² + (z − z0)²), m = 4aρ / ((ρ + a)² + (z − z0)²),
    in units where 4πε₀ = 1, and 1/distance for a point charge (a = 0).
    K is evaluated through the complementary parameter 1 − m, which stays
    accurate next to the ring where m → 1.
    """
    dz = z - z0
    den = (rho + a) ** 2 + dz**2
    return (2.0 / np.pi) * ellipkm1(((rho - a) ** 2 + dz**2) / den) / np.sqrt(den)


def _imaged_potential(rho, z, a, z0):
    """:func:`_ring_potential` of a ring and its negative image in z = 0."""
    return _ring_potential(rho, z, a, z0) - _ring_potential(rho, z, a, -z0)


class _Meridian:
    """
    The meridian of an axisymmetric protrusion, from its apex down to z = 0.

    A chain of circular arcs and straight segments, each given exactly, so
    that a point and its outward normal are known at any arc length.

    Parameters
    ----------
    pieces : list of tuple
        ``("arc", zc, R, α0, α1)``: the circle of radius R centred on the
        axis at height zc, from polar angle α0 to α1 (from +z);
        ``("line", (ρ0, z0), (ρ1, z1))``: a straight segment.
    """

    def __init__(self, pieces):
        self.pieces = [p for p in pieces if self._length(p) > 0.0]
        self.lengths = np.array([self._length(p) for p in self.pieces])
        self.starts = np.concatenate(([0.0], np.cumsum(self.lengths)[:-1]))
        self.total = float(self.lengths.sum())

    @staticmethod
    def _length(p):
        if p[0] == "arc":
            return p[2] * (p[4] - p[3])
        (r0, z0), (r1, z1) = p[1], p[2]
        return math.hypot(r1 - r0, z1 - z0)

    def at(self, s):
        """(ρ, z, n_ρ, n_z) at arc length *s* from the apex."""
        i = min(
            int(np.searchsorted(self.starts, s, side="right")) - 1, len(self.pieces) - 1
        )
        p, u = self.pieces[i], s - self.starts[i]
        if p[0] == "arc":
            _, zc, R, a0, _ = p
            al = a0 + u / R
            return R * math.sin(al), zc + R * math.cos(al), math.sin(al), math.cos(al)
        (r0, z0), (r1, z1) = p[1], p[2]
        L = self.lengths[i]
        tr, tz = (r1 - r0) / L, (z1 - z0) / L
        return r0 + u * tr, z0 + u * tz, -tz, tr


def _collocation(meridian, r_tip):
    """
    Arc lengths of the collocation points, apex first.

    The spacing follows the local length scale: a quarter of the tip radius
    at the apex, then the distance from the apex or the radius of the
    body, whichever is smaller, the height above the plane near the base,
    where the protrusion meets its image, and the distance to a joint
    between pieces, where the curvature jumps.  It grows by at most
    ``_CSM_GROWTH`` from one point to the next.
    """
    rt = 0.25 * r_tip
    out, dx_prev = [0.0], _CSM_SPACING * rt
    while True:
        x = out[-1]
        rho, z, _, _ = meridian.at(x)
        dx = _CSM_SPACING * max(rt, min(max(rho, rt), x + rt))
        joint = min((abs(x - j) for j in meridian.starts[1:]), default=math.inf)
        dx = min(
            dx,
            _CSM_GROWTH * dx_prev,
            _CSM_SPACING * max(z, 1e-3 * rt),
            _CSM_SPACING * max(joint, rt),
        )
        dx_prev = dx
        if x + dx >= meridian.total - 0.5 * dx:
            break
        out.append(x + dx)
        if len(out) > _CSM_MAX_POINTS:
            raise ValueError(
                f"the protrusion needs more than {_CSM_MAX_POINTS} charges; it is "
                f"too slender for the charge simulation"
            )
    return np.array(out)


def _csm_solve(meridian, r_tip):
    """
    Charges that make the protrusion an equipotential at zero in E₀ = 1.

    One ring charge per collocation point, set inward along the normal by
    ``_CSM_ASSIGNMENT`` local spacings (a point charge on the axis at the
    apex), each with its negative image in the plane z = 0.  Their strengths
    solve Σ_j q_j Φ_j(ρ_i, z_i) = z_i, which cancels the applied potential
    −z at every collocation point.

    The solution is checked halfway between collocation points, where it is
    least accurate: a surface potential above ``_CSM_RESIDUAL_LIMIT`` times
    the applied potential it should cancel there is an error rather than a
    result.

    Returns
    -------
    a, z0, q : ndarray
        Radius, height and charge of each ring (4πε₀ = 1, E₀ = 1).
    """
    s = _collocation(meridian, r_tip)
    pts = np.array([meridian.at(x) for x in s])
    ds = np.diff(np.append(s, meridian.total))
    ds = np.maximum(ds, np.concatenate(([ds[0]], ds[:-1])))
    rho, z, nr, nz = pts.T
    a = np.maximum(rho - _CSM_ASSIGNMENT * ds * nr, 0.0)
    a[0] = 0.0
    z0 = z - _CSM_ASSIGNMENT * ds * nz
    P = _imaged_potential(rho[:, None], z[:, None], a[None, :], z0[None, :])
    q = np.linalg.solve(P, z)

    mid = np.array([meridian.at(x) for x in 0.5 * (s[1:] + s[:-1])])
    phi = (
        -mid[:, 1]
        + _imaged_potential(
            mid[:, 0][:, None], mid[:, 1][:, None], a[None, :], z0[None, :]
        )
        @ q
    )
    residual = float(np.max(np.abs(phi) / np.maximum(mid[:, 1], r_tip)))
    if not residual < _CSM_RESIDUAL_LIMIT:
        raise ValueError(
            f"the charge simulation leaves a surface potential of {residual:.2g} "
            f"of the applied one, above {_CSM_RESIDUAL_LIMIT:g}; the shape is "
            f"outside the range it resolves"
        )
    return a, z0, q


def _csm_axial_field(s, height, a, z0, q):
    """E/E₀ on the axis at *s* above the tip, from the charges of :func:`_csm_solve`."""
    z = height + s
    dp, dm = z - z0, z + z0
    a2 = a * a
    return 1.0 + float(q @ (dp / (a2 + dp * dp) ** 1.5 - dm / (a2 + dm * dm) ** 1.5))


# ── The shapes ───────────────────────────────────────────────────────────────


@dataclasses.dataclass(eq=False)
class Protrusion:
    """
    An axisymmetric protrusion on an electrode.

    Attributes
    ----------
    height : float
        Height h of the tip above the electrode surface, in metres.
    """

    height: float

    shape = ""

    def enhancement(self, s: float) -> float:
        """E/E₀ on the axis at distance *s* (m) above the tip."""
        raise NotImplementedError

    @property
    def tip_radius(self) -> float:
        """Radius of curvature of the tip, in metres."""
        raise NotImplementedError

    @property
    def beta(self) -> float:
        """Field-enhancement factor at the tip, E(tip)/E₀."""
        return self.enhancement(0.0)

    @property
    def label(self) -> str:
        """Human-readable description for labels and file headers."""
        raise NotImplementedError


@dataclasses.dataclass(eq=False)
class Spheroid(Protrusion):
    """
    Half-spheroid of height h and base radius b: a needle for h > b, a
    hemisphere for h = b, a flat bump for h < b.  Closed form,
    :func:`_spheroid_axial_field`.
    """

    base_radius: float = 0.0
    shape = "spheroid"

    def enhancement(self, s):
        return _spheroid_axial_field(s, self.height, self.base_radius)

    @property
    def tip_radius(self):
        return self.base_radius**2 / self.height

    @property
    def label(self):
        return (
            f"spheroid protrusion h = {self.height * 1e3:.4g} mm, "
            f"b = {self.base_radius * 1e3:.4g} mm"
        )


@dataclasses.dataclass(eq=False)
class _SimulatedProtrusion(Protrusion):
    """A protrusion solved by the charge simulation method on construction."""

    def __post_init__(self):
        # Solved in units of the tip radius, so that the conditioning does
        # not depend on the absolute size.
        r = self.tip_radius
        self._r = r
        self._h = self.height / r
        self._a, self._z0, self._q = _csm_solve(self._meridian(), 1.0)

    def _meridian(self) -> _Meridian:
        """The meridian in units of the tip radius."""
        raise NotImplementedError

    def enhancement(self, s):
        return _csm_axial_field(s / self._r, self._h, self._a, self._z0, self._q)


@dataclasses.dataclass(eq=False)
class Rod(_SimulatedProtrusion):
    """
    Cylinder of radius R with a hemispherical cap, total height h ≥ R: a
    whisker, fibre or wire end.  Charge simulation, :func:`_csm_solve`.
    """

    radius: float = 0.0
    shape = "rod"

    @property
    def tip_radius(self):
        return self.radius

    def _meridian(self):
        h = self.height / self.radius
        return _Meridian(
            [
                ("arc", h - 1.0, 1.0, 0.0, math.pi / 2),
                ("line", (1.0, h - 1.0), (1.0, 0.0)),
            ]
        )

    @property
    def label(self):
        return (
            f"rod protrusion h = {self.height * 1e3:.4g} mm, "
            f"R = {self.radius * 1e3:.4g} mm"
        )


@dataclasses.dataclass(eq=False)
class Cone(_SimulatedProtrusion):
    """
    Cone of half-angle θ whose apex is rounded by a sphere of radius r
    tangent to it, total height h: a burr or a sharp asperity.  Charge
    simulation, :func:`_csm_solve`.
    """

    radius: float = 0.0
    half_angle: float = 0.0  # degrees
    shape = "cone"

    @property
    def tip_radius(self):
        return self.radius

    @property
    def base_radius(self) -> float:
        """Radius of the cone where it meets the electrode, in metres."""
        th = math.radians(self.half_angle)
        return (self.height - self.radius + self.radius / math.sin(th)) * math.tan(th)

    def _meridian(self):
        th = math.radians(self.half_angle)
        h = self.height / self.radius
        zc = h - 1.0
        rho1, z1 = math.cos(th), zc + math.sin(th)
        return _Meridian(
            [
                ("arc", zc, 1.0, 0.0, math.pi / 2 - th),
                ("line", (rho1, z1), (self.base_radius / self.radius, 0.0)),
            ]
        )

    @property
    def label(self):
        return (
            f"cone protrusion h = {self.height * 1e3:.4g} mm, "
            f"r = {self.radius * 1e3:.4g} mm, {self.half_angle:.4g} deg"
        )


#: Shape keyword → (class, argument names, help) for ``--protrusion``.
SHAPES = {
    "spheroid": (Spheroid, ("H_mm", "B_mm"), "height and base radius"),
    "cone": (Cone, ("H_mm", "R_mm", "ANGLE_deg"), "height, tip radius, half-angle"),
    "rod": (Rod, ("H_mm", "R_mm"), "height and radius"),
}

#: Limits of the range the shapes are validated over.
ROD_MAX_ASPECT = 1000.0
CONE_MAX_ASPECT = 1e5
CONE_ANGLE_RANGE = (5.0, 80.0)


def parse_protrusion_spec(tokens) -> Protrusion:
    """
    Build a :class:`Protrusion` from ``--protrusion`` tokens.

    ``spheroid H B``, ``cone H R ANGLE`` or ``rod H R``, lengths in mm and
    the cone half-angle in degrees.

    Raises
    ------
    ValueError
        On an unknown shape, a wrong number of values, or a shape outside
        the validated range.
    """
    shape = tokens[0]
    if shape not in SHAPES:
        raise ValueError(
            f"unknown shape {shape!r}; use one of {', '.join(sorted(SHAPES))}"
        )
    cls, names, what = SHAPES[shape]
    if len(tokens) != len(names) + 1:
        raise ValueError(f"{shape} requires {what} ({' '.join(names)})")
    try:
        vals = [float(t) for t in tokens[1:]]
    except ValueError:
        raise ValueError(f"{shape}: values must be numbers, got {tokens[1:]}")
    if not all(v > 0.0 for v in vals):
        raise ValueError(f"{shape}: {what} must be positive")
    h, r = vals[0] * 1e-3, vals[1] * 1e-3
    if shape == "spheroid":
        return Spheroid(h, base_radius=r)
    if shape == "rod":
        if h < r:
            raise ValueError("rod: the height must be at least the radius")
        if h / r > ROD_MAX_ASPECT:
            raise ValueError(f"rod: h/R above {ROD_MAX_ASPECT:g} is not resolved")
        return Rod(h, radius=r)
    th = vals[2]
    lo, hi = CONE_ANGLE_RANGE
    if not lo <= th <= hi:
        raise ValueError(f"cone: half-angle must be between {lo:g} and {hi:g} deg")
    if h <= r * (1.0 - math.sin(math.radians(th))):
        raise ValueError("cone: the height does not reach past the rounded tip")
    if h / r > CONE_MAX_ASPECT:
        raise ValueError(f"cone: h/r above {CONE_MAX_ASPECT:g} is not resolved")
    return Cone(h, radius=r, half_angle=th)
