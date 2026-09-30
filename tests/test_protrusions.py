# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Protrusions: the spheroid closed form and the charge-simulated shapes."""

import math

import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.special import lpmv

import incept1d.protrusions as prot
from incept1d.protrusions import Cone, Rod, parse_protrusion_spec


class TestSpheroidClosedForm:
    """
    The on-axis field of a half-spheroid on a plane, against closed forms.

    Independently of :func:`incept1d.protrusions._spheroid_axial_field`, which sums
    a series or evaluates Legendre functions depending on the shape, the
    potential of a grounded prolate spheroid in a uniform field is built here
    from the distances to its foci: ξ = (r₊ + r₋)/2c and η = (r₋ − r₊)/2c are
    prolate spheroidal coordinates, and φ = E₀(−cξη + h η Q₁(ξ)/Q₁(ξ₀)).
    """

    H, B = 2.0, 0.7  # prolate, arbitrary units

    @classmethod
    def _potential(cls, rho, z):
        c = np.sqrt(cls.H**2 - cls.B**2)
        rp, rm = np.hypot(rho, z - c), np.hypot(rho, z + c)
        xi, eta = (rp + rm) / (2 * c), (rm - rp) / (2 * c)

        def Q1(x):
            return x * np.arctanh(1.0 / x) - 1.0

        return -c * xi * eta + cls.H * eta * Q1(xi) / Q1(cls.H / c)  # E0 = 1

    def test_potential_is_harmonic(self):
        """F40: ∇²φ = 0 in cylindrical coordinates outside the spheroid."""
        e = 1e-4
        phi = self._potential
        for rho, z in [(0.5, 2.5), (1.5, 1.0), (3.0, 4.0)]:
            d2r = (phi(rho + e, z) - 2 * phi(rho, z) + phi(rho - e, z)) / e**2
            d1r = (phi(rho + e, z) - phi(rho - e, z)) / (2 * e)
            d2z = (phi(rho, z + e) - 2 * phi(rho, z) + phi(rho, z - e)) / e**2
            assert abs(d2r + d1r / rho + d2z) < 1e-5 * (abs(d2r) + abs(d2z) + 1)

    def test_spheroid_and_plane_are_grounded(self):
        """F41: φ = 0 on the surface ρ²/b² + z²/h² = 1 and on the plane z = 0."""
        for th in np.linspace(0.05, np.pi / 2, 7):
            rho, z = self.B * np.sin(th), self.H * np.cos(th)
            assert abs(self._potential(rho, z)) < 1e-9
        assert abs(self._potential(2.0, 0.0)) < 1e-12

    def test_axial_field_is_minus_the_gradient(self):
        """F42: the closed form is −∂φ/∂z on the axis."""
        from incept1d.protrusions import _spheroid_axial_field

        e = 1e-6
        for s in [0.01, 0.3, 1.0, 5.0, 50.0]:
            z = self.H + s
            dphi = (self._potential(0.0, z + e) - self._potential(0.0, z - e)) / (2 * e)
            assert _spheroid_axial_field(s, self.H, self.B) == pytest.approx(
                -dphi, rel=1e-6
            )

    @pytest.mark.parametrize("aspect", [1.05, 1.5, 3.0, 10.0, 1e3])
    def test_prolate_enhancement(self, aspect):
        """F43: β = 1/((ξ₀² − 1)(ξ₀ artanh(1/ξ₀) − 1)), ξ₀ = h/c."""
        from incept1d.protrusions import spheroid_enhancement

        x0 = aspect / np.sqrt(aspect**2 - 1.0)
        beta = 1.0 / ((x0**2 - 1.0) * (x0 * np.arctanh(1.0 / x0) - 1.0))
        assert spheroid_enhancement(aspect, 1.0) == pytest.approx(beta, rel=1e-7)

    @pytest.mark.parametrize("aspect", [0.01, 0.3, 0.7, 0.95])
    def test_oblate_enhancement(self, aspect):
        """F44: β = 1/((1 + ζ₀²)(1 − ζ₀ arccot ζ₀)), ζ₀ = h/c."""
        from incept1d.protrusions import spheroid_enhancement

        z0 = aspect / np.sqrt(1.0 - aspect**2)
        beta = 1.0 / ((1.0 + z0**2) * (1.0 - z0 * np.arctan(1.0 / z0)))
        assert spheroid_enhancement(aspect, 1.0) == pytest.approx(beta, rel=1e-9)

    def test_limits(self):
        """F45: hemisphere 3, flat disc → 1, needle → (h/b)²/(ln(2h/b) − 1)."""
        from incept1d.protrusions import spheroid_enhancement

        assert spheroid_enhancement(1.0, 1.0) == pytest.approx(3.0, rel=1e-14)
        assert spheroid_enhancement(1e-6, 1.0) == pytest.approx(1.0, abs=1e-5)
        a = 1e4
        assert spheroid_enhancement(a, 1.0) == pytest.approx(
            a**2 / (np.log(2 * a) - 1.0), rel=1e-3
        )

    @pytest.mark.parametrize("b", [0.5, np.sqrt(0.75), 1.0, np.sqrt(1.25), 2.0])
    def test_continuous_across_the_series_switch(self, b):
        """F46: the series and closed-form branches meet without a jump."""
        from incept1d.protrusions import _spheroid_axial_field

        for s in [0.0, 0.2, 1.0]:
            lo = _spheroid_axial_field(s, 1.0, b * (1 - 1e-10))
            hi = _spheroid_axial_field(s, 1.0, b * (1 + 1e-10))
            assert lo == pytest.approx(hi, rel=1e-8)


class TestChargeSimulation:
    """The cone and the rod, which have no closed form, against what is known."""

    def test_a_hemispherical_rod_is_the_hemisphere(self):
        """P1: a rod of height R is a hemisphere, E/E₀ = 1 + 2R³/z³."""
        R = 1e-3
        rod = Rod(R, radius=R)
        for s in [0.0, 0.1 * R, R, 10 * R]:
            z = R + s
            assert rod.enhancement(s) == pytest.approx(1 + 2 * R**3 / z**3, rel=1e-3)

    @pytest.mark.parametrize("aspect", [5.0, 10.0, 100.0, 300.0])
    def test_rod_enhancement_follows_the_published_fit(self, aspect):
        """
        P2: β = 1.2 (h/R + 2.15)^0.9 within its few-percent accuracy.

        The fit is [EdgcombeValdre2001]'s to their own finite-element
        solutions for a hemisphere-capped cylinder on a plane.
        """
        beta = Rod(aspect * 1e-6, radius=1e-6).beta
        assert beta == pytest.approx(1.2 * (aspect + 2.15) ** 0.9, rel=0.05)

    @pytest.mark.parametrize("theta", [10.0, 30.0, 50.0])
    def test_cone_field_follows_the_sharp_cone_power_law(self, theta):
        """
        P3: between the tip radius and the height, E ∝ s^(ν−1).

        Near the apex of a sharp grounded cone of half-angle θ the potential
        is R^ν P_ν(cos α), with ν the root of P_ν(−cos θ) = 0 that vanishes
        on the cone surface.
        """
        nu = brentq(lambda v: lpmv(0, v, -math.cos(math.radians(theta))), 1e-6, 1.0)
        r = 1e-6
        cone = Cone(1e5 * r, radius=r, half_angle=theta)
        g1, g2 = cone.enhancement(100 * r), cone.enhancement(1000 * r)
        assert math.log10(g2 / g1) == pytest.approx(nu - 1.0, abs=5e-3)

    @pytest.mark.parametrize(
        "make",
        [
            lambda: Rod(50e-6, radius=1e-6),
            lambda: Cone(10e-6, radius=1e-6, half_angle=10.0),
            lambda: Cone(10e-6, radius=1e-6, half_angle=60.0),
            lambda: Cone(100e-6, radius=1e-6, half_angle=80.0),
        ],
        ids=["rod", "cone10", "cone60", "cone80"],
    )
    def test_converged_in_the_spacing(self, make, monkeypatch):
        """P4: halving the collocation spacing moves the field by < 1e-3."""
        coarse = make()
        monkeypatch.setattr(prot, "_CSM_SPACING", 0.5 * prot._CSM_SPACING)
        fine = make()
        h = coarse.height
        for s in [0.0, 0.01 * h, 0.1 * h, h]:
            assert coarse.enhancement(s) == pytest.approx(fine.enhancement(s), rel=1e-3)

    def test_the_enhancement_is_scale_free(self):
        """P5: E/E₀ depends on shape alone, so s scales with the protrusion."""
        a = Cone(10e-6, radius=1e-6, half_angle=30.0)
        b = Cone(10e-3, radius=1e-3, half_angle=30.0)
        for x in [0.0, 0.5, 5.0]:
            assert a.enhancement(x * 1e-6) == pytest.approx(
                b.enhancement(x * 1e-3), rel=1e-12
            )

    def test_an_unresolved_shape_is_an_error(self, monkeypatch):
        """P6: a residual above the limit is refused, not returned."""
        monkeypatch.setattr(prot, "_CSM_RESIDUAL_LIMIT", 1e-12)
        with pytest.raises(ValueError, match="surface potential"):
            Rod(10e-6, radius=1e-6)

    def test_the_cone_meets_the_plane_at_its_base_radius(self):
        """P7: the meridian ends on z = 0, at the radius the class reports."""
        cone = Cone(10e-6, radius=1e-6, half_angle=30.0)
        m = cone._meridian()
        rho, z, _, _ = m.at(m.total)
        assert z == pytest.approx(0.0, abs=1e-12)
        assert rho * 1e-6 == pytest.approx(cone.base_radius, rel=1e-12)

    def test_sharper_is_stronger(self):
        """P8: β grows as the cone narrows and as the rod lengthens."""
        cones = [Cone(20e-6, radius=1e-6, half_angle=t).beta for t in (60, 30, 10)]
        rods = [Rod(h * 1e-6, radius=1e-6).beta for h in (2, 10, 50)]
        assert np.all(np.diff(cones) > 0) and np.all(np.diff(rods) > 0)


class TestParsing:
    def test_shapes_and_units(self):
        """P9: lengths in mm, the cone angle in degrees."""
        c = parse_protrusion_spec(["cone", "2", "0.01", "15"])
        assert (c.height, c.radius, c.half_angle) == pytest.approx((2e-3, 1e-5, 15))
        r = parse_protrusion_spec(["rod", "1", "0.05"])
        assert (r.height, r.radius) == pytest.approx((1e-3, 5e-5))

    @pytest.mark.parametrize(
        "spec, match",
        [
            (["rod", "2000", "1"], "not resolved"),
            (["cone", "0.4", "1", "30"], "rounded tip"),
            (["cone", "1", "0.1", "2"], "half-angle"),
        ],
    )
    def test_outside_the_validated_range(self, spec, match):
        """P10."""
        with pytest.raises(ValueError, match=match):
            parse_protrusion_spec(spec)


class TestSlenderSpheroid:
    @pytest.mark.parametrize("aspect", [1e4, 1e6, 1e8, 1e12])
    def test_needle_keeps_its_digits(self, aspect):
        """
        P11: β of a slender needle, from x0 − 1 = b²/(c(h + c)) and
        artanh(1/x0) = ½ log1p(2/(x0 − 1)); atanh(1/x0) lost every digit.
        """
        from incept1d.protrusions import spheroid_enhancement

        h, b = aspect, 1.0
        c = math.sqrt((h - b) * (h + b))
        xm1 = b * b / (c * (h + c))
        x0 = 1.0 + xm1
        beta = 1.0 / (xm1 * (x0 + 1.0) * (x0 * 0.5 * math.log1p(2.0 / xm1) - 1.0))
        assert spheroid_enhancement(h, b) == pytest.approx(beta, rel=1e-12)

    @pytest.mark.parametrize(
        "spec, match",
        [
            (["spheroid", "inf", "1"], "finite"),
            (["rod", "1", "nan"], "finite"),
            (["spheroid", "1e7", "1"], "validated range"),
        ],
    )
    def test_unphysical_sizes_are_refused(self, spec, match):
        """P12."""
        with pytest.raises(ValueError, match=match):
            parse_protrusion_spec(spec)
