# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Ionization integrals across the gap."""

import numpy as np
import pytest

from incept1d.fields import FieldDistribution
from incept1d.ionization import aed_integral, eig_integral

_trapz = getattr(np, "trapezoid", None) or np.trapz


class TestUniform:
    def test_matches_the_closed_form(self, toy, toy_raw):
        """X1: a uniform field makes the integral a product."""
        EN, p, T, d = 400.0, 1.0, 293.0, 10e-3
        expected = max(0.0, toy_raw.alpha(EN, p, T) - toy_raw.eta(EN, p, T)) * d
        got = aed_integral(EN, p, d, toy, T, FieldDistribution("uniform"), 200)
        assert got == pytest.approx(expected, rel=1e-12)

    def test_is_clamped_at_zero(self, toy, toy_raw):
        """Below the alpha = eta crossover there is no net ionization."""
        EN, p, T, d = 30.0, 1.0, 293.0, 10e-3
        assert toy_raw.alpha(EN, p, T) < toy_raw.eta(EN, p, T)
        assert aed_integral(EN, p, d, toy, T, FieldDistribution("uniform"), 200) == 0.0


class TestNonUniform:
    def test_monotone_profile_matches_fine_quadrature(self, toy, toy_raw):
        """X2: sphere-plane, integrated over the active region only."""
        EN, p, T, d = 400.0, 1.0, 293.0, 10e-3
        fd = FieldDistribution("sphere-plane", 50e-3)
        f = fd.build(d)
        xi = np.linspace(0.0, 1.0, 40001)
        net = np.array(
            [toy_raw.alpha(EN * f(x), p, T) - toy_raw.eta(EN * f(x), p, T) for x in xi]
        )
        expected = float(_trapz(np.maximum(net, 0.0), xi)) * d
        got = aed_integral(EN, p, d, toy, T, fd, 4000)
        assert got == pytest.approx(expected, rel=1e-3)

    def test_symmetric_profile_captures_both_active_regions(self, toy, toy_raw):
        """
        X3: a sphere-sphere gap is active near *both* electrodes.

        The quadrature must not assume a single region starting at xi = 0, as
        it would for a monotone profile.  Comparing against a fine reference
        integral over the whole gap is the check that both ends contribute.
        """
        # R = 5 mm, d = 20 mm puts the alpha = eta crossover (48.3 Td for the
        # toy defaults) between f(0.5) = 0.545 and f(0) = 2.68, so the gap is
        # active at both electrodes and dead in the middle.
        EN, p, T, d = 40.0, 1.0, 293.0, 20e-3
        fd = FieldDistribution("sphere-sphere", 5e-3)
        assert not fd.is_monotone_decreasing
        f = fd.build(d)
        xi = np.linspace(0.0, 1.0, 40001)
        net = np.array(
            [toy_raw.alpha(EN * f(x), p, T) - toy_raw.eta(EN * f(x), p, T) for x in xi]
        )
        active = net > 0.0
        # The configuration must genuinely have a dead zone in the middle,
        # otherwise the test would pass for the wrong reason.
        assert active[0] and active[-1] and not active[len(active) // 2]
        expected = float(_trapz(np.maximum(net, 0.0), xi)) * d
        got = aed_integral(EN, p, d, toy, T, fd, 4000)
        assert got == pytest.approx(expected, rel=1e-2)

    def test_fieldline_profile_is_supported(self, toy, tmp_path):
        """X3b: a tabulated line takes the same general quadrature path."""
        d = 10e-3
        s = np.linspace(0.0, d, 201)
        E = 1.0 + np.cos(2 * np.pi * s / d)  # two maxima, like sphere-sphere
        p = tmp_path / "line.dat"
        p.write_text("\n".join(f"{a} {b}" for a, b in zip(s, E + 0.1)))
        fd = FieldDistribution.from_fieldline(str(p))
        assert not fd.is_monotone_decreasing
        val = aed_integral(400.0, 1.0, d, toy, 293.0, fd, 2000)
        assert np.isfinite(val) and val > 0.0


class TestEigenvalueIntegral:
    def test_equals_the_alpha_eta_integral_without_detachment(self, toy_path):
        """X4: with delta = 0 the apparent coefficient is exactly alpha - eta."""
        from incept1d.mechanism import load_mechanism

        mod = load_mechanism(toy_path, {"delta_0": 0.0})
        EN, p, T, d = 400.0, 1.0, 293.0, 10e-3
        fd = FieldDistribution("uniform")
        assert eig_integral(EN, p, d, mod, T, fd, 200) == pytest.approx(
            aed_integral(EN, p, d, mod, T, fd, 200), rel=1e-9
        )

    def test_detachment_raises_the_apparent_integral(self, toy, toy_path):
        """With delta > 0, lambda_+ > alpha - eta, so the integral is larger."""
        from incept1d.mechanism import load_mechanism

        nodet = load_mechanism(toy_path, {"delta_0": 0.0})
        EN, p, T, d = 400.0, 1.0, 293.0, 10e-3
        fd = FieldDistribution("uniform")
        assert eig_integral(EN, p, d, toy, T, fd, 200) > eig_integral(
            EN, p, d, nodet, T, fd, 200
        )


class TestRobustness:
    @pytest.mark.parametrize("fn", [aed_integral, eig_integral])
    def test_extreme_input_returns_nan_not_an_exception(self, toy, fn):
        """X5: the CLI sweeps blindly, so overflow must not raise."""
        val = fn(1e30, 1.0, 1e6, toy, 293.0, FieldDistribution("uniform"), 50)
        assert np.isnan(val) or np.isfinite(val)
