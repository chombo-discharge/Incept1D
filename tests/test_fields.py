# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Gap geometry: normalisation, symmetry, and the tabulated field-line loader."""

import dataclasses

import numpy as np
import pytest

from incept1d.fields import (
    FieldDistribution,
    load_fieldline,
    parse_field_spec,
    parse_field_unit,
)
from incept1d.protrusions import Cone, Rod, Spheroid

_trapz = getattr(np, "trapezoid", None) or np.trapz


def _integral(f, n=20001):
    xi = np.linspace(0.0, 1.0, n)
    return float(_trapz([f(x) for x in xi], xi))


def _write_line(path, rows, header=None, sep=" "):
    with open(path, "w") as fh:
        if header:
            fh.write(header + "\n")
        for r in rows:
            fh.write(sep.join(f"{v:.12e}" for v in r) + "\n")
    return str(path)


@pytest.fixture
def sphere_plane_samples():
    """The analytic sphere-plane profile sampled along a straight line."""
    d, R = 20e-3, 50e-3
    f = FieldDistribution("sphere-plane", R).build(d)
    xi = np.linspace(0.0, 1.0, 401)
    E = 1.5e6 * np.array([f(x) for x in xi])
    return d, f, xi, E


class TestNormalisation:
    @pytest.mark.parametrize(
        "fd",
        [
            FieldDistribution("uniform"),
            FieldDistribution("sphere-plane", 50e-3),
            FieldDistribution("sphere-plane", 5e-3),
            FieldDistribution("sphere-sphere", 50e-3),
            FieldDistribution("sphere-sphere", 5e-3),
            FieldDistribution("coaxial", coax_a=1e-3, coax_b=10e-3),
            FieldDistribution("coaxial", coax_a=0.5e-3, coax_b=25e-3),
            FieldDistribution("hyperboloid-plane", tip_R=20e-3),
            FieldDistribution("hyperboloid-plane", tip_R=5e-3),
        ],
    )
    @pytest.mark.parametrize("d", [1e-3, 20e-3, 100e-3])
    def test_profile_integrates_to_one(self, fd, d):
        """F1: every geometry satisfies the normalisation the solver assumes."""
        assert _integral(fd.build(d)) == pytest.approx(1.0, rel=1e-6)

    def test_uniform_is_exactly_one(self):
        """F2."""
        f = FieldDistribution("uniform").build(10e-3)
        assert all(f(x) == 1.0 for x in np.linspace(0.0, 1.0, 11))


class TestShape:
    def test_sphere_plane_decreases_monotonically(self):
        """F3: the field falls from the sphere surface to the plane."""
        f = FieldDistribution("sphere-plane", 50e-3).build(20e-3)
        v = np.array([f(x) for x in np.linspace(0.0, 1.0, 200)])
        assert np.all(np.diff(v) < 0.0)
        assert v[0] > 1.0 > v[-1]

    def test_sphere_sphere_is_symmetric_about_midgap(self):
        """F3b: f(xi) = f(1 - xi)."""
        f = FieldDistribution("sphere-sphere", 50e-3).build(20e-3)
        for x in np.linspace(0.0, 0.5, 25):
            assert f(x) == pytest.approx(f(1.0 - x), rel=1e-10)

    @pytest.mark.parametrize(
        "fd, symmetric, monotone",
        [
            (FieldDistribution("uniform"), True, True),
            (FieldDistribution("sphere-plane", 50e-3), False, True),
            (FieldDistribution("sphere-sphere", 50e-3), True, False),
            (FieldDistribution("coaxial", coax_a=1e-3, coax_b=10e-3), False, True),
            (FieldDistribution("hyperboloid-plane", tip_R=1e-3), False, True),
        ],
    )
    def test_geometry_flags(self, fd, symmetric, monotone):
        """F4: the flags the solver and the quadrature branch on."""
        assert fd.is_symmetric is symmetric
        assert fd.is_monotone_decreasing is monotone

    def test_small_gap_tends_to_uniform(self):
        """F5: d/R -> 0 recovers the uniform field."""
        f = FieldDistribution("sphere-plane", 1.0).build(1e-4)
        v = np.array([f(x) for x in np.linspace(0.0, 1.0, 50)])
        assert np.allclose(v, 1.0, atol=1e-3)


class TestHyperboloidPlane:
    """The profile is the axis of a Laplace solution with the right electrodes.

    Independently of the closed form in :mod:`incept1d.fields`, build the
    potential from the distances to the two foci: with foci at z = ±a,
    η = (r₋ − r₊)/(2a) is a prolate spheroidal coordinate, its level sets are
    confocal hyperboloids (η = 0 the plane z = 0), and artanh η is harmonic.
    """

    D, R = 10e-3, 0.5e-3

    @classmethod
    def _potential(cls, rho, z):
        a = np.sqrt(cls.D * (cls.D + cls.R))
        eta = (np.hypot(rho, z + a) - np.hypot(rho, z - a)) / (2.0 * a)
        return np.arctanh(eta) / np.arctanh(cls.D / a)  # U = 1 on the tip

    def test_potential_is_harmonic(self):
        """F30: ∇²φ = 0 in cylindrical coordinates, off axis in the gap."""
        h = 1e-6
        for rho, z in [(0.3e-3, 5e-3), (2e-3, 8e-3), (4e-3, 1e-3)]:
            phi = self._potential
            d2r = (phi(rho + h, z) - 2 * phi(rho, z) + phi(rho - h, z)) / h**2
            d1r = (phi(rho + h, z) - phi(rho - h, z)) / (2 * h)
            d2z = (phi(rho, z + h) - 2 * phi(rho, z) + phi(rho, z - h)) / h**2
            scale = abs(d2z) + abs(d2r)
            assert abs(d2r + d1r / rho + d2z) < 1e-4 * scale

    def test_electrodes_are_the_plane_and_a_tip_of_radius_r(self):
        """F31: φ = 0 on z = 0, φ = 1 at z = d, vertex curvature radius r."""
        assert np.allclose(self._potential(np.linspace(0, 0.05, 11), 0.0), 0.0)
        assert self._potential(0.0, self.D) == pytest.approx(1.0)
        # The tip surface φ = 1 near the axis: z(ρ) ≈ d − ρ²/(2r).
        rho = 1e-6
        a = np.sqrt(self.D * (self.D + self.R))
        eta0 = self.D / a
        # On the hyperboloid z²/eta0² − ρ²/(1 − eta0²) = a², solved for z.
        z = eta0 * np.sqrt(a**2 + rho**2 / (1.0 - eta0**2))
        assert self._potential(rho, z) == pytest.approx(1.0, abs=1e-12)
        assert rho**2 / (2.0 * (z - self.D)) == pytest.approx(self.R, rel=1e-4)

    def test_profile_is_the_axial_field(self):
        """F32: f(ξ) = −∂φ/∂z · d/U at z = d(1 − ξ)."""
        f = FieldDistribution("hyperboloid-plane", tip_R=self.R).build(self.D)
        h = 1e-8
        for xi in np.linspace(0.0, 1.0, 11):
            z = self.D * (1.0 - xi)
            dphi = (self._potential(0.0, z + h) - self._potential(0.0, z - h)) / (2 * h)
            assert f(xi) == pytest.approx(dphi * self.D, rel=1e-6)

    @pytest.mark.parametrize("r", [1e-3, 1e-4, 1e-5])
    def test_sharp_tip_integrates_to_one(self, r):
        """F33: normalisation where the trapezoid rule of F1 cannot resolve the tip."""
        from scipy.integrate import quad

        d = 0.1
        f = FieldDistribution("hyperboloid-plane", tip_R=r).build(d)
        val = quad(f, 0.0, 1.0, points=[r / d], limit=200, epsabs=0, epsrel=1e-12)
        assert val[0] == pytest.approx(1.0, rel=1e-10)

    def test_sharp_tip_recovers_the_classic_field(self):
        """F34: E_tip → 2U / (r ln(4d/r)) as r/d → 0."""
        d, r = 10e-3, 1e-6
        f = FieldDistribution("hyperboloid-plane", tip_R=r).build(d)
        assert f(0.0) / d == pytest.approx(2.0 / (r * np.log(4 * d / r)), rel=1e-3)

    def test_blunt_tip_tends_to_uniform(self):
        """F35: r/d → ∞ recovers the plane-parallel gap."""
        f = FieldDistribution("hyperboloid-plane", tip_R=1e3).build(1e-3)
        assert np.allclose([f(x) for x in np.linspace(0, 1, 11)], 1.0, atol=1e-5)
        f = FieldDistribution("hyperboloid-plane", tip_R=1e12).build(1e-3)
        assert f(0.3) == 1.0

    def test_profile_depends_on_d_over_r_only(self):
        """F36: scaling tip radius and gap together leaves f unchanged."""
        f = FieldDistribution("hyperboloid-plane", tip_R=self.R).build(self.D)
        g = FieldDistribution("hyperboloid-plane", tip_R=7 * self.R).build(7 * self.D)
        for x in np.linspace(0.0, 1.0, 11):
            assert f(x) == pytest.approx(g(x), rel=1e-12)


class TestCoaxial:
    A, B = 1e-3, 10e-3

    def _f(self, a=A, b=B):
        return FieldDistribution("coaxial", coax_a=a, coax_b=b).build(b - a)

    def test_field_falls_as_one_over_r(self):
        """F12: Gauss's law, E·r is constant across the gap."""
        f = self._f()
        r = lambda x: self.A + x * (self.B - self.A)  # noqa: E731
        er = [f(x) * r(x) for x in np.linspace(0.0, 1.0, 41)]
        assert np.allclose(er, er[0], rtol=1e-12)

    def test_surface_field_matches_the_exact_solution(self):
        """F13: E(a) = U / (a ln(b/a)), i.e. f(0) = (b − a) / (a ln(b/a))."""
        a, b = self.A, self.B
        assert self._f()(0.0) == pytest.approx((b - a) / (a * np.log(b / a)))
        assert self._f()(1.0) == pytest.approx((b - a) / (b * np.log(b / a)))

    def test_profile_depends_on_the_radius_ratio_only(self):
        """F14: scaling both radii leaves f unchanged, and so does d."""
        fd = FieldDistribution("coaxial", coax_a=self.A, coax_b=self.B)
        g = FieldDistribution("coaxial", coax_a=7 * self.A, coax_b=7 * self.B)
        for x in np.linspace(0.0, 1.0, 11):
            assert fd.build(1e-3)(x) == pytest.approx(g.build(5e-2)(x), rel=1e-14)

    def test_thin_annulus_tends_to_uniform(self):
        """F15: b/a -> 1 recovers the plane-parallel gap."""
        f = self._f(a=1.0, b=1.0 + 1e-4)
        assert np.allclose([f(x) for x in np.linspace(0, 1, 11)], 1.0, atol=1e-4)
        assert self._f(a=1.0, b=1.0 + 1e-14)(0.3) == 1.0

    def test_gap_is_fixed_by_the_radii(self):
        fd = FieldDistribution("coaxial", coax_a=self.A, coax_b=self.B)
        assert fd.fixed_gap_length == pytest.approx(self.B - self.A)
        assert FieldDistribution("sphere-plane", 50e-3).fixed_gap_length is None

    def test_tabulated_copy_gives_the_same_profile(self, tmp_path):
        """F16: the analytic profile and a sampled field line of it agree."""
        f = self._f()
        s = np.linspace(0.0, self.B - self.A, 4001)
        path = _write_line(tmp_path / "coax.dat", [(si, f(si / s[-1])) for si in s])
        g = FieldDistribution.from_fieldline(path).build(s[-1])
        for x in np.linspace(0.0, 1.0, 21):
            assert g(x) == pytest.approx(f(x), rel=1e-4)


class TestFieldLineLayouts:
    def test_two_four_and_six_column_layouts_agree(
        self, tmp_path, sphere_plane_samples
    ):
        """F6: the same physical line, written three ways, gives one profile."""
        d, _, xi, E = sphere_plane_samples
        s = xi * d
        p2 = _write_line(tmp_path / "c2.dat", zip(s, E))
        p4 = _write_line(
            tmp_path / "c4.dat", zip(s, np.zeros_like(s), np.zeros_like(s), E)
        )
        # 6-column: same |E| carried as a vector along +x
        p6 = _write_line(
            tmp_path / "c6.dat",
            zip(
                s,
                np.zeros_like(s),
                np.zeros_like(s),
                E,
                np.zeros_like(s),
                np.zeros_like(s),
            ),
        )
        fds = [FieldDistribution.from_fieldline(p) for p in (p2, p4, p6)]
        assert all(fd.fieldline_length == pytest.approx(d, rel=1e-12) for fd in fds)
        probe = np.linspace(0.0, 1.0, 101)
        ref = [fds[0].build(d)(x) for x in probe]
        for fd in fds[1:]:
            assert [fd.build(d)(x) for x in probe] == pytest.approx(ref, rel=1e-12)

    def test_field_vector_magnitude_is_used(self, tmp_path):
        """F6b: a 6-column file uses |E|, so the vector direction is irrelevant."""
        s = np.linspace(0.0, 1.0, 21)
        E = 1.0 + s
        rows_pos = list(zip(s, 0 * s, 0 * s, E, 0 * s, 0 * s))
        rows_neg = list(zip(s, 0 * s, 0 * s, -E, 0 * s, 0 * s))
        a = FieldDistribution.from_fieldline(_write_line(tmp_path / "p.dat", rows_pos))
        b = FieldDistribution.from_fieldline(_write_line(tmp_path / "n.dat", rows_neg))
        assert a.fieldline_f == pytest.approx(b.fieldline_f)

    @pytest.mark.parametrize(
        "unit, scale", [("m", 1.0), ("cm", 1e-2), ("mm", 1e-3), ("um", 1e-6)]
    )
    def test_length_units(self, tmp_path, unit, scale):
        """F7: the unit scales the arc length and leaves f(xi) untouched."""
        s = np.linspace(0.0, 10.0, 11)
        E = 2.0 - 0.05 * s
        fd = FieldDistribution.from_fieldline(
            _write_line(tmp_path / f"{unit}.dat", zip(s, E)), unit
        )
        assert fd.fieldline_length == pytest.approx(10.0 * scale, rel=1e-12)
        assert _integral(fd.build(fd.fieldline_length)) == pytest.approx(1.0, rel=1e-6)

    def test_arc_length_of_a_curved_line(self, tmp_path):
        """F8: arc length is the cumulative Euclidean distance, not |x| extent."""
        n, R = 2001, 1.0
        th = np.linspace(0.0, np.pi, n)
        rows = zip(R * np.cos(th), R * np.sin(th), np.zeros(n), np.ones(n))
        fd = FieldDistribution.from_fieldline(_write_line(tmp_path / "arc.dat", rows))
        assert fd.fieldline_length == pytest.approx(np.pi * R, rel=1e-5)

    def test_round_trip_reproduces_the_analytic_profile(
        self, tmp_path, sphere_plane_samples
    ):
        """
        F9: sample the exact sphere-plane profile, write it, load it back.

        This is the end-to-end check that arc length, normalisation and
        interpolation compose correctly.
        """
        d, f_exact, xi, E = sphere_plane_samples
        fd = FieldDistribution.from_fieldline(
            _write_line(
                tmp_path / "sp.csv", zip(xi * d * 1e3, E), header="s_mm,E", sep=","
            )
        )
        f = fd.build(d)
        err = max(abs(f(x) - f_exact(x)) / f_exact(x) for x in np.linspace(0, 1, 500))
        assert err < 1e-5
        assert _integral(f) == pytest.approx(1.0, rel=1e-6)

    def test_profile_is_independent_of_gap_length(self, tmp_path):
        """A tabulated line describes a shape; build(d) must not rescale it."""
        s = np.linspace(0.0, 1.0, 11)
        fd = FieldDistribution.from_fieldline(
            _write_line(tmp_path / "x.dat", zip(s, 2.0 - s))
        )
        a = [fd.build(1e-3)(x) for x in np.linspace(0, 1, 21)]
        b = [fd.build(1.0)(x) for x in np.linspace(0, 1, 21)]
        assert a == pytest.approx(b)


class TestFieldLineRobustness:
    def test_comments_headers_and_separators(self, tmp_path):
        """F10: CSV, whitespace, '#' comments and a text header all parse."""
        p = tmp_path / "messy.csv"
        p.write_text(
            "# a comment\n"
            "s_mm,E_Vpm\n"
            "0.0, 2.0\n"
            "\n"
            "1.0; 1.5   # trailing comment\n"
            "2.0 1.0\n"
        )
        s, E, _ = load_fieldline(str(p), "mm")
        assert len(s) == 3
        assert E == pytest.approx([2.0, 1.5, 1.0])

    def test_duplicate_positions_are_dropped(self, tmp_path):
        p = _write_line(tmp_path / "dup.dat", [(0.0, 1.0), (0.0, 1.0), (1.0, 2.0)])
        s, E, _ = load_fieldline(p)
        assert len(s) == 2

    @pytest.mark.parametrize(
        "rows, match",
        [
            ([(0.0, 1.0), (0.0, 1.0)], "zero length"),
            ([(0.0, 0.0), (1.0, 0.0)], "finite and not all zero"),
            ([(0.0, 1.0)], "at least two points"),
        ],
    )
    def test_degenerate_inputs_are_rejected(self, tmp_path, rows, match):
        """F10b: every failure names what is wrong with the file."""
        p = _write_line(tmp_path / "bad.dat", rows)
        with pytest.raises(ValueError, match=match):
            load_fieldline(p)

    def test_a_line_traced_backwards_reads_the_same(self, tmp_path):
        """
        F10c: the direction a line was traced in is not the reader's business.

        A two-column export often carries a signed axis coordinate rather
        than an arc length, and whether it runs up or down depends on which
        electrode the trace started from.  Both give the same line.
        """
        fwd = _write_line(tmp_path / "f.dat", [(0.0, 3.0), (1.0, 2.0), (2.0, 1.0)])
        back = _write_line(tmp_path / "b.dat", [(-5.0, 3.0), (-6.0, 2.0), (-7.0, 1.0)])
        s_f, e_f, _ = load_fieldline(fwd)
        s_b, e_b, _ = load_fieldline(back)
        assert s_f == pytest.approx(s_b)
        assert e_f == pytest.approx(e_b)
        assert s_f[-1] == pytest.approx(2.0)

    def test_an_offset_coordinate_is_measured_from_the_first_row(self, tmp_path):
        """F10d: only distances along the line matter, not where it sits."""
        p = _write_line(tmp_path / "off.dat", [(-0.055, 2.0), (-0.065, 1.0)])
        s, _, _ = load_fieldline(p, "m")
        assert s[0] == 0.0
        assert s[-1] == pytest.approx(0.010)

    def test_two_column_arc_length_is_still_taken_as_given(self, tmp_path):
        """F10e: a column that already is an arc length must not change."""
        p = _write_line(tmp_path / "arc.dat", [(0.0, 3.0), (0.5, 2.0), (2.0, 1.0)])
        s, _, _ = load_fieldline(p)
        assert s == pytest.approx([0.0, 0.5, 2.0])

    def test_a_column_running_both_ways_is_rejected(self, tmp_path):
        """
        F10f: this is the one reading the file cannot settle.

        An increasing column reads the same as an arc length or as a
        coordinate, and a decreasing one can only be a coordinate.  A column
        that does both is neither, and is most often one coordinate of a
        curved line, where the distance travelled is understated.
        """
        p = _write_line(
            tmp_path / "zig.dat", [(0.0, 3.0), (2.0, 2.0), (1.0, 2.0), (3.0, 1.0)]
        )
        with pytest.raises(ValueError, match="both up and down"):
            load_fieldline(p)

    def test_how_the_position_column_was_read_is_reported(self, tmp_path):
        """F10g: the caller must be able to say which reading was used."""
        up = _write_line(tmp_path / "up.dat", [(0.0, 3.0), (1.0, 1.0)])
        down = _write_line(tmp_path / "down.dat", [(-0.05, 3.0), (-0.06, 1.0)])
        xyz = _write_line(
            tmp_path / "xyz.dat", [(0.0, 0.0, 0.0, 3.0), (1.0, 0.0, 0.0, 1.0)]
        )
        assert "coordinate" in load_fieldline(up)[2]
        assert load_fieldline(down)[2] == "decreasing coordinate"
        assert "cumulative distance" in load_fieldline(xyz)[2]

    def test_wrong_column_count_is_rejected(self, tmp_path):
        p = _write_line(tmp_path / "c3.dat", [(0.0, 1.0, 2.0), (1.0, 1.0, 2.0)])
        with pytest.raises(ValueError, match="expected 2"):
            load_fieldline(p)

    def test_ragged_file_is_rejected(self, tmp_path):
        p = tmp_path / "ragged.dat"
        p.write_text("0.0 1.0\n1.0 1.0 3.0\n")
        with pytest.raises(ValueError, match="Inconsistent number of columns"):
            load_fieldline(str(p))

    def test_empty_file_is_rejected(self, tmp_path):
        p = tmp_path / "empty.dat"
        p.write_text("# only a comment\n")
        with pytest.raises(ValueError, match="No numeric data"):
            load_fieldline(str(p))

    def test_unknown_unit_is_rejected(self, tmp_path):
        p = _write_line(tmp_path / "u.dat", [(0.0, 1.0), (1.0, 1.0)])
        with pytest.raises(ValueError, match="Unknown length unit"):
            load_fieldline(p, "furlong")

    def test_missing_file(self):
        with pytest.raises(OSError):
            load_fieldline("/nonexistent/field/line.dat")


class TestParseFieldSpec:
    def test_uniform(self):
        assert parse_field_spec(["uniform"]).field_type == "uniform"

    @pytest.mark.parametrize("kind", ["sphere-plane", "sphere-sphere"])
    def test_sphere_radius_is_millimetres(self, kind):
        """F11: the CLI takes mm, the object stores metres."""
        fd = parse_field_spec([kind, "50"])
        assert fd.field_type == kind
        assert fd.sphere_R == pytest.approx(50e-3)

    def test_tip_radius_is_millimetres(self):
        fd = parse_field_spec(["hyperboloid-plane", "0.1"])
        assert fd.field_type == "hyperboloid-plane"
        assert fd.tip_R == pytest.approx(1e-4)
        assert fd.sphere_R is None

    def test_coaxial_radii_are_millimetres(self):
        fd = parse_field_spec(["coaxial", "1", "10"])
        assert fd.field_type == "coaxial"
        assert fd.coax_a == pytest.approx(1e-3)
        assert fd.coax_b == pytest.approx(10e-3)

    def test_fieldline(self, tmp_path):
        p = _write_line(tmp_path / "f.dat", [(0.0, 2.0), (1.0, 1.0)])
        fd = parse_field_spec(["fieldline", p, "mm"])
        assert fd.field_type == "fieldline"
        assert fd.fieldline_length == pytest.approx(1e-3)

    @pytest.mark.parametrize(
        "spec, match",
        [
            (["banana"], "Unknown field type"),
            (["sphere-plane"], "requires a sphere radius"),
            (["fieldline"], "requires a data file"),
            (["coaxial", "1"], "inner and an outer radius"),
            (["coaxial", "10", "1"], "inner radius < outer radius"),
            (["coaxial", "0", "1"], "inner radius < outer radius"),
            (["coaxial", "a", "1"], "must be numbers"),
            (["hyperboloid-plane"], "requires a tip radius"),
            (["hyperboloid-plane", "0"], "must be positive"),
            (["hyperboloid-plane", "x"], "must be a number"),
        ],
    )
    def test_errors(self, spec, match):
        with pytest.raises(ValueError, match=match):
            parse_field_spec(spec)

    def test_label_mentions_the_geometry(self, tmp_path):
        assert "uniform" in FieldDistribution("uniform").label
        assert "50" in FieldDistribution("sphere-plane", 50e-3).label
        lbl = FieldDistribution("coaxial", coax_a=1e-3, coax_b=10e-3).label
        assert "coaxial" in lbl and "10 mm" in lbl
        lbl = FieldDistribution("hyperboloid-plane", tip_R=0.2e-3).label
        assert "hyperboloid" in lbl and "0.2 mm" in lbl


class TestFieldLineScaling:
    """The file's units must not reach the answer."""

    @staticmethod
    def _write(tmp_path, name, s_col, e_col):
        path = tmp_path / name
        np.savetxt(path, np.c_[s_col, e_col])
        return str(path)

    @pytest.fixture
    def two_unit_systems(self, tmp_path):
        """The same 20 mm, 2:1 linear ramp at 100 kV, written two ways."""
        s_mm = np.linspace(0.0, 20.0, 41)
        e_kv_mm = np.linspace(20.0 / 3.0, 10.0 / 3.0, 41)
        si = self._write(tmp_path, "si.csv", s_mm * 1e-3, e_kv_mm * 1e6)
        eng = self._write(tmp_path, "eng.csv", s_mm, e_kv_mm)
        return si, eng

    def test_profile_is_independent_of_the_units(self, two_unit_systems):
        """F20: the normalised profile is what the solver sees, and it is units-free."""
        si, eng = two_unit_systems
        a = parse_field_spec(["fieldline", si, "m"])
        b = parse_field_spec(["fieldline", eng, "mm"])
        assert a.fieldline_length == pytest.approx(b.fieldline_length)
        assert a.fieldline_f == pytest.approx(b.fieldline_f)
        xi = np.linspace(0.0, 1.0, 17)
        fa, fb = a.build(a.fieldline_length), b.build(b.fieldline_length)
        assert [fa(x) for x in xi] == pytest.approx([fb(x) for x in xi])

    def test_field_integral_is_not_claimed_to_be_a_voltage(self, two_unit_systems):
        """F21: the integral carries the file's own units, and differs between them.

        It is kept as a units check, which is only useful if it is reported
        as what it is.  The length column is converted to metres either way,
        so the two differ by exactly the field-unit ratio: 1 kV/mm is
        1e6 V/m, and 1e5 V vs 0.1 kV*m/mm is that factor.
        """
        si, eng = two_unit_systems
        a = parse_field_spec(["fieldline", si, "m"])
        b = parse_field_spec(["fieldline", eng, "mm"])
        assert a.fieldline_integral == pytest.approx(1e5, rel=1e-6)
        assert b.fieldline_integral == pytest.approx(0.1, rel=1e-6)
        assert a.fieldline_integral / b.fieldline_integral == pytest.approx(1e6)
        assert a.fieldline_applied_voltage is None
        assert b.fieldline_applied_voltage is None

    def test_declared_excitation_survives_both_unit_systems(self, two_unit_systems):
        """F22: --fieldline-voltage is what a reported ratio may be divided by."""
        si, eng = two_unit_systems
        for path, unit in ((si, "m"), (eng, "mm")):
            fd = parse_field_spec(["fieldline", path, unit], applied_voltage_kv=100.0)
            assert fd.fieldline_applied_voltage == pytest.approx(1e5)

    def test_declared_excitation_is_rejected_for_other_geometries(self):
        """F23: for every other geometry the voltage is a result, not an input."""
        for spec in (["uniform"], ["sphere-plane", "50"], ["sphere-sphere", "50"]):
            with pytest.raises(ValueError, match="only to"):
                parse_field_spec(spec, applied_voltage_kv=100.0)

    def test_a_nonpositive_excitation_is_rejected(self, two_unit_systems):
        """F24: a ratio against zero or a negative voltage is meaningless."""
        si, _ = two_unit_systems
        for bad in (0.0, -1.0):
            with pytest.raises(ValueError, match="positive"):
                parse_field_spec(["fieldline", si, "m"], applied_voltage_kv=bad)


class TestPolarityLabel:
    """The label names the electrode that is the anode in positive polarity."""

    @pytest.mark.parametrize(
        "field, expected",
        [
            (FieldDistribution("uniform"), "positive"),
            (FieldDistribution("sphere-plane", 5e-3), "sphere=positive"),
            (FieldDistribution("sphere-sphere", 5e-3), "sphere=positive"),
            (FieldDistribution("coaxial", coax_a=1e-3, coax_b=10e-3), "inner=positive"),
            (FieldDistribution("hyperboloid-plane", tip_R=1e-4), "tip=positive"),
        ],
    )
    def test_labels(self, field, expected):
        assert field.polarity_label("positive") == expected
        assert field.polarity_label("negative") == expected.replace(
            "positive", "negative"
        )


# ── Protrusions ──────────────────────────────────────────────────────────────


def _with_protrusion(fd, h, b):
    fd.protrusion = Spheroid(h, base_radius=b)
    return fd


def _quad_integral(f):
    from scipy.integrate import quad

    edges = [0.0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    return sum(
        quad(f, a, b, epsabs=0.0, epsrel=1e-10, limit=500)[0]
        for a, b in zip(edges[:-1], edges[1:])
    )


class TestProtrusion:
    H = 0.5e-3

    def test_hemisphere_on_a_uniform_field(self):
        """F47: f = (1 + 2h³/z³)·d / (D − h³/D²), z = h + ξd, D = d + h."""
        h, d = self.H, 10e-3
        D = d + h
        f = _with_protrusion(FieldDistribution("uniform"), h, h).build(d)
        for x in np.linspace(0.0, 1.0, 11):
            z = h + x * d
            exact = (1 + 2 * h**3 / z**3) * d / (D - h**3 / D**2)
            assert f(x) == pytest.approx(exact, rel=1e-9)

    @pytest.mark.parametrize(
        "fd",
        [
            FieldDistribution("uniform"),
            FieldDistribution("sphere-plane", 50e-3),
            FieldDistribution("sphere-sphere", 50e-3),
            FieldDistribution("hyperboloid-plane", tip_R=5e-3),
            FieldDistribution("coaxial", coax_a=5e-3, coax_b=20e-3),
        ],
    )
    @pytest.mark.parametrize(
        "pr",
        [
            Spheroid(0.5e-3, base_radius=0.5e-3),
            Spheroid(1e-3, base_radius=0.05e-3),
            Spheroid(0.2e-3, base_radius=1e-3),
            Cone(1e-3, radius=0.02e-3, half_angle=20.0),
            Rod(1e-3, radius=0.05e-3),
        ],
        ids=["hemisphere", "needle", "bump", "cone", "rod"],
    )
    def test_normalised_and_decreasing(self, fd, pr):
        """F48: ∫f = 1 (the voltage is kept), and f falls away from the tip."""
        fd = dataclasses.replace(fd, protrusion=pr)
        f = fd.build(fd.fixed_gap_length or 10e-3)
        assert _quad_integral(f) == pytest.approx(1.0, rel=1e-8)
        if fd.field_type != "sphere-sphere":
            assert fd.is_monotone_decreasing
            v = np.array([f(x) for x in np.linspace(0.0, 1.0, 400)])
            assert np.all(np.diff(v) < 0.0)

    @pytest.mark.parametrize(
        "fd",
        [
            FieldDistribution("uniform"),
            FieldDistribution("sphere-plane", 50e-3),
            FieldDistribution("hyperboloid-plane", tip_R=5e-3),
        ],
    )
    def test_a_vanishing_protrusion_leaves_the_background(self, fd):
        """
        F49: h → 0 gives back the electrodes alone, away from the tip.

        At the tip itself the enhancement is β however small the protrusion
        is; it is the region it covers that shrinks.
        """
        d = 10e-3
        ref = fd.build(d)
        f = _with_protrusion(dataclasses.replace(fd), 1e-9, 1e-9).build(d)
        for x in np.linspace(0.01, 1.0, 11):
            assert f(x) == pytest.approx(ref(x), rel=1e-6)

    def test_gap_is_measured_from_the_tip(self):
        """F50: fixed geometries lose h; the background is read from the tip."""
        a, b, h = 5e-3, 20e-3, self.H
        fd = _with_protrusion(FieldDistribution("coaxial", coax_a=a, coax_b=b), h, h)
        assert fd.fixed_gap_length == pytest.approx(b - a - h)
        # Far from the tip the enhancement is gone, so f follows 1/r.
        f = fd.build(fd.fixed_gap_length)
        r = lambda x: a + h + x * (b - a - h)  # noqa: E731
        assert f(0.9) / f(1.0) == pytest.approx(r(1.0) / r(0.9), rel=1e-4)

    def test_flags_and_labels(self):
        """F51: a protrusion breaks the symmetry and is named in the labels."""
        fd = _with_protrusion(FieldDistribution("uniform"), self.H, self.H)
        assert not fd.is_symmetric
        assert fd.polarity_label("negative") == "protrusion=negative"
        assert "spheroid protrusion h = 0.5 mm" in fd.label
        sp = _with_protrusion(FieldDistribution("sphere-plane", 5e-3), self.H, self.H)
        assert sp.polarity_label("positive") == "sphere=positive"


class TestProtrusionParsing:
    def test_spec_is_parsed_in_mm(self):
        """F52."""
        fd = parse_field_spec(["uniform"], protrusion=["spheroid", "1", "0.2"])
        pr = fd.protrusion
        assert (pr.height, pr.base_radius) == pytest.approx((1e-3, 0.2e-3))

    @pytest.mark.parametrize(
        "spec, match",
        [
            (["pyramid", "1", "1"], "unknown shape"),
            (["spheroid", "1"], "height and base radius"),
            (["spheroid", "1", "x"], "must be numbers"),
            (["rod", "1", "2"], "at least the radius"),
            (["cone", "1", "0.1", "85"], "half-angle"),
            (["cone", "1", "0.1"], "tip radius, half-angle"),
            (["spheroid", "0", "1"], "positive"),
        ],
    )
    def test_bad_specs_are_rejected(self, spec, match):
        """F53."""
        with pytest.raises(ValueError, match=match):
            parse_field_spec(["uniform"], protrusion=spec)

    def test_a_protrusion_filling_the_gap_is_rejected(self):
        """F54: coaxial b − a = 1 mm cannot hold a 2 mm protrusion."""
        with pytest.raises(ValueError, match="leaves no gap"):
            parse_field_spec(["coaxial", "1", "2"], protrusion=["spheroid", "2", "1"])

    def test_notes_on_scale(self, capsys):
        """F55: h against the electrode radius is reported, not refused."""
        parse_field_spec(["sphere-plane", "1"], protrusion=["spheroid", "0.5", "0.5"])
        assert "of the sphere radius" in capsys.readouterr().err
        parse_field_spec(["sphere-plane", "100"], protrusion=["spheroid", "0.5", "0.5"])
        assert capsys.readouterr().err == ""

    def test_gap_note(self):
        """F56: the gap check needs the gap, and is last in the list."""
        from incept1d.fields import protrusion_notes

        fd = _with_protrusion(FieldDistribution("uniform"), 1e-3, 1e-3)
        assert protrusion_notes(fd) == []
        assert "smallest gap" in protrusion_notes(fd, 2e-3)[-1]
        assert protrusion_notes(fd, 1.0) == []


class TestReverseField:
    @staticmethod
    def _line(tmp_path):
        s = np.linspace(0.0, 10.0, 11)
        return _write_line(tmp_path / "up.dat", zip(s, 1.0 + 0.3 * s))

    @pytest.mark.parametrize(
        "spec, electrode",
        [
            (["sphere-plane", "5"], "plane"),
            (["hyperboloid-plane", "0.5"], "plane"),
            (["coaxial", "1", "10"], "outer"),
        ],
    )
    def test_reversed_profile_is_the_mirror_image(self, spec, electrode):
        """F57: f_rev(ξ) = f(1 − ξ), and ξ = 0 is named in the labels."""
        fwd = parse_field_spec(spec)
        rev = parse_field_spec(spec, reverse=True)
        d = fwd.fixed_gap_length or 10e-3
        f, g = fwd.build(d), rev.build(d)
        for x in np.linspace(0.0, 1.0, 21):
            assert g(x) == pytest.approx(f(1.0 - x), rel=1e-12)
        assert rev.polarity_label("negative") == f"{electrode}=negative"
        assert f"xi = 0 at the {electrode}" in rev.label
        assert fwd.is_monotone_decreasing and not rev.is_monotone_decreasing
        assert rev.fixed_gap_length == fwd.fixed_gap_length

    def test_a_field_line_is_read_from_its_last_point(self, tmp_path):
        """F58: for a table, reversing makes the last row ξ = 0."""
        path = self._line(tmp_path)
        fwd = parse_field_spec(["fieldline", path])
        rev = parse_field_spec(["fieldline", path], reverse=True)
        L = fwd.fieldline_length
        f, g = fwd.build(L), rev.build(L)
        for x in np.linspace(0.0, 1.0, 21):
            assert g(x) == pytest.approx(f(1.0 - x), rel=1e-12)
        assert rev.polarity_label("positive") == "end=positive"

    @pytest.mark.parametrize("spec", [["uniform"], ["sphere-sphere", "5"]])
    def test_a_symmetric_gap_is_unchanged(self, spec, capsys):
        """F59: nothing to reverse, and the command says so."""
        fd = parse_field_spec(spec, reverse=True)
        assert not fd.reversed
        assert "changes nothing" in capsys.readouterr().err

    def test_a_protrusion_on_the_plane(self):
        """F60: reversed sphere-plane puts the protrusion on the plane."""
        fd = parse_field_spec(
            ["sphere-plane", "50"], reverse=True, protrusion=["spheroid", "0.1", "0.1"]
        )
        bg = parse_field_spec(["sphere-plane", "50"]).build(10e-3 + 0.1e-3)
        f = fd.build(10e-3)
        # Far from the tip, the plane side of the background, renormalised.
        ratio = f(0.5) / bg(1.0 - (0.1e-3 + 0.5 * 10e-3) / 10.1e-3)
        assert f(0.9) / bg(1.0 - (0.1e-3 + 0.9 * 10e-3) / 10.1e-3) == pytest.approx(
            ratio, rel=1e-3
        )
        assert f(0.0) > 2.5 * f(0.5)
        assert fd.polarity_label("positive") == "plane=positive"

    def test_protrusion_at_the_weak_end_is_pointed_out(self, tmp_path, capsys):
        """F61: a line rising towards ξ = 1 suggests --reverse-field."""
        path = self._line(tmp_path)
        parse_field_spec(["fieldline", path], protrusion=["spheroid", "1", "1"])
        assert "--reverse-field" in capsys.readouterr().err
        parse_field_spec(
            ["fieldline", path], reverse=True, protrusion=["spheroid", "1", "1"]
        )
        assert "--reverse-field" not in capsys.readouterr().err

    def test_the_scale_note_follows_the_electrode(self, capsys):
        """F62: reversed coaxial compares h with the outer radius."""
        parse_field_spec(["coaxial", "0.5", "5"], protrusion=["spheroid", "0.2", "0.2"])
        assert "inner radius" in capsys.readouterr().err
        parse_field_spec(
            ["coaxial", "0.5", "5"], reverse=True, protrusion=["spheroid", "0.2", "0.2"]
        )
        assert capsys.readouterr().err == ""


def test_a_protrusion_on_a_plate_is_not_uniform():
    """F63: the constant-field shortcut must not swallow the protrusion."""
    assert FieldDistribution("uniform").is_uniform
    fd = _with_protrusion(FieldDistribution("uniform"), 1e-3, 1e-3)
    assert fd.field_type == "uniform" and not fd.is_uniform


class TestFieldUnit:
    @pytest.mark.parametrize(
        "text, value",
        [
            ("V/m", 1.0),
            ("kV/mm", 1e6),
            ("kV/cm", 1e5),
            ("V/mm", 1e3),
            ("MV/m", 1e6),
            ("mV/um", 1e3),
            ("1e3*V/m", 1e3),
            ("2.5*kV/cm", 2.5e5),
        ],
    )
    def test_units(self, text, value):
        """F70: the value of one unit of the field column in V/m."""
        assert parse_field_unit(text) == pytest.approx(value, rel=1e-15)

    @pytest.mark.parametrize(
        "text",
        [
            "kV",
            "kV/ft",
            "kv/mm",
            "x*V/m",
            "-1*V/m",
            "V/m/s",
            "0*V/m",
            "*V/m",
            "1e400*V/m",
        ],
    )
    def test_rejected(self, text):
        """F71."""
        with pytest.raises(ValueError):
            parse_field_unit(text)

    @staticmethod
    def _line(tmp_path, scale):
        # 20 mm, |E| from 2 to 1 kV/mm: ∫|E| ds = 30 kV, written in V/m / scale.
        s = np.linspace(0.0, 20.0, 41)
        return _write_line(tmp_path / "l.dat", zip(s, (2.0 - 0.05 * s) * 1e6 / scale))

    @pytest.mark.parametrize(
        "unit, scale", [("V/m", 1.0), ("kV/mm", 1e6), ("kV/cm", 1e5)]
    )
    def test_the_integral_is_the_excitation(self, tmp_path, unit, scale):
        """F72: any declared unit gives the same 30 kV, and the same profile."""
        path = self._line(tmp_path, scale)
        fd = parse_field_spec(["fieldline", path, "mm", unit])
        assert fd.fieldline_line_voltage == pytest.approx(30e3, rel=1e-12)
        assert fd.fieldline_applied_voltage == pytest.approx(30e3, rel=1e-12)
        assert "∫|E| ds" in fd.fieldline_voltage_source
        ref = parse_field_spec(["fieldline", path, "mm"]).build(20e-3)
        f = fd.build(20e-3)
        assert all(f(x) == ref(x) for x in np.linspace(0, 1, 11))

    def test_order_and_default_length(self, tmp_path):
        """F73: the field unit may come first, and the length defaults to m."""
        path = self._line(tmp_path, 1e6)
        a = parse_field_spec(["fieldline", path, "kV/mm", "mm"])
        assert a.fieldline_line_voltage == pytest.approx(30e3)
        b = parse_field_spec(["fieldline", path, "kV/mm"])  # read as metres
        assert b.fieldline_line_voltage == pytest.approx(30e6)

    def test_declared_voltage_wins_and_a_mismatch_is_noted(self, tmp_path, capsys):
        """F74: --fieldline-voltage overrides the integral; 2 % is the tolerance."""
        path = self._line(tmp_path, 1e6)
        fd = parse_field_spec(
            ["fieldline", path, "mm", "kV/mm"], applied_voltage_kv=30.3
        )
        assert fd.fieldline_applied_voltage == pytest.approx(30.3e3)
        assert capsys.readouterr().err == ""
        parse_field_spec(["fieldline", path, "mm", "kV/mm"], applied_voltage_kv=100)
        assert "the declared excitation is used" in capsys.readouterr().err

    @pytest.mark.parametrize(
        "extra", [["mm", "cm"], ["kV/mm", "V/m"], ["mm", "kV/mm", "x"]]
    )
    def test_malformed_specs(self, tmp_path, extra):
        """F75."""
        path = self._line(tmp_path, 1e6)
        with pytest.raises(ValueError):
            parse_field_spec(["fieldline", path, *extra])


class TestAuditFixes:
    @pytest.mark.parametrize(
        "spec",
        [["uniform", "kV/mm"], ["sphere-plane", "5", "x"], ["sphere-plane", "-1"]],
    )
    def test_extra_or_bad_tokens_are_refused(self, spec):
        """F80: tokens an analytic geometry has no use for are not dropped."""
        with pytest.raises(ValueError):
            parse_field_spec(spec)

    def test_a_wide_protrusion_is_noted(self, capsys):
        """F81: the footprint, not only the height, must be small."""
        parse_field_spec(["sphere-plane", "5"], protrusion=["spheroid", "0.2", "4"])
        assert "base radius" in capsys.readouterr().err

    def test_sphere_sphere_names_the_protrusion(self):
        """F82: 'sphere=' would not say which sphere."""
        fd = parse_field_spec(["sphere-sphere", "50"], protrusion=["rod", "1", "0.1"])
        assert fd.polarity_label("negative") == "protrusion=negative"

    def test_fixed_gap_wording(self, tmp_path):
        """F83: with a protrusion the fixed gap is L - h, and scaling says so."""
        from incept1d.fields import describe_fixed_gap, fixed_gap_warning

        s = np.linspace(0.0, 20.0, 21)
        path = _write_line(tmp_path / "l.dat", zip(s, 2.0 - 0.05 * s))
        fd = parse_field_spec(
            ["fieldline", path, "mm"], protrusion=["spheroid", "1", "1"]
        )
        assert "d = L - h = 20 - 1 = 19 mm" in describe_fixed_gap(fd)
        assert fixed_gap_warning(fd, [19.0]) is None
        assert "protrusion at its given size" in fixed_gap_warning(fd, [20.0])
        assert describe_fixed_gap(FieldDistribution("uniform")) is None

    def test_header_lines(self, tmp_path):
        """F84: one set of geometry lines for every command's results file."""
        from incept1d.fields import field_header_lines

        s = np.linspace(0.0, 20.0, 21)
        path = _write_line(tmp_path / "l.dat", zip(s, 2.0 - 0.05 * s))
        fd = parse_field_spec(
            ["fieldline", path, "mm", "kV/mm"],
            reverse=True,
            protrusion=["rod", "1", "0.05"],
        )
        keys = [line.split(":")[0] for line in field_header_lines(fd)]
        assert keys == [
            "Protrusion",
            "Reversed",
            "Polarity",
            "Field line",
            "Arc length",
            "∫|E| ds",
            "U_applied",
        ]
