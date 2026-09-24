# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Gap geometry: normalisation, symmetry, and the tabulated field-line loader."""

import numpy as np
import pytest

from incept1d.fields import FieldDistribution, load_fieldline, parse_field_spec

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
        ],
    )
    def test_labels(self, field, expected):
        assert field.polarity_label("positive") == expected
        assert field.polarity_label("negative") == expected.replace(
            "positive", "negative"
        )
