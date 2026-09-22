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


class TestFieldLineLayouts:
    def test_two_four_and_six_column_layouts_agree(self, tmp_path, sphere_plane_samples):
        """F6: the same physical line, written three ways, gives one profile."""
        d, _, xi, E = sphere_plane_samples
        s = xi * d
        p2 = _write_line(tmp_path / "c2.dat", zip(s, E))
        p4 = _write_line(tmp_path / "c4.dat", zip(s, np.zeros_like(s), np.zeros_like(s), E))
        # 6-column: same |E| carried as a vector along +x
        p6 = _write_line(
            tmp_path / "c6.dat",
            zip(s, np.zeros_like(s), np.zeros_like(s), E, np.zeros_like(s), np.zeros_like(s)),
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

    @pytest.mark.parametrize("unit, scale", [("m", 1.0), ("cm", 1e-2), ("mm", 1e-3), ("um", 1e-6)])
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

    def test_round_trip_reproduces_the_analytic_profile(self, tmp_path, sphere_plane_samples):
        """
        F9: sample the exact sphere-plane profile, write it, load it back.

        This is the end-to-end check that arc length, normalisation and
        interpolation compose correctly.
        """
        d, f_exact, xi, E = sphere_plane_samples
        fd = FieldDistribution.from_fieldline(
            _write_line(tmp_path / "sp.csv", zip(xi * d * 1e3, E), header="s_mm,E", sep=",")
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
        s, E = load_fieldline(str(p), "mm")
        assert len(s) == 3
        assert E == pytest.approx([2.0, 1.5, 1.0])

    def test_duplicate_positions_are_dropped(self, tmp_path):
        p = _write_line(tmp_path / "dup.dat", [(0.0, 1.0), (0.0, 1.0), (1.0, 2.0)])
        s, E = load_fieldline(p)
        assert len(s) == 2

    @pytest.mark.parametrize(
        "rows, match",
        [
            ([(0.0, 1.0), (-1.0, 1.0)], "monotonically increasing"),
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
        ],
    )
    def test_errors(self, spec, match):
        with pytest.raises(ValueError, match=match):
            parse_field_spec(spec)

    def test_label_mentions_the_geometry(self, tmp_path):
        assert "uniform" in FieldDistribution("uniform").label
        assert "50" in FieldDistribution("sphere-plane", 50e-3).label
