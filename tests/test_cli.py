# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Command-line wiring.

These are cheap but disproportionately valuable: a stale call signature once
broke the ionization-integral script for an unknown length of time, because
nothing imported it and nothing invoked it.
"""

import importlib
import os
import subprocess
import sys
import numpy as np

import pytest

from incept1d.cli import COMMANDS, build_parser, main

TOY = os.path.join(os.path.dirname(__file__), "toy", "toy_mechanism.py")


def _run(*args):
    """Invoke the console entry point in-process and capture the exit code."""
    try:
        return main(list(args)) or 0
    except SystemExit as exc:  # argparse exits this way
        return exc.code or 0


class TestWiring:
    def test_every_command_module_implements_the_contract(self):
        """C2: HELP / DESCRIPTION / add_arguments / run on every subcommand."""
        assert COMMANDS, "no subcommands registered"
        for name, modname in COMMANDS.items():
            mod = importlib.import_module(modname)
            assert isinstance(mod.HELP, str) and mod.HELP
            assert isinstance(mod.DESCRIPTION, str) and mod.DESCRIPTION
            assert callable(mod.add_arguments)
            assert callable(mod.run)

    def test_parser_builds_with_every_subcommand(self):
        """C1: a broken add_arguments would raise here."""
        parser = build_parser()
        actions = [a for a in parser._actions if hasattr(a, "choices") and a.choices]
        registered = set()
        for a in actions:
            registered |= set(a.choices)
        assert set(COMMANDS) <= registered

    @pytest.mark.parametrize("cmd", sorted(COMMANDS))
    def test_subcommand_help_exits_cleanly(self, cmd, capsys):
        """C1b: `incept1d <cmd> --help` must not crash."""
        assert _run(cmd, "--help") == 0
        assert capsys.readouterr().out

    def test_top_level_help(self, capsys):
        assert _run("--help") == 0
        assert "incept1d" in capsys.readouterr().out


class TestPdiv:
    def test_end_to_end_on_the_toy_mechanism(self, capsys):
        """C3: a real solve through the CLI, no plotting."""
        rc = _run(
            "pdiv",
            TOY,
            "--p",
            "1",
            "--pd-min",
            "1",
            "--pd-max",
            "50",
            "--pd-num",
            "4",
            "--no-plot",
        )
        out = capsys.readouterr().out
        assert rc == 0
        assert "pd (bar·mm)" in out and "E/N (Td)" in out

    def test_write_to_file_round_trip(self, tmp_path, capsys):
        """
        C4: the header must describe the data.

        Column count and the declared column names have to agree, or a reader
        (including our own pgfplots figures) silently plots the wrong series.
        """
        out_file = tmp_path / "curve.dat"
        rc = _run(
            "pdiv",
            TOY,
            "--p",
            "1",
            "--pd-min",
            "1",
            "--pd-max",
            "50",
            "--pd-num",
            "5",
            "--no-plot",
            "--write-to-file",
            str(out_file),
        )
        capsys.readouterr()
        assert rc == 0 and out_file.exists()

        lines = out_file.read_text().splitlines()
        declared = [ln for ln in lines if ln.startswith("# Column ")]
        data = [ln for ln in lines if ln.strip() and not ln.startswith("#")]
        assert declared and data
        assert len(data[0].split()) == len(declared)

    def test_metadata_header_is_present(self, tmp_path, capsys):
        out_file = tmp_path / "curve.dat"
        _run(
            "pdiv",
            TOY,
            "--p",
            "1",
            "--pd-min",
            "1",
            "--pd-max",
            "10",
            "--pd-num",
            "2",
            "--no-plot",
            "--write-to-file",
            str(out_file),
        )
        capsys.readouterr()
        text = out_file.read_text()
        for field in (
            "# Date:",
            "# Git:",
            "# Command:",
            "# Mechanism:",
            "# Temperature:",
            "# Field type:",
        ):
            assert field in text

    def test_fixed_distance_header_records_the_pressures(self, tmp_path, capsys):
        """
        C5: in fixed-d mode the pressures vary along the sweep.

        The header used to print an empty '# Pressures:' because it reported
        the (unused) --p list rather than the pressures actually solved.
        """
        out_file = tmp_path / "fixed_d.dat"
        _run(
            "pdiv",
            TOY,
            "--d",
            "10",
            "--pd-min",
            "1",
            "--pd-max",
            "50",
            "--pd-num",
            "4",
            "--no-plot",
            "--write-to-file",
            str(out_file),
        )
        capsys.readouterr()
        line = next(
            ln
            for ln in out_file.read_text().splitlines()
            if ln.startswith("# Pressures:")
        )
        assert line.split(":", 1)[1].strip(), "header declares no pressures"

    @staticmethod
    def _fieldline(tmp_path):
        """A 20 mm line whose |E| falls by half, in unstated field units."""
        line = tmp_path / "line.dat"
        line.write_text("\n".join(f"{s} {2.0 - 0.05 * s}" for s in range(0, 21)))
        return str(line)

    def _run_fieldline(self, line, *extra):
        return _run(
            "pdiv",
            TOY,
            "--field",
            "fieldline",
            line,
            "mm",
            *extra,
            "--pd-min",
            "1",
            "--pd-max",
            "20",
            "--pd-num",
            "3",
            "--no-plot",
        )

    def test_fieldline_through_the_cli(self, tmp_path, capsys):
        """C7: --field fieldline is accepted and solves."""
        rc = self._run_fieldline(self._fieldline(tmp_path))
        out = capsys.readouterr().out
        assert rc == 0
        assert "E/N (Td)" in out

    def test_no_ratio_is_reported_without_a_declared_excitation(self, tmp_path, capsys):
        """
        C7b: an undeclared excitation must produce no ratio, not a wrong one.

        The field units of the file are unknown, so the line integral is not
        a voltage and nothing can be divided by it.  Reporting the ratio
        anyway made the same line in kV/mm disagree with itself in V/m by
        six orders of magnitude.
        """
        rc = self._run_fieldline(self._fieldline(tmp_path))
        out = capsys.readouterr().out
        assert rc == 0
        assert "U*/U_applied" not in out
        assert "U*/U_file" not in out

    def test_declared_excitation_is_reported_as_a_ratio(self, tmp_path, capsys):
        """C7c: --fieldline-voltage brings the ratio column back, correctly."""
        rc = self._run_fieldline(
            self._fieldline(tmp_path), "--fieldline-voltage", "100"
        )
        out = capsys.readouterr().out
        assert rc == 0
        assert "U*/U_applied" in out, "the scale factor column is missing"

    def test_a_field_line_sweeps_pressure_at_fixed_geometry(self, tmp_path, capsys):
        """
        C7e: a tabulated line is one geometry at one size.

        Its arc length fixes the gap, so pd must vary through p alone.  A
        swept d would silently rescale the whole electrode arrangement,
        which is not what importing a line means.
        """
        out = tmp_path / "sim.dat"
        rc = self._run_fieldline(self._fieldline(tmp_path), "--write-to-file", str(out))
        assert rc == 0
        assert "using d = L" in capsys.readouterr().out

        data = np.genfromtxt(out)
        # With d pinned at L, the sweep column is the pressure itself.
        sweep, p_col, d_col = data[:, 0], data[:, 1], data[:, 2]
        assert np.allclose(d_col, d_col[0]), "the gap length must not vary"
        assert d_col[0] == pytest.approx(20.0), "the gap must be the arc length"
        assert sweep.max() > sweep.min(), "the pressure must be what varies"
        assert np.allclose(sweep, p_col), "the sweep column must be the pressure"

        header = out.read_text()
        assert "# Column 1: p_bar\n" in header
        assert "pd_bar_mm" not in header, "p*d says nothing the pressure does not"

    def test_fixed_pressure_is_refused_for_a_field_line(self, tmp_path, capsys):
        """
        C7f: --p would make the gap the swept variable.

        Every point of such a sweep is a differently sized copy of the
        imported arrangement, which is not a question a tabulated line can
        answer, so it is refused rather than warned about.
        """
        rc = self._run_fieldline(self._fieldline(tmp_path), "--p", "1.0")
        assert rc != 0
        err = capsys.readouterr().err
        assert "--p cannot be used with --field fieldline" in err
        assert "--pd-num 1" in err, "the error must say how to ask for one pressure"

    def _run_coaxial(self, *extra):
        return _run(
            "pdiv",
            TOY,
            "--field",
            "coaxial",
            "1",
            "10",
            *extra,
            "--pd-min",
            "5",
            "--pd-max",
            "50",
            "--pd-num",
            "3",
            "--no-plot",
        )

    def test_coaxial_sweeps_pressure_at_fixed_geometry(self, tmp_path, capsys):
        """C7g: the radii fix the gap at b − a, so pd varies through p."""
        out = tmp_path / "sim.dat"
        rc = self._run_coaxial("--write-to-file", str(out))
        stdout = capsys.readouterr().out
        assert rc == 0
        assert "using d = L = 9 mm" in stdout
        assert "inner=positive" in stdout and "inner=negative" in stdout

        header = out.read_text()
        assert "# Radii:       a = 1 mm, b = 10 mm" in header
        names = [
            ln.split(":", 1)[1].strip()
            for ln in header.splitlines()
            if ln.startswith("# Column ")
        ]
        data = np.genfromtxt(out)
        col = {n: data[:, i] for i, n in enumerate(names)}
        d_cols = [c for n, c in col.items() if n.startswith("d_mm[")]
        p_cols = [c for n, c in col.items() if n.startswith("p_bar[")]
        assert len(d_cols) == len(p_cols) == 2, "one curve per polarity"

        # Every curve, both polarities: the gap is b - a at every point, the
        # pressure is the swept coordinate, and p * (b - a) reproduces the
        # requested pd grid -- so all of the pd variation is pressure.
        sweep = col["p_bar"]
        pd_requested = np.logspace(np.log10(5.0), np.log10(50.0), 3)
        for d_col, p_col in zip(d_cols, p_cols):
            assert np.all(d_col == 9.0), "the gap must be b - a"
            assert np.array_equal(p_col, sweep), "p must be the sweep"
            assert np.allclose(p_col * d_col, pd_requested, rtol=1e-6)  # %.7g
        assert sweep.max() / sweep.min() == pytest.approx(10.0), "p spans pd"

    def test_fixed_pressure_is_refused_for_coaxial(self, capsys):
        """C7h: as for a field line, --p would resize the arrangement."""
        rc = self._run_coaxial("--p", "1.0")
        assert rc != 0
        assert "--p cannot be used with --field coaxial" in capsys.readouterr().err

    def test_declared_excitation_is_refused_for_a_uniform_gap(self, tmp_path, capsys):
        """C7d: there the voltage is a result, so accepting it would mislead."""
        rc = _run("pdiv", TOY, "--fieldline-voltage", "100", "--no-plot")
        assert rc != 0


class TestArgumentValidation:
    @pytest.mark.parametrize(
        "args",
        [
            ("pdiv", TOY, "--streamer-criterion", "0"),
            ("pdiv", TOY, "--dx", "0"),
            ("pdiv", TOY, "--field", "banana"),
            ("pdiv", TOY, "--field", "sphere-plane"),
            ("field", "--field", "uniform"),  # --d is required here
        ],
    )
    def test_bad_arguments_exit_with_an_error(self, args):
        """C6: argparse must reject these, not crash or silently proceed."""
        with pytest.raises(SystemExit) as exc:
            main(list(args))
        assert exc.value.code != 0


class TestConsoleScript:
    def test_installed_entry_point_runs(self):
        """The packaging actually produces a working `incept1d` command."""
        r = subprocess.run(
            [sys.executable, "-m", "incept1d.cli"],
            capture_output=True,
            text=True,
            env={**os.environ, "MPLBACKEND": "Agg"},
        )
        # No subcommand given -> argparse error, but the module must import.
        assert "Traceback" not in r.stderr or "the following arguments" in r.stderr
