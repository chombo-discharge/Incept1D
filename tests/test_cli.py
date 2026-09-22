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
        rc = _run("pdiv", TOY, "--p", "1", "--pd-min", "1", "--pd-max", "50",
                  "--pd-num", "4", "--no-plot")
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
        rc = _run("pdiv", TOY, "--p", "1", "--pd-min", "1", "--pd-max", "50",
                  "--pd-num", "5", "--no-plot", "--write-to-file", str(out_file))
        capsys.readouterr()
        assert rc == 0 and out_file.exists()

        lines = out_file.read_text().splitlines()
        declared = [ln for ln in lines if ln.startswith("# Column ")]
        data = [ln for ln in lines if ln.strip() and not ln.startswith("#")]
        assert declared and data
        assert len(data[0].split()) == len(declared)

    def test_metadata_header_is_present(self, tmp_path, capsys):
        out_file = tmp_path / "curve.dat"
        _run("pdiv", TOY, "--p", "1", "--pd-min", "1", "--pd-max", "10",
             "--pd-num", "2", "--no-plot", "--write-to-file", str(out_file))
        capsys.readouterr()
        text = out_file.read_text()
        for field in ("# Date:", "# Git:", "# Command:", "# Mechanism:",
                      "# Temperature:", "# Field type:"):
            assert field in text

    def test_fixed_distance_header_records_the_pressures(self, tmp_path, capsys):
        """
        C5: in fixed-d mode the pressures vary along the sweep.

        The header used to print an empty '# Pressures:' because it reported
        the (unused) --p list rather than the pressures actually solved.
        """
        out_file = tmp_path / "fixed_d.dat"
        _run("pdiv", TOY, "--d", "10", "--pd-min", "1", "--pd-max", "50",
             "--pd-num", "4", "--no-plot", "--write-to-file", str(out_file))
        capsys.readouterr()
        line = next(
            ln for ln in out_file.read_text().splitlines()
            if ln.startswith("# Pressures:")
        )
        assert line.split(":", 1)[1].strip(), "header declares no pressures"

    def test_fieldline_through_the_cli(self, tmp_path, capsys):
        """C7: --field fieldline is accepted and reports the input excitation."""
        line = tmp_path / "line.dat"
        line.write_text("\n".join(f"{s} {2.0 - 0.05 * s}" for s in range(0, 21)))
        rc = _run("pdiv", TOY, "--field", "fieldline", str(line), "mm",
                  "--pd-min", "1", "--pd-max", "20", "--pd-num", "3", "--no-plot")
        out = capsys.readouterr().out
        assert rc == 0
        assert "U*/U_file" in out, "the scale factor column is missing"


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
            capture_output=True, text=True,
            env={**os.environ, "MPLBACKEND": "Agg"},
        )
        # No subcommand given -> argparse error, but the module must import.
        assert "Traceback" not in r.stderr or "the following arguments" in r.stderr
