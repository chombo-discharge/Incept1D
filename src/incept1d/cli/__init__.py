# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
``incept1d`` command-line interface.

Each subcommand lives in its own module of this package and exposes

* ``HELP`` — one-line summary shown in ``incept1d --help``,
* ``DESCRIPTION`` — longer text shown in ``incept1d <command> --help``,
* ``add_arguments(parser)`` — registers the command's arguments,
* ``run(args, parser)`` — executes the command.

The physics lives in the top-level ``incept1d`` modules; the CLI modules only
parse arguments, drive the solvers, print tables and plot.
"""

import argparse
import importlib

from incept1d import __version__

#: Subcommand name → module (imported lazily so ``--help`` stays fast).
COMMANDS = {
    "pdiv": "incept1d.cli.pdiv",
    "eigenvalues": "incept1d.cli.eigenvalues",
    "ionization": "incept1d.cli.ionization",
    "growth": "incept1d.cli.growth",
    "field": "incept1d.cli.field",
    "chombo": "incept1d.cli.chombo",
}


def build_parser():
    """Return the top-level ``incept1d`` argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="incept1d",
        description=(
            "Incept1D — discharge inception in a 1-D drift-reaction model: "
            "inception curves (PDIV vs. p·d), local eigenvalues, ionization integrals, temporal "
            "growth rates and transport tables for 3-D solvers."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)
    for name, modname in COMMANDS.items():
        mod = importlib.import_module(modname)
        p = sub.add_parser(name, help=mod.HELP, description=mod.DESCRIPTION)
        mod.add_arguments(p)
        p.set_defaults(_run=mod.run, _parser=p)
    return parser


def main(argv=None):
    """Entry point of the ``incept1d`` console script."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args._run(args, args._parser)
