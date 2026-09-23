# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Tabulate the closed-form Paschen curve for the gases in ``gases.json``.

    U = B (pd) / [ln(A pd) - ln(ln(1 + 1/gamma))]

evaluated from the very coefficients the mechanism uses, so the figure in the
documentation compares the solver against its own analytic limit rather than
against an external table.

Usage::

    python3 examples/paschen/closed_form.py OUTPUT.dat
"""

import sys

import numpy as np

from incept1d.mechanism import load_mechanism

GASES = ["helium", "argon", "air"]
MECHANISM = "mechanisms/paschen/paschen.py"


def main(out_path, n=60):
    # Sparse on purpose: these are plotted as markers over the solver's
    # continuous curve, so they must read as discrete samples.
    pd = np.logspace(np.log10(5e-6), np.log10(1e-1), n)  # bar m
    columns = [pd * 1e3]
    for gas in GASES:
        raw = load_mechanism(MECHANISM, {"gas": gas})._mod
        columns.append(np.array([raw.paschen_voltage(x) / 1e3 for x in pd]))

    with open(out_path, "w") as fh:
        fh.write("# Closed-form Paschen curve of the Townsend model\n")
        fh.write("#   U = B pd / [ln(A pd) - ln(ln(1 + 1/gamma))]\n")
        fh.write(
            "# Columns: pd (bar mm)  " + "  ".join(f"U_{g} (kV)" for g in GASES) + "\n"
        )
        for row in np.column_stack(columns):
            fh.write(
                "  ".join("nan" if np.isnan(v) else f"{v:.8e}" for v in row) + "\n"
            )
    print(f"Results written to: {out_path}")


if __name__ == "__main__":
    main(sys.argv[1])
