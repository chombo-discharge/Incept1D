# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Tabulate the closed-form Paschen curve for the gases in ``gases.json``.

    U = B (pd) / [ln(A pd) - ln(ln(1 + 1/gamma))]

evaluated from the very coefficients the configurations give the mechanism,
so the figure in the documentation compares the solver against its own
analytic limit rather than against an external table.

Usage::

    python3 examples/paschen/closed_form.py OUTPUT.dat
"""

import json
import sys

import numpy as np

CONFIGS = "mechanisms/paschen/gases.json"
#: [cm^-1 Torr^-1] -> [m^-1 bar^-1]
PER_CM_TORR = 100.0 * 750.061682704


def paschen_voltage(pd, A, B, gamma):
    """Breakdown voltage [V] at pd [bar m]; NaN below the left asymptote."""
    denom = np.log(A * PER_CM_TORR * pd) - np.log(np.log1p(1.0 / gamma))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denom > 0.0, B * PER_CM_TORR * pd / denom, np.nan)


def main(out_path, n=60):
    # Sparse on purpose: these are plotted as markers over the solver's
    # continuous curve, so they must read as discrete samples.
    pd = np.logspace(np.log10(5e-6), np.log10(1e-1), n)  # bar m
    with open(CONFIGS) as fh:
        gases = json.load(fh)["configurations"]
    columns = [pd * 1e3]
    for g in gases:
        columns.append(paschen_voltage(pd, g["A"], g["B"], g["gamma0"]) / 1e3)

    with open(out_path, "w") as fh:
        fh.write("# Closed-form Paschen curve of the Townsend model\n")
        fh.write("#   U = B pd / [ln(A pd) - ln(ln(1 + 1/gamma))]\n")
        fh.write(
            "# Columns: pd (bar mm)  "
            + "  ".join(f"U_{g['label'].lower()} (kV)" for g in gases)
            + "\n"
        )
        for row in np.column_stack(columns):
            fh.write(
                "  ".join("nan" if np.isnan(v) else f"{v:.8e}" for v in row) + "\n"
            )
    print(f"Results written to: {out_path}")


if __name__ == "__main__":
    main(sys.argv[1])
