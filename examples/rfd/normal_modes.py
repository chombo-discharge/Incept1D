# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Tabulate the self-sustained RFD field by normal modes, for both radii.

Evaluates ``normal_mode_threshold`` of ``mechanisms/air/rfd/air_rfd.py``,
the steady-state solution of the model by its normal modes, with the
attenuation lengths of each configuration in ``radius.json``, so that the
figure in the documentation compares the solver against its own model
rather than against an external table.

Usage::

    python3 examples/rfd/normal_modes.py OUTPUT.dat
"""

import sys

import numpy as np

from incept1d.mechanism import load_mechanism, read_json_configs

MECHANISM = "mechanisms/air/rfd/air_rfd.py"
CONFIGS = "mechanisms/air/rfd/radius.json"


def main(out_path, n=25):
    # Sparse on purpose: plotted as markers over the solver's curve.
    L = np.logspace(np.log10(25.0), np.log10(5e3), n)  # m
    labels, columns = [], [L]
    for cfg in read_json_configs([CONFIGS]):
        raw = load_mechanism(MECHANISM, cfg)._mod
        labels.append(cfg["label"])
        columns.append(np.array([raw.normal_mode_threshold(x) / 1e3 for x in L]))

    with open(out_path, "w") as fh:
        fh.write("# Self-sustained RFD field by normal modes, sea-level density\n")
        fh.write(
            "# Columns: L (m)  " + "  ".join(f"E_s (kV/m) [{s}]" for s in labels) + "\n"
        )
        for row in np.column_stack(columns):
            fh.write(
                "  ".join("nan" if np.isnan(v) else f"{v:.8e}" for v in row) + "\n"
            )
    print(f"Results written to: {out_path}")


if __name__ == "__main__":
    main(sys.argv[1])
