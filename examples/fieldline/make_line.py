# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Write the same synthetic field line in two unit conventions.

The line is a 20 mm path along which |E| falls linearly by a factor of two,
as it would along the axis of a mildly divergent gap, and it is scaled to an
excitation of 100 kV:

    <|E|> = 100 kV / 20 mm = 5 kV/mm,   |E| : 6.667 -> 3.333 kV/mm

A linear ramp is used because the integral is exact,

    integral of |E| ds = L * <|E|> = 100 kV,

so anything the solver reports about the scale of the input can be checked
by hand.  The point of writing it twice is that the two files describe the
*same* line: one in SI (m, V/m), one in the units an engineering field
solver tends to export (mm, kV/mm).  Only the shape enters the solve, so
every physical result must be identical.

Usage:  python3 examples/fieldline/make_line.py OUTDIR
"""

import os
import sys

import numpy as np

#: Arc length of the line.
LENGTH_MM = 20.0

#: Excitation the line is scaled to.
VOLTAGE_KV = 100.0

#: Ratio of the field at the two ends.
RATIO = 2.0

#: Number of tabulated points.
N_POINTS = 41


def build():
    """Return (s_mm, E_kV_per_mm) for the linear ramp."""
    s_mm = np.linspace(0.0, LENGTH_MM, N_POINTS)
    mean = VOLTAGE_KV / LENGTH_MM
    # A linear ramp with this ratio and this mean.
    hi = 2.0 * mean * RATIO / (1.0 + RATIO)
    lo = 2.0 * mean / (1.0 + RATIO)
    return s_mm, np.linspace(hi, lo, N_POINTS)


def main(outdir):
    s_mm, E_kV_per_mm = build()

    np.savetxt(
        os.path.join(outdir, "line_si.csv"),
        np.c_[s_mm * 1e-3, E_kV_per_mm * 1e6],
        header="s_m            E_V_per_m",
        fmt="%.9e",
    )
    np.savetxt(
        os.path.join(outdir, "line_engineering.csv"),
        np.c_[s_mm, E_kV_per_mm],
        header="s_mm           E_kV_per_mm",
        fmt="%.9e",
    )

    integral = np.trapezoid(E_kV_per_mm, s_mm)
    print(f"Line length      {LENGTH_MM:g} mm, {N_POINTS} points")
    print(f"|E| range        {E_kV_per_mm[0]:.4f} -> {E_kV_per_mm[-1]:.4f} kV/mm")
    print(f"∫|E| ds          {integral:.6f} kV   (exact: {VOLTAGE_KV:g} kV)")
    print(f"Written to       {outdir}/line_si.csv, {outdir}/line_engineering.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__))
