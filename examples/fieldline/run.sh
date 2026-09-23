#!/usr/bin/env bash
# Inception along a tabulated field line, and what the file's units do to it
# (Docs: Examples -> "A tabulated field line").
#
# The same 20 mm line is written twice by make_line.py -- once in SI (m, V/m)
# and once in engineering units (mm, kV/mm) -- and solved twice.  Only the
# shape of the profile enters the solve, so both runs must agree exactly.
# --fieldline-voltage declares the excitation the line was computed at, which
# is what makes U*/U_applied meaningful in either unit system.
#
# Usage, from the repository root:
#
#     bash examples/fieldline/run.sh
#     OUT=/tmp/fieldline bash examples/fieldline/run.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

MECH=mechanisms/air/air_pancheshnyi.py
OUT=${OUT:-examples/fieldline}
PD_NUM=${PD_NUM:-40}
mkdir -p "$OUT"

python3 examples/fieldline/make_line.py "$OUT"

# The gap length is the arc length of the line, so the sweep is in pressure.
incept1d pdiv "$MECH" \
    --field fieldline "$OUT/line_si.csv" m \
    --fieldline-voltage 100 \
    --pd-min 1 --pd-max 100 --pd-num "$PD_NUM" \
    --no-plot \
    --write-to-file "$OUT/sim_si.dat"

incept1d pdiv "$MECH" \
    --field fieldline "$OUT/line_engineering.csv" mm \
    --fieldline-voltage 100 \
    --pd-min 1 --pd-max 100 --pd-num "$PD_NUM" \
    --no-plot \
    --write-to-file "$OUT/sim_engineering.dat"
