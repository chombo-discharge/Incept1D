#!/usr/bin/env bash
# Compute the classical Paschen curves of helium, argon and air with Incept1D
# (Docs: Examples -> "The classical Paschen curve").
#
# The mechanism is the textbook Townsend model -- alpha/p = A exp(-Bp/E), no
# attachment, no photoionization, ion-induced secondary emission only -- so the
# result has a closed form and this run is a validation rather than a
# prediction.  See mechanisms/paschen/paschen.py.
#
# Writes two files into $OUT (default: this directory):
#
#     sim.dat           the solver's inception curves
#     closed_form.dat   the analytic curve from the same coefficients
#
# Nothing third-party is involved, so the whole comparison ships with the
# repository.
#
# Usage, from the repository root:
#
#     bash examples/paschen/run.sh
#     OUT=/tmp/paschen bash examples/paschen/run.sh
#
# PD_NUM (default 120) sets the number of pd points.
set -euo pipefail
cd "$(dirname "$0")/../.."

MECH=mechanisms/paschen/paschen.py
CFG=mechanisms/paschen/gases.json
OUT=${OUT:-examples/paschen}
PD_NUM=${PD_NUM:-120}
mkdir -p "$OUT"

incept1d pdiv "$MECH" "$CFG" \
    --p 1.0 \
    --pd-min 3e-3 --pd-max 1e2 \
    --pd-num "$PD_NUM" \
    --no-plot \
    --write-to-file "$OUT/sim.dat"

python3 examples/paschen/closed_form.py "$OUT/closed_form.dat"
