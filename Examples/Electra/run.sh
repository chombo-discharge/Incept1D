#!/usr/bin/env bash
# Compute the Incept1D curve for the comparison against the Dakin et al.
# (ELECTRA No. 32) compilation of breakdown voltages in dry air
# (Docs: Examples -> "Quasi-uniform gaps: the ELECTRA Paschen curve").
#
# A sphere-sphere gap with R = 1000 mm at p = 1 bar is used as a nearly
# uniform field; the gap length is swept so that pd covers 5e-3 to 500 bar mm.
# Writes Sim.dat into $OUT (default: this directory).
#
# Usage, from the repository root:
#
#     bash Examples/Electra/run.sh
#     OUT=/tmp/electra PD_NUM=40 bash Examples/Electra/run.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

OUT=${OUT:-Examples/Electra}
PD_NUM=${PD_NUM:-100}
mkdir -p "$OUT"

incept1d pdiv mechanisms/Air/Air_Pancheshnyi.py mechanisms/Air/NoDetachment.json \
    --pd-min 5E-3 --pd-max 500 --pd-num "$PD_NUM" \
    --field sphere-sphere 1000 \
    --streamer-criterion 18 --no-plot \
    --dx 5 25 0.05 \
    --write-to-file "$OUT/Sim.dat"
