#!/usr/bin/env bash
# Compute the Incept1D curves for the IEC 60052 sphere-gap comparison
# (Docs: Examples -> "Sphere gaps: the IEC 60052 standard").
#
# One run per sphere diameter D (cm).  The sphere radius passed to --field is
# D/2 in mm, and the pd range covers the tabulated gap spacings at 1 bar.
# Each run writes sim_<D>cm.dat into $OUT (default: this directory).
#
# Usage, from the repository root:
#
#     bash examples/iec60052/run.sh            # all diameters
#     bash examples/iec60052/run.sh 10 50      # only D = 10 cm and 50 cm
#     OUT=/tmp/iec bash examples/iec60052/run.sh 200
#
# The full set takes tens of minutes.  PD_NUM (default 100) sets the number
# of pd points per curve.
set -euo pipefail
cd "$(dirname "$0")/../.."

MECH=mechanisms/air/air_pancheshnyi.py
CFG=mechanisms/air/nodetachment.json
OUT=${OUT:-examples/IEC60052}
PD_NUM=${PD_NUM:-100}
mkdir -p "$OUT"

# pd_max (bar mm) for each diameter, matching the extent of the IEC tables.
declare -A PDMAX=( [10]=55 [15]=90 [25]=150 [50]=300 [75]=400 [100]=550 [150]=800 [200]=1000 )

DIAMS=("$@")
if [ ${#DIAMS[@]} -eq 0 ]; then DIAMS=(10 15 25 50 75 100 150 200); fi

for D in "${DIAMS[@]}"; do
    R_mm=$(python3 -c "print($D * 10 / 2)")
    incept1d pdiv "$MECH" "$CFG" \
        --pd-min 1 --pd-max "${PDMAX[$D]}" --pd-num "$PD_NUM" \
        --field sphere-sphere "$R_mm" \
        --streamer-criterion 18 --no-plot \
        --write-to-file "$OUT/sim_${D}cm.dat"
done
