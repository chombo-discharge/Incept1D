# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Compare the self-sustained RFD field with the paper's own results, point by point.

Færder et al. publish their code and results on Zenodo
(doi:10.5281/zenodo.22983026, CC-BY 4.0).  Its folder
``selfsustainend-rfd-threshold-by-nme`` holds E_s(L) from their normal-mode
expansion (wide region) and from their time-dependent fluid model (both
radii).  This script reads those three files from the unpacked archive,
solves the same thresholds with Incept1D at the same L, and writes both to
one table.  The paper's data is theirs and is not part of this repository;
see examples/rfd/README.md.

Usage::

    python3 examples/rfd/paper_comparison.py \\
        PATH/TO/selfsustainend-rfd-threshold-by-nme examples/rfd/paper.dat
"""

import functools
import os
import sys

import numpy as np

from incept1d.constants import kB
from incept1d.fields import FieldDistribution
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import load_mechanism, read_json_configs
from incept1d.parallel import parallel_map, physical_cores
from incept1d.solver import riccati_criterion

MECHANISM = "mechanisms/air/rfd/air_rfd.py"
CONFIGS = "mechanisms/air/rfd/radius.json"
#: Sea level, n~ = 1, as in the paper: the authors' density 2.688e25 m^-3.
T = 273.15
P = 2.688e25 * kB * T / 1e5

#: The archive's files: normal modes in V/m, time-dependent model in units of E_t.
NME_WIDE = "eselfsus-vs-L-from-determinant-photons-and-positrons.npz"
NUM_WIDE = (
    "self-sustained-feedback-threshold-modeled-L=10-12000_dz=Ld200_logLscale"
    "_cone90deg_Rcloud=5000m.npz"
)
NUM_NARROW = (
    "self-sustained-feedback-threshold-modeled-L=10-12000_dz=Ld200_logLscale"
    "_cone90deg_Rcloud=100m.npz"
)
E_T = 284.0  # kV/m, the unit of the time-dependent results


def _read(folder, name):
    data = np.load(os.path.join(folder, name), allow_pickle=False)
    return np.asarray(data["Lcloud_data"], float), np.asarray(
        data["e_self_sus_data"], float
    )


def _incept1d(L):
    """Incept1D's E_s [kV/m] at each L [m], for the R >> L and R = 100 m cases."""
    cfgs = {c["label"]: c for c in read_json_configs([CONFIGS])}
    det = functools.partial(riccati_criterion, field_dist=FieldDistribution("uniform"))
    cases = [load_mechanism(MECHANISM, cfgs[k]) for k in ("R >> L", "R = 100 m")]

    def task(i, previous):
        mod = cases[i // len(L)]
        roots = find_all_breakdown_EN(
            P * L[i % len(L)], mod, P, T, first_only=True, det_fn=det
        )
        return mod._mod.E_field(roots[0], P, T) / 1e3 if roots else np.nan

    out = np.full(2 * len(L), np.nan)

    def store(i, value):
        out[i] = value

    parallel_map(task, len(out), physical_cores(), on_result=store)
    return out[: len(L)], out[len(L) :]


def main(folder, out_path):
    L, nme_wide = _read(folder, NME_WIDE)
    L_w, num_wide = _read(folder, NUM_WIDE)
    L_n, num_narrow = _read(folder, NUM_NARROW)
    if not (np.array_equal(L, L_w) and np.array_equal(L, L_n)):
        raise ValueError("the three files do not share their L values")
    ours_wide, ours_narrow = _incept1d(L)

    with open(out_path, "w") as fh:
        fh.write(
            "# Self-sustained RFD field, sea-level density: Incept1D and "
            "Faerder et al. (2026), doi:10.5281/zenodo.22983026, CC-BY 4.0\n"
        )
        fh.write(
            "# Columns: L (m)  E_s (kV/m): Incept1D [R >> L]  "
            "Incept1D [R = 100 m]  paper normal modes [R >> L]  "
            "paper time-dependent [R >> L]  paper time-dependent [R = 100 m]\n"
        )
        rows = zip(
            L, ours_wide, ours_narrow, nme_wide / 1e3, num_wide * E_T, num_narrow * E_T
        )
        for row in rows:
            fh.write("  ".join(f"{v:.10e}" for v in row) + "\n")
    print(f"Results written to: {out_path}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
