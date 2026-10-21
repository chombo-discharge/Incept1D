# Relativistic feedback discharge

The field at which relativistic feedback in a thundercloud becomes
self-sustaining, against the length `L` of the high-field region. The
mechanism `mechanisms/air/rfd/air_rfd.py` is the relativistic fluid model
of Færder et al.: runaway electrons, positrons, and one forward and one
backward photon group, with every coefficient fitted to GEANT4 simulations
in that paper. Nothing is emitted at either end of the region, so the
feedback closes entirely in the gas.

The two horizontal radii of the paper — a region much wider than it is
long, and one 100 m across — are the two configurations in
`mechanisms/air/rfd/radius.json`, so one command produces both curves. At
fixed pressure, `pd` is `L`: the range below is `L` = 25 m to 5 km at sea
level.

## Run it

From the repository root:

```bash
incept1d pdiv mechanisms/air/rfd/air_rfd.py mechanisms/air/rfd/radius.json \
    --p 1.01325 --T 293.15 \
    --pd-min 2.5e4 --pd-max 5e6 --pd-num 60 \
    --write-to-file examples/rfd/sim.dat
```

and, for the normal-mode solution of the same model:

```bash
python3 examples/rfd/normal_modes.py examples/rfd/normal_modes.dat
```

## What you get

`sim.dat` with the self-sustained field `E_Vm` in column 6 for the wide
region and column 16 for the narrow one, and the gap `d_mm` in columns 3
and 13. `normal_modes.dat` has `L` and the normal-mode field for each
radius. The two agree to solver tolerance.

The field falls from conventional breakdown, 3 MV/m near `L` = 23 m,
toward 305 kV/m at kilometre scale, where positrons start to make runaway
electrons; a narrow region needs a slightly higher field, because photons
escape it sideways.

`make -C docs/figures rfd` runs both and builds the documentation figure
(Examples → *Relativistic feedback discharge*).

## Comparison with the paper

The figure above checks the solver against its own model. To compare with
the paper's published results, get them from the authors: their code and
data are on Zenodo, [doi:10.5281/zenodo.22983026](https://doi.org/10.5281/zenodo.22983026)
(CC-BY 4.0), and are not part of this repository. Unpack the archive and
run

```bash
python3 examples/rfd/paper_comparison.py \
    PATH/TO/Faerder2026b_openresearch/selfsustainend-rfd-threshold-by-nme \
    examples/rfd/paper.dat
```

which reads their normal-mode `E_s(L)` (wide region) and their
time-dependent fluid-model results (both radii), solves the same thresholds
with Incept1D at the same `L`, and writes both side by side. Then
`make -C docs/figures rfdpaper` builds the comparison figure, and the test
suite checks Incept1D against their normal modes. `paper.dat` is listed in
`.gitignore`: please keep it that way.

## Things to try

* **Altitude.** Lower `--p`: every coefficient scales with density, so the
  field at threshold scales with it and `E/N` against `pd` does not change.
* **Another radius.** Add a configuration to `radius.json` with your own
  `photon_attenuation_fwd` and `photon_attenuation_bck` (in m at sea
  level).
* **A non-uniform field.** Add `--field sphere-plane 1e5`, a sphere of
  100 m radius: the coefficients are then applied locally along the field.
