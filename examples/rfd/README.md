# Relativistic feedback discharge

The field at which relativistic feedback in a thundercloud becomes
self-sustaining, against the length `L` of the high-field region. The
mechanism `mechanisms/air/rfd/air_rfd.py` is the relativistic fluid model
of Færder et al.: runaway electrons, positrons, and one forward and one
backward photon group, with every coefficient fitted to GEANT4 simulations
in that paper. Nothing is emitted at either end of the region, so the
feedback closes entirely in the gas.

The region has no electrodes: its two ends are open, nothing enters, and
whatever leaves is lost, which Incept1D expresses as a gap with every
emission yield zero. Its horizontal radius `R` enters only through the
photon attenuation lengths, fitted by the authors for `R` = 5 km (which the
paper calls R ≫ L) and 100 m; the two are the configurations in
`mechanisms/air/rfd/radius.json`, so one command produces both curves. The
paper's results are all at sea level, which is the authors' sea-level
density 2.688e25 m⁻³: `--p 1.01371 --T 273.15`. At fixed pressure, `pd` is
`L`: the range below is `L` = 10 m to 12 km, as in the paper.

## Run it

From the repository root:

```bash
incept1d pdiv mechanisms/air/rfd/air_rfd.py mechanisms/air/rfd/radius.json \
    --p 1.01371 --T 273.15 \
    --pd-min 1e4 --pd-max 1.2e7 --pd-num 80 \
    --write-to-file examples/rfd/sim.dat
```

## What you get

`sim.dat` with the self-sustained field `E_Vm` in column 6 for `R` = 5 km
and column 16 for `R` = 100 m, and the gap `d_mm` in columns 3
and 13.

The field falls from conventional breakdown, 3 MV/m near `L` = 23 m,
toward 305 kV/m at kilometre scale, where positrons start to make runaway
electrons; the narrower region needs a slightly higher field, because photons
escape it sideways.

`make -C docs/figures rfd` runs the same command and builds the
documentation figure (Examples → *Relativistic feedback discharge*).

## Comparison with the paper

The paper's published results are not part of this repository. The
authors' code and data are on Zenodo,
[doi:10.5281/zenodo.22983026](https://doi.org/10.5281/zenodo.22983026)
(CC-BY 4.0). Unpack the archive and run

```bash
python3 examples/rfd/paper_comparison.py \
    PATH/TO/Faerder2026b_openresearch/selfsustainend-rfd-threshold-by-nme \
    examples/rfd/paper.dat
```

which reads their normal-mode `E_s(L)` (`R` = 5 km) and their
time-dependent fluid-model results (both radii), solves the same thresholds
with Incept1D at the same `L`, and writes both side by side. With
`paper.dat` in place, `make -C docs/figures rfd` adds their results and the
difference to the figure, and the test suite checks Incept1D against their
normal modes. `paper.dat` is listed in `.gitignore`: please keep it that
way.

## Things to try

* **Altitude.** Lower `--p`. The model is the authors': the avalanche
  thresholds scale with the relative density, and every other coefficient
  is its sea-level function of the field, times the relative density. The
  paper publishes sea-level results only, so nothing here is checked
  against it.
* **Another radius.** Add a configuration to `radius.json` with your own
  `photon_attenuation_fwd` and `photon_attenuation_bck` (in m at sea
  level).
* **A non-uniform field.** Add `--field sphere-plane 1e5`, a sphere of
  100 m radius: the coefficients are then applied locally along the field.
