# Quasi-uniform gaps: the ELECTRA Paschen curve

Dakin *et al.* (ELECTRA No. 32) compiled breakdown measurements in dry air
from several laboratories and gap types into one Paschen curve, spanning
`pd` from about 3e-3 to 400 bar·mm. This example compares the inception
criterion against that compilation over five decades of `pd`, from the
Paschen minimum to the regime where negative-ion detachment dominates.

The compilation mixes geometries, so a single nearly uniform one stands in
for all of them: a sphere-sphere gap with `R = 1000` mm at 1 bar, sweeping
the gap length. The full air scheme, the no-detachment variant and the
streamer criterion `∫α dx = 18` are computed together.

## Run it

From the repository root:

```bash
incept1d pdiv mechanisms/air/air_pancheshnyi.py mechanisms/air/nodetachment.json \
    --pd-min 5E-3 --pd-max 500 --pd-num 100 \
    --field sphere-sphere 1000 \
    --streamer-criterion 18 \
    --dx 5 25 0.05 \
    --write-to-file examples/electra/sim.dat
```

At small `pd` the propagator is cheap, so the adaptive grid budget is
reduced (`--dx 5 25 0.05`) to keep the run short.

## What you get

A figure of `U*` and `(E/N)*` against `pd`, and `sim.dat` with the `EN_Td`
and `U_kV` columns of each configuration and of the streamer criterion.
`make -C docs/figures electra` builds the documentation figure from it.

## Reference data is not included

The breakdown voltages of

> T. W. Dakin et al., *Breakdown of gases in uniform fields: Paschen curves for
> nitrogen, air and sulphur hexafluoride*, ELECTRA No. 32, pp. 61-82

are copyrighted by CIGRE and are **not distributed with this repository**. The
comparison figure in the documentation is therefore built from the computed
curve alone.

If you have access to the article, you can reproduce the full comparison
locally. Create `tablec2_air.dat` here with two whitespace-separated columns
and any number of `#` comment lines:

```
# Columns: p*d (bar.mm)   breakdown voltage (kV, crest value)
2.6e-3    0.520
3.0e-3    0.440
```

`docs/figures/electra.tex` guards every reference overlay with
`\IfFileExists`, so the figure picks the file up automatically on the next
`make -C docs figures`.

This filename is listed in `.gitignore`: please keep it that way.
