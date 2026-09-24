# Sphere gaps: the IEC 60052 standard

IEC 60052 tabulates the disruptive-discharge voltage of standard sphere
gaps — two equal spheres of diameter `D` between 2 cm and 200 cm, one of
them earthed — against gap spacing, for voltage calibration. The field
between two spheres is known in closed form, so these tables test the
inception criterion in a *non-uniform* field over a wide range of `pd`
without any numerical field calculation.

This example computes the inception curves of dry air for the standard
diameters with `--field sphere-sphere R` (radius in mm). The pressure is
held at 1 bar and the gap spacing is swept, so `pd` in bar·mm equals the
spacing in mm. The full air scheme and a no-detachment variant are both
computed (`mechanisms/air/pancheshnyi/nodetachment.json` holds the two configurations),
and the streamer criterion `∫α dx = 18` is solved alongside for comparison.

## Run it

From the repository root, for the largest sphere, `D = 200` cm:

```bash
incept1d pdiv mechanisms/air/pancheshnyi/air_pancheshnyi.py mechanisms/air/pancheshnyi/nodetachment.json \
    --pd-min 1 --pd-max 1000 --pd-num 100 \
    --field sphere-sphere 1000 \
    --streamer-criterion 18 \
    --write-to-file examples/iec60052/sim_200cm.dat
```

This takes a few minutes. The sphere-sphere gap is symmetric, so both
polarities are written but are identical.

## What you get

A figure of `U*` and `(E/N)*` against `pd`, one curve per configuration and
one for the streamer criterion, and `sim_200cm.dat`. The column layout is
listed in the file header; Examples → *Sphere gaps: the IEC 60052 standard*
in the documentation explains which columns the figure plots.

## The other diameters

The documentation figure uses all eight diameters, each with a `pd` range
matched to its table. Change `--field sphere-sphere`, `--pd-max` and the
output name:

| `D` (cm) | `--field sphere-sphere` | `--pd-max` |
|---------:|------------------------:|-----------:|
| 10       | 50                      | 55         |
| 15       | 75                      | 90         |
| 25       | 125                     | 150        |
| 50       | 250                     | 300        |
| 75       | 375                     | 400        |
| 100      | 500                     | 550        |
| 150      | 750                     | 800        |
| 200      | 1000                    | 1000       |

`make -C docs/figures iec60052` runs the whole set and builds the figure.

## Reference data is not included

The disruptive-discharge voltages tabulated in

> IEC 60052:2002, *Voltage measurement by means of standard air gaps*

are copyrighted by the International Electrotechnical Commission and are **not
distributed with this repository**. The comparison figure in the documentation
is therefore built from the computed curves alone.

If you hold a copy of the standard, you can reproduce the full comparison
locally. Create one file per diameter here, named `iec60052_<D>cm.dat`, with
two whitespace-separated columns and any number of `#` comment lines:

```
# Sphere diameter D = 10 cm, one sphere earthed, dry air at 20 C and 1013 mbar.
# Columns: p*d (bar.mm, = gap spacing in mm at 1 bar)   voltage (kV, peak)
   5.0   16.8
   6.0   19.9
```

`docs/figures/iec60052.tex` guards every reference overlay with
`\IfFileExists`, so the figure picks the files up automatically on the next
`make -C docs figures` and ignores the ones you do not provide.

These filenames are listed in `.gitignore`: please keep it that way.
