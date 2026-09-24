# A tabulated field line

Any geometry beyond spheres, planes and cylinders needs an electrostatic
solve, and what comes back is a table of `|E|` along a field line. This
example runs the inception criterion along such a line with
`--field fieldline FILE UNIT`, and shows what the units of that file do and
do not affect.

The line is synthetic so that every number can be checked by hand: 20 mm
long, with `|E|` falling linearly by a factor of two, scaled to an
excitation of 100 kV. `make_line.py` writes it twice, as `line_si.csv`
(m, V/m) and `line_engineering.csv` (mm, kV/mm). Only the *shape* of the
profile enters the solve, so both must give identical results. The length
unit is declared on the command line; the field unit cannot be, and
`--fieldline-voltage` declares the excitation the line was computed at so
that `U*/U_applied` can be reported.

## Run it

From the repository root, first write the two files:

```bash
python3 examples/fieldline/make_line.py examples/fieldline
```

then solve along the line in engineering units:

```bash
incept1d pdiv mechanisms/air/pancheshnyi/air_pancheshnyi.py \
    --field fieldline examples/fieldline/line_engineering.csv mm \
    --fieldline-voltage 100 \
    --pd-min 1 --pd-max 100 --pd-num 40 \
    --write-to-file examples/fieldline/sim_engineering.dat
```

The arc length fixes the gap at `d = L = 20` mm, so this is a **pressure
sweep** (0.05 to 5 bar); `--p` is refused, and the results are reported
against `p`.

## What you get

A figure of `U*` and `(E/N)*` against `p`, one curve per polarity
(`start=positive`: the first tabulated point is the anode), a table with a
`U*/U_applied` column, and `sim_engineering.dat`.

The check: run the same command on `line_si.csv m`, writing
`sim_si.dat`. The two output files must agree to the last digit.

## Things to try

* **Check a file before solving.**
  `incept1d field --field fieldline examples/fieldline/line_engineering.csv mm --fieldline-voltage 100`
  plots the profile and says which field unit the line integral implies.
* **Your own line.** Export `s |E|`, `x y z |E|` or `x y z Ex Ey Ez` from
  any field solver; see Python modules → *Gap geometry* for the layouts.
* **Reverse the line.** Exporting from the other electrode swaps the
  polarity labels, which confirms which end is which.
