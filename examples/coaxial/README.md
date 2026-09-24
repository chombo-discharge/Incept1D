# Coaxial cylinders

A thin wire inside a grounded cylinder is the textbook corona geometry. The
field between the conductors is exact, `E(r) = U / (r ln(b/a))`, and it is
strongly non-uniform near the wire: ionization is confined to a thin layer
around the wire while most of the gap is attaching.

This example computes the inception voltage of dry air in such an
arrangement with `--field coaxial A B` (radii in mm). The radii fix the gap
at `d = b − a`, so the sweep is a **pressure sweep** at fixed geometry, and
the results are reported against `p`. Both polarities are computed:
`inner=positive` has the inner conductor as the anode.

## Run it

From the repository root, for a 1 mm wire inside a 10 mm cylinder, from 0.1
to 10 bar:

```bash
incept1d pdiv mechanisms/air/air_pancheshnyi.py \
    --field coaxial 1 10 \
    --pd-min 0.9 --pd-max 90 --pd-num 30 \
    --dx 10 400 0.01 \
    --write-to-file examples/coaxial/sim_a1mm.dat
```

`--pd-min`/`--pd-max` are asked for in p·d, so they are the pressure range
multiplied by `b − a = 9` mm. `incept1d pdiv` confirms the gap on startup
(`using d = L = 9 mm`), refuses `--p`, and reports every result against `p`.

`--dx 10 400 0.01` is finer than the default integration grid. Near a thin
wire the field falls off over a distance comparable to its radius, and the
default grid is up to 6 % off for `a = 0.25` mm at 10 bar; this setting is
within 0.2 % of a much finer one for every radius here.

## What you get

A figure of `U*` and `(E/N)*` against `p` for both polarities, a summary
table on stdout, and `sim_a1mm.dat` with the full curves. Column 1 is `p`
in bar; each polarity then has `p_bar`, `d_mm`, `EN_Td`, `U_kV`, `E_Vm`.
`EN_Td` and `E_Vm` are the *mean* field `U/d`. The field that decides
inception is the one at the wire, `E(a)/N = f(0) · EN_Td` with
`f(0) = (b − a) / (a ln(b/a))`; `docs/figures/coaxial.tex` plots it against
`p·a`.

## The documentation figure

The figure in the documentation (Examples → *Coaxial cylinders*) uses three
wires in the same 10 mm cylinder, each from 0.1 to 10 bar:

| `--field`          | `--pd-min` | `--pd-max` |
|--------------------|-----------:|-----------:|
| `coaxial 0.25 10`  | 0.975      | 97.5       |
| `coaxial 1 10`     | 0.9        | 90         |
| `coaxial 4 10`     | 0.6        | 60         |

`make -C docs/figures coaxial` runs all three and builds the figure.

## Things to try

* **Another wire.** Change `A`, and the pd range with it. Where the wires
  overlap in `p·a`, their wire fields `E(a)/N` coincide to within a few
  per cent, whatever `b/a` is.
* **Another outer radius.** The wire region decides inception, so the
  wire field should barely move: from `b = 10` to 25 mm it changes by less
  than 2 % for `a = 1` and 4 mm. A very large `b/a` with a thin wire at
  high pressure puts most of a long gap deep in attachment, which is where the
  determinant formulation runs out of range (Numerics → *The determinant*):
  with `coaxial 0.25 25` the positive curve develops a kink above about
  1 bar and then loses its root.
* **Look at the profile first.** `incept1d field --field coaxial 1 10`
  plots `f(ξ)` and prints `f(0)`, the ratio of the wire field to `U/d`.
