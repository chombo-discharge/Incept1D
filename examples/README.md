# Examples

Reference data and reproduction scripts for the worked examples in the
documentation (`docs/source/examples/`).

| Directory   | Reference                                                    | Script   |
|-------------|--------------------------------------------------------------|----------|
| `paschen/`  | closed-form Paschen law, computed alongside (self-contained) | `run.sh` |
| `iec60052/` | IEC 60052 standard sphere-gap voltages (not included)        | `run.sh` |
| `electra/`  | Dakin et al., ELECTRA 32, Table C2 (not included)            | `run.sh` |

`paschen/` needs nothing external: it compares the solver against its own
analytic limit.

**The IEC 60052 and ELECTRA reference tables are not part of this repository.** They are copyrighted
by the IEC and by CIGRE respectively. Each directory's `README.md` explains
where to obtain them and what to name the file if you want the comparison
overlay; the figures build without them.

Each `run.sh` is run from the repository root and writes the Incept1D
results (`sim*.dat`) next to the reference tables (or into `$OUT`).  The
`sim*.dat` files are build products and are not committed.

The figures in the documentation are built from these scripts and the
pgfplots sources in `docs/figures/` — see `docs/figures/Makefile` and
`make figures` in `docs/`.
