# Examples

Reference data and reproduction scripts for the worked examples in the
documentation (`Docs/source/Examples/`).

| Directory   | Reference data                                               | Script   |
|-------------|--------------------------------------------------------------|----------|
| `IEC60052/` | IEC 60052 standard sphere-gap voltages, `IEC60052_<D>cm.dat` | `run.sh` |
| `Electra/`  | Dakin et al., ELECTRA 32, Table C2 (`TableC2_Air.dat`)       | `run.sh` |

Each `run.sh` is run from the repository root and writes the Incept1D
results (`Sim*.dat`) next to the reference tables (or into `$OUT`).  The
`Sim*.dat` files are build products and are not committed.

The figures in the documentation are built from these scripts and the
pgfplots sources in `Docs/figures/` — see `Docs/figures/Makefile` and
`make figures` in `Docs/`.
