# Examples

Worked examples for the documentation (`docs/source/examples/`). Each
directory has a `README.md` that says what the example shows and gives the
command to run it, from the repository root. Change the command and run it
again: that is the point.

| Directory    | What it shows                                              | Reference                                 |
|--------------|------------------------------------------------------------|-------------------------------------------|
| `paschen/`   | The solver against the closed-form Paschen law             | computed alongside (self-contained)       |
| `fieldline/` | Inception along a tabulated field line, and its units      | none needed (self-contained)              |
| `coaxial/`   | Corona inception on a wire in a cylinder                   | none needed (self-contained)              |
| `iec60052/`  | Standard sphere gaps                                       | IEC 60052 tables (not included)           |
| `electra/`   | Quasi-uniform gaps over five decades of `pd`               | Dakin et al., ELECTRA 32 (not included)   |

**The IEC 60052 and ELECTRA reference tables are not part of this
repository.** They are copyrighted by the IEC and by CIGRE respectively.
Each of those READMEs explains where to obtain them and what to name the
file if you want the comparison overlay; the figures build without them.

The commands write their results (`sim*.dat`) into the example's own
directory. Those files are outputs and are not committed.

The figures in the documentation are built from the same commands (run
headless, into `docs/build/figures/`) and the pgfplots sources in
`docs/figures/` — see `docs/figures/Makefile` and `make figures` in `docs/`.
