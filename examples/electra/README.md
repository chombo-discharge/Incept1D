# ELECTRA quasi-uniform gaps

`run.sh` computes the Incept1D inception curve for a quasi-uniform air gap and
writes `sim.dat`.

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
