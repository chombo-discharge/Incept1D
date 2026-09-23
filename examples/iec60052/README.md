# IEC 60052 sphere gaps

`run.sh` computes Incept1D inception curves for the standard sphere gaps of
IEC 60052 and writes `sim_<D>cm.dat`, one file per sphere diameter.

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
