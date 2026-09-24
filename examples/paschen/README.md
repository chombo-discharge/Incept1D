# The classical Paschen curve

A verification, not a prediction. The mechanism
`mechanisms/paschen/paschen.py` is the textbook Townsend model — one
ionizing reaction with `α/p = A exp(−Bp/E)`, no attachment, no
photoionization, and electrons released from the cathode only by arriving
positive ions. That reduced model has a closed-form inception condition,
Paschen's law,

    U = B (pd) / [ ln(A pd) − ln(ln(1 + 1/γ)) ],

so running the full solver on it checks the solver against algebra. The
coefficients for helium, argon and air are three configurations in
`mechanisms/paschen/gases.json`, so one command produces all three curves.
Nothing third-party is involved.

## Run it

From the repository root:

```bash
incept1d pdiv mechanisms/paschen/paschen.py mechanisms/paschen/gases.json \
    --p 1.0 \
    --pd-min 3e-3 --pd-max 1e2 --pd-num 120 \
    --write-to-file examples/paschen/sim.dat
```

and, for the closed form from the very same coefficients:

```bash
python3 examples/paschen/closed_form.py examples/paschen/closed_form.dat
```

## What you get

A figure of `U*` and `(E/N)*` against `pd` for the three gases, and
`sim.dat` with one column group per gas (`U_kV` is column 5 for helium, 15
for argon, 25 for air). `closed_form.dat` has `pd` and the analytic voltage
for each gas. The two agree to plotting accuracy; the solver curves stop at
the left-branch asymptote, where no breakdown is possible at any voltage.

`make -C docs/figures paschen` runs both and builds the documentation
figure (Examples → *The classical Paschen curve*).

## Things to try

* **Another γ.** Edit `gamma0` in `gases.json`: the minimum moves as
  `(pd)_min ∝ ln(1 + 1/γ)`, in the solver and in the closed form alike.
* **A non-uniform gap.** Add `--field sphere-plane 50`. There is no closed
  form any more, and the two polarities separate.
