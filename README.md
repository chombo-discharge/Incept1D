# Incept1D

[![Documentation](https://github.com/chombo-discharge/Incept1D/actions/workflows/docs.yml/badge.svg)](https://chombo-discharge.github.io/Incept1D/)

`Incept1D` computes the **inception (breakdown) condition** of a one-dimensional
drift-reaction model of a gas discharge gap, including negative-ion transport
and detachment, ion conversion, two-stream photoionization, and cathode
secondary emission.

Rather than asking whether an electron avalanche reaches a critical size, the
code solves the full boundary-value problem for the augmented species/photon
system across the gap and evaluates the determinant criterion

$$\det \mathbf{Q}(\lambda) = 0,$$

whose root in $E/N$ at fixed $p\cdot d$ is the inception field. The criterion
accounts for both electrode feedback (ion- and photon-induced secondary
emission) and volume feedback (photoionization, detachment), so it remains
valid where the classical Townsend and streamer criteria do not.

📖 **[Documentation](https://chombo-discharge.github.io/Incept1D/)** — theory,
numerics, module reference and worked examples.

## What it does

Given a plasma-chemistry *mechanism file*, `Incept1D` can

- compute **inception curves** (generalized Paschen curves): inception voltage
  and reduced field vs. $p\cdot d$, for uniform, sphere-plane, sphere-sphere or
  **tabulated field-line** geometries, including curved field lines exported
  from a 3-D electrostatic solver;
- inspect the **local eigenvalues** of the reaction-transport matrix
  $\mathbf{R}\mathbf{V}^{-1}$ vs. $E/N$;
- compare the full criterion against the classical **ionization integral** and
  the streamer criterion;
- compute the **temporal growth rate** $\lambda$ of the discharge above the
  inception voltage;
- export transport and rate-coefficient tables for the 3-D plasma solver
  [chombo-discharge](https://github.com/chombo-discharge/chombo-discharge).

## Installation

Requires Python ≥ 3.10. An editable install is recommended, so the `incept1d`
command always runs the sources in your checkout:

```bash
git clone https://github.com/chombo-discharge/Incept1D.git
cd Incept1D
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

Optional extras: `pip install -e '.[docs]'` for Sphinx, `pip install -e '.[dev]'`
for pytest, pre-commit, black and flake8.

## Quick start

```bash
# Inception curve for dry air at 1 bar, uniform gap
incept1d pdiv mechanisms/air/air_pancheshnyi.py --p 1 --pd-min 1e-2 --pd-max 1e3
```

```
     pd (bar·mm)     p (bar)    d (mm)      E/N (Td)        U (kV)         E (V/m)
      1.0000e-02           1      0.01     1521.4349        0.3761      3.7610e+07
      1.5999e-02           1     0.016      979.5071        0.3874      2.4213e+07
      2.5595e-02           1    0.0256      700.6243        0.4433      1.7319e+07
      ...
      3.0888e+02           1     308.9       98.9278      755.3766      2.4455e+06
```

The minimum of that curve is the Paschen minimum — here 387 V at
$p\cdot d \approx 0.016$ bar·mm for dry air.

Every tool is a subcommand of `incept1d` (`incept1d --help` lists them all):

| Command | Purpose |
|---|---|
| `incept1d pdiv` | Inception curve PDIV($p\cdot d$) — the main command |
| `incept1d eigenvalues` | Eigenvalues of the local transport matrix $\mathbf{R}\mathbf{V}^{-1}$ |
| `incept1d ionization` | Ionization integrals vs. applied voltage |
| `incept1d growth` | Temporal growth rate $\lambda$ above inception |
| `incept1d field` | Plot the normalised field profile $f(\xi)$ of a gap |
| `incept1d chombo` | Transport/rate tables for `chombo-discharge` |

A gap geometry is chosen with `--field`, and JSON configuration files may follow
the mechanism to overlay sensitivity variants in one figure:

```bash
incept1d pdiv mechanisms/air/air_pancheshnyi.py mechanisms/air/databases.json \
    --field sphere-plane 50 --d 20
```

See the [quick start](https://chombo-discharge.github.io/Incept1D/introduction/quickstart.html)
for a guided walk-through.

## Repository layout

```
src/incept1d/        the Python package (solver, inception curves, CLI)
mechanisms/air/      dry-air mechanism family, swarm data, configurations
examples/            scripts reproducing the worked examples
docs/                Sphinx documentation sources
```

## Documentation

The documentation is built with Sphinx and deployed to GitHub Pages on every
push to `main`. To build it locally:

```bash
pip install -e '.[docs]'
make -C docs html          # -> docs/build/html/index.html
```

Figures are **built, not shipped**: `make -C docs figures` runs the solver over
the worked examples and compiles the pgfplots sources. Use `PD_NUM=10 make -C docs html`
for a fast, coarse build.

## Contributing

Pull requests are welcome. Please run `pre-commit install` once per clone; the
hooks run `black`, `flake8` and a Sphinx build before every commit. Changes that
affect the physics or numerics should cite the corresponding equation in
`docs/source/theory/` and show before/after inception curves — see the
[contributing guide](https://chombo-discharge.github.io/Incept1D/maintenance/contributing.html)
and the pull request template.

## Licensing

`Incept1D` is licensed under **GPL-3.0-or-later** and is
[REUSE](https://reuse.software) compliant: every file declares its copyright
holder and licence, either through an inline SPDX header or an entry in
[`REUSE.toml`](REUSE.toml). Run `reuse lint` to check.

The electron swarm data under `mechanisms/air/*.txt` is **not** covered by the
project licence. It is retrieved from the [LXCat](https://www.lxcat.net)
open-access database and redistributed verbatim, with its original headers
intact, under [`LicenseRef-LXCat`](LICENSES/LicenseRef-LXCat.txt); copyright
rests with the contributing databases (IST-Lisbon, Phelps, Biagi, Morgan,
TRINITI, Viehland). If you use it, cite the database named in the header of
the file you used — that is the reference format its authors ask for.

The IEC 60052 sphere-gap voltages and the CIGRE ELECTRA Table C2 breakdown
curve used by the worked examples are copyrighted by their publishers and are
**not distributed here**. The example figures are built from Incept1D's own
computed curves; if you hold a copy of either reference you can drop it into
the example directory to overlay it — see `examples/iec60052/README.md` and
`examples/electra/README.md`.

## Authors

Developed at [SINTEF Energy Research](https://www.sintef.no/en/sintef-energy/).
