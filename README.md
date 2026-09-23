# Incept1D

[![CI](https://github.com/chombo-discharge/Incept1D/actions/workflows/ci.yml/badge.svg)](https://github.com/chombo-discharge/Incept1D/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-github.io-blue)](https://chombo-discharge.github.io/Incept1D/)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue)](LICENSES/)
[![REUSE compliant](https://img.shields.io/badge/REUSE-compliant-brightgreen)](https://reuse.software)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21916821.svg)](https://doi.org/10.5281/zenodo.21916821)

**At what voltage does this gap break down?**

`Incept1D` answers that for a one-dimensional drift-reaction model of a gas
discharge gap, and it does so without assuming the answer is an avalanche of a
particular size. It solves the boundary-value problem for the coupled
species/photon system across the gap and finds the root of

$$\det \mathbf{Q}(\lambda) = 0 .$$

That criterion carries both electrode feedback — ion- and photon-induced
secondary emission — and volume feedback — photoionization, detachment, ion
conversion. It therefore stays valid in the regimes where the classical
Townsend and streamer criteria quietly stop being right: electronegative gases,
long gaps where negative-ion detachment matters, and short gaps where the
cathode rather than the avalanche decides the onset.

📖 **[Documentation](https://chombo-discharge.github.io/Incept1D/)** — theory,
numerics, configuration reference and worked examples.

## Features

|  | |
|---|---|
| **Any gas** | The chemistry lives in a *mechanism file* outside the package, so a new gas is a new file rather than a patch to the solver. Dry air ships with the code. |
| **Any gap** | Uniform, sphere-plane and sphere-sphere gaps are analytic. A **tabulated field line** from any electrostatic solver can be used directly, curvature included. |
| **Sensitivity in one run** | JSON configuration files overlay variants — cross-section databases, cathode yields, individual reactions scaled or switched off — side by side in one figure. |
| **Beyond the threshold** | Above the inception voltage, the temporal growth rate $\lambda$ says how *fast* the discharge grows, not merely that it does. |
| **Checked against algebra** | A reduced limit of the model has a closed-form solution, and the test suite asserts the solver reproduces it to a relative error below $10^{-9}$. |
| **Hands off to 3-D** | Transport and rate tables export to [chombo-discharge](https://github.com/chombo-discharge/chombo-discharge), so the same chemistry drives the 3-D simulation. |

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
      1.0000e-02           1      0.01     1521.4350        0.3761      3.7610e+07
      1.3434e-02           1  0.013434     1136.0110        0.3773      2.8082e+07
      1.8047e-02           1 0.0180472      892.1123        0.3980      2.2053e+07
      ...
      1.0000e+03           1      1000       98.4322     2433.2480      2.4332e+06
```

The minimum of that curve is the Paschen minimum — 373 V at
$p\cdot d = 0.0113$ bar·mm for this dry-air scheme.

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

### Using your own field line

Export $|E|$ along a field line from any electrostatic solver — as `s |E|`,
`x y z |E|` or `x y z Ex Ey Ez` — and hand it over. Only the *shape* of the
profile enters the solve, so the field units of the file never reach the
answer, and the arc length fixes the gap:

```bash
incept1d pdiv mechanisms/air/air_pancheshnyi.py \
    --field fieldline line.csv mm --fieldline-voltage 100
```

`--fieldline-voltage` declares the excitation the line was computed at, which
is what lets the result be reported as a factor on it — a reported
`U*/U_applied` of 0.37 means the gap breaks down at 37 % of the voltage you
solved for. Check a file
before spending a sweep on it with `incept1d field --field fieldline line.csv mm`,
which plots the profile and identifies the field units.

See the [quick start](https://chombo-discharge.github.io/Incept1D/introduction/quickstart.html)
for a guided walk-through.

## Repository layout

```
src/incept1d/        the Python package (solver, inception curves, CLI)
mechanisms/air/      dry-air mechanism family, swarm data, configurations
examples/            scripts reproducing the worked examples
tests/               the test suite
docs/                Sphinx documentation sources
```

## Tests

```bash
pip install -e '.[dev]'
pytest                    # the whole suite
pytest -m "not slow"      # skips the tests that load real swarm data
```

The suite is built in three layers, strongest first: verification against
closed-form solutions, unit tests of each module, and invariants that hold
whatever the numbers are. Invariants are preferred to pinned values — a
wrong-but-self-consistent curve passes a pinned value, which is how a real
root-finding bug survived here once.

## Documentation

Built with Sphinx and deployed to GitHub Pages on every push to `main`. To
build it locally:

```bash
pip install -e '.[docs]'
make -C docs html          # -> docs/build/html/index.html
```

Figures are **built, not shipped**: `make -C docs figures` runs the solver over
the worked examples and compiles the pgfplots sources. Use
`PD_NUM=10 make -C docs html` for a fast, coarse build.

## Contributing

Pull requests are welcome. Please run `pre-commit install` once per clone; the
hooks run `black`, `flake8`, `reuse lint` and a Sphinx build before every
commit. Changes that affect the physics or numerics should cite the
corresponding equation in `docs/source/theory/` and show before/after inception
curves — see the
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

Published standards, journal tables and datasets used by the worked examples
are copyrighted by their publishers and are **not distributed here**. Those
examples commit the *calculation* and document the expected input file, so a
reader who holds a copy can supply it and reproduce the comparison — see
`examples/iec60052/README.md` and `examples/electra/README.md`.

## Authors

Developed at [SINTEF Energy Research](https://www.sintef.no/en/sintef-energy/).
