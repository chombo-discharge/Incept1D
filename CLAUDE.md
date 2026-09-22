# CLAUDE.md

Guidance for Claude Code (or any future contributor) working in this repository.

## What this program does

Incept1D computes the **inception (breakdown) condition** for a 1-D
drift-reaction model of a discharge gap (electrons + positive/negative ions +
two-stream photoionization), and derived quantities (inception curves PDIV(p·d),
ionization-integral curves, temporal growth rates, transport-coefficient
tables for 3-D simulation codes).

**Before touching any of the physics code, read the Theory chapter of the
documentation, `Docs/source/Theory/`.**  Those pages are the single source
of truth for the equations the code implements (the manuscript they were
derived from is not part of the repository and must not be referenced
from the docs):

- `Theory/Overview.rst`, eq. `eq_drift_reaction`: the governing drift-reaction PDE.
- `Theory/Transport.rst`, eq. `eq_augmented_ode`: the augmented first-order
  ODE `∂_x θ = A_aug θ` with block matrix `A_aug = [[A,B,B],[C,-D,0],[-C,0,D]]`
  — this is exactly `incept1d.solver._build_A_aug`; eq. `eq_theta_soln` is the
  propagator `M(d)`.
- `Theory/Photoionization.rst`, eq. `eq_two_stream`: the two-stream photon
  transport supplying the `B`, `C`, `D` blocks.
- `Theory/SecondaryEmission.rst`, eqs. `eq_see_condition` / `eq_Q0`: the
  cathode boundary conditions.
- `Theory/InceptionCriterion.rst`, eqs. `eq_Qd` / `eq_Q_system` /
  `eq_det_criterion`: `Q(λ) θ_0 = 0` and `det Q(λ=0) = 0` — this is
  `incept1d.solver._assemble_det_Q` / `incept1d.solver.inception_det`, whose root in
  `E/N` at fixed `p·d` is what `incept1d.inception.compute_inception_curve` scans for.
  The `3×3` reduced model (`eq_generalized_paschen`, `eq_standard_paschen`)
  is the closed-form sanity check for the attachment/detachment physics.
- `Theory/AirScheme.rst`, table `tab_reactions`: the reaction list
  implemented in `mechanisms/Air/Air_Pancheshnyi.py` / `mechanisms/Air/Air_2body.py`.
- `Docs/source/Numerics/`: how the propagator, determinant and root
  finding are actually implemented.

Code docstrings still cite "manuscript, eq. NNN" in places; those numbers
refer to an external LaTeX source and have drifted. When code and docs
disagree, the docs equations (by label, not number) are the intended
behaviour.

## Layout and module map

The code is a Python package, `incept1d`, in `src/incept1d/` (src layout;
`pyproject.toml` at the root). Install once per clone with
`pip install -e .` (add `[dev]` / `[docs]` extras as needed); this provides
the `incept1d` console command. Everything is driven through subcommands:

```
incept1d pdiv | eigenvalues | ionization | growth | field | chombo  [options]
```

Library (physics) and CLI (argparse, printing, plotting) are separate: each
`incept1d/cli/<command>.py` exposes `HELP`, `DESCRIPTION`,
`add_arguments(parser)` and `run(args, parser)`, and is registered in
`incept1d/cli/__init__.py:COMMANDS`. Physics never lives in `cli/`.

| Module | Role |
|---|---|
| `constants.py` | Physical constants (`kB`, `Q`, `c_light`) from `scipy.constants`. Everything else imports from here instead of hardcoding constants. |
| `reactions.py` | Declarative reaction-string parser (`"e + N2 -> 2e + N2+"`) that assembles the reaction-rate matrix `R` (`build_R` / the pre-compiled fast path `compile_reactions` + `build_R_from_compiled`). Used by mechanism files, not by the solvers directly. |
| `fields.py` | Gap-geometry abstraction (`FieldDistribution`): uniform / sphere-plane / sphere-sphere / tabulated field-line profiles `f(ξ)`, `ξ∈[0,1]`, normalised so `∫f dξ = 1`. Shared `--field` CLI parsing (`add_field_argument` / `parse_field_spec`). |
| `mechanism.py` | `load_mechanism` execs a mechanism file + optional JSON config (`read_json_configs`) into a `Mechanism` object; `REQUIRED_ATTRS` is the mechanism interface. |
| `solver.py` | Core solver: `_build_A_aug` (augmented ODE matrix), `midpoint_propagator` / `magnus2_propagator` / adaptive stepping (`parse_dx_spec`, `DX_*_DEFAULT`), `_assemble_det_Q`, and `inception_det` which evaluates `det Q(λ)` for given `E/N`, `p·d`, geometry. |
| `inception.py` | `find_all_breakdown_EN` finds the `E/N` roots of `det Q = 0`; `compute_inception_curve` tracks them (branches) over a `p·d` sweep → the inception curve PDIV(p·d). CLI: `cli/pdiv.py`. |
| `eigenvalues.py` | Diagnostic: eigenvalues of the *local* transport matrix `A = R V⁻¹` vs `E/N` (no gap integration); `max_real_eigenvalue` is shared with `ionization`/`growth`. CLI: `cli/eigenvalues.py`. |
| `ionization.py` | Ionization integrals `∫max(α−η,0)dx` and `∫max(Re λ_max(RV⁻¹),0)dx` across the gap, for comparison against the full `det Q` criterion. CLI: `cli/ionization.py`. |
| `growth.py` | For voltages above the inception voltage `V*`, solves `det Q(λ,E/N)=0` for the temporal growth rate `λ>0`. CLI: `cli/growth.py`. |
| `chombo.py` | Transport-/rate-coefficient tables from a mechanism file for the external 3-D `chombo-discharge` solver; has its own `load_raw_mechanism` because it needs the raw `REACTIONS` list. CLI: `cli/chombo.py`. |
| `output.py` | `write_metadata_header` — the date / git revision / command-line block at the top of every `--write-to-file` output. Use it; do not re-implement the git lookup. |
| `cli/` | `incept1d` entry point (`cli/__init__.py`) and one module per subcommand (`cli/field.py` is the standalone field-profile plotter). |

### Mechanism files (`mechanisms/Air/Air_Pancheshnyi.py`, `mechanisms/Air/Air_2body.py`)

Mechanism files live under `mechanisms/<family>/` and are **data, not part
of the package**: they are `exec`'d by `incept1d.mechanism.load_mechanism`
via `importlib`, after optionally injecting override variables (`_GAMMA0`,
`_EREF`, `BOLSIG_FILE`, ...) from a companion `Config.py` in the same
directory (see `mechanisms/Air/Config.py`). They import the package normally
(`from incept1d.constants import kB, Q`), so the package must be installed.
A mechanism module must expose the fixed interface
(`incept1d.mechanism.REQUIRED_ATTRS`):

- `SPECIES` (ordered list of tracked species names), `ELECTRON_INDEX`
- `get_R(EN, p, T, multipliers=None)` — reaction-rate matrix `R`
- `get_V(EN, p, T)` — diagonal drift-velocity matrix
- `get_B/get_C/get_kappa` — photoionization coupling and absorption (may
  return zero-width arrays if the mechanism has no photon groups)
- `get_Pi_e / get_Pi_plus / get_Pi_minus` — row-selection matrices for
  electrons / positive ions / negative ions
- `get_gamma_plus`, `get_gamma_plus_with`, `get_gamma_Psi` — secondary
  electron emission efficiencies (ion- and photon-induced)

When writing or editing a mechanism file, follow the coordinate convention
documented at the top of `mechanisms/Air/Air_Pancheshnyi.py` (cathode at
`x=0`, anode at `x=d`, sign convention for `V`). `Config.py` implements the
`pre_exec_vars()` / `post_exec_init()` / `mechanism_params()` protocol that
`load_mechanism` expects — copy that pattern for a new mechanism family
rather than inventing a new config mechanism.

### Typical call graph

```
mechanism file (+ Config.py, *.json) ──► incept1d.mechanism.load_mechanism ──► Mechanism
                                                                                  │
                     ┌──────────────────────┬─────────────────────┬───────────────┤
                     ▼                      ▼                     ▼               ▼
        solver.inception_det       eigenvalues.*         ionization.*        growth.*
        inception.* (PDIV curve)
                     ▲                      ▲                     ▲               ▲
             cli/pdiv.py           cli/eigenvalues.py     cli/ionization.py   cli/growth.py
```

`reactions.py` and `fields.py` sit underneath everything (used by mechanism
files and by the solvers/CLIs respectively); `constants.py` sits under all
of them.

## Working conventions

- **Formatting/linting**: `black` (line length 88) and `flake8` are wired up
  via `.pre-commit-config.yaml` / `pyproject.toml` / `.flake8`. Run
  `pre-commit install` once per clone; `pre-commit run --all-files` to check
  everything. Note: the existing `.py` files predate `black` and use manual
  column alignment in places — `black` will reformat any file it touches, so
  expect a real diff the first time a given file is committed.
- **Docs**: Sphinx sources live in `Docs/source/`, organised as one
  directory per chapter (`Introduction/`, `Theory/`, `Numerics/`,
  `Modules/`, `Examples/`, `Maintenance/`).  Build with `make html` from
  `Docs/` (or `python3 -m sphinx -W -b html Docs/source Docs/build/html`
  from the repo root; prefer `python3 -m sphinx` over the bare
  `sphinx-build`, which on at least one dev machine resolves to a pipx shim
  without numpy).  The build runs with `-W`: any warning fails it, and a
  `dummy` build is a pre-commit hook.  Figures are **built, not shipped**:
  `Docs/figures/Makefile` runs `Examples/*/run.sh` and `mechanisms/Air/Zheleznyak.py`,
  compiles the pgfplots `.tex` sources in `Docs/figures/`, and drops
  PDF/PNG into the git-ignored `Docs/source/figures/`; `make html` triggers
  it.  The full IEC computation takes tens of minutes the first time
  (`PD_NUM=30 make figures` for a quick check).  Do not add pre-rendered
  figures or reference the manuscript.  New public functions get
  NumPy-style docstrings (autodoc); docs pages `literalinclude` the
  functions that implement an equation rather than re-typing them; cite
  literature with `[Key]_` and add the entry to `ZZReferences.rst`
  (unreferenced citations fail the build).
- **Packaging**: `pyproject.toml` (setuptools, src layout) is the single
  source of dependencies and of the `incept1d` entry point. No
  `sys.path.insert` hacks anywhere: package modules use absolute
  `incept1d.*` imports, mechanism files import the installed package, tests
  and Sphinx import the installed (editable) package. Adding a subcommand =
  new `cli/<name>.py` + one line in `cli/__init__.py:COMMANDS` + a docs
  page. Shared CLI options (`--field`, `--dx`) come from
  `fields.add_field_argument` / `solver.parse_dx_spec`.
- **Physics-affecting changes**: if you change a rate coefficient, a
  boundary condition, or the augmented-matrix assembly, cite the
  corresponding equation label / table in `Docs/source/Theory/` in the
  commit message or docstring, update the theory page if the model itself
  changed, and check the closed-form limit with `mechanisms/Air/Paschen.json`.
- **Examples**: `Examples/<Name>/` holds reference data with a provenance
  header and a `run.sh` that reproduces the calculation; the corresponding
  docs page lives in `Docs/source/Examples/` and its figure source in
  `Docs/figures/<Name>.tex`.  Keep the three in sync (column indices in the
  `.tex` follow the `--write-to-file` header layout).
