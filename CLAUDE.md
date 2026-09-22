# CLAUDE.md

Guidance for Claude Code (or any future contributor) working in this repository.

## What this program does

Incept1D computes the **inception (breakdown) condition** for a 1-D
drift-reaction model of a discharge gap (electrons + positive/negative ions +
two-stream photoionization), and derived quantities (Paschen curves,
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
  — this is exactly `Inception._build_A_aug`; eq. `eq_theta_soln` is the
  propagator `M(d)`.
- `Theory/Photoionization.rst`, eq. `eq_two_stream`: the two-stream photon
  transport supplying the `B`, `C`, `D` blocks.
- `Theory/SecondaryEmission.rst`, eqs. `eq_see_condition` / `eq_Q0`: the
  cathode boundary conditions.
- `Theory/InceptionCriterion.rst`, eqs. `eq_Qd` / `eq_Q_system` /
  `eq_det_criterion`: `Q(λ) θ_0 = 0` and `det Q(λ=0) = 0` — this is
  `Inception._assemble_det_Q` / `Inception.inception_det`, whose root in
  `E/N` at fixed `p·d` is what `Inception.compute_paschen_curve` scans for.
  The `3×3` reduced model (`eq_generalized_paschen`, `eq_standard_paschen`)
  is the closed-form sanity check for the attachment/detachment physics.
- `Theory/AirScheme.rst`, table `tab_reactions`: the reaction list
  implemented in `Air/Air_Pancheshnyi.py` / `Air/Air_2body.py`.
- `Docs/source/Numerics/`: how the propagator, determinant and root
  finding are actually implemented.

Code docstrings still cite "manuscript, eq. NNN" in places; those numbers
refer to an external LaTeX source and have drifted. When code and docs
disagree, the docs equations (by label, not number) are the intended
behaviour.

## Module map

All modules below live at the repo root and import each other directly
(no package/`src` layout, no `__init__.py`). Run them from the repo root.

| Module | Role |
|---|---|
| `Constants.py` | Physical constants (`kB`, `Q`, `c_light`) from `scipy.constants`. Everything else imports from here instead of hardcoding constants. |
| `Reactions.py` | Declarative reaction-string parser (`"e + N2 -> 2e + N2+"`) that assembles the reaction-rate matrix `R` (`build_R` / the pre-compiled fast path `compile_reactions` + `build_R_from_compiled`). Used by mechanism files, not by the solvers directly. |
| `FieldDistributions.py` | Gap-geometry abstraction (`FieldDistribution`): uniform / sphere-plane / sphere-sphere field profiles `f(ξ)`, `ξ∈[0,1]`, normalised so `∫f dξ = 1`. Shared `--field` CLI parsing used by `Inception.py`, `IonizationIntegral.py`, `Lambda.py`. |
| `Inception.py` | Core solver. `load_mechanism` loads a mechanism file + optional JSON config into a `Mechanism` object; `inception_det` integrates the augmented ODE across the gap and evaluates `det Q(λ)`; `compute_paschen_curve` finds and tracks all `E/N` roots (branches) over a `p·d` sweep. Also a CLI (`main`) that plots Paschen curves. |
| `Eigenvalues.py` | Diagnostic: eigenvalues of the *local* transport matrix `A = R V⁻¹` vs `E/N` (no gap integration). A positive real eigenvalue means locally growing charge density. Imports `load_mechanism` from `Inception.py`. |
| `IonizationIntegral.py` | Plots the classical ionization integral `∫max(α−η,0)dx` alongside `∫max(Re λ_max(RV⁻¹),0)dx` vs. applied voltage, for comparison against the full `det Q` inception criterion. |
| `Lambda.py` | For voltages above the inception voltage `V*`, solves `det Q(λ,E/N)=0` for the temporal growth rate `λ>0` (discharge growth rate above threshold). |
| `CreateChomboDischargeData.py` | Exports transport-coefficient / rate-coefficient tables from a mechanism file for use by the external 3-D `chombo-discharge` plasma solver (the "3D plasma simulations" mentioned in the paper). Not part of the inception solve itself. |

### Mechanism files (e.g. `Air/Air_Pancheshnyi.py`, `Air/Air_2body.py`)

These are **not imported as Python packages** — they are `exec`'d by
`Inception.load_mechanism` via `importlib`, after optionally injecting
override variables (`_GAMMA0`, `_EREF`, `BOLSIG_FILE`, ...) from a companion
`Config.py` in the same directory (see `Air/Config.py`). A mechanism module
must expose the fixed interface (`Inception._REQUIRED_ATTRS`):

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
documented at the top of `Air/Air_Pancheshnyi.py` (cathode at `x=0`, anode at
`x=d`, sign convention for `V`). `Config.py` implements the
`pre_exec_vars()` / `post_exec_init()` / `mechanism_params()` protocol that
`load_mechanism` expects — copy that pattern for a new mechanism family
rather than inventing a new config mechanism.

### Typical call graph

```
mechanism.py (+ Config.py) ──► Inception.load_mechanism ──► Mechanism
                                                                │
                        ┌───────────────────────┬──────────────┼───────────────────────┐
                        ▼                       ▼              ▼                       ▼
                Inception.main          Eigenvalues.main  IonizationIntegral.main  Lambda.main
             (Paschen curves)         (local eigenvalues)   (ionization integral)   (growth rate)
```

`Reactions.py` and `FieldDistributions.py` sit underneath everything (used by
mechanism files and by the solvers/CLIs respectively); `Constants.py` sits
under all of them.

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
  `Docs/figures/Makefile` runs `Examples/*/run.sh` and `Air/Zheleznyak.py`,
  compiles the pgfplots `.tex` sources in `Docs/figures/`, and drops
  PDF/PNG into the git-ignored `Docs/source/figures/`; `make html` triggers
  it.  The full IEC computation takes tens of minutes the first time
  (`PD_NUM=30 make figures` for a quick check).  Do not add pre-rendered
  figures or reference the manuscript.  New public functions get
  NumPy-style docstrings (autodoc); docs pages `literalinclude` the
  functions that implement an equation rather than re-typing them; cite
  literature with `[Key]_` and add the entry to `ZZReferences.rst`
  (unreferenced citations fail the build).
- **No build system yet**: there is no `setup.py`/`pyproject` package
  metadata beyond tool config — scripts are run directly with
  `python Inception.py ...` from the repo root, and mechanism files resolve
  their own relative imports via `sys.path.insert(0, ...)` at the top of the
  file. Preserve that pattern rather than introducing package-relative
  imports.
- **Physics-affecting changes**: if you change a rate coefficient, a
  boundary condition, or the augmented-matrix assembly, cite the
  corresponding equation label / table in `Docs/source/Theory/` in the
  commit message or docstring, update the theory page if the model itself
  changed, and check the closed-form limit with `Air/Paschen.json`.
- **Examples**: `Examples/<Name>/` holds reference data with a provenance
  header and a `run.sh` that reproduces the calculation; the corresponding
  docs page lives in `Docs/source/Examples/` and its figure source in
  `Docs/figures/<Name>.tex`.  Keep the three in sync (column indices in the
  `.tex` follow the `--write-to-file` header layout).
