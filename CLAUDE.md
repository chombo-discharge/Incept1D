# CLAUDE.md

Guidance for Claude Code (or any future contributor) working in this repository.

## What this program does

Incept1D computes the **inception (breakdown) condition** for a 1-D
drift-reaction model of a discharge gap (electrons + positive/negative ions +
two-stream photoionization), and derived quantities (Paschen curves,
ionization-integral curves, temporal growth rates, transport-coefficient
tables for 3-D simulation codes).

**Before touching any of the physics code, read `concepts.tex`, section
"Methods" → "Theoretical model" (`\subsection{Theoretical model}`, starting
around line 184, `\label{sec:model}`).** That section is the single source of
truth for the equations the code implements:

- Eq. "drift_reaction" / "two_stream": the governing PDEs (species drift +
  reactions + two-stream photon transport).
- Eq. "augmented_ode": the augmented first-order ODE system `∂_x θ = A_aug θ`
  with block matrix `A_aug = [[A,B,B],[C,-D,0],[-C,0,D]]` — this is exactly
  `Inception._build_A_aug`.
- Eq. "theta_soln" / the `Q_0`, `Q_d` block systems: the boundary-value
  problem `Q(λ) θ_0 = 0` — this is exactly `Inception._assemble_det_Q`.
- Eq. "det_criterion": `det Q(λ=0) = 0` is the inception threshold — this is
  `Inception.inception_det`, and its root in `E/N` at fixed `p·d` is what
  `Inception.compute_paschen_curve` / `Inception.main` scan for.
- The `3×3` reduced model (§ "Standard Paschen law") and
  `eq:generalized_paschen` are the closed-form sanity check for the
  attachment/detachment physics — useful for validating changes to the
  eigenvalue/propagator code against a case with a known analytic answer.
- Table "reactions" (§ "A minimal scheme for dry air") is the reaction list
  implemented in `Air/Air_Pancheshnyi.py` / `Air/Air_2body.py`.

Whenever code and manuscript disagree, treat the manuscript equation numbers
cited in the module docstrings (e.g. "manuscript, eq. 333") as the intended
behavior, and check whether the equation numbers just drifted (LaTeX
renumbers on edits) before assuming the code is wrong.

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
- **Docs**: Sphinx sources live in `Docs/`. Build locally with
  `python3 -m sphinx -b html Docs Docs/_build/html` (also runs as a
  pre-commit hook so broken autodoc imports / RST are caught before commit).
  Prefer `python3 -m sphinx` over the bare `sphinx-build` command: on at
  least one dev machine a `sphinx-build` shim from an unrelated pipx venv
  (missing numpy/scipy/matplotlib) shadows the one with the right
  dependencies earlier in `PATH`. New
  functions/classes intended for the public API should get a NumPy-style
  docstring (the existing modules are already documented this way) so
  `autodoc` picks them up; docs source files should `literalinclude` the
  handful of functions that directly implement a manuscript equation rather
  than re-typing them.
- **No build system yet**: there is no `setup.py`/`pyproject` package
  metadata beyond tool config — scripts are run directly with
  `python Inception.py ...` from the repo root, and mechanism files resolve
  their own relative imports via `sys.path.insert(0, ...)` at the top of the
  file. Preserve that pattern rather than introducing package-relative
  imports.
- **Physics-affecting changes**: if you change a rate coefficient, a
  boundary condition, or the augmented-matrix assembly, cite the
  corresponding equation/table in `concepts.tex` in the commit message or
  docstring, the same way the existing code does.
