# CLAUDE.md

Guidance for Claude Code (or any future contributor) working in this repository.

## What this program does

Incept1D computes the **inception (breakdown) condition** for a 1-D
drift-reaction model of a discharge gap (electrons + positive/negative ions +
two-stream photoionization), and derived quantities (inception curves PDIV(p·d),
ionization-integral curves, temporal growth rates, transport-coefficient
tables for 3-D simulation codes).

**Before touching any of the physics code, read the Theory chapter of the
documentation, `docs/source/theory/`.**  Those pages are the single source
of truth for the equations the code implements (the manuscript they were
derived from is not part of the repository and must not be referenced
from the docs).  The chapter builds the criterion in order — what is
transported, then the photons, then the electrodes, then the assembly:

- `theory/overview.rst`, eq. `eq_drift_reaction`: the governing drift-reaction PDE.
- `theory/transport.rst`, eq. `eq_flux_ode`: the charged-species equation in
  flux form, and the eigenvalues of `R V⁻¹`.
- `theory/photoionization.rst`, eq. `eq_two_stream`: the two-stream photon
  transport supplying the `B`, `C`, `D` blocks.
- `theory/secondaryemission.rst`, eqs. `eq_see_condition` / `eq_Q0`: the
  electrode boundary conditions and the row-selection operators.
- `theory/augmented.rst`, eq. `eq_augmented_ode`: the augmented first-order
  ODE `∂_x θ = A_aug θ` with block matrix `A_aug = [[A,B,B],[C,-D,0],[-C,0,D]]`
  — this is exactly `incept1d.solver._build_A_aug`; eq. `eq_theta_soln` is the
  propagator `M(d)`.
- `theory/inceptioncriterion.rst`, eqs. `eq_Qd` / `eq_Q_system` /
  `eq_det_criterion`: `Q(λ) θ_0 = 0` and `det Q(λ=0) = 0` — this is
  `incept1d.solver.inception_det` (`--criterion detq`).  The default,
  `incept1d.solver.riccati_criterion`, evaluates the same condition without
  forming `M` (`numerics/riccati.rst`, eqs. `eq_riccati` /
  `eq_riccati_criterion`: a reflection operator carried from the anode,
  `g = 1 − loop gain`).  Their root in `E/N` at fixed `p·d` is what
  `incept1d.inception.compute_inception_curve` scans for.
  The `3×3` reduced model (`eq_generalized_paschen`, `eq_standard_paschen`)
  is the closed-form sanity check for the attachment/detachment physics, and
  is what `tests/closed_form.py` transcribes.
- `docs/source/numerics/`: how the propagator, the two criteria
  (`numerics/riccati.rst`, `numerics/determinant.rst`) and root finding
  are actually implemented.  Every photon group is propagated explicitly:
  folding thick groups into `A` removed photon feedback and must not come
  back.

The reaction schemes are inputs, not part of the derivation: they live in
`docs/source/configuration/examples/` (`air.rst`, table `tab_reactions`, is
the scheme implemented by `mechanisms/air/pancheshnyi/air_pancheshnyi.py`).

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
| `fields.py` | Gap-geometry abstraction (`FieldDistribution`): uniform / sphere-plane / sphere-sphere / coaxial / tabulated field-line profiles `f(ξ)`, `ξ∈[0,1]`, normalised so `∫f dξ = 1`. Shared `--field` CLI parsing (`add_field_argument` / `parse_field_spec`). |
| `mechanism.py` | `load_mechanism` execs a mechanism file + optional JSON config (`read_json_configs`) into a `Mechanism` object; `REQUIRED_ATTRS` is the mechanism interface. |
| `solver.py` | Core solver: `_build_A_aug` (augmented ODE matrix, every photon group explicit), `midpoint_propagator` / `magnus2_propagator` / adaptive stepping from a field-following start (`_field_following_edges`, `parse_dx_spec`, `DX_*_DEFAULT`), the two criteria — `riccati_criterion` (default; `_riccati_g`, adding–doubling in `_slab_scattering` / `_star`, pole test `_spectral_radius`) and `inception_det` (`det Q(λ)`, compound-matrix fallback `_det_Q_compound`) — and `CRITERIA` / `add_criterion_argument` for `--criterion`. |
| `inception.py` | `find_all_breakdown_EN` finds the `E/N` roots of the criterion (bottom-up scan, floor check below 10 Td, + → − rule for Riccati); `compute_inception_curve` tracks them (branches) over a `p·d` sweep, in parallel blocks with `jobs` → the inception curve PDIV(p·d). CLI: `cli/pdiv.py`. |
| `eigenvalues.py` | Diagnostic: eigenvalues of the *local* transport matrix `A = R V⁻¹` vs `E/N` (no gap integration); `max_real_eigenvalue` is shared with `ionization`/`growth`. CLI: `cli/eigenvalues.py`. |
| `ionization.py` | Ionization integrals `∫max(α−η,0)dx` and `∫max(Re λ_max(RV⁻¹),0)dx` across the gap, for comparison against the full `det Q` criterion. CLI: `cli/ionization.py`. |
| `growth.py` | For voltages above the inception voltage `V*`, finds the temporal growth rate `λ>0` as the largest real root of the criterion in `λ` (`find_lambda_for_voltage`; `solve_voltage` / `voltage_sweep` let many sweeps share a worker pool). CLI: `cli/growth.py` (`--pressure`/`--distance` ranges, `--jobs`, `--verify`). |
| `chombo.py` | Transport-/rate-coefficient tables from a mechanism file for the external 3-D `chombo-discharge` solver; has its own `load_raw_mechanism` because it needs the raw `REACTIONS` list. CLI: `cli/chombo.py`. |
| `parallel.py` | `physical_cores` (default `--jobs`: physical cores only, 1 if undeterminable) and `parallel_map` (forked workers, contiguous warm-started blocks, results streamed back in completion order); `limit_blas_threads` runs in `cli/__init__.py` before NumPy loads. |
| `output.py` | `write_metadata_header` — the date / git revision / command-line block at the top of every `--write-to-file` output. Use it; do not re-implement the git lookup. |
| `cli/` | `incept1d` entry point (`cli/__init__.py`) and one module per subcommand (`cli/field.py` is the standalone field-profile plotter). |

### Mechanism files (`mechanisms/air/pancheshnyi/air_pancheshnyi.py`, `mechanisms/air/2body/air_2body.py`)

Mechanism files live one per directory under `mechanisms/<gas>/<scheme>/`
and are **data, not part of the package**: they are `exec`'d by
`incept1d.mechanism.load_mechanism` via `importlib`, after optionally
injecting override variables (`_GAMMA0`, `_EREF`, `BOLSIG_FILE`, ...) from
the companion `config.py` in the same directory. Each directory also holds
the JSON configurations written for that mechanism: a configuration belongs
to one mechanism, and the loader rejects reaction multipliers that name no
reaction of it. What the mechanisms of one gas share lives in
`mechanisms/<gas>/` — for air, `air_config.py` (the configuration class each
`config.py` subclasses, restricting the accepted keys), `zheleznyak.py`, and
the LXCat data in `lxcat/` — and is loaded by path with
`incept1d.mechanism.load_helper`, never via `sys.path`. They import the package normally
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
documented at the top of `mechanisms/air/pancheshnyi/air_pancheshnyi.py` (cathode at
`x=0`, anode at `x=d`, sign convention for `V`). `config.py` implements the
`pre_exec_vars()` / `post_exec_init()` / `mechanism_params()` protocol that
`load_mechanism` expects — copy that pattern for a new mechanism rather
than inventing a new config mechanism.

### Typical call graph

```
mechanism file (+ config.py, *.json) ──► incept1d.mechanism.load_mechanism ──► Mechanism
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

These mirror `docs/source/maintenance/`; that chapter is the contributor-facing
version of the same rules, so change both together.

- **Before committing**: `pre-commit run --all-files` and `python3 -m pytest`
  must pass. Run `pre-commit install` once per clone. Hooks: `black`
  (line length 88), `flake8`, `reuse lint`, and a Sphinx `dummy` build. If
  the Sphinx hook fails oddly right after moving or renaming a docs file,
  clear its cache: `rm -rf /tmp/incept1d-docs-precommit`.

- **Tests**: `tests/`, run with `pytest`. `-m "not slow"` skips the ones that
  load real swarm data. Three layers, strongest first: verification against
  closed forms (`tests/toy/` plus `tests/closed_form.py`, which must stay an
  independent transcription of the documented algebra — never import solver
  code into it), unit tests, then invariants. **Prefer an invariant to a
  pinned value**: a wrong-but-self-consistent curve passes a pin, which is
  exactly how the warm-start bug survived. Add tests with the change, not
  after.

- **Packaging**: `pyproject.toml` (setuptools, src layout) is the single
  source of dependencies and of the `incept1d` entry point. No
  `sys.path.insert` hacks anywhere: package modules use absolute
  `incept1d.*` imports, mechanism files import the installed package, tests
  and Sphinx import the installed (editable) package. Adding a subcommand =
  new `cli/<name>.py` + one line in `cli/__init__.py:COMMANDS` + a docs
  page. Shared CLI options (`--field`, `--dx`) come from
  `fields.add_field_argument` / `solver.parse_dx_spec`.

- **CI**: one workflow, `.github/workflows/ci.yml`, because `needs:` cannot
  reach across workflow files. Jobs `reuse`, `tests`, `rst`, `docs` all feed
  the aggregate `CI-passed`, which is the single context branch protection
  on `main` requires; Pages deployment is gated behind it. A new check must
  be added to `needs:` of `CI-passed` or it gates nothing. `main` takes
  squash merges only, and admin bypass is enabled.

- **Licensing (REUSE)**: the project is GPL-3.0-or-later and must stay REUSE
  compliant — `reuse lint` is a hook and a CI job. New `.py` files get the
  SPDX header; everything else is declared in `REUSE.toml`, whose existing
  globs usually already cover it.
  **Never add content the project has no right to redistribute.** Published
  standards, journal tables, figures and datasets are copyrighted, and a
  small extract used for validation is still redistribution. Commit the
  *calculation* instead, document the expected input file, and let a reader
  supply the data — two worked examples are built that way. For data that
  *is* redistributable: keep it verbatim with its original header, annotate
  it in `REUSE.toml` with its **real** rights holder (never SINTEF by
  default), and add a `LICENSES/LicenseRef-*.txt` if its terms are not an
  SPDX licence. Mislabelling third-party data as project-owned is worse than
  leaving it undeclared, because `reuse lint` then passes.

- **Docs**: Sphinx sources in `docs/source/`, one directory per chapter
  (`introduction/`, `theory/`, `numerics/`, `configuration/`, `modules/`,
  `examples/`, `maintenance/`). Build with `make html` from `docs/`; prefer
  `python3 -m sphinx` over bare `sphinx-build`, which on at least one dev
  machine resolves to a pipx shim without numpy. The build runs with `-W`,
  so any warning fails it. New pages go in the right `toctree` in
  `index.rst`. New public functions get NumPy-style docstrings (autodoc).
  Cite literature with `[Key]_` and add the entry to `zzreferences.rst` —
  unreferenced citations fail the build.

- **Docs style**: bullet and numbered list items begin with a capital letter.
  The generic chapters (Theory, Numerics, and the interface pages of
  Configuration files) must not name a specific mechanism, gas or data file;
  implementation detail belongs in `configuration/examples/`. Subscripts
  follow one scheme — `e`, `+`, `-` for species, `Ψ^±` for photon streams,
  `j` for photon groups — so `Π_-` is negative ions and `Π_{Ψ^-}` is
  backward photons. User documentation does not show implementation code:
  the `modules/` pages describe what a command does, what it takes and what
  it writes, and leave the source to `automodule`. Where an example page
  does show a file, `literalinclude` a whole object or a whole file, never
  `:lines:`, and pass `:dedent:` for a method.

- **Command pages** in `docs/source/modules/` follow one shape: an opening
  paragraph saying what the command answers, then Inputs, Examples, Outputs
  and API reference.

- **Figures** are built, not shipped: `docs/figures/Makefile` runs the
  examples' `incept1d` commands and the mechanism helpers, compiles the pgfplots
  `.tex` sources, and drops PDF/PNG into the git-ignored
  `docs/source/figures/`; `make html` triggers it. `make figures` builds
  only the self-contained figures and takes seconds. Figures that compare
  against published reference data are opt-in (`OPTIONAL_FIGURES`), because
  that data is not in the repository and those runs dominate the cost.
  Do not add pre-rendered figures or reference the manuscript.

- **Physics-affecting changes**: if you change a rate coefficient, a
  boundary condition, or the augmented-matrix assembly, cite the
  corresponding equation label in `docs/source/theory/` in the commit
  message or docstring, update the theory page if the model itself changed,
  and check a closed-form limit. Include before/after inception curves in
  the PR, and confirm that cases the change should *not* affect are
  unchanged.

- **Examples**: `examples/<name>/` holds a `README.md`, not a script: what
  the example shows and represents, and one main `incept1d` command (the
  first ```` ```bash ```` block starting with `incept1d`) that a reader can
  paste and change. The docs page in `docs/source/examples/<name>.rst`
  shows that same command, and `docs/figures/Makefile` repeats it headless
  (plus any extra cases the figure needs) — `tests/test_docs.py` checks
  that the command parses and that the docs page matches the README. The
  figure source is `docs/figures/<name>.tex`; column indices follow the
  `--write-to-file` header layout. Where a
  comparison needs data we cannot ship, the example documents the expected
  filename and column layout instead, and the figure is opt-in.
