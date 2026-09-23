# Summary

### Background

### Solution

### Side-effects

### Alternative solutions

### PR-review checklist

Mandatory:

- [ ] I have run `pre-commit run --all-files` (black, flake8, Sphinx) and it passed.
- [ ] I have run the test suite and made sure that it passed.
- [ ] I have added all relevant user documentation to Sphinx, and `make -C docs html` builds without warnings.
- [ ] I have added NumPy-style docstrings to all new public functions and classes.
- [ ] I have added appropriate labels to this PR.
- [ ] I have added/revised proper licensing and copyright information.
- [ ] I have run a PR review using `@claude review`.

If this PR changes the physics or the numerics (rate coefficients, boundary
conditions, the augmented matrix, propagators, root finding):

- [ ] I have cited the corresponding equation label or table in `docs/source/theory/` in the commit message or docstring, and updated the theory page if the model itself changed.
- [ ] I have checked the closed-form limit with `mechanisms/air/paschen.json` and it still reproduces the standard Paschen law.
- [ ] I have included before/after inception curves (or equivalent evidence) for a case the change is expected to affect.
- [ ] I have confirmed that cases the change is *not* expected to affect are unchanged.

If this PR adds or modifies a mechanism:

- [ ] I have given the swarm-data source and a literature reference for every rate.
- [ ] I have added a configuration file that reduces the scheme to a limit with a known answer.
