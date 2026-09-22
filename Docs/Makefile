# Minimal makefile for the Incept1D Sphinx documentation.
#
# Sources live in source/, build products go to build/<builder>/.
#
# `python3 -m sphinx` is used instead of the bare `sphinx-build` executable
# on purpose: on some machines a `sphinx-build` shim from an unrelated
# virtualenv (without numpy/scipy/matplotlib) shadows the one with the
# right dependencies earlier in PATH.  Override with e.g.
# `make html SPHINXBUILD="sphinx-build"` if needed.

SPHINXOPTS    ?= -W --keep-going
SPHINXBUILD   ?= python3 -m sphinx
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: help clean dummy figures Makefile

help:
	@echo "Incept1D documentation.  Targets:"
	@echo "  make figures    compute example data and build the figures (slow the first time)"
	@echo "  make html       figures + HTML into build/html"
	@echo "  make latexpdf   figures + PDF into build/latex"
	@echo "  make dummy      full Sphinx check without output (pre-commit / CI)"
	@echo "  make clean      remove build/ and the generated figures"
	@echo
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf $(BUILDDIR)

# Full pipeline (autodoc, cross-references, RST parsing) but no output files.
# Used by the pre-commit hook and CI as a fast correctness check.
dummy:
	@$(SPHINXBUILD) -b dummy "$(SOURCEDIR)" "$(BUILDDIR)/dummy" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets (html, latexpdf, linkcheck,
# ...) to Sphinx using the "make mode" option, after making sure the
# figures exist.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile figures
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
