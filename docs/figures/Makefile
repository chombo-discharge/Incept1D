# Build the figures of the Incept1D documentation from the reference data and
# the solver, so that nothing pre-rendered is shipped.
#
#   make            build every figure (PDF + PNG) into ../source/figures/
#   make <Figure>   build one figure, e.g. `make Electra`
#   make clean      remove build products
#
# Pipeline, per figure:
#   1. data:   examples/*/run.sh (incept1d pdiv) or mechanisms/air/zheleznyak.py
#              write .dat files into build/figures/<name>/, next to a copy of
#              any reference tables the user has placed in examples/ (those are
#              third-party and not shipped -- see the READMEs there);
#   2. pgfplots: <name>.tex is compiled with pdflatex in that directory;
#   3. the PDF (LaTeX builder) and a 150 dpi PNG (HTML builder) are copied to
#      docs/source/figures/, which is ignored by git.
#
# The incept1d pdiv runs are the expensive step (tens of minutes for the full
# IEC set).  They are only re-run when the solver, the mechanism, or the
# run script changes; set PD_NUM=<n> to use fewer pd points for a quick build.
#
# Requirements: the incept1d package installed (pip install -e .), pdflatex with pgfplots and
# siunitx (TeX Live: texlive-latex-extra texlive-pictures texlive-science),
# pdftoppm (poppler-utils).

ROOT   := $(abspath ../..)
BUILD  := $(abspath ../build/figures)
OUT    := $(abspath ../source/figures)
FIGDIR := $(abspath .)

export MPLBACKEND := Agg
export PD_NUM

FIGURES := iec60052 electra zheleznyakfit
DIAMS   := 10 15 25 50 75 100 150 200

# Anything that changes the physics or numerics invalidates the computed data.
SOLVER_DEPS := $(wildcard $(ROOT)/src/incept1d/*.py) $(ROOT)/src/incept1d/cli/pdiv.py \
               $(wildcard $(ROOT)/mechanisms/air/*.py) $(wildcard $(ROOT)/mechanisms/air/*.json) \
               $(wildcard $(ROOT)/mechanisms/air/*.txt)

.PHONY: all clean data $(FIGURES)
all: $(foreach f,$(FIGURES),$(OUT)/$(f).pdf $(OUT)/$(f).png)

# ---- 1. data -----------------------------------------------------------------

IEC_SIM := $(foreach D,$(DIAMS),$(BUILD)/iec60052/sim_$(D)cm.dat)
# Reference tables are OPTIONAL: the IEC 60052 and CIGRE ELECTRA values are
# copyrighted by their publishers and are not part of this repository.  Both
# lists come from $(wildcard ...) and are simply empty when the user has not
# supplied the files, in which case the figures show the computed curves alone
# (each overlay in the .tex is guarded by \IfFileExists).
IEC_REF     := $(patsubst $(ROOT)/examples/iec60052/%,$(BUILD)/iec60052/%,$(wildcard $(ROOT)/examples/iec60052/iec60052_*.dat))
ELECTRA_REF := $(patsubst $(ROOT)/examples/electra/%,$(BUILD)/electra/%,$(wildcard $(ROOT)/examples/electra/tablec2_air.dat))

$(BUILD)/iec60052/sim_%cm.dat: $(ROOT)/examples/iec60052/run.sh $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	OUT=$(dir $@) bash $< $*

$(BUILD)/iec60052/iec60052_%.dat: $(ROOT)/examples/iec60052/iec60052_%.dat
	@mkdir -p $(dir $@) && cp $< $@

$(BUILD)/electra/sim.dat: $(ROOT)/examples/electra/run.sh $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	OUT=$(dir $@) bash $<

$(BUILD)/electra/tablec2_air.dat: $(ROOT)/examples/electra/tablec2_air.dat
	@mkdir -p $(dir $@) && cp $< $@

$(BUILD)/zheleznyakfit/zheleznyak_absorption_curve.dat: $(ROOT)/mechanisms/air/zheleznyak.py
	@mkdir -p $(dir $@)
	cd $(dir $@) && python3 $< --ngroups 8 --threshold 1E-6 --write-to-file zheleznyak

data: $(IEC_SIM) $(IEC_REF) $(BUILD)/electra/sim.dat $(ELECTRA_REF) \
      $(BUILD)/zheleznyakfit/zheleznyak_absorption_curve.dat

# ---- 2. pgfplots + 3. outputs ---------------------------------------------------

DATA_iec60052     := $(IEC_SIM) $(IEC_REF)
DATA_electra      := $(BUILD)/electra/sim.dat $(ELECTRA_REF)
DATA_zheleznyakfit := $(BUILD)/zheleznyakfit/zheleznyak_absorption_curve.dat

# One rule set per figure: tex + data -> PDF in the build dir -> PDF/PNG in OUT.
define FIGURE_RULES
$(BUILD)/$(1)/$(1).pdf: $(FIGDIR)/$(1).tex $(FIGDIR)/preamble.tex $(DATA_$(1))
	cd $(BUILD)/$(1) && TEXINPUTS=$(FIGDIR): pdflatex -interaction=nonstopmode -halt-on-error $(FIGDIR)/$(1).tex > pdflatex.log
	@test "$$$$(pdfinfo $$@ | awk '/^Pages:/ {print $$$$2}')" = 1 || { echo "ERROR: $$@ has more than one page - the plot overflowed the paper size in $(1).tex"; exit 1; }

$(OUT)/$(1).pdf: $(BUILD)/$(1)/$(1).pdf
	@mkdir -p $(OUT) && cp $$< $$@

$(OUT)/$(1).png: $(BUILD)/$(1)/$(1).pdf
	@mkdir -p $(OUT) && pdftoppm -png -r 150 -singlefile $$< $(OUT)/$(1)

# Convenience target: `make zheleznyakfit` builds one figure.
$(1): $(OUT)/$(1).pdf $(OUT)/$(1).png
endef
$(foreach f,$(FIGURES),$(eval $(call FIGURE_RULES,$(f))))

clean:
	rm -rf $(BUILD) $(OUT)
