# Build the figures of the Incept1D documentation from the reference data and
# the solver, so that nothing pre-rendered is shipped.
#
#   make            build every figure (PDF + PNG) into ../source/figures/
#   make <Figure>   build one figure, e.g. `make Electra`
#   make clean      remove build products
#
# Pipeline, per figure:
#   1. data:   Examples/*/run.sh (incept1d pdiv) or mechanisms/Air/Zheleznyak.py write
#              .dat files into build/figures/<name>/, next to a copy of the
#              reference tables from Examples/;
#   2. pgfplots: <name>.tex is compiled with pdflatex in that directory;
#   3. the PDF (LaTeX builder) and a 150 dpi PNG (HTML builder) are copied to
#      Docs/source/figures/, which is ignored by git.
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

FIGURES := IEC60052 Electra ZheleznyakFit
DIAMS   := 10 15 25 50 75 100 150 200

# Anything that changes the physics or numerics invalidates the computed data.
SOLVER_DEPS := $(wildcard $(ROOT)/src/incept1d/*.py) $(ROOT)/src/incept1d/cli/pdiv.py \
               $(wildcard $(ROOT)/mechanisms/Air/*.py) $(wildcard $(ROOT)/mechanisms/Air/*.json) \
               $(wildcard $(ROOT)/mechanisms/Air/*.txt)

.PHONY: all clean data $(FIGURES)
all: $(foreach f,$(FIGURES),$(OUT)/$(f).pdf $(OUT)/$(f).png)

# ---- 1. data -----------------------------------------------------------------

IEC_SIM := $(foreach D,$(DIAMS),$(BUILD)/IEC60052/Sim_$(D)cm.dat)
IEC_REF := $(patsubst $(ROOT)/Examples/IEC60052/%,$(BUILD)/IEC60052/%,$(wildcard $(ROOT)/Examples/IEC60052/IEC60052_*.dat))

$(BUILD)/IEC60052/Sim_%cm.dat: $(ROOT)/Examples/IEC60052/run.sh $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	OUT=$(dir $@) bash $< $*

$(BUILD)/IEC60052/IEC60052_%.dat: $(ROOT)/Examples/IEC60052/IEC60052_%.dat
	@mkdir -p $(dir $@) && cp $< $@

$(BUILD)/Electra/Sim.dat: $(ROOT)/Examples/Electra/run.sh $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	OUT=$(dir $@) bash $<

$(BUILD)/Electra/TableC2_Air.dat: $(ROOT)/Examples/Electra/TableC2_Air.dat
	@mkdir -p $(dir $@) && cp $< $@

$(BUILD)/ZheleznyakFit/zheleznyak_absorption_curve.dat: $(ROOT)/mechanisms/Air/Zheleznyak.py
	@mkdir -p $(dir $@)
	cd $(dir $@) && python3 $< --ngroups 8 --threshold 1E-6 --write-to-file zheleznyak

data: $(IEC_SIM) $(IEC_REF) $(BUILD)/Electra/Sim.dat $(BUILD)/Electra/TableC2_Air.dat \
      $(BUILD)/ZheleznyakFit/zheleznyak_absorption_curve.dat

# ---- 2. pgfplots + 3. outputs ---------------------------------------------------

DATA_IEC60052     := $(IEC_SIM) $(IEC_REF)
DATA_Electra      := $(BUILD)/Electra/Sim.dat $(BUILD)/Electra/TableC2_Air.dat
DATA_ZheleznyakFit := $(BUILD)/ZheleznyakFit/zheleznyak_absorption_curve.dat

# One rule set per figure: tex + data -> PDF in the build dir -> PDF/PNG in OUT.
define FIGURE_RULES
$(BUILD)/$(1)/$(1).pdf: $(FIGDIR)/$(1).tex $(FIGDIR)/preamble.tex $(DATA_$(1))
	cd $(BUILD)/$(1) && TEXINPUTS=$(FIGDIR): pdflatex -interaction=nonstopmode -halt-on-error $(FIGDIR)/$(1).tex > pdflatex.log
	@test "$$$$(pdfinfo $$@ | awk '/^Pages:/ {print $$$$2}')" = 1 || { echo "ERROR: $$@ has more than one page - the plot overflowed the paper size in $(1).tex"; exit 1; }

$(OUT)/$(1).pdf: $(BUILD)/$(1)/$(1).pdf
	@mkdir -p $(OUT) && cp $$< $$@

$(OUT)/$(1).png: $(BUILD)/$(1)/$(1).pdf
	@mkdir -p $(OUT) && pdftoppm -png -r 150 -singlefile $$< $(OUT)/$(1)

# Convenience target: `make ZheleznyakFit` builds one figure.
$(1): $(OUT)/$(1).pdf $(OUT)/$(1).png
endef
$(foreach f,$(FIGURES),$(eval $(call FIGURE_RULES,$(f))))

clean:
	rm -rf $(BUILD) $(OUT)
