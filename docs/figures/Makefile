# Build the figures of the Incept1D documentation from the reference data and
# the solver, so that nothing pre-rendered is shipped.
#
#   make            build every figure (PDF + PNG) into ../source/figures/
#   make <Figure>   build one figure, e.g. `make Electra`
#   make clean      remove build products
#
# Pipeline, per figure:
#   1. data:   incept1d pdiv (the commands below, which are the examples'
#              README commands with --no-plot and the output redirected) or
#              mechanisms/air/zheleznyak.py write .dat files into
#              build/figures/<name>/, next to a copy of
#              any reference tables the user has placed in examples/ (those are
#              third-party and not shipped -- see the READMEs there);
#   2. pgfplots: <name>.tex is compiled with pdflatex in that directory;
#   3. the PDF (LaTeX builder) and a 150 dpi PNG (HTML builder) are copied to
#      docs/source/figures/, which is ignored by git.
#
# The incept1d pdiv runs are the expensive step (tens of minutes for the full
# IEC set).  They are only re-run when the solver, a mechanism, or this
# Makefile changes; set PD_NUM=<n> to use fewer pd points for a quick build.
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

# Figures built by default: self-contained and cheap.  `paschen` verifies the
# solver against its own closed form; `zheleznyakfit` shows the quality of the
# two-stream fit; `coaxial` shows corona inception on a wire in a cylinder.
# Together they take about a minute.
FIGURES := paschen zheleznyakfit coaxial

# Opt-in figures.  These compare against IEC 60052 and CIGRE ELECTRA reference
# tables that are copyrighted and therefore not part of this repository (see
# examples/*/README.md), so by default they would plot the computed curves
# with nothing to compare them against.  They are also by far the most
# expensive: the eight IEC sphere diameters alone are ~75 % of the total
# figure cost.  Build one explicitly after supplying the data:
#
#     make -C docs/figures iec60052
#
OPTIONAL_FIGURES := iec60052 electra

ALL_FIGURES := $(FIGURES) $(OPTIONAL_FIGURES)
DIAMS   := 10 15 25 50 75 100 150 200

# Anything that changes the physics or numerics invalidates the computed data.
SOLVER_DEPS := $(wildcard $(ROOT)/src/incept1d/*.py) $(ROOT)/src/incept1d/cli/pdiv.py \
               $(wildcard $(ROOT)/mechanisms/air/*.py) $(wildcard $(ROOT)/mechanisms/air/*.json) \
               $(wildcard $(ROOT)/mechanisms/air/*.txt) \
               $(wildcard $(ROOT)/mechanisms/paschen/*.py) $(wildcard $(ROOT)/mechanisms/paschen/*.json) \
               $(FIGDIR)/Makefile

# Every run below is an example's README command, run from the repository
# root, headless, into the build directory.  Keep the two in step: the README
# is what a reader runs, this is what the figure shows.
PDIV    := cd $(ROOT) && incept1d pdiv
AIR     := mechanisms/air/air_pancheshnyi.py
NODET   := mechanisms/air/nodetachment.json

.PHONY: all clean data $(ALL_FIGURES)
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

# IEC 60052: sphere radius (mm) and pd_max (bar mm) per diameter (cm), matching
# the extent of the IEC tables.  See examples/iec60052/README.md.
IEC_R_10 := 50
IEC_R_15 := 75
IEC_R_25 := 125
IEC_R_50 := 250
IEC_R_75 := 375
IEC_R_100 := 500
IEC_R_150 := 750
IEC_R_200 := 1000
IEC_PDMAX_10 := 55
IEC_PDMAX_15 := 90
IEC_PDMAX_25 := 150
IEC_PDMAX_50 := 300
IEC_PDMAX_75 := 400
IEC_PDMAX_100 := 550
IEC_PDMAX_150 := 800
IEC_PDMAX_200 := 1000

$(BUILD)/iec60052/sim_%cm.dat: $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	$(PDIV) $(AIR) $(NODET) \
	    --pd-min 1 --pd-max $(IEC_PDMAX_$*) --pd-num $(or $(PD_NUM),100) \
	    --field sphere-sphere $(IEC_R_$*) \
	    --streamer-criterion 18 --no-plot \
	    --write-to-file $@

$(BUILD)/iec60052/iec60052_%.dat: $(ROOT)/examples/iec60052/iec60052_%.dat
	@mkdir -p $(dir $@) && cp $< $@

$(BUILD)/electra/sim.dat: $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	$(PDIV) $(AIR) $(NODET) \
	    --pd-min 5E-3 --pd-max 500 --pd-num $(or $(PD_NUM),100) \
	    --field sphere-sphere 1000 \
	    --streamer-criterion 18 --no-plot \
	    --dx 5 25 0.05 \
	    --write-to-file $@

$(BUILD)/electra/tablec2_air.dat: $(ROOT)/examples/electra/tablec2_air.dat
	@mkdir -p $(dir $@) && cp $< $@

$(BUILD)/paschen/sim.dat: $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	$(PDIV) mechanisms/paschen/paschen.py mechanisms/paschen/gases.json \
	    --p 1.0 \
	    --pd-min 3e-3 --pd-max 1e2 --pd-num $(or $(PD_NUM),120) \
	    --no-plot \
	    --write-to-file $@

$(BUILD)/paschen/closed_form.dat: $(ROOT)/examples/paschen/closed_form.py $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	cd $(ROOT) && python3 examples/paschen/closed_form.py $@

# Coaxial: three wires in a 10 mm cylinder, each from 0.1 to 10 bar, so the
# pd range (bar mm) is the pressure range times b - a.  The table in
# examples/coaxial/README.md lists the same numbers.
COAX_RADII := 0.25 1 4
COAX_SIM   := $(foreach A,$(COAX_RADII),$(BUILD)/coaxial/sim_a$(A)mm.dat)
COAX_PDMIN_0.25 := 0.975
COAX_PDMAX_0.25 := 97.5
COAX_PDMIN_1 := 0.9
COAX_PDMAX_1 := 90
COAX_PDMIN_4 := 0.6
COAX_PDMAX_4 := 60

$(BUILD)/coaxial/sim_a%mm.dat: $(SOLVER_DEPS)
	@mkdir -p $(dir $@)
	$(PDIV) $(AIR) \
	    --field coaxial $* 10 \
	    --pd-min $(COAX_PDMIN_$*) --pd-max $(COAX_PDMAX_$*) \
	    --pd-num $(or $(PD_NUM),30) \
	    --dx 10 400 0.01 \
	    --no-plot \
	    --write-to-file $@

$(BUILD)/zheleznyakfit/zheleznyak_absorption_curve.dat: $(ROOT)/mechanisms/air/zheleznyak.py
	@mkdir -p $(dir $@)
	cd $(dir $@) && python3 $< --ngroups 8 --threshold 1E-6 --write-to-file zheleznyak

# `data` covers the default figures only; the opt-in ones pull their own.
data: $(BUILD)/paschen/sim.dat $(BUILD)/paschen/closed_form.dat \
      $(BUILD)/zheleznyakfit/zheleznyak_absorption_curve.dat $(COAX_SIM)

# ---- 2. pgfplots + 3. outputs ---------------------------------------------------

DATA_iec60052     := $(IEC_SIM) $(IEC_REF)
DATA_electra      := $(BUILD)/electra/sim.dat $(ELECTRA_REF)
DATA_zheleznyakfit := $(BUILD)/zheleznyakfit/zheleznyak_absorption_curve.dat
DATA_paschen      := $(BUILD)/paschen/sim.dat $(BUILD)/paschen/closed_form.dat
DATA_coaxial      := $(COAX_SIM)

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
$(foreach f,$(ALL_FIGURES),$(eval $(call FIGURE_RULES,$(f))))

clean:
	rm -rf $(BUILD) $(OUT)
