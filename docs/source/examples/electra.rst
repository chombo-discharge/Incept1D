.. _Chap:Examples:Electra:

Quasi-uniform gaps: the ELECTRA Paschen curve
=============================================

.. contents:: On this page
   :local:
   :depth: 1

Dakin *et al.* [Dakin1974]_ compiled breakdown measurements in dry air
from several laboratories, obtained in plane-plane, sphere-plane and
sphere-sphere gaps, into a single Paschen curve spanning :math:`pd` from
about :math:`3\times10^{-3}` to 400 bar·mm (Table C2 of that report).
This example compares the inception criterion against that compilation
over five decades of :math:`pd`, from the Paschen minimum to the regime
where negative-ion detachment dominates.

Reference data (not included)
-----------------------------

.. important::

   The ELECTRA breakdown voltages are copyrighted by CIGRE and are **not
   distributed with this repository**.  The figure below is therefore built
   from the computed curve alone.

If you have access to the article you can reproduce the full comparison
locally.  Place ``tablec2_air.dat`` in ``examples/electra/``, with columns
:math:`pd` (bar·mm) and breakdown voltage (kV, crest value) at 20 °C; ``#``
comment lines are ignored.  The overlays in ``docs/figures/electra.tex`` are
guarded by ``\IfFileExists``, so the next ``make -C docs figures`` picks the
file up automatically.  See ``examples/electra/README.md``.

Running the calculation
-----------------------

The compilation mixes geometries, so a single nearly uniform geometry is
used for the calculation: a sphere-sphere gap with :math:`R = 1000` mm at
:math:`p = 1` bar, sweeping the gap length.  As in the previous example
the no-detachment variant and the streamer criterion are computed
alongside.  At small :math:`pd` the propagator is cheap, so the adaptive
grid budget is reduced (``--dx 5 25 0.05``) to keep the run short:

.. literalinclude:: ../../../examples/electra/run.sh
   :language: bash
   :start-at: incept1d pdiv

Reading the output
------------------

``examples/electra/sim.dat`` contains, per configuration, the reduced
field ``EN_Td`` and voltage ``U_kV`` columns.  The voltage panel of the
figure plots ``U_kV[Baseline]``, ``U_kV[No detachment]`` and
``U_kV[Streamer (C=18.0), Baseline]`` against column 1; the field panel
plots the corresponding ``EN_Td`` columns.  For the measured data the mean
reduced field is recovered from the voltage as
:math:`E/N = U/(pd\,N_0)` with :math:`N_0 = 2.47\times10^{25}`
m\ :sup:`-3` bar\ :sup:`-1` at 20 °C.

Result
------

As for the sphere-gap example, the figure is built during the documentation
build from the output of ``run.sh`` and the pgfplots source
``docs/figures/electra.tex``.

As for the sphere gaps, the figure is **not built with the documentation**:
without ``tablec2_air.dat`` there is nothing to compare against.  Supply the
file as described above and build it explicitly:

.. code-block:: bash

   make -C docs/figures electra

which writes ``docs/source/figures/electra.pdf`` and ``.png``: a) breakdown
voltage and b) mean reduced field vs. :math:`pd`, with the ELECTRA values as
symbols and the inception criterion as a line, the no-detachment calculation
dashed and the streamer criterion dotted; the insets show
:math:`pd > 100` bar·mm on a linear scale.

The calculation follows the measured curve across the whole range: it
slightly under-predicts the data by about 5 % at :math:`pd > 1` bar·mm and
over-predicts near the Paschen minimum.  The streamer criterion, by
contrast, over-predicts everywhere — by a factor of two at
:math:`pd = 2\times10^{-2}` bar·mm, where the inception criterion agrees
with the data to within a few percent — because a fixed number of
e-foldings is a poor description of inception in short gaps, where
secondary emission from the cathode, not avalanche size, controls the
onset.  At large :math:`pd` the no-detachment curve departs from the data
in the same way as in the sphere-gap example: the measured mean field
keeps falling with :math:`pd` (panel b) because detachment lowers the
effective critical field below :math:`\alpha = \eta`.

Variations to try
-----------------

* ``--all-branches --plot-separate-branches`` reveals additional roots of
  the determinant at low :math:`pd`; the lowest branch is the physical
  inception voltage.
* ``mechanisms/air/paschen.json`` reduces the chemistry to the textbook limit
  (:eq:`eq_standard_paschen`); the difference to the baseline isolates
  the combined effect of ion conversion, detachment and photon feedback.
* ``mechanisms/air/see.json`` shows the sensitivity to the cathode yield and to the
  photon cone angle, which matters most near the Paschen minimum.
