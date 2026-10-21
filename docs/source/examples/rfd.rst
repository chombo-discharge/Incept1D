.. _Chap:Examples:RFD:

Relativistic feedback discharge
===============================

.. contents:: On this page
   :local:
   :depth: 1

How strong must the field in a thundercloud be for relativistic feedback to
sustain itself?  Relativistic runaway electrons avalanche along the field
and emit photons forward and backward; forward photons make positrons, and
positrons and backward photons run back and make new runaway electrons.
Above a field :math:`E_\mathrm{s}(L)`, set by the length :math:`L` of the
high-field region, this loop sustains itself with nothing coming in.

The mechanism ``mechanisms/air/rfd/air_rfd.py`` is the relativistic fluid
model of [Faerder2026b]_, described in :ref:`Chap:RFDScheme`.  Its
self-sustained condition is the inception condition of ``Incept1D`` with
every electrode yield zero, so ``incept1d pdiv`` computes
:math:`E_\mathrm{s}(L)` directly.

Open ends, not electrodes
-------------------------

There are no electrodes here.  The gap is the high-field region itself,
:math:`0 \le x \le L`, with the runaway electrons drifting toward
:math:`+x`; outside it the field is too weak for them to multiply, so
whatever leaves the region is lost and nothing comes in.  Those are exactly
the boundary conditions of :ref:`Chap:SecondaryEmission` with every yield
zero: at :math:`x = 0` no runaway electrons and no forward photons enter,
at :math:`x = L` no positrons and no backward photons.  The output keeps
the names of a gap: ``d`` is :math:`L`, :math:`pd` is the density times
:math:`L`, and the voltage is :math:`E_\mathrm{s} L`, the potential drop
across the region.

The model is one-dimensional, as the paper's is.  The horizontal radius
:math:`R` of the region enters only through the photon attenuation lengths,
which the authors fitted with GEANT4 to include the photons that escape
sideways, for :math:`R` = 5 km and 100 m.  The paper calls the 5 km case
:math:`R \gg L`; that holds for regions up to about a kilometre long, and
the same attenuation lengths are used, by the paper and here, over the
whole range of :math:`L`.

Running the calculation
-----------------------

The two radii of the high-field region in the paper are two configurations,
so one invocation produces both curves.  The paper's results are all at sea
level, which is the authors' sea-level density
:math:`2.688\times10^{25}` m\ :sup:`-3`: ``--p 1.01371 --T 273.15``.  At
fixed pressure :math:`pd` is :math:`L`; the range is :math:`L` = 10 m to
12 km, as in the paper:

.. code-block:: bash

   incept1d pdiv mechanisms/air/rfd/air_rfd.py mechanisms/air/rfd/radius.json \
       --p 1.01371 --T 273.15 \
       --pd-min 1e4 --pd-max 1.2e7 --pd-num 80 \
       --write-to-file examples/rfd/sim.dat

The test suite checks the solver against the normal-mode solution of the
same model (:eq:`eq_rfd_normal_modes`), transcribed independently, to a
relative error below :math:`10^{-8}` with both criteria.

Comparison with the paper
-------------------------

The paper's published results are not distributed with this repository.
The authors' code and data are on Zenodo, `doi:10.5281/zenodo.22983026
<https://doi.org/10.5281/zenodo.22983026>`_ (CC-BY 4.0).  With the archive
unpacked,

.. code-block:: bash

   python3 examples/rfd/paper_comparison.py \
       PATH/TO/Faerder2026b_openresearch/selfsustainend-rfd-threshold-by-nme \
       examples/rfd/paper.dat

reads two kinds of their results, solves the same thresholds with
``Incept1D`` at the same :math:`L`, and writes both to one table.  With
``paper.dat`` in place, the figure below adds their results and a panel
with the difference; without it, the figure shows the ``Incept1D`` curves
alone, as in the published documentation.

* **Normal modes** (:math:`R` = 5 km): the steady-state threshold of the same
  model, found on a 100 V/m grid of fields.  ``Incept1D`` agrees with it to
  within that grid, 0.05 kV/m or :math:`1.5\times10^{-4}`, for every
  :math:`L` from 10 m to 5 km; beyond that their minimum of the determinant
  is too flat to resolve as finely, and the difference stays below
  :math:`6\times10^{-4}`.  This is the like-for-like comparison, and the test
  suite asserts it when ``paper.dat`` is present.
* **Time-dependent fluid model** (both radii): the authors' simulation, with
  the threshold found by bisection on whether the runaway density keeps
  growing.  It lies 0.1–2 % below ``Incept1D``, as it does below their own
  normal modes.  Beyond about 5 km it falls to the avalanche threshold
  :math:`E_\mathrm{t}`, which their normal modes do not; their growth tests
  judge a finite simulation, in which a single avalanche across so long a
  region already grows strongly, and the steady-state threshold cannot
  lie below the 305.4 kV/m where positrons start to feed back.

Result
------

.. _fig_rfd:

.. figure:: ../figures/rfd.*
   :width: 100%
   :align: center

   Self-sustained field of relativistic feedback against the length of the
   high-field region, at sea-level density.  Lines are the ``Incept1D``
   solution for both radii.  Dashed: the runaway avalanche threshold
   :math:`E_\mathrm{t}`; dotted: conventional breakdown, 3 MV/m.  When the
   reader has supplied the paper's results, circles are their normal modes,
   dotted lines their time-dependent model, and the lower panel the
   difference from ``Incept1D`` at the same :math:`L`.

The field reaches 3 MV/m between 22 and 24 m, where the paper places it.
The curve has two limits:

* **Short regions** need conventional breakdown fields: the loop has no
  room to close, and :math:`E_\mathrm{s}` reaches 3 MV/m near
  :math:`L = 23` m.
* **Long regions** approach 305.4 kV/m, the field where positrons start to
  make runaway electrons by Bhabha scattering, slightly above the avalanche
  threshold :math:`E_\mathrm{t} = 284` kV/m.

A region of 100 m radius needs a slightly higher field than one of 5 km,
because photons escape it sideways; the difference is a few per cent.

Altitude
--------

Lowering ``--p`` lowers the relative density :math:`\tilde{n}`.  The model
is then the authors' (:ref:`Chap:RFDScheme`): the avalanche thresholds scale
with :math:`\tilde{n}`, and every other coefficient is its sea-level
function of the field times :math:`\tilde{n}`.  The paper publishes
sea-level results only, so nothing off sea level is checked against it.
