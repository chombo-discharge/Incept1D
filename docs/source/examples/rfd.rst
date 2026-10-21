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

Running the calculation
-----------------------

The two radii of the high-field region in the paper are two configurations,
so one invocation produces both curves.  At fixed pressure :math:`pd` is
:math:`L`; the range is :math:`L` = 25 m to 5 km at sea level:

.. code-block:: bash

   incept1d pdiv mechanisms/air/rfd/air_rfd.py mechanisms/air/rfd/radius.json \
       --p 1.01325 --T 293.15 \
       --pd-min 2.5e4 --pd-max 5e6 --pd-num 60 \
       --write-to-file examples/rfd/sim.dat

``examples/rfd/normal_modes.py`` tabulates the normal-mode solution of the
same model (:eq:`eq_rfd_normal_modes`) from the same coefficients, so the
comparison checks the solver against the model rather than the model against
anything else.

Result
------

.. _fig_rfd:

.. figure:: ../figures/rfd.*
   :width: 100%
   :align: center

   Self-sustained field of relativistic feedback against the length of the
   high-field region, at sea-level density.  Lines are the ``Incept1D``
   solution; circles are the normal-mode solution of the same model.
   Dashed: the runaway avalanche threshold :math:`E_\mathrm{t}`; dotted:
   conventional breakdown, 3 MV/m.

The two agree to solver tolerance; the test suite asserts a relative error
below :math:`10^{-8}` with both criteria, and that the field reaches 3 MV/m
between 22 and 24 m, where the paper places it.  The curve has two limits:

* **Short regions** need conventional breakdown fields: the loop has no
  room to close, and :math:`E_\mathrm{s}` reaches 3 MV/m near
  :math:`L = 23` m.
* **Long regions** approach 305.4 kV/m, the field where positrons start to
  make runaway electrons by Bhabha scattering, slightly above the avalanche
  threshold :math:`E_\mathrm{t} = 284` kV/m.

A region 100 m across needs a slightly higher field than a wide one,
because photons escape it sideways; the difference is a few per cent.

Comparison with the paper
-------------------------

The figure above checks the solver against its own model.  The paper's
published results are a separate check: the authors' code and data are on
Zenodo, `doi:10.5281/zenodo.22983026
<https://doi.org/10.5281/zenodo.22983026>`_ (CC-BY 4.0), and are not
distributed with this repository.  With the archive unpacked,

.. code-block:: bash

   python3 examples/rfd/paper_comparison.py \
       PATH/TO/Faerder2026b_openresearch/selfsustainend-rfd-threshold-by-nme \
       examples/rfd/paper.dat

reads two kinds of their results, solves the same thresholds with
``Incept1D`` at the same :math:`L`, and writes both to one table:

* **Normal modes** (wide region): the steady-state threshold of the same
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

The comparison figure is not built with the documentation, because the data
is not in the repository; with ``paper.dat`` in place,

.. code-block:: bash

   make -C docs/figures rfdpaper

writes ``docs/source/figures/rfdpaper.pdf`` and ``.png``.

Altitude
--------

Every coefficient scales with the air density, so at altitude the field at
threshold scales with the density too, and :math:`E/N` against :math:`pd` is
the same curve.  Lowering ``--p`` shows it.
