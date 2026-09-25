.. Incept1D documentation master file.

Welcome to the ``Incept1D`` user documentation
***********************************************

.. important::

   ``Incept1D`` decides whether a gas discharge **starts**.  It solves a
   one-dimensional drift-reaction model of a discharge gap — electrons,
   positive and negative ions, two-stream photoionization and cathode
   secondary emission — and returns the voltage at which the gap breaks
   down.  The code is hosted at `GitHub
   <https://github.com/chombo-discharge/Incept1D>`_ together with the source
   files for this documentation.

.. attention::

   If you publish results obtained with ``Incept1D``, you must cite the paper
   that presents the model — see :ref:`Chap:Citing`.

What it does
============

.. list-table::
   :header-rows: 1
   :widths: 26 46 28

   * - Command
     - Answers
     - See
   * - ``incept1d pdiv``
     - **Inception curves.**  The inception voltage and reduced field against
       :math:`pd` — a generalized Paschen curve.
     - :ref:`Chap:Inception`
   * - ``incept1d growth``
     - **Growth rate.**  How fast the discharge grows once the voltage exceeds
       inception.
     - :ref:`Chap:Lambda`
   * - ``incept1d ionization``
     - **Ionization integrals.**  The classical and apparent integrals, for
       comparison against the full criterion.
     - :ref:`Chap:IonizationIntegral`
   * - ``incept1d eigenvalues``
     - **Local growth modes.**  Eigenvalues of the reaction-transport matrix
       :math:`\bm{R}\bm{V}^{-1}` against :math:`E/N`.
     - :ref:`Chap:Eigenvalues`
   * - ``incept1d field``
     - **Gap geometry.**  The normalised field profile of a gap, before you
       spend time on a sweep.
     - :ref:`Chap:FieldDistributions`
   * - ``incept1d chombo``
     - **Transport tables.**  Export for the 3-D plasma solver
       `chombo-discharge <https://github.com/chombo-discharge/chombo-discharge>`_.
     - :ref:`Chap:ThirdParty`

Features
========

**Any gas.**  The chemistry lives in a *mechanism file* outside the package,
so a new gas is a new file rather than a patch to the solver
(:ref:`Chap:NewMechanisms`).  Dry air is shipped, with several ready-made
variants.

**Any gap geometry.**  Uniform, sphere-plane, sphere-sphere and coaxial gaps
are analytic.  A **tabulated field line** from an external electrostatic solver
is read straight from file, curvature included (:ref:`Chap:FieldLines`).

**A criterion that does not assume the answer.**  Rather than asking whether
an avalanche reaches a critical size, the solver evaluates
:math:`\det\bm{Q}(\lambda) = 0` for the whole coupled system, so electrode
feedback and volume feedback are on the same footing
(:ref:`Chap:InceptionCriterion`).

**Electronegative gases done properly.**  Attachment, detachment and ion
conversion are carried as species, so the inception field is set by the
apparent ionization coefficient :math:`\lambda_{\max}` and not by
:math:`\alpha = \eta` alone.

**Photon feedback along the field line.**  A multigroup two-stream model
carries the ionizing photons, so photoionization in the gas and the
photoelectric effect at the cathode are both part of the system rather than
lumped into an effective yield (:ref:`Chap:Photoionization`).

**Answers you can check.**  Every result carries a header recording the git
commit and the command that produced it, and the solver is verified against
closed-form limits rather than against stored output
(:ref:`Chap:Examples:Paschen`, :ref:`Chap:TestSuite`).

Getting started
===============

.. code-block:: bash

   git clone https://github.com/chombo-discharge/Incept1D.git
   cd Incept1D && pip install -e .
   incept1d pdiv mechanisms/paschen/paschen.py mechanisms/paschen/gases.json --p 1

:ref:`Chap:Installation` covers the install in full, and
:ref:`Chap:QuickStart` walks through two complete calculations.

How this documentation is organised
===================================

* **Introduction** — What you need, where to get the code, how to run a
  first calculation, and how to cite it.
* **Theory** — The drift-reaction model, two-stream photoionization,
  secondary emission and the determinant inception criterion.  Every
  equation is cross-referenced to the function that implements it.
* **Numerics** — How the propagator is built, how the determinant is
  evaluated robustly, how roots are found and tracked, and what a
  calculation costs.
* **Configuration files** — The mechanism module, the ``config.py``
  protocol and the JSON format: what a gas needs and the interfaces it must
  satisfy.
* **Python modules** — What each module and command does, and its
  command-line interface.
* **Examples** — Worked calculations: the classical Paschen curve, and
  comparisons against the IEC 60052 sphere-gap standard and the Dakin
  *et al.* (ELECTRA) breakdown compilation.
* **Maintenance** — The test suite, continuous integration, and how to
  contribute.

The API reference is generated from the NumPy-style docstrings in the code.
The version at `chombo-discharge.github.io/Incept1D
<https://chombo-discharge.github.io/Incept1D/>`_ is rebuilt by continuous
integration from every commit to ``main``; see :ref:`Chap:Infrastructure`
for how to build it locally.

Notation
========

Throughout the theory pages, bold upright symbols (:math:`\bm{R}`,
:math:`\bm{V}`, :math:`\bm{Q}`, :math:`\bm{\Pi}`) are matrices and arrow
symbols (:math:`\vec{n}`, :math:`\vec{\theta}`, :math:`\vec{\gamma}`) are
column vectors.  Subscripts name the species or stream a quantity belongs to,
and are used consistently: :math:`\mathrm{e}` for electrons, :math:`+` and
:math:`-` for positive and negative ions, and :math:`\Psi^+` or
:math:`\Psi^-` for the forward and backward photon streams — so
:math:`\bm{\Pi}_-` selects negative ions while :math:`\bm{\Pi}_{\Psi^-}`
selects backward photons.  A subscript :math:`j` always indexes a photon
group.  Reduced fields
:math:`E/N` are given in Townsend
(:math:`1\,\mathrm{Td} = 10^{-21}\,\mathrm{V\,m^2}`), pressures in bar, gap
distances in mm, and the product :math:`pd` in bar·mm.  Inside the code all
quantities are SI except where a function docstring explicitly says
otherwise; the command-line interfaces accept bar, mm, kV and Td.

This documentation was built from commit |commit|.

.. This is for getting rid of the TOC in html view.
.. raw:: html

   <style>
   section#introduction,
   section#theory,
   section#numerics,
   section#configuration-files,
   section#python-modules,
   section#examples,
   section#maintenance,
   section#references {
	 display:none;
   }
   </style>

.. only:: latex

   .. toctree::
      :caption: Contents

Introduction
************

.. toctree::
   :maxdepth: 3
   :caption: Introduction
   :hidden:

   introduction/prerequisites
   introduction/obtaining
   introduction/installation
   introduction/quickstart
   introduction/citing

Theory
******

.. toctree::
   :maxdepth: 3
   :caption: Theory
   :hidden:

   theory/overview
   theory/transport
   theory/photoionization
   theory/secondaryemission
   theory/augmented
   theory/inceptioncriterion

Numerics
********

.. toctree::
   :maxdepth: 3
   :caption: Numerics
   :hidden:

   numerics/overview
   numerics/propagator
   numerics/riccati
   numerics/determinant
   numerics/rootfinding
   numerics/cost

Configuration files
*******************

.. toctree::
   :maxdepth: 3
   :caption: Configuration files
   :hidden:

   configuration/overview
   configuration/mechanismfile
   configuration/configpy
   configuration/jsonfiles
   configuration/reactions
   configuration/examples/index

Python modules
**************

.. toctree::
   :maxdepth: 3
   :caption: Python modules
   :hidden:

   modules/overview
   modules/inception
   modules/eigenvalues
   modules/ionizationintegral
   modules/lambda
   modules/fielddistributions
   modules/reactions
   modules/constants
   modules/thirdparty

Examples
********

.. toctree::
   :maxdepth: 3
   :caption: Examples
   :hidden:

   examples/paschen
   examples/fieldline
   examples/coaxial
   examples/iec60052
   examples/electra

Maintenance
***********

.. toctree::
   :maxdepth: 3
   :caption: Maintenance
   :hidden:

   maintenance/testsuite
   maintenance/infrastructure
   maintenance/contributing

References
**********

.. toctree::
   :maxdepth: 3
   :caption: References
   :hidden:

   zzreferences
