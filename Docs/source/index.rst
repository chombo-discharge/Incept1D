.. Incept1D documentation master file.

Welcome to the ``Incept1D`` user documentation
***********************************************

.. important::

   ``Incept1D`` computes the inception (breakdown) condition of a
   one-dimensional drift-reaction model of a gas discharge gap, including
   negative-ion transport and detachment, ion conversion, two-stream
   photoionization, and cathode secondary emission.  The code is hosted at
   `GitHub <https://github.com/chombo-discharge/Incept1D>`_ together with the source
   files for this documentation.

``Incept1D`` is a small collection of Python scripts, not a library.  Given a
plasma-chemistry *mechanism file* it can

* compute generalized Paschen curves (breakdown voltage and reduced field
  vs. :math:`pd`) for uniform, sphere-plane, sphere-sphere, or tabulated
  field-line geometries (:ref:`Chap:Inception`);
* inspect the local eigenvalues of the reaction-transport matrix
  :math:`\bm{R}\bm{V}^{-1}` vs. :math:`E/N` (:ref:`Chap:Eigenvalues`);
* compare the full inception criterion against the classical ionization
  integral and the streamer criterion (:ref:`Chap:IonizationIntegral`);
* compute the temporal growth rate of the discharge above the inception
  voltage (:ref:`Chap:Lambda`);
* export transport and rate-coefficient tables for the 3-D plasma solver
  `chombo-discharge <https://github.com/chombo-discharge/chombo-discharge>`_
  (:ref:`Chap:CreateChomboDischargeData`).

The underlying theory is derived in :ref:`Chap:TheoryOverview` and the
following pages, and the documentation maps every equation onto the
function that implements it.

This documentation was built from commit |commit|.

.. This is for getting rid of the TOC in html view.
.. raw:: html

   <style>
   section#introduction,
   section#theory,
   section#numerics,
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

   Introduction/Documentation
   Introduction/Prerequisites
   Introduction/Obtaining
   Introduction/Installation
   Introduction/QuickStart

Theory
******

.. toctree::
   :maxdepth: 3
   :caption: Theory
   :hidden:

   Theory/Overview
   Theory/Transport
   Theory/Photoionization
   Theory/SecondaryEmission
   Theory/InceptionCriterion
   Theory/AirScheme

Numerics
********

.. toctree::
   :maxdepth: 3
   :caption: Numerics
   :hidden:

   Numerics/Overview
   Numerics/Propagator
   Numerics/Determinant
   Numerics/RootFinding
   Numerics/Cost

Python modules
**************

.. toctree::
   :maxdepth: 3
   :caption: Python modules
   :hidden:

   Modules/Overview
   Modules/Inception
   Modules/Eigenvalues
   Modules/IonizationIntegral
   Modules/Lambda
   Modules/CreateChomboDischargeData
   Modules/FieldDistributions
   Modules/Reactions
   Modules/Constants
   Modules/Configuration
   Modules/ModifyingReactions
   Modules/NewMechanisms
   Modules/FieldLines

Examples
********

.. toctree::
   :maxdepth: 3
   :caption: Examples
   :hidden:

   Examples/IEC60052
   Examples/Electra

Maintenance
***********

.. toctree::
   :maxdepth: 3
   :caption: Maintenance
   :hidden:

   Maintenance/TestSuite
   Maintenance/Infrastructure
   Maintenance/Contributing

References
**********

.. toctree::
   :maxdepth: 3
   :caption: References
   :hidden:

   ZZReferences
