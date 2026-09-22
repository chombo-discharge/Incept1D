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

``Incept1D`` is a small Python package with a command-line front end.  Given a
plasma-chemistry *mechanism file* it can

* compute inception curves, i.e. generalized Paschen curves (inception voltage and reduced field
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

   introduction/documentation
   introduction/prerequisites
   introduction/obtaining
   introduction/installation
   introduction/quickstart

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
   theory/inceptioncriterion
   theory/airscheme

Numerics
********

.. toctree::
   :maxdepth: 3
   :caption: Numerics
   :hidden:

   numerics/overview
   numerics/propagator
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
   modules/createchombodischargedata
   modules/fielddistributions
   modules/reactions
   modules/constants
   modules/fieldlines

Examples
********

.. toctree::
   :maxdepth: 3
   :caption: Examples
   :hidden:

   examples/paschen
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
