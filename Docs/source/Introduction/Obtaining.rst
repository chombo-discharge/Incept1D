.. _Chap:Obtaining:

How to obtain
=============

The source code is hosted on GitHub at |repo|.  Clone it with

.. code-block:: bash

   git clone git@github.com:rmrsk/Incept1D.git

or, over HTTPS,

.. code-block:: bash

   git clone https://github.com/rmrsk/Incept1D.git

The repository contains everything needed to run the examples in this
documentation: the solver scripts, the dry-air mechanism files with their
swarm data, a set of ready-made configuration files, and the
documentation sources.

Repository layout
-----------------

.. code-block:: text

   Incept1D/
   ├── Inception.py                  Core solver and Paschen-curve CLI
   ├── Eigenvalues.py                Local eigenvalues of R V^-1
   ├── IonizationIntegral.py         Ionization integrals vs. voltage
   ├── Lambda.py                     Temporal growth rate above threshold
   ├── CreateChomboDischargeData.py  Transport tables for chombo-discharge
   ├── FieldDistributions.py         Gap geometry / field profiles
   ├── Reactions.py                  Reaction-string parser, R-matrix assembly
   ├── Constants.py                  Physical constants
   ├── Air/                          Dry-air mechanism family
   │   ├── Air_Pancheshnyi.py        Three-body attachment scheme (reference)
   │   ├── Air_2body.py              Explicit Bloch-Bradbury two-body scheme
   │   ├── Config.py                 Configuration protocol for the Air family
   │   ├── Zheleznyak.py             Two-stream fit to the Zheleznyak model
   │   ├── *.json                    Ready-made configuration sets
   │   └── *.txt                     BOLSIG+ swarm data, LXCat ion mobilities
   ├── Examples/                     Reference data for the worked examples
   ├── Docs/                         This documentation (Sphinx)
   └── tests/                        Test suite (pytest; see Maintenance)

Versioning
----------

There are no tagged releases yet; the ``main`` branch is the reference
version and is the one exercised by continuous integration
(:ref:`Chap:Infrastructure`).  Every data file written by the scripts records
the git commit it was produced with in its header, so results can always be
traced back to a specific version of the code.
