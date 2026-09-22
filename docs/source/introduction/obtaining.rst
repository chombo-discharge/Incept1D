.. _Chap:Obtaining:

How to obtain
=============

The source code is hosted on GitHub at |repo|.  Clone it with

.. code-block:: bash

   git clone git@github.com:chombo-discharge/Incept1D.git

or, over HTTPS,

.. code-block:: bash

   git clone https://github.com/chombo-discharge/Incept1D.git

The repository contains everything needed to run the examples in this
documentation: the ``incept1d`` package, the dry-air mechanism files with their
swarm data, a set of ready-made configuration files, and the
documentation sources.

Repository layout
-----------------

.. code-block:: text

   Incept1D/
   ├── pyproject.toml                Package metadata, dependencies, `incept1d` entry point
   ├── src/incept1d/                 The Python package
   │   ├── solver.py                 Core solver: augmented ODE, propagators, det Q
   │   ├── inception.py              E/N roots of det Q and inception-curve tracking
   │   ├── mechanism.py              Mechanism / configuration loading
   │   ├── eigenvalues.py            Local eigenvalues of R V^-1
   │   ├── ionization.py             Ionization integrals across the gap
   │   ├── growth.py                 Temporal growth rate above threshold
   │   ├── chombo.py                 Transport tables for chombo-discharge
   │   ├── fields.py                 Gap geometry / field profiles
   │   ├── reactions.py              Reaction-string parser, R-matrix assembly
   │   ├── constants.py              Physical constants
   │   ├── output.py                 Metadata header of the result files
   │   └── cli/                      `incept1d <command>` front ends
   ├── mechanisms/air/               Dry-air mechanism family
   │   ├── air_pancheshnyi.py        Three-body attachment scheme (reference)
   │   ├── air_2body.py              Explicit Bloch-Bradbury two-body scheme
   │   ├── config.py                 Configuration protocol for the Air family
   │   ├── zheleznyak.py             Two-stream fit to the Zheleznyak model
   │   ├── *.json                    Ready-made configuration sets
   │   └── *.txt                     BOLSIG+ swarm data, LXCat ion mobilities
   ├── examples/                     Reference data for the worked examples
   ├── docs/                         This documentation (Sphinx)
   └── tests/                        Test suite (pytest; see Maintenance)

Versioning
----------

There are no tagged releases yet; the ``main`` branch is the reference
version and is the one exercised by continuous integration
(:ref:`Chap:Infrastructure`).  Every data file written by the commands records
the git commit it was produced with in its header, so results can always be
traced back to a specific version of the code.
