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

Versioning
----------

There are no tagged releases yet; the ``main`` branch is the reference
version and is the one exercised by continuous integration
(:ref:`Chap:Infrastructure`).  Every data file written by the commands records
the git commit it was produced with in its header, so results can always be
traced back to a specific version of the code.
