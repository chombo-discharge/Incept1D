.. _Chap:QuickStart:

Quick start
===========

This page walks through the three most common calculations.  All commands
are run from the repository root.

An inception curve in a uniform field
-------------------------------------

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py --p 1.0 --pd-min 1e-2 --pd-max 1e3

This loads the dry-air mechanism ``mechanisms/air/air_pancheshnyi.py`` with its default
(baseline) configuration, sweeps :math:`pd` logarithmically from
:math:`10^{-2}` to :math:`10^{3}` bar·mm at a fixed pressure of 1 bar, and for
every :math:`pd` finds the reduced field :math:`E/N` at which the inception
determinant :math:`\det\bm{Q}(\lambda{=}0)` vanishes (see
:ref:`Chap:InceptionCriterion`).  A two-panel figure is shown: the breakdown
voltage :math:`U(pd)` and the corresponding reduced field :math:`E/N(pd)`.
A table with the same numbers is printed to the terminal.

To write the curve to a file instead of (or in addition to) plotting it:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py --p 1.0 --no-plot \
       --write-to-file pdiv_air_1bar.dat

The output is a tab-separated table with a self-describing header that
records the date, git commit, full command line, and one column group per
configuration and polarity.

Comparing configurations
------------------------

The physics of a mechanism can be varied without touching the Python file by
passing one or more JSON configuration files (see :ref:`Chap:Configuration`).
For example, ``mechanisms/air/nodetachment.json`` defines two configurations, the
baseline and one with both electron-detachment reactions switched off:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py mechanisms/air/nodetachment.json \
       --p 1.0 --pd-min 1 --pd-max 1e3

Both configurations are computed and plotted in the same figure, giving
the "with / without detachment" comparison directly.

A non-uniform field
-------------------

Sphere-plane and sphere-sphere gaps are supported analytically.  For a
sphere-sphere gap with sphere radius 50 mm, additionally overlaying the
streamer criterion :math:`\int\max(\alpha-\eta,0)\,dx = 18`:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py mechanisms/air/nodetachment.json \
       --pd-min 1 --pd-max 55 --pd-num 100 \
       --field sphere-sphere 50 --streamer-criterion 18

For sphere-plane gaps (which are not symmetric) both polarities are computed
and reported separately.  Tabulated field lines from an external
electrostatic solver are supported through ``--field fieldline FILE``; see
:ref:`Chap:FieldLines`.

Where to go next
----------------

* :ref:`Chap:ModulesOverview` — what every command does and its options.
* :ref:`Chap:Configuration` — the JSON configuration format.
* :ref:`Chap:Examples:IEC60052` and :ref:`Chap:Examples:Electra` — full
  comparisons against reference breakdown data.
