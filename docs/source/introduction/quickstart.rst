.. _Chap:QuickStart:

Quick start
===========

Two complete calculations, from opposite ends of the range: one whose answer
is known in closed form, and one that is a genuine engineering case.  All
commands are run from the repository root.

The classical Paschen curve
---------------------------

Start here, because the answer can be checked by hand.  ``mechanisms/paschen``
is the textbook Townsend model — one ionizing reaction, no attachment, no
photoionization, secondary emission at the cathode — for helium, argon and
air:

.. code-block:: bash

   incept1d pdiv mechanisms/paschen/paschen.py mechanisms/paschen/gases.json --p 1.0

``gases.json`` holds three configurations, one per gas, so all three curves
are computed and plotted together.  A two-panel figure appears — the
inception voltage :math:`U(pd)` and the reduced field :math:`E/N(pd)` — and
the same numbers are printed as a table.

Each curve should show the familiar Paschen shape: a steep left branch, a
minimum, and a slowly rising right branch.  The minima are the check:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Gas
     - :math:`U_\mathrm{min}`
     - :math:`(pd)_\mathrm{min}`
   * - Helium
     - 142 V
     - 0.056 bar·mm
   * - Argon
     - 188 V
     - 0.014 bar·mm
   * - Air
     - 305 V
     - 0.011 bar·mm

Those are not fitted values: for this model both follow from
:math:`\ln(1 + \gamma^{-1})` in closed form, and
:ref:`Chap:Examples:Paschen` compares the whole curve against the analytic
result.  If your run reproduces them, the installation is sound.

Dry air in a sphere gap
-----------------------

Now a real case.  Dry air is electronegative, so attachment, detachment and
ion conversion all matter, and a realistic gap is not uniform.  The following
sweeps a sphere-sphere gap of 50 mm sphere radius and overlays the streamer
criterion :math:`\int\max(\alpha-\eta,0)\,dx = 18`:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py \
       --pd-min 1 --pd-max 55 --pd-num 100 \
       --field sphere-sphere 50 --streamer-criterion 18

Three things are worth looking at in the result.

* **The field profile is not flat.**  ``--field`` normalises the on-axis
  profile so that :math:`\int_0^1 f\,d\xi = 1`, which makes :math:`E/N` in the
  output the *gap-average* reduced field; the field at the sphere surface is
  higher.  Use ``incept1d field --field sphere-sphere 50 --d 20`` to see the
  profile itself.
* **The streamer criterion is a different curve.**  It is the classical
  avalanche-size condition, plotted alongside the inception criterion rather
  than instead of it.  Where the two separate, the discharge is not
  streamer-limited — that comparison is the subject of
  :ref:`Chap:IonizationIntegral`.
* **Polarity.**  A sphere-sphere gap is symmetric, so both polarities give
  the same answer and the curves coincide.  Repeat with ``--field
  sphere-plane 50`` and they separate, because the high-field electrode is
  then either the anode or the cathode.

To keep the numbers instead of the figure:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py \
       --pd-min 1 --pd-max 55 --pd-num 100 \
       --field sphere-sphere 50 --streamer-criterion 18 \
       --no-plot --write-to-file air_sphere50.dat

The output is a tab-separated table with a self-describing header recording
the date, the git commit, the full command line, and one column group per
configuration and polarity.

Varying the physics
-------------------

Neither example edited a Python file.  A mechanism's parameters — which
cross-section database to use, which reactions to scale or switch off, the
secondary-emission yield — are changed by passing JSON configuration files
after the mechanism, and each becomes its own curve in the same figure:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py mechanisms/air/nodetachment.json \
       --p 1.0 --pd-min 1 --pd-max 1e3

``nodetachment.json`` holds the baseline and a variant with both
electron-detachment reactions switched off, so this produces the
"with / without detachment" comparison directly.  See
:ref:`Chap:Configuration` for the format and
:ref:`Chap:ConfigurationOverview` for the files a mechanism is made of.

Where to go next
----------------

* :ref:`Chap:ModulesOverview` — What every command does and its options.
* :ref:`Chap:FieldLines` — Using a tabulated field line from an external
  electrostatic solver, curvature included.
* :ref:`Chap:Examples:Paschen` — The closed-form verification in full.
* :ref:`Chap:Examples:IEC60052` and :ref:`Chap:Examples:Electra` —
  Comparisons against reference breakdown data.
