.. _Chap:Examples:FieldLine:

A tabulated field line
======================

.. contents:: On this page
   :local:
   :depth: 1

The analytic geometries cover a sphere against a plane and a sphere against
a sphere.  Everything else — a rod, an edge, a contaminated surface, a
triple junction — needs an electrostatic solve, and what comes back is a
table of :math:`|E|` along a field line.  This example runs the inception
criterion along such a line, and shows what the units of that file do and
do not affect.

The line
--------

The line used here is synthetic, so that every number on this page can be
checked by hand: 20 mm long, with :math:`|E|` falling linearly by a factor
of two from one end to the other, scaled to an excitation of 100 kV.

.. math::

   \langle|E|\rangle = \frac{100\ \mathrm{kV}}{20\ \mathrm{mm}}
                     = 5\ \mathrm{kV/mm},
   \qquad
   |E| : 6.667 \rightarrow 3.333\ \mathrm{kV/mm}.

A linear ramp is convenient because the line integral is exact,
:math:`\int|E|\,ds = L\,\langle|E|\rangle = 100` kV, which is the number the
commands report back.  ``make_line.py`` writes it twice:

.. list-table::
   :header-rows: 1
   :widths: 34 32 34

   * - File
     - Length column
     - Field column
   * - ``line_si.csv``
     - metres
     - V/m
   * - ``line_engineering.csv``
     - millimetres
     - kV/mm

These describe the *same* line.  The second is the form an engineering
field solver tends to export, and the pair exists to make one point
concrete: the file's units must not reach the answer.

What the units do and do not affect
-----------------------------------

Only the **shape** of the profile enters the solve.
:meth:`~incept1d.fields.FieldDistribution.from_fieldline` normalises the
tabulated field so that :math:`\int_0^1 f\,d\xi = 1`, which divides out the
absolute scale entirely.  The gap length is the arc length, and that *is*
unit-dependent — which is why the length unit is declared on the command
line (``m``, ``cm``, ``mm``, ``um``).

The field column is a different matter: the file does not say what unit it
is in, and nothing in the table reveals it.  A profile is therefore enough
to compute the inception voltage :math:`U^*`, but not enough to say how far
the *supplied* excitation is from inception.  That ratio needs the
excitation itself, which is what ``--fieldline-voltage`` declares:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py \
       --field fieldline examples/fieldline/line_engineering.csv mm \
       --fieldline-voltage 100 \
       --pd-min 1 --pd-max 100 --pd-num 40

The sweep is asked for in :math:`pd`, as for any other geometry, and comes
back in :math:`p`.

Because the declared voltage is divided into a computed :math:`U^*`, and
:math:`U^*` comes from the normalised shape, the resulting ratio is exact
whatever the field column happens to be in.  Without the flag the ratio is
not reported at all, rather than reported wrongly.

``∫|E| ds`` is still printed, in the file's own units, and
``incept1d field`` turns it into a statement you can check.  Its ratio to
the declared excitation *is* the field unit expressed in V/m, so an
unlabelled column can be identified rather than guessed:

.. code-block:: console

   $ incept1d field --field fieldline examples/fieldline/line_engineering.csv mm \
         --fieldline-voltage 100
   Geometry:  field line line_engineering.csv, L = 20 mm,  d = 20.0 mm
   ∫|E| ds  = 0.1   (in the field units of the file × m)
   U_applied = 100 kV   (--fieldline-voltage)
              consistent with |E| tabulated in kV/mm or MV/m (×1e+06 V/m); the solve is unaffected either way

The SI file reports ``consistent with |E| tabulated in V/m`` instead.  A
ratio matching no common unit is the useful failure: it means the length
unit is wrong, the wrong column was exported, or the line does not span the
whole gap — the last of which is a modelling mistake rather than a units
one, and the only one of the three that changes the answer.

Running the calculation
-----------------------

``examples/fieldline/make_line.py`` writes the two files, and the solve
along the engineering-unit copy is

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py \
       --field fieldline examples/fieldline/line_engineering.csv mm \
       --fieldline-voltage 100 \
       --pd-min 1 --pd-max 100 --pd-num 40 \
       --write-to-file examples/fieldline/sim_engineering.dat

Running the same command on ``line_si.csv m`` gives the SI companion.

This is a **pressure sweep**, and for a tabulated field line it is the only
sweep that makes sense — which is also why the results come back against
:math:`p` in bar rather than :math:`pd`: with :math:`d` pinned at
:math:`L`, the two differ by a constant.  The line describes one geometry at one size: its
arc length fixes the gap at :math:`d = L = 20` mm, and :math:`pd` is varied
by varying :math:`p` alone — here from 0.05 to 5 bar across the requested
1 to 100 bar·mm.  Neither ``--p`` nor ``--d`` is given, which is what
selects this mode; ``incept1d pdiv`` says so on startup:

.. code-block:: text

   --field fieldline: no --p/--d given, using d = L = 20 mm (arc length of the tabulated line).

Fixing ``--p`` would make :math:`d` the swept variable instead, so every
point would be a differently sized copy of the imported arrangement.  There
is no useful reading of that sweep, so it is refused:

.. code-block:: text

   incept1d pdiv: error: --p cannot be used with --field fieldline: it would sweep the
   gap length, and every d ≠ L = 20 mm is the electrode arrangement at a different
   size.  Omit --p for a pressure sweep at d = L, or for the single point p = 1 bar
   use --pd-min 20 --pd-max 20 --pd-num 1.

Passing a ``--d`` other than :math:`L` is allowed — :math:`f(\xi)` is
invariant under a geometric rescaling, so it solves the same arrangement at
a different size — but it says so, since that is a modelling decision
rather than a grid one.

The line is not symmetric — :math:`|E|` is highest at the first tabulated
point — so both polarities are computed: ``start=positive`` means the first
point is the anode.

Result
------

Each run opens a two-panel figure — the inception voltage :math:`U^*` and
the reduced field :math:`(E/N)^*` against the pressure, one curve per
polarity — and prints the table behind it.  Add ``--no-plot``, or set
``MPLBACKEND=Agg``, for a headless run; the tables are written either way.

The two runs are compared directly:

.. code-block:: text

           p (bar)      E/N (Td)        U (kV)         E (V/m)    U*/U_applied
        5.0000e-02      164.8633        4.0754      2.0377e+05          0.0408   start=positive
        5.0000e-02      164.3114        4.0618      2.0309e+05          0.0406   start=negative

Both files give this same table — the two output files are identical to the
last digit, which is the property worth having and what the test suite
asserts.  Had the field column been read as a physical quantity rather than
a shape, the ``kV/mm`` file would have disagreed with the SI one by a factor
of :math:`10^6`.

The small split between the polarities is real: with a monotonically
falling profile the high-field end is either the anode or the cathode, and
secondary emission happens at the cathode.  It is the same effect as the
sphere-plane polarity difference, and it disappears for a uniform gap.

At :math:`pd = 1` bar·mm the inception voltage is about 4 % of the declared
100 kV, so the supplied excitation is far above inception — the ratio is
what tells you that at a glance, in whatever units the file arrived in.

Variations to try
-----------------

* **Your own line.**  Export :math:`|E|` along a field line from any
  electrostatic solver as ``s |E|``, ``x y z |E|`` or ``x y z Ex Ey Ez``;
  see :ref:`Chap:FieldLines` for the accepted layouts.
* **Check it first.**  ``incept1d field --field fieldline FILE mm`` plots
  the normalised profile and prints the line integral, which catches a
  units or column mistake before a sweep is spent on it.
* **Reverse the line.**  Exporting from the other electrode swaps which
  polarity is which, and is a quick way to confirm that the polarity
  labelling matches your electrode arrangement.
