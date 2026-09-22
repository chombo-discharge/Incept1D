.. _Chap:FieldLines:

Using a precomputed field line
==============================

.. contents:: On this page
   :local:
   :depth: 1

For electrode geometries beyond spheres and planes, the field along a
field line can be computed with an external electrostatic solver (a finite
element package, ``chombo-discharge``, …), exported as a table, and used
directly with ``--field fieldline FILE [UNIT]``.  The 1-D model is then
integrated along that line exactly as it is along the axis of a sphere gap.

File format
-----------

A plain numeric table, whitespace- or comma-separated; header lines and
lines starting with ``#`` are skipped.  The layout is inferred from the
number of numeric columns:

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Columns
     - Interpretation
   * - 2
     - ``s  |E|`` — arc length and field magnitude
   * - 4
     - ``x  y  z  |E|`` — position and field magnitude
   * - 6
     - ``x  y  z  Ex  Ey  Ez`` — position and field vector

For the 4- and 6-column layouts the arc length is the cumulative
point-to-point Euclidean distance.  Rows are used in file order: the
**first row defines** :math:`\xi = 0` and the last :math:`\xi = 1`.  Only
:math:`|E|` enters the model, so the direction of the field vector and the
absolute field units are irrelevant — the profile is normalised.  ``UNIT``
(``m``, ``cm``, ``mm``, ``um``; default ``m``) is the unit of the length
columns.

Example (``line.csv``, lengths in mm):

.. code-block:: text

   # s_mm   E_kV_per_mm
   0.0      12.4
   0.5      9.8
   1.0      7.9
   ...
   20.0     3.1

How the line is used
--------------------

:meth:`FieldDistributions.FieldDistribution.from_fieldline` parametrises
the profile by normalised arc length :math:`\xi = s/L` and normalises it
so that :math:`\int_0^1 f\,d\xi = 1`.  Consequently

* the gap length is the arc length, :math:`d = L`;
* the reference field is the mean field along the line, :math:`E_\mathrm{ref}
  = U/L` with :math:`U = \int|E|\,ds` the voltage drop along the line;
* every voltage reported by the commands is this line integral.

Polarity
........

The file does not say which end is the anode.  ``incept1d pdiv`` evaluates
both polarities and labels them ``start=positive`` (the first tabulated
point is the anode) and ``start=negative`` (the first point is the
cathode).  Use the one that matches your electrode arrangement, or export
the line in the direction that makes ``start`` the electrode of interest.

Sweeps
......

* With neither ``--p`` nor ``--d``, ``incept1d pdiv`` uses :math:`d = L`
  and sweeps the pressure — the natural sweep for a fixed geometry.
* With ``--p``, the gap length :math:`d` is swept.  Since :math:`f(\xi)`
  is fixed, :math:`d \ne L` corresponds to the *same* electrode arrangement
  scaled geometrically by :math:`d/L`.

Resolution
..........

The profile is linearly interpolated between tabulated points.  Export
enough points to resolve the field near the high-field electrode; the
adaptive integration grid (``--dx``) refines the propagator, not the
profile.  Check the imported profile visually first:

.. code-block:: bash

   incept1d field --field fieldline line.csv mm

Limitations
-----------

* Flux-tube divergence — the spreading of neighbouring field lines — is
  neglected, as for the sphere geometries.  The model is a drift-reaction
  balance *along* the line.
* The two-stream photon model tracks photons along the line only
  (:ref:`Chap:Photoionization`); in a strongly curved field it cannot
  account for photons that cross between field lines.
* Only ``|E|`` matters, so a line that reverses direction is treated as if
  the field had not reversed.  Export a line that runs from one electrode
  to the other.

Full example
------------

.. code-block:: bash

   # Field along a line exported from a 3-D electrostatic solver, mm units,
   # sweeping pressure from 0.1 to 10 bar at the line's own length
   incept1d pdiv mechanisms/Air/Air_Pancheshnyi.py --field fieldline line.csv mm \
       --pd-min 2 --pd-max 200 --pd-num 40 --streamer-criterion 18 \
       --write-to-file line_pdiv.dat
