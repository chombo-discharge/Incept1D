.. _Chap:FieldDistributions:

Gap geometry — ``incept1d.fields``
==================================

.. contents:: On this page
   :local:
   :depth: 1

``incept1d field`` abstracts the electrode geometry away from the
solvers.  Everything downstream only ever needs a normalised field profile
:math:`f(\xi)` on :math:`\xi \in [0, 1]` with

.. math::

   \int_0^1 f(\xi)\,d\xi = 1, \qquad E(x) = E_\mathrm{ref}\,f(x/d),
   \qquad E_\mathrm{ref} = \frac{U}{d},

so that the "reference" reduced field :math:`E_\mathrm{ref}/N` used on the
axes of every plot and output file is always the *mean* field
:math:`U/d`, whatever the geometry.  :math:`\xi = 0` is the high-field
electrode (sphere surface, inner conductor, or either plate for a uniform
field) and :math:`\xi = 1` the low-field side.

Geometries
----------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - ``--field``
     - Profile
   * - ``uniform``
     - :math:`f \equiv 1`.  Symmetric; polarity irrelevant.
   * - ``sphere-plane R``
     - Sphere of radius :math:`R` (mm) above a grounded plane.  The on-axis
       field is computed *exactly* from the bispherical image-charge series
       (:func:`~incept1d.fields._sphere_sphere_axial_field` with the plane as the
       mirror sphere).  Not symmetric: ``incept1d pdiv`` computes both polarities
       (``sphere=positive`` / ``sphere=negative``).
   * - ``sphere-sphere R``
     - Two equal spheres of radius :math:`R` (mm) at :math:`\pm U/2`.
       Exact bispherical series.  Symmetric.
   * - ``coaxial A B``
     - Coaxial cylinders with inner radius :math:`a` and outer radius
       :math:`b` (mm).  Exact radial field, see `Coaxial cylinders`_.  Not
       symmetric: ``incept1d pdiv`` computes both polarities
       (``inner=positive`` / ``inner=negative``).
   * - ``fieldline FILE [UNIT]``
     - Tabulated :math:`|E|` along an arbitrary (curved) field line from
       an external electrostatic solver.  See below.

Since :math:`f(\xi)` is invariant under a geometric rescaling of the
electrode arrangement, :meth:`~incept1d.fields.FieldDistribution.build` returns the same
profile for any :math:`d`: a :math:`pd` sweep at fixed pressure corresponds
to scaled copies of the same geometry (sphere radius *and* gap scale
together), while a sweep at fixed :math:`d` varies only the pressure.  For a
tabulated field line the arc length fixes :math:`d = L`, and for coaxial
cylinders the radii fix :math:`d = b - a`, so only fixed-:math:`d` sweeps are
meaningful there.

Coaxial cylinders
-----------------

Between coaxial cylinders at potential difference :math:`U` the field is
radial and exact,

.. math::

   E(r) = \frac{U}{r\,\ln(b/a)}, \qquad a \le r \le b .

With :math:`r = a + \xi\,(b - a)` and :math:`d = b - a`, the normalised
profile is

.. math::

   f(\xi) = \frac{\rho - 1}{\big(1 + \xi(\rho - 1)\big)\ln\rho},
   \qquad \rho = \frac{b}{a},

which integrates to one exactly and depends on the radius ratio alone.
The field at the inner conductor, where inception starts, is

.. math::

   E(a) = \frac{U}{a\,\ln(b/a)} = f(0)\,\frac{U}{d}.

The radial direction is a field line, so the 1-D model is integrated along a
radius; the growth of the cross-section with :math:`r` is neglected, as the
spreading of flux tubes is for every other geometry.

The radii fix the gap at :math:`d = b - a`, and a :math:`pd` sweep is a
**pressure sweep** at that gap, exactly as for a tabulated field line (see
`Sweeps`_): ``incept1d pdiv`` without ``--p``/``--d`` uses :math:`d = b - a`
and reports against :math:`p`, ``--p`` is refused, and a ``--d`` other than
:math:`b - a` solves the arrangement scaled by :math:`d/(b-a)` and says so.
:ref:`Chap:Examples:Coaxial` works through a wire in a cylinder.

The ``FieldDistribution`` class
-------------------------------

:class:`~incept1d.fields.FieldDistribution` is a small dataclass holding the
geometry type and, where applicable, the sphere radius or the tabulated
profile.  :meth:`~incept1d.fields.FieldDistribution.build` turns that static
description into the callable :math:`f(\xi)` for a specific gap length.


``add_field_argument`` / ``parse_field_spec`` implement the shared
``--field SPEC`` syntax used by every subcommand, so geometry parsing is
written and tested once.

.. _Chap:FieldLines:

Tabulated field lines
---------------------

For electrode geometries beyond spheres and planes, the field along a field
line can be computed with an external electrostatic solver, exported as a
table, and used directly with ``--field fieldline FILE [UNIT]``.  The 1-D
model is then integrated along that line exactly as it is along the axis of
a sphere gap.

File format
...........

A plain numeric table, whitespace- or comma-separated; header lines and
lines starting with ``#`` are skipped.  The layout is inferred from the
number of numeric columns:

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Columns
     - Interpretation
   * - 2
     - ``s  |E|`` — position along the line and field magnitude
   * - 4
     - ``x  y  z  |E|`` — position and field magnitude
   * - 6
     - ``x  y  z  Ex  Ey  Ez`` — position and field vector

The arc length is accumulated from the positions in every layout: the
point-to-point Euclidean distance for the 4- and 6-column forms, and
:math:`|\Delta s|` for the 2-column one.  The first column of a 2-column
file therefore does not have to be an arc length already — a signed axis
coordinate works, and so does one that runs downwards, which is what an
export gives when the line was traced from the other electrode.  A column
that *is* an increasing arc length is unchanged by this.

Rows are used in file order: the **first row defines** :math:`\xi = 0` and
the last :math:`\xi = 1`, whichever way the coordinate runs.  Only
:math:`|E|` enters the model, so the direction of the field vector and the
absolute field units are irrelevant — the profile is normalised.  ``UNIT``
(``m``, ``cm``, ``mm``, ``um``; default ``m``) is the unit of the length
columns, and getting it wrong is what the check under
`Declaring the excitation`_ is for.

Example (``line.csv``, lengths in mm):

.. code-block:: text

   # s_mm   E_kV_per_mm
   0.0      12.4
   0.5      9.8
   1.0      7.9
   ...
   20.0     3.1

Declaring the excitation
........................

The length unit is declared on the command line, but the **field** unit is
not, and the file does not reveal it.  That is enough for the solve, which
sees only the normalised shape, but not for reporting how far the supplied
excitation is from inception.  ``--fieldline-voltage U_KV`` declares that
excitation, and the commands then report :math:`U^*/U_\mathrm{applied}`,
the factor by which the excitation must be scaled to reach inception.  The
ratio is exact whatever the field column is in, because :math:`U^*` comes
from the shape alone.  Without the flag no ratio is reported.

The line integral :math:`\int|E|\,ds` is printed in the file's own units as
a check: if the file is in V/m it must equal the declared excitation, and a
disagreement means a wrong length unit, a wrong column, or a line that does
not span the whole gap.  :ref:`Chap:Examples:FieldLine` works this through.

How the line is used
....................

:meth:`~incept1d.fields.FieldDistribution.from_fieldline` parametrises
the profile by normalised arc length :math:`\xi = s/L` and normalises it
so that :math:`\int_0^1 f\,d\xi = 1`.  Consequently

* The gap length is the arc length, :math:`d = L`;
* The reference field is the mean field along the line, :math:`E_\mathrm{ref}
  = U/L` with :math:`U = \int|E|\,ds` the voltage drop along the line;
* Every voltage reported by the commands is this line integral.

Polarity
~~~~~~~~

The file does not say which end is the anode.  ``incept1d pdiv`` evaluates
both polarities and labels them ``start=positive`` (the first tabulated
point is the anode) and ``start=negative`` (the first point is the
cathode).  Use the one that matches your electrode arrangement, or export
the line in the direction that makes ``start`` the electrode of interest.

Sweeps
~~~~~~

A tabulated line describes one geometry at one size.  Its arc length fixes
the gap, :math:`d = L`, so a :math:`pd` sweep along a field line is a
**pressure sweep**: the geometry is held and :math:`pd` varies through
:math:`p` alone.  That is what ``incept1d pdiv`` does when neither ``--p``
nor ``--d`` is given, and it is the only sweep that keeps the imported
geometry intact.

``--p`` is therefore **refused** for a field line: fixing the pressure
makes the gap the swept variable, so every point of the sweep would be a
differently sized copy of the imported arrangement.  For a single point at
one pressure, ask for it directly — :math:`pd = p\,L` with ``--pd-num 1``,
which the error message spells out.

Since :math:`d` is pinned, :math:`pd` is the pressure times a constant and
says nothing the pressure does not.  The results are therefore reported
against :math:`p` in bar rather than :math:`pd` in bar·mm: that is the
x-axis of the figures and the first column of ``--write-to-file``, and the
gap appears once in each curve's label instead of as a column that never
changes.  The sweep is still *requested* in :math:`pd` through
``--pd-min`` / ``--pd-max``.

The gap length can still be set explicitly with ``--d``.  That is allowed,
because :math:`f(\xi)` is invariant under a geometric rescaling, so
:math:`d \ne L` solves the *same* arrangement scaled by :math:`d/L` — every
electrode dimension, not only the gap.  It rescales the geometry once and
on purpose rather than across a sweep, so it is reported on stderr rather
than refused.

Resolution
~~~~~~~~~~

The profile is linearly interpolated between tabulated points.  Export
enough points to resolve the field near the high-field electrode; the
adaptive integration grid (``--dx``) refines the propagator, not the
profile.  Check the imported profile visually first:

.. code-block:: bash

   incept1d field --field fieldline line.csv mm

Limitations
...........

* Flux-tube divergence — the spreading of neighbouring field lines — is
  neglected, as for the sphere geometries.  The model is a drift-reaction
  balance *along* the line.
* The two-stream photon model tracks photons along the line only
  (:ref:`Chap:Photoionization`); in a strongly curved field it cannot
  account for photons that cross between field lines.
* Only ``|E|`` matters, so a line that reverses direction is treated as if
  the field had not reversed.  Export a line that runs from one electrode
  to the other.

Standalone use
--------------

The module doubles as a diagnostic: it plots the normalised profile and the
initial integration grid for a geometry, which is a quick way to check a
field-line file before using it in a solve.

.. code-block:: bash

   incept1d field --field sphere-plane 50 --d 20
   incept1d field --field coaxial 1 10
   incept1d field --field fieldline line.csv mm

API reference
-------------

.. automodule:: incept1d.fields
   :members:
   :undoc-members:
   :private-members: _sphere_sphere_axial_field, _coaxial_field, _compute_n_steps
