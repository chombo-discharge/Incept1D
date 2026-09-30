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
electrode (sphere surface, hyperboloid tip, inner conductor, or either plate
for a uniform field) and :math:`\xi = 1` the low-field side, unless the
profile is reversed (see `Reversing the profile`_).  The electrode at
:math:`\xi = 0` can carry a protrusion (``--protrusion``, see
`Protrusions`_).

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
       mirror sphere).  Not symmetric: ``incept1d pdiv`` and ``incept1d growth`` compute both polarities
       (``sphere=positive`` / ``sphere=negative``).
   * - ``sphere-sphere R``
     - Two equal spheres of radius :math:`R` (mm) at :math:`\pm U/2`.
       Exact bispherical series.  Symmetric.
   * - ``hyperboloid-plane R``
     - Hyperboloid of revolution with tip radius of curvature :math:`r`
       (mm) above a grounded plane, the usual model of a point electrode.
       Exact on-axis field, see `Hyperboloid-plane`_.  Not symmetric:
       ``incept1d pdiv`` and ``incept1d growth`` compute both polarities
       (``tip=positive`` / ``tip=negative``).
   * - ``coaxial A B``
     - Coaxial cylinders with inner radius :math:`a` and outer radius
       :math:`b` (mm).  Exact radial field, see `Coaxial cylinders`_.  Not
       symmetric: ``incept1d pdiv`` and ``incept1d growth`` compute both
       polarities (``inner=positive`` / ``inner=negative``).
   * - ``fieldline FILE [LENGTH] [FIELD]``
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

Reversing the profile
---------------------

``--reverse-field`` reads the profile from the other electrode,
:math:`f(\xi) \to f(1 - \xi)`, for any geometry.  The physics is
unchanged; what changes is which electrode is at :math:`\xi = 0`, and so
which one the polarity labels name and a protrusion sits on:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - ``--field``
     - :math:`\xi = 0` as given
     - :math:`\xi = 0` reversed
   * - ``sphere-plane``
     - ``sphere``
     - ``plane``
   * - ``hyperboloid-plane``
     - ``tip``
     - ``plane``
   * - ``coaxial``
     - ``inner`` conductor
     - ``outer`` conductor
   * - ``fieldline``
     - ``start`` (first row)
     - ``end`` (last row)

So ``plane=negative`` on a reversed sphere-plane gap is the same case as
``sphere=positive`` without it, and gives the same inception field; the
option exists to put a protrusion on the plane or the outer conductor, and
to take a field line from whichever end it was exported.  A uniform or
sphere-sphere gap is symmetric, and is left unchanged with a note.  The
reversal is recorded in the labels and in the ``--write-to-file`` header.

Hyperboloid-plane
-----------------

A point electrode is usually modelled as a hyperboloid of revolution above a
plane, because both surfaces are coordinate surfaces of prolate spheroidal
coordinates and the potential is then known in closed form [Coelho1971]_.
With the tip at distance :math:`d` from the plane and tip radius of
curvature :math:`r`, the common foci lie on the axis at

.. math::

   z = \pm a, \qquad a = \sqrt{d\,(d + r)},

and the potential depends on the coordinate
:math:`\eta = (r_- - r_+)/2a` alone, :math:`r_\pm` being the distances to
the two foci: :math:`\phi = U \operatorname{artanh}\eta \,/\,
\operatorname{artanh}\eta_0`, with :math:`\eta = 0` the plane and
:math:`\eta_0 = d/a` the tip.  On the axis :math:`\eta = z/a`, so the field
at height :math:`z` above the plane is

.. math::

   E(z) = \frac{U}{\operatorname{artanh} k}\,\frac{a}{a^2 - z^2},
   \qquad k = \frac{d}{a} = \sqrt{\frac{d}{d + r}} .

With :math:`z = d(1 - \xi)` the normalised profile is

.. math::

   f(\xi) = \frac{k}{\big(1 - k^2 (1 - \xi)^2\big)\operatorname{artanh} k},

which integrates to one exactly and depends on :math:`d/r` alone.  The field
at the tip,

.. math::

   E(d) = \frac{2 U a}{r\,d\,\ln\dfrac{1 + k}{1 - k}}
   \;\longrightarrow\; \frac{2U}{r \ln(4d/r)} \quad (r \ll d),

is the classic point-plane estimate, and :math:`f \to 1` as
:math:`r/d \to \infty`.  As for the sphere-plane gap, :math:`r` is held
fixed while a :math:`pd` sweep at fixed pressure varies :math:`d`, so
:math:`d/r` changes along the sweep.

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

Protrusions
-----------

A surface defect — a burr, a particle, a scratch edge — enhances the field
locally and can move inception well below the value for a smooth
electrode.  ``--protrusion`` puts a conducting, axisymmetric protrusion
of height :math:`h` on the :math:`\xi = 0` electrode of any ``--field``
geometry, on the axis.  Three shapes are available (lengths in mm):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ``--protrusion``
     - Shape
   * - ``spheroid H B``
     - Half-spheroid of base radius :math:`b`: a needle (:math:`h > b`),
       a hemispherical boss (:math:`h = b`) or a flat bump
       (:math:`h < b`).  Closed form, see `Spheroid`_.
   * - ``cone H R ANGLE``
     - Cone of half-angle :math:`\theta` (degrees, 5 to 80) whose apex
       is rounded by a sphere of radius :math:`r` tangent to it: a burr or
       a sharp asperity.  See `Cone and rod`_.
   * - ``rod H R``
     - Cylinder of radius :math:`R` with a hemispherical cap
       (:math:`h \ge R`): a whisker, a fibre or a wire end.  See
       `Cone and rod`_.

Near its tip every smooth protrusion looks like a paraboloid of the tip
radius of curvature, and the field falls from :math:`\beta E_0` over a
distance of that order whatever the shape.  The shapes differ in how far
the enhancement reaches beyond that, out to distances of order :math:`h`
— and since the criterion integrates the ionization along the gap, that
tail can matter as much as :math:`\beta`.  A spheroid is the natural
choice when the tip radius and height are what is known; a cone when the
opening angle is; a rod for a slender body of constant thickness.

Spheroid
........

A half-spheroid on a grounded plane in a uniform field :math:`E_0` is,
with its mirror image, a whole spheroid in a uniform field, and the
potential is known in closed form in spheroidal coordinates with focal
distance :math:`c = \sqrt{|h^2 - b^2|}`.  On the axis, at height
:math:`z \ge h` above the plane,

.. math::
   :label: eq_spheroid_axial

   \frac{E(z)}{E_0} = 1 + \frac{N(z)}{D}, \qquad
   D = \sum_{k \ge 1} \frac{\tau^{k-1}}{2k + 1}, \qquad
   N = w^3 \sum_{k \ge 1} \frac{2k}{2k + 1}\,(\tau w^2)^{k-1},

with :math:`\tau = (h^2 - b^2)/h^2` and :math:`w = h/z`.  The one expression
holds for prolate (:math:`\tau > 0`), spherical (:math:`\tau = 0`) and
oblate (:math:`\tau < 0`) shapes; where the series converge slowly they are
summed in closed form through the Legendre function of the second kind,
:math:`Q_1(x) = x \operatorname{artanh}(1/x) - 1` (prolate) or
:math:`1 - x \operatorname{arccot} x` (oblate).  At the tip this is the
field-enhancement factor,

.. math::

   \beta = \frac{1}{(\xi_0^2 - 1)\,\big(\xi_0 \operatorname{artanh}(1/\xi_0) - 1\big)},
   \qquad \xi_0 = \frac{h}{c},

which is :math:`3` for a hemisphere — where :eq:`eq_spheroid_axial` is
:math:`1 + 2h^3/z^3` —, grows as :math:`(h/b)^2 / (\ln(2h/b) - 1)` for a
needle, and tends to :math:`1` for a flat disc.  The enhancement is
confined to a few tip radii of curvature :math:`b^2/h` and a few heights
:math:`h` from the tip.  ``incept1d field`` prints :math:`\beta` and
the tip radius :math:`b^2/h` for every shape.

Cone and rod
............

Neither shape has a closed form, so its enhancement is computed once, when
the options are read, with the charge simulation method
[Singer1974]_.  Ring charges :math:`q_j` of radius :math:`a_j` at height
:math:`z_j` are placed just inside the surface, each with its negative
image in the plane; a ring of unit charge has the potential

.. math::

   \Phi(\rho, z; a, z_j) = \frac{2}{\pi}\,
   \frac{K(m)}{\sqrt{(\rho + a)^2 + (z - z_j)^2}}, \qquad
   m = \frac{4 a \rho}{(\rho + a)^2 + (z - z_j)^2},

in units where :math:`4\pi\varepsilon_0 = 1`, with :math:`K` the complete
elliptic integral of the first kind (a point charge at the apex, where
:math:`a = 0`).  The strengths make the surface an equipotential at one
collocation point per charge,

.. math::

   \sum_j q_j \big[\Phi(\rho_i, z_i; a_j, z_j) - \Phi(\rho_i, z_i; a_j, -z_j)\big]
   = E_0 z_i,

which cancels the applied potential :math:`-E_0 z` there, and the axial
field follows from the charges in closed form.  The collocation points are
spaced by a fraction of the local length scale — the tip radius at the
apex, the body radius along the flank, the height above the plane near
the base, where the protrusion meets its image — so that the number of
charges grows only logarithmically with :math:`h` for a cone and linearly
with :math:`h/R` for a rod.

The solution is checked halfway between collocation points, where it is
least accurate: a surface potential above :math:`2\times10^{-3}` of the
applied one is an error, not a result.  Over the accepted range —
half-angles from 5° to 80° and :math:`h/r \le 10^5` for a cone,
:math:`h/R \le 1000` for a rod — halving the spacing changes the axial
field by less than :math:`10^{-3}`, and the method reproduces what is known
in closed form:

* A rod of height :math:`R` is a hemisphere, and matches
  :math:`1 + 2R^3/z^3`;
* Between the tip radius and the height, the field of a cone falls as
  :math:`s^{\nu - 1}`, with :math:`\nu` the root of
  :math:`P_\nu(-\cos\theta) = 0` that describes the field near the
  apex of a sharp cone;
* The rod's :math:`\beta` agrees with the fit
  :math:`1.2\,(h/R + 2.15)^{0.9}` of [EdgcombeValdre2001]_ within its
  few-percent accuracy.

The gap and the voltage
.......................

With a protrusion, :math:`d` is the gap **from its tip** to the other
electrode, so the electrode surface is at :math:`D = d + h`.  The profile
is the background profile of the electrodes alone, built for :math:`D` and
read from the tip onwards, times the enhancement
:eq:`eq_spheroid_axial`, and normalised again over the tip-to-electrode
path,

.. math::

   f(\xi) = \frac{f_\mathrm{bg}\big((h + \xi d)/D\big)\,
   E(h + \xi d)/E_0}{\int_0^1 f_\mathrm{bg}\big((h + \xi' d)/D\big)\,
   E(h + \xi' d)/E_0 \, d\xi'} .

The renormalisation is what keeps :math:`\int E\,dx = U`: the protrusion
concentrates the same voltage into a stronger field at its tip and a
slightly weaker one elsewhere, and the reference field :math:`U/d` of every
output stays the mean field along the tip-to-electrode path.  Where the
geometry fixes the gap (coaxial cylinders, a tabulated field line), the gap
is that length less :math:`h`.

Validity
........

The profile is exact for a protrusion on a plate in a uniform field (to
the accuracy of the charge simulation for a cone or a rod), but for one
thing: the closed form assumes nothing above the protrusion, so
the other electrode is an equipotential only up to the protrusion's
dipole field, a relative error of order :math:`(h/d)^3`.  On a curved
electrode the enhancement is superposed on the background, which needs
:math:`h` small against the length over which the background varies: the
sphere radius :math:`R`, the tip radius :math:`r`, the inner radius
:math:`a` (outer radius :math:`b` when reversed), or, for a tabulated
line, :math:`f/|df/ds|` near :math:`\xi = 0`.  On a reversed plane the
field varies over the gap, which the gap check covers.
The commands print a ``note:`` on stderr when :math:`h` exceeds a tenth of
any of these (the gap check uses the smallest gap of the sweep), and
compute anyway.

A protrusion is also where the 1-D model is most stretched.  The ionization
region shrinks towards the size of the tip, over which the real field lines
fan out and photons escape sideways; both are neglected (see
`Limitations`_), and increasingly so as :math:`h/b` grows.

For a tabulated field line, the protrusion sits on the first point.  If the
line is stronger at its far end, the command says so and suggests
``--reverse-field`` (see `Reversing the profile`_).  A line exported from a model
that already contains the protrusion must not be given ``--protrusion``
too: that counts the enhancement twice.

The ``FieldDistribution`` class
-------------------------------

:class:`~incept1d.fields.FieldDistribution` is a small dataclass holding the
geometry type and, where applicable, the sphere radius or the tabulated
profile, and the protrusion if there is one.
:meth:`~incept1d.fields.FieldDistribution.build` turns that static
description into the callable :math:`f(\xi)` for a specific gap length.


``add_field_argument`` / ``parse_field_spec`` implement the shared
``--field SPEC`` syntax used by every subcommand, so geometry parsing is
written and tested once.

.. _Chap:FieldLines:

Tabulated field lines
---------------------

For electrode geometries beyond spheres and planes, the field along a field
line can be computed with an external electrostatic solver, exported as a
table, and used directly with ``--field fieldline FILE [LENGTH] [FIELD]``.  The 1-D
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
absolute field units are irrelevant to the solve — the profile is
normalised.  ``LENGTH`` (``m``, ``cm``, ``mm``, ``um``; default ``m``) is
the unit of the length columns, and getting it wrong is what the check
under `Declaring the excitation`_ is for.  ``FIELD`` is the optional unit
of the field column, ``[MULT*]VOLT/LENGTH`` with ``VOLT`` one of ``mV``,
``V``, ``kV``, ``MV`` — for instance ``kV/mm``, ``V/m``, ``kV/cm``, or
``1e3*V/m`` for a column in thousands of V/m.  The two tokens may come in
either order, since only a field unit contains ``/``.

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

The solve sees only the normalised shape, so it needs no field unit.
Reporting how far the supplied excitation is from inception does: the
commands then report :math:`U^*/U_\mathrm{applied}`, the factor by which the
excitation must be scaled to reach inception, and that needs
:math:`U_\mathrm{applied}`.  There are two ways to give it:

* Declare the field unit, ``--field fieldline FILE mm kV/mm``.  The line
  integral :math:`\int|E|\,ds` is then a voltage, and it is the excitation;
* Declare the excitation, ``--fieldline-voltage U_KV``.  With a field unit
  as well, this takes precedence, and a line integral more than 2 % away
  from it is noted on stderr.

The ratio is exact whichever way the voltage is given, because :math:`U^*`
comes from the shape alone.  With neither, no ratio is reported, rather
than a wrong one.

Without a field unit the line integral is printed in the file's own units
as a check: if the file is in V/m it must equal the declared excitation,
and a disagreement means a wrong length unit, a wrong column, or a line
that does not span the whole gap.  :ref:`Chap:Examples:FieldLine` works
this through.

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

The file does not say which end is the anode.  ``incept1d pdiv`` and
``incept1d growth`` evaluate both polarities and labels them ``start=positive`` (the first tabulated
point is the anode) and ``start=negative`` (the first point is the
cathode).  Use the one that matches your electrode arrangement.

``--reverse-field`` reads the rows last to first, which makes the last
tabulated point :math:`\xi = 0`, labelled ``end``, without editing the
file (see `Reversing the profile`_).

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

The profile is linearly interpolated between tabulated points.  The
integration grid follows the table: it is refined wherever the tabulated
field changes by more than 20 % across a segment, including at a narrow
peak between the electrodes, so a feature the table resolves is resolved
by the solve too, and a protrusion added with ``--protrusion`` is refined
in the same way.  What the grid cannot do is add detail the table lacks:
export enough points to resolve the field near the high-field electrode
and at any local peak.  Check the imported profile visually first:

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
   incept1d field --field hyperboloid-plane 0.1 --d 20
   incept1d field --field coaxial 1 10
   incept1d field --field fieldline line.csv mm
   incept1d field --field uniform --d 10 --protrusion spheroid 0.5 0.1
   incept1d field --field uniform --d 10 --protrusion cone 0.5 0.01 20
   incept1d field --field sphere-plane 50 --d 20 --reverse-field --protrusion rod 0.5 0.02

API reference
-------------

.. automodule:: incept1d.fields
   :members:
   :undoc-members:
   :private-members: _sphere_sphere_axial_field, _hyperboloid_plane_field,
      _coaxial_field

.. automodule:: incept1d.protrusions
   :members:
   :private-members: _spheroid_axial_field, _csm_solve, _ring_potential
