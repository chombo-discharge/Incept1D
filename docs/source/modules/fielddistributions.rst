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
electrode (sphere surface, or either plate for a uniform field) and
:math:`\xi = 1` the low-field side.

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
       (:func:`_sphere_sphere_axial_field` with the plane as the mirror
       sphere).  Not symmetric: ``incept1d pdiv`` computes both polarities
       (``sphere=positive`` / ``sphere=negative``).
   * - ``sphere-sphere R``
     - Two equal spheres of radius :math:`R` (mm) at :math:`\pm U/2`.
       Exact bispherical series.  Symmetric.
   * - ``fieldline FILE [UNIT]``
     - Tabulated :math:`|E|` along an arbitrary (curved) field line from
       an external electrostatic solver.  See :ref:`Chap:FieldLines`.

Since :math:`f(\xi)` is invariant under a geometric rescaling of the
electrode arrangement, :meth:`FieldDistribution.build` returns the same
profile for any :math:`d`: a :math:`pd` sweep at fixed pressure corresponds
to scaled copies of the same geometry (sphere radius *and* gap scale
together), while a sweep at fixed :math:`d` varies only the pressure.  For a
tabulated field line the arc length fixes :math:`d = L`, so only fixed-
:math:`d` sweeps are meaningful.

The ``FieldDistribution`` class
-------------------------------

:class:`~incept1d.fields.FieldDistribution` is a small dataclass holding the
geometry type
and, where applicable, the sphere radius or the tabulated profile;
:meth:`FieldDistribution.build` turns that static description into the
callable :math:`f(\xi)` for a specific gap length.

.. literalinclude:: ../../../src/incept1d/fields.py
   :language: python
   :pyobject: FieldDistribution.build

``add_field_argument`` / ``parse_field_spec`` implement the shared
``--field SPEC`` syntax used by every subcommand, so geometry parsing is
written and tested once.

Standalone use
--------------

The module doubles as a diagnostic: it plots the normalised profile and the
initial integration grid for a geometry, which is a quick way to check a
field-line file before using it in a solve.

.. code-block:: bash

   incept1d field --field sphere-plane 50 --d 20
   incept1d field --field fieldline line.csv mm

API reference
-------------

.. automodule:: incept1d.fields
   :members:
   :undoc-members:
   :private-members: _sphere_sphere_axial_field, _compute_n_steps
