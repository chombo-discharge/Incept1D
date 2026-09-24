.. _Chap:Examples:Coaxial:

Coaxial cylinders
=================

.. contents:: On this page
   :local:
   :depth: 1

A thin wire inside a grounded cylinder is the textbook corona geometry.  The
field is exact and strongly non-uniform near the wire, so ionization is
confined to a thin layer around it while most of the gap is attaching.  This
example computes the inception voltage of such an arrangement in dry air
with ``--field coaxial``, for three wire radii across two decades of
pressure, and shows that the result is governed by the field at the wire
and the product :math:`pa`.

The geometry
------------

``--field coaxial A B`` is a pair of coaxial cylinders with inner radius
:math:`a` and outer radius :math:`b`, in mm.  The field is radial,

.. math::

   E(r) = \frac{U}{r\,\ln(b/a)},

so the highest field is at the surface of the inner conductor,
:math:`E(a) = U / (a\ln(b/a))`, and the 1-D model is integrated along a
radius from :math:`r = a` (:math:`\xi = 0`) to :math:`r = b`
(:math:`\xi = 1`).  See :ref:`Chap:FieldDistributions` for the
normalised profile.

The radii fix the gap at :math:`d = b - a`.  As for a tabulated field line,
the :math:`pd` sweep is therefore a sweep in pressure at that fixed gap, and
the results are reported against :math:`p`.  The arrangement is not
symmetric, so both polarities are computed: ``inner=positive`` has the inner
conductor as the anode, ``inner=negative`` as the cathode.

Running the calculation
-----------------------

For a 1 mm wire inside a 10 mm cylinder, from 0.1 to 10 bar:

.. code-block:: bash

   incept1d pdiv mechanisms/air/pancheshnyi/air_pancheshnyi.py \
       --field coaxial 1 10 \
       --pd-min 0.9 --pd-max 90 --pd-num 30 \
       --dx 10 400 0.01 \
       --write-to-file examples/coaxial/sim_a1mm.dat

The sweep is *requested* in :math:`pd`, so the range is the pressure range
multiplied by :math:`b - a = 9` mm.  The figure below repeats this for three
inner radii, 0.25, 1 and 4 mm, in the same cylinder;
``examples/coaxial/README.md`` lists the :math:`pd` range of each.
``incept1d pdiv`` confirms the gap on startup:

.. code-block:: text

   --field coaxial: no --p/--d given, using d = L = 9 mm (outer minus inner radius).

The ``--dx`` setting is finer than the default.  Near a thin wire the field
falls off over a distance comparable to :math:`a`, and at high pressure the
default grid resolves that too coarsely: for :math:`a = 0.25` mm at 10 bar it
is 6 % off, while ``--dx 10 400 0.01`` is within 0.2 % of a much finer grid
for every radius here.  See :ref:`Chap:Numerics:Propagator`.

Each output file carries the radii and the polarity convention in its
header:

.. code-block:: text

   # Radii:       a = 1 mm, b = 10 mm
   # Polarity:    inner=positive → inner conductor is anode (+),  inner=negative → inner conductor is cathode (−)

Result
------

.. figure:: ../figures/coaxial.*
   :width: 90%
   :align: center

   Reduced field at the surface of the inner conductor at inception,
   :math:`E(a)/N`, against :math:`pa`, for three inner radii in a cylinder
   of radius :math:`b = 10` mm.  Inner conductor positive (solid) and
   negative (dashed).

The output files give the *mean* reduced field :math:`U^*/(Nd)`; the figure
multiplies it by :math:`f(0) = (b-a)/(a\ln(b/a))` to get the field at the
wire, and plots it against :math:`pa`.  Three things stand out.

* **The radii collapse onto one curve.**  Each inner radius covers its own
  range of :math:`pa`, but where they overlap they agree to within 5 %,
  although their radius ratios :math:`b/a` differ by a factor of sixteen.
  Inception is decided in the ionizing layer around the wire, whose extent
  in units of the mean free path is set by :math:`pa`; the rest of the gap
  hardly matters.
* **The wire field falls towards the α = η field.**  A thin wire at low
  pressure needs about 770 Td at its surface, because the ionizing layer is
  only a few ionization lengths thick.  As :math:`pa` grows the layer
  thickens and the required field drops, to about 135 Td at
  :math:`pa = 40` bar·mm — approaching, and always staying above, the field
  at which ionization and attachment balance in this mechanism
  (114 Td at 1 bar, 121 Td at 10 bar).
* **The polarities are close.**  A positive inner conductor needs a
  slightly higher field than a negative one, by 0.4 to 4 %.

Variations to try
-----------------

* **Change the outer radius.**  If the wire region decides inception,
  :math:`b` should barely matter.  With :math:`b = 25` mm instead of 10 mm
  the curves for :math:`a = 1` and 4 mm move by less than 2 %.  Keep an
  eye on the thinnest wire at high pressure, though.  A large radius ratio
  puts most of a long gap deep in attachment, which is where the
  determinant formulation runs out of range
  (:ref:`Chap:Numerics:Determinant`): with :math:`b = 25` mm and
  :math:`a = 0.25` mm the positive-polarity curve develops a kink above
  about 1 bar and then loses its root.
* **Look at the profile first.**  ``incept1d field --field coaxial 1 10``
  plots :math:`f(\xi)` and the integration grid, and prints :math:`f(0)`,
  the ratio of the surface field to the mean field :math:`U/d`.
* **Compare with the ionization integral.**  ``incept1d ionization`` with
  the same ``--field`` evaluates :math:`\int\max(\alpha-\eta,0)\,dx`
  along the radius, the quantity a streamer or avalanche-size criterion
  would use.
