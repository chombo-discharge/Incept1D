.. _Chap:ExampleConfigurations:

Example configurations
======================

The preceding pages specify the interfaces.  These pages work them through on
real mechanisms, from one end of the range to the other.

:ref:`Chap:PaschenMechanism` is the smallest mechanism that does anything:
two species, one reaction, no data files at all, and a closed-form answer to
check against.  Read it to see the whole interface at once, and copy it when
starting a new gas.

:ref:`Chap:AirScheme` is a complete one: six species, a declarative reaction
table, tabulated swarm data, multigroup photoionization, field-dependent
cathode yields, and a family of ready-made configurations built on top.  Read
it to see what each part of the interface looks like when the chemistry is
not trivial.

:ref:`Chap:RFDScheme` is not a gas discharge at all: runaway electrons,
positrons and photons in a thundercloud, whose feedback closes in the gas
with nothing emitted at either end.  Read it to see one-way photon groups
(``PHOTON_DIRECTIONS``) and how far the interface stretches.

.. toctree::
   :maxdepth: 2

   paschen
   air
   rfd
