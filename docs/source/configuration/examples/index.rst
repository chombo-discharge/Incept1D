.. _Chap:ExampleConfigurations:

Example configurations
======================

The preceding pages specify the formats.  These pages work them through on
real mechanisms, from one end of the range to the other.

:ref:`Chap:PaschenMechanism` is the smallest mechanism file that does
anything: two species, one reaction, no data files at all, and a closed-form
answer to check against.  Read it to see a whole mechanism file at once, and
copy it when starting a new gas.

:ref:`Chap:AirScheme` is a complete one: six species, a reaction table,
tabulated swarm data, multigroup photoionization, field-dependent cathode
yields, and a family of ready-made configurations built on top.  Read it to
see what each block of a mechanism file looks like when the chemistry is not
trivial.

:ref:`Chap:RFDScheme` is not a gas discharge at all: runaway electrons,
positrons and photons in a thundercloud, whose feedback closes in the gas
with nothing emitted at either end.  It is a Python mechanism, because its
coefficients switch form at thresholds and its photon coupling is not a
source times an absorption; read it to see one-way photon groups
(``PHOTON_DIRECTIONS``) and how far the Python interface stretches.

.. toctree::
   :maxdepth: 2

   paschen
   air
   rfd
