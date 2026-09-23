.. _Chap:ExampleConfigurations:

Example configurations
======================

The preceding pages specify the interfaces.  These two work them through on
real mechanisms, from opposite ends of the range.

:ref:`Chap:PaschenMechanism` is the smallest mechanism that does anything:
two species, one reaction, no data files at all, and a closed-form answer to
check against.  Read it to see the whole interface at once, and copy it when
starting a new gas.

:ref:`Chap:AirScheme` is a complete one: six species, a declarative reaction
table, tabulated swarm data, multigroup photoionization, field-dependent
cathode yields, and a family of ready-made configurations built on top.  Read
it to see what each part of the interface looks like when the chemistry is
not trivial.

.. toctree::
   :maxdepth: 2

   paschen
   air
