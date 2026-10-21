.. _Chap:Configuration:

Configuration files
===================

.. contents:: On this page
   :local:
   :depth: 1

Every solver command accepts configuration files after the mechanism.  A
configuration changes the *parameters* of a mechanism — which data it reads,
how strongly a reaction contributes, what the electrode surface does —
without editing the mechanism, and several configurations can be run and
plotted side by side in one invocation.  This is how a sensitivity study is
set up: one file, one curve per entry.

File format
-----------

A file contains either a single configuration object, a bare list of
objects, or, the usual form, a list under a ``"configurations"`` key:

.. code-block:: json

   {
     "configurations": [
       {
         "label": "Baseline"
       },
       {
         "label": "No detachment",
         "reaction_multipliers": {
           "B- + M -> e + B + M": 0.0
         }
       },
       {
         "label": "High yield",
         "gamma0": 1e-2
       }
     ]
   }

Configurations are independent — a later entry does **not** inherit from an
earlier one — so every key a configuration needs must appear in it.  Keys
beginning with ``_``, such as ``"_comment"``, are ignored, which is the
only way to write a comment in JSON.

Several files may be given at once, and every configuration in every file
is run.  That makes a file a reusable set rather than a whole study: one
file holding the cross-section variants and another holding the surface
variants can be combined on the command line without editing either.

What a configuration can set
----------------------------

For a mechanism file (:ref:`Chap:MechanismFiles`), a configuration may
contain exactly these keys:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Key
     - Meaning
   * - ``label``
     - The string shown in legends, terminal tables and output-file headers.
       There is one curve per configuration, so it should be short and
       distinguishing.
   * - any parameter
     - Replaces the value the mechanism file gives that parameter: a number
       for a number, a file name for a file.
   * - ``reaction_multipliers``
     - ``{reaction string: factor}``.  ``0.0`` disables a reaction, ``2.0``
       doubles it.  See :ref:`Chap:ModifyingReactions`.
   * - ``positive``, ``negative``
     - Parameters for one polarity only; see below.

Anything else is an error, and so is a multiplier for a reaction the
mechanism does not have.  A configuration therefore belongs to one
mechanism: one written for another, or one with a misspelt parameter, fails
when it is loaded instead of solving a different problem than it describes.
What can be configured is exactly the mechanism's ``parameters`` block, so
the mechanism file is also the list of every key its configurations may use,
with its default.

A relative file name in a configuration is looked for beside the
configuration file first, and then beside the mechanism file, so a
configuration can bring its own data:

.. code-block:: json

   {
     "label": "My cross sections",
     "cross_sections": "my_swarm_data.txt"
   }

Per-polarity overrides
----------------------

In a gap that is not symmetric, the cathode is a different electrode in the
two polarities, so the surface parameters need not be the same.  A nested
``positive`` or ``negative`` object overrides parameters for that polarity
alone:

.. code-block:: json

   {
     "label": "Steel sphere, brass plane",
     "positive": { "gamma0": 5e-4 },
     "negative": { "gamma0": 2e-3 }
   }

Any parameter may be overridden this way; nothing else may.  The polarity
is the sign of the electrode at :math:`\xi = 0` of the field profile
(:ref:`Chap:FieldDistributions`); in a symmetric gap with no overrides the
two polarities are the same problem and are solved once.

Python mechanisms
-----------------

A configuration for a Python mechanism (:ref:`Chap:NewMechanisms`) has the
same file format, but what its keys mean is decided by the mechanism
directory's ``config.py`` (:ref:`Chap:ConfigPy`), which should refuse keys
it does not know.  The conventional names — ``reaction_multipliers``,
``gamma0``, ``positive``, ``negative`` — mean the same there.
