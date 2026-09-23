.. _Chap:Configuration:

JSON configuration files
========================

.. contents:: On this page
   :local:
   :depth: 1

Every solver command accepts JSON configuration files after the mechanism
path.  A configuration changes the *parameters* of a mechanism — which data
it reads, how strongly a reaction contributes, what the electrode surface
does — without editing the Python file, and several configurations can be
run and plotted side by side in one invocation.  This is how a sensitivity
study is set up: one file, one curve per entry.

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
           "A- + M -> e + A + M": 0.0
         }
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

What the keys mean
------------------

A configuration file is not validated against a fixed schema.  It is handed
to the mechanism directory's ``config.py``, which decides what each key
means (:ref:`Chap:ConfigPy`), and unknown keys are ignored so that one file
can be shared between mechanisms understanding different subsets of it.

One key is universal:

``label``
   The string shown in legends, terminal tables and output-file headers.
   There is one curve per configuration, so it should be short and
   distinguishing.

The rest are the mechanism's own.  A set of conventional names has
nevertheless established itself, because these map directly onto the
parameters :class:`~incept1d.mechanism.Mechanism` applies:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Key
     - Conventional meaning
   * - ``reaction_multipliers``
     - Dictionary ``{reaction string: factor}``.  ``0.0`` disables a
       reaction, ``2.0`` doubles it.  See :ref:`Chap:ModifyingReactions`.
   * - ``xi_photo``, ``xi_emit``
     - Scale factors on the photoionization coupling and the photoemission
       yields (:ref:`Chap:Photoionization`,
       :ref:`Chap:SecondaryEmission`).
   * - ``gamma0``, ``gamma1``, ``eref``, ``beta``
     - Parameters of the ion-induced secondary-emission yield.
   * - ``positive``, ``negative``
     - Nested objects overriding any of the above for one polarity only.

A new mechanism should use these names for these meanings rather than
inventing its own, so that a reader who knows one mechanism can read the
configurations of another.  Anything genuinely specific to a mechanism —
which data file to read, which gas to select — is named by that mechanism
and documented with it.

Per-polarity overrides
----------------------

In a gap that is not symmetric, the cathode is a different electrode in the
two polarities, so the surface parameters need not be the same.  A nested
``positive`` or ``negative`` object overrides the top-level values for that
polarity alone:

.. code-block:: json

   {
     "label": "Steel sphere, brass plane",
     "positive": { "gamma0": 5e-4 },
     "negative": { "gamma0": 2e-3 }
   }

This works because :meth:`~incept1d.mechanism.Mechanism.resolve` is called
once per polarity, after the module has been executed.  A parameter read
while the module executes — a data file name, say — is fixed for the whole
run and therefore cannot be overridden this way;
:ref:`Chap:ConfigPy` explains which parameters fall on which side of that
line.

How configurations are applied
------------------------------

A configuration file is inert on its own: the mechanism directory's
``config.py`` is what turns it into overrides.  That interface is described
in :ref:`Chap:ConfigPy`, and the two worked examples,
:ref:`Chap:PaschenMechanism` and :ref:`Chap:AirScheme`, show the keys each
of those mechanisms actually defines.
