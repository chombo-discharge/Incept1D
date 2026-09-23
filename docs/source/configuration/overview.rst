.. _Chap:ConfigurationOverview:

Overview of the configuration files
===================================

.. contents:: On this page
   :local:
   :depth: 1

``Incept1D`` ships no chemistry of its own.  Everything the solver knows about
a gas — which species exist, how fast they are made and destroyed, how quickly
they drift, how many electrons a returning ion releases from the cathode —
comes from files outside the package, under ``mechanisms/``.  This chapter
describes those files and the interfaces they must satisfy.

The unit is a **mechanism directory**.  One directory holds one gas, and the
variants of that gas are configurations within it:

.. code-block:: text

   mechanisms/air/
   ├── air_pancheshnyi.py    the mechanism module      (required)
   ├── config.py             the configuration class   (optional)
   ├── baseline.json         configuration files       (optional)
   ├── nodetachment.json
   ├── paschen.json
   ├── lisbon.txt            data the module reads     (optional)
   └── o2m_mobility.txt

A run names the module, and optionally one or more configuration files:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py mechanisms/air/databases.json

Which files are needed
----------------------

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - File
     - Required
     - Role
   * - ``<gas>_<scheme>.py``
     - **yes**
     - The mechanism module.  Defines the species list and the functions the
       solver calls to build :math:`\bm{R}`, :math:`\bm{V}`, the photon
       coupling and the secondary-emission yields.  Must expose the interface
       in :ref:`Chap:NewMechanisms`.
   * - ``config.py``
     - No
     - A ``Config`` class that turns a JSON dictionary into overrides for the
       module.  Needed only if the mechanism is to accept configuration files
       at all; see :ref:`Chap:ConfigPy`.
   * - ``*.json``
     - No
     - Named parameter sets — cross-section database, rate multipliers,
       photoionization and SEE settings.  See :ref:`Chap:Configuration`.
   * - data files
     - No
     - Whatever the module chooses to read: BOLSIG+ swarm output, LXCat
       mobility tables, fitted absorption curves.  The package never opens
       them; the module does, relative to its own location.

A mechanism with no ``config.py`` and no JSON is perfectly valid — it is then
a fixed gas with no adjustable parameters.  Passing a configuration file to
such a mechanism is an error rather than a silent no-op, so a mistyped path
cannot quietly leave the defaults in place.

How the three fit together
--------------------------

:func:`incept1d.mechanism.load_mechanism` is the only place these files meet.
Given a module path and a raw dictionary, it takes five steps:

1. Look for ``config.py`` next to the module.  If present, build
   ``Config.from_dict(raw, mech_dir)``; if absent and a dictionary was
   supplied, raise.
2. Inject ``config.pre_exec_vars()`` into the module namespace **before** the
   module body runs, so module-level constants can be overridden.
3. Execute the module and check it against
   :data:`incept1d.mechanism.REQUIRED_ATTRS`, naming anything missing.
4. Call ``config.post_exec_init(mod)`` for work that needs the executed
   module — fitting a photoionization model, for instance.
5. Wrap the result in a :class:`~incept1d.mechanism.Mechanism` with
   ``config.mechanism_params()``, which bakes rate multipliers, scale factors
   and per-polarity overrides into the accessors.

.. code-block:: text

   raw JSON dict ──► Config.from_dict ──► Config
                                           │
                    pre_exec_vars() ───────┤ injected before exec
                                           │
   <gas>_<scheme>.py ──► exec ──► module ──┤ post_exec_init(mod)
                                           │
                    mechanism_params() ────┴──► Mechanism ──► the solver

Everything downstream sees only the :class:`~incept1d.mechanism.Mechanism`.
By the time a solver calls ``mod.get_R(EN, p, T)`` the multipliers have
already been applied, so no other part of the code has to know that
configurations exist.

.. important::

   Mechanism modules are **executed, not imported as a package**.  They are
   ordinary Python files run by ``importlib`` with override variables injected
   into their namespace.  That is what makes ``pre_exec_vars`` possible, and
   it is why a mechanism resolves its data files against ``__file__`` rather
   than relying on the import system.

Where to go next
----------------

* :ref:`Chap:NewMechanisms` — The module interface, function by function,
  and how to write one for a new gas.
* :ref:`Chap:ConfigPy` — The ``Config`` protocol.
* :ref:`Chap:Configuration` — The JSON format and the keys the air family
  understands.
* :ref:`Chap:ModifyingReactions` — The declarative reaction strings and how
  multipliers select them.

Two complete mechanisms are documented as worked examples:

* :ref:`Chap:PaschenMechanism` — The smallest mechanism that does
  anything: two species, one reaction, no data files, and an answer that can
  be checked in closed form.
* :ref:`Chap:AirScheme` — The reference dry-air scheme: six species, a
  declarative reaction table, multigroup photoionization and a family of
  ready-made configurations.
