.. _Chap:ConfigurationOverview:

Overview of the configuration files
===================================

.. contents:: On this page
   :local:
   :depth: 1

``Incept1D`` builds no chemistry of its own.  Everything the solver knows
about a gas — which species exist, how fast they are made and destroyed, how
quickly they drift, how many electrons a returning ion releases from the
cathode — is supplied to it, from files outside the package under
``mechanisms/``.  The solver never derives a rate, fits a coefficient or
assumes a reaction; it reads what a mechanism gives it and solves the system
that results.  This chapter describes those files.

Two kinds of file
-----------------

A **mechanism file** describes a gas: a JSON file listing its species,
reactions, rate coefficients, photoionization model and electrode yields,
each coefficient a formula, and the parameters that may be changed with
their defaults (:ref:`Chap:MechanismFiles`).  It is checked when it is
loaded, so a misspelt name or a reaction that does not conserve charge is
reported before anything is solved.

A **configuration file** sets that mechanism's parameters — another
cross-section database, a different cathode yield, a reaction switched off
— and names the result (:ref:`Chap:Configuration`).  A file may hold many
configurations, one curve each.

A run names the mechanism, and optionally configuration files:

.. code-block:: bash

   incept1d pdiv gas.json studies.json --p 1.0

Mechanisms live one per directory, with the configurations written for
them; mechanisms of the same gas sit side by side under a gas directory,
which also holds the data they share:

.. code-block:: text

   mechanisms/<gas>/
   ├── <scheme>/
   │   ├── <gas>_<scheme>.json   the mechanism file
   │   ├── baseline.json         configuration files
   │   └── sensitivity.json
   └── data/                     tables the mechanisms read

Python mechanisms
-----------------

Some models do not fit the format: a rate that changes form at a threshold,
a photon coupling that is not a source times an absorption, a model that is
not a list of reactions at all.  Such a mechanism is a Python module
instead, exposing the matrices the solver needs as functions
(:ref:`Chap:NewMechanisms`), with an optional ``config.py`` that turns
configuration files into its parameters (:ref:`Chap:ConfigPy`).  The solver
sees both kinds the same way, and a run names a ``.py`` mechanism exactly as
it names a ``.json`` one.  Use it when the mechanism file cannot say what
the model needs, and not otherwise: a Python mechanism is executed, so
running someone else's means trusting their code, and nothing in it is
checked beyond the shapes the solver receives.

Where to go next
----------------

* :ref:`Chap:MechanismFiles` — The mechanism file, block by block.
* :ref:`Chap:Configuration` — Configuration files: parameters, reaction
  multipliers, per-polarity overrides.
* :ref:`Chap:ModifyingReactions` — Scaling, changing, adding and removing
  reactions.
* :ref:`Chap:NewMechanisms` and :ref:`Chap:ConfigPy` — Python mechanisms.

Three complete mechanisms are documented as worked examples
(:ref:`Chap:ExampleConfigurations`): the classical Townsend model, a
six-species scheme for dry air with photoionization, and relativistic
feedback in a thundercloud, which is a Python mechanism.

API reference
-------------

:func:`incept1d.mechanism.load_mechanism` reads either kind of mechanism
with one configuration and returns what the solver uses; every command
starts here.

.. automodule:: incept1d.mechanism
   :members:
   :undoc-members:
