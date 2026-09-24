.. _Chap:ConfigPy:

The ``config.py`` protocol
==========================

.. contents:: On this page
   :local:
   :depth: 1

Why a mechanism needs one
-------------------------

A mechanism module is **executed once**, by ``importlib``, and everything it
computes at module level — reading a swarm-data table, fitting a photon
decomposition, building the reaction list — happens during that single
execution.  Nothing downstream can reach back in afterwards and change it.

That is a problem as soon as you want to ask a question of the form "what if
this were different?", which is most of what the code is for:

* Which cross-section database should the swarm data come from?
* What is the cathode yield, given that it is uncertain by orders of
  magnitude?
* How much does the answer move if detachment is switched off?
* Does the positive polarity of this gap want different surface parameters
  from the negative one?

None of these are new physics, so none of them should require a new mechanism
file.  But the first two must be settled *before* the module runs, because the
module reads them at import.

``config.py`` exists to bridge that gap.  It is a small adapter that turns an
inert JSON dictionary (:ref:`Chap:Configuration`) into two distinct things: a
set of values injected into the module namespace **before** it executes, and a
set of parameters applied **afterwards**, every time an accessor is called.
With it, one mechanism file supports an unlimited family of variants, each a
few lines of JSON, and several of them can be solved and plotted side by side
in a single run.

A mechanism with no ``config.py`` is still perfectly valid — it is then a
fixed gas with no adjustable parameters.  Passing a configuration file to such
a mechanism raises ``ImportError`` rather than silently ignoring it, so a
mistyped path cannot quietly leave the defaults in place.

Where the file lives
--------------------

:func:`incept1d.mechanism.load_mechanism` looks for a file called exactly
``config.py`` in the same directory as the mechanism module, and expects it to
define a class called ``Config``.  It is executed, not imported, so it must
not rely on being part of a package.  Code that the configuration classes of
several mechanisms share belongs in a separate file, loaded with
:func:`incept1d.mechanism.load_helper`; the ``config.py`` next to each
mechanism then only states what differs, such as which keys it accepts.

The interface
-------------

``Config`` must provide one class method, three instance methods and one
attribute.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Member
     - Contract
   * - ``from_dict(d, cfg_dir) -> Config``
     - Class method.  Builds an instance from the raw JSON dictionary.
       ``cfg_dir`` is the directory of the JSON file, so that relative paths
       inside it (a cross-section file, say) can be resolved.  **Unknown keys
       must be ignored**, so that one configuration file can be shared by
       mechanisms that understand different subsets of it.
   * - ``pre_exec_vars() -> dict``
     - Names injected into the module namespace *before* the module body
       runs.  This is the only way to influence module-level constants,
       because the body executes exactly once.  Returning an empty dict is
       fine.
   * - ``post_exec_init(mod) -> None``
     - Called after the module body has executed and passed the interface
       check.  For work that needs the finished module — fitting a photon
       decomposition, for instance, which cannot run until the rate functions
       exist.  May be a no-op.
   * - ``mechanism_params() -> dict``
     - Keyword arguments forwarded to :class:`~incept1d.mechanism.Mechanism`.
       These are applied at *call* time inside the accessors, so they can vary
       per polarity.  The accepted keys are fixed by ``Mechanism``:
       ``reaction_multipliers``, ``xi_photo``, ``xi_emit``, ``gamma0``,
       ``gamma1``, ``eref``, ``beta``, ``pos_override``, ``neg_override``.
   * - ``label``
     - Attribute.  The string shown in tables, legends and file headers.  One
       curve per configuration, so make it short and distinguishing.

Choosing between the two override paths
---------------------------------------

The commonest mistake when writing a new ``config.py`` is putting a parameter
in the wrong one, so the rule is worth stating plainly.

**Use** ``pre_exec_vars`` **when the module needs the value while it is being
executed.**  Which data file to read is the archetype: the module opens it at
module level, so the name has to be present beforehand.  The module picks such
values up with ``globals().get``, which returns the injected value if there is
one and the default otherwise:

.. code-block:: python

   CROSS_SECTIONS = globals().get(
       "CROSS_SECTIONS", os.path.join(_HERE, "default_swarm_data.txt")
   )

**Use** ``mechanism_params`` **when the value only affects what the accessors
return.**  Rate multipliers, the photoionization and photoemission scale
factors and the per-polarity surface parameters all belong here.  They cost
nothing to change, they do not require re-executing the module, and —
importantly — they are the only ones that can differ between the positive and
negative polarity of the same run, because
:meth:`incept1d.mechanism.Mechanism.resolve` applies them per polarity.

If a parameter could plausibly go either way, prefer ``mechanism_params``: a
value read at import is fixed for the whole run, while one applied in the
accessor can still be varied afterwards.

A minimal implementation
------------------------

Nothing in the protocol requires a dataclass, and a mechanism with a single
tunable parameter needs only this:

.. code-block:: python

   from dataclasses import dataclass

   @dataclass
   class Config:
       label: str = ""
       gamma0: float = 1e-2

       def pre_exec_vars(self):
           return {"_GAMMA0": self.gamma0}

       def post_exec_init(self, mod):
           pass

       def mechanism_params(self):
           return {}

       @classmethod
       def from_dict(cls, d, cfg_dir=None):
           known = {f.name for f in cls.__dataclass_fields__.values()}
           return cls(**{k: v for k, v in d.items() if k in known})

The ``from_dict`` body above is the idiom worth copying: filtering against the
dataclass fields is what makes unknown keys harmless, which is what lets one
JSON file be shared between mechanisms that understand different subsets of
it.

The test suite exercises exactly this shape in ``tests/toy/config.py``, so it
is a working starting point rather than an illustration.

Worked examples
---------------

:ref:`Chap:PaschenMechanism` shows the protocol at its simplest: one key
selects the gas, which has to be known before the module body runs because it
fixes the coefficients at module level.

:ref:`Chap:AirScheme` shows it carrying a realistic load: a choice of
cross-section database, four cathode-yield parameters, per-reaction
multipliers, photon scale factors, and separate surface parameters for the two
polarities of an asymmetric gap.
