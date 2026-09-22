.. _Chap:ConfigPy:

The ``config.py`` protocol
==========================

.. contents:: On this page
   :local:
   :depth: 1

A mechanism directory may contain a ``config.py`` defining a class called
``Config``.  It is the adapter between an inert JSON dictionary
(:ref:`Chap:Configuration`) and the mechanism module
(:ref:`Chap:NewMechanisms`), and it is what makes a mechanism configurable at
all.  Without it the module is a fixed gas: passing a configuration file to a
mechanism that has no ``config.py`` raises ``ImportError`` rather than
silently ignoring the file.

:func:`incept1d.mechanism.load_mechanism` discovers the file by name, in the
same directory as the mechanism module.  It is executed, not imported, so it
must not rely on being part of a package.

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
     - Names injected into the module namespace *before* the module body runs.
       This is the only way to influence module-level constants, because the
       body executes exactly once.  Returning an empty dict is fine.
   * - ``post_exec_init(mod) -> None``
     - Called after the module body has executed and passed the interface
       check.  For work that needs the finished module — the air family uses
       it to run ``init_photoionization(ngroups, cone_angle)``.  May be a
       no-op.
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
in the wrong one, so the rule is worth stating plainly:

**Use** ``pre_exec_vars`` **when the module needs the value while it is being
executed.**  Which cross-section file to read is the archetype: the module
opens it at module level, so the name has to be present beforehand.  The
module picks such values up with ``globals().get``:

.. code-block:: python

   BOLSIG_FILE = globals().get("BOLSIG_FILE", os.path.join(_HERE, "lisbon.txt"))

**Use** ``mechanism_params`` **when the value only affects what the accessors
return.**  Rate multipliers, the photoionization and photoemission scale
factors and the per-polarity secondary-emission overrides all belong here.
They cost nothing to change, they do not require re-executing the module, and
— importantly — they are the only ones that can differ between the positive
and negative polarity of the same run, because
:meth:`incept1d.mechanism.Mechanism.resolve` applies them per polarity.

Worked example
--------------

The air family's implementation is the reference.  ``from_dict`` filters the
raw dictionary against the dataclass fields, which is what makes unknown keys
harmless:

.. literalinclude:: ../../../mechanisms/air/config.py
   :language: python
   :pyobject: Config.from_dict

``pre_exec_vars`` maps configuration names onto the module-level constants the
mechanism reads at execution time:

.. literalinclude:: ../../../mechanisms/air/config.py
   :language: python
   :pyobject: Config.pre_exec_vars

and ``mechanism_params`` hands the run-time scalings to ``Mechanism``:

.. literalinclude:: ../../../mechanisms/air/config.py
   :language: python
   :pyobject: Config.mechanism_params

A minimal version
-----------------

Nothing in the protocol requires a dataclass or the air family's key set.  A
mechanism with a single tunable parameter needs only this:

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

The test suite exercises exactly this shape in ``tests/toy/config.py``, so it
is a working starting point rather than an illustration.
