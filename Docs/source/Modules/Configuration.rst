.. _Chap:Configuration:

Configuration files
===================

.. contents:: On this page
   :local:
   :depth: 1

Every solver script accepts JSON configuration files after the mechanism
path.  A configuration changes the *parameters* of a mechanism — cross
sections, rate multipliers, photoionization and SEE settings — without
editing the Python file, and several configurations can be run and plotted
side by side in one invocation.  This is how sensitivity studies are set
up.

File format
-----------

A file contains either a single configuration object, a bare list of
objects, or (the usual form) a list under a ``"configurations"`` key:

.. code-block:: json

   {
     "configurations": [
       {
         "label": "Baseline",
         "cross_sections": "Lisbon.txt"
       },
       {
         "label": "No detachment",
         "cross_sections": "Lisbon.txt",
         "reaction_multipliers": {
           "O2- + O2 -> e + O2 + O2": 0.0,
           "O- + N2 -> e + N2O":      0.0
         }
       }
     ]
   }

Configurations are independent — later entries do **not** inherit from
earlier ones — so repeat every key you need in each entry.  Keys beginning
with ``_`` (e.g. ``"_comment"``) are ignored.  Relative paths in
``cross_sections`` are resolved first against the directory of the JSON
file, then against the mechanism's directory.

Keys for the ``Air`` mechanism family
-------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Key
     - Default
     - Meaning
   * - ``label``
     - ``"Baseline"``
     - Name used in legends, terminal tables and output-file headers.
   * - ``cross_sections``
     - ``Lisbon.txt``
     - BOLSIG+ output file for electron transport and rates
       :math:`k_1, k_2, k_3` (:ref:`Chap:AirScheme`).  Shipped:
       ``Lisbon.txt``, ``Phelps.txt``, ``Biagi.txt``, ``Trinity.txt``,
       ``Morgan.txt``.
   * - ``reaction_multipliers``
     - ``{}``
     - Dictionary ``{reaction string: factor}``.  ``0.0`` disables a
       reaction, ``2.0`` doubles it.  See :ref:`Chap:ModifyingReactions`.
   * - ``ngroups``
     - ``3``
     - Number of two-stream photon groups in the Zheleznyak fit
       (:ref:`Chap:Photoionization`).
   * - ``cone_angle``
     - ``45.0``
     - Half-opening angle :math:`\theta_\mathrm{cone}` in degrees;
       :math:`\Delta\Omega/4\pi = (1 - \cos\theta)/2`.  ``0`` disables
       photon feedback.
   * - ``xi_photo``
     - ``1.0``
     - Scale factor on the photoionization coupling :math:`\bm{B}`.
   * - ``xi_emit``
     - ``1.0``
     - Scale factor on the photoemission yields :math:`\vec{\gamma}_\Psi`.
   * - ``gamma0``, ``gamma1``, ``eref``, ``beta``
     - ``1e-3``, ``0``, ``1.7e9``, ``1``
     - Ion-induced SEE yield
       :math:`\gamma = \gamma_0 + \gamma_1\exp(-E_\mathrm{ref}/\beta E)`
       (:ref:`Chap:SecondaryEmission`).  :math:`E_\mathrm{ref}` in V/m.
   * - ``positive``, ``negative``
     - —
     - Nested objects with any of ``xi_photo``, ``xi_emit``, ``gamma0``,
       ``gamma1``, ``eref``, ``beta`` that override the top-level values for
       one polarity only (the *cathode* differs between polarities in a
       non-symmetric gap).

Example with polarity-specific cathode yields:

.. code-block:: json

   {
     "label": "Steel sphere, brass plane",
     "positive": { "gamma0": 5e-4 },
     "negative": { "gamma0": 2e-3 }
   }

Shipped configuration sets
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - File
     - Contents
   * - ``Air/Baseline.json``
     - The reference configuration, all multipliers 1.
   * - ``Air/NoDetachment.json``
     - Baseline plus a variant with reactions 5 and 6 (collisional and
       associative detachment) switched off — the "no detachment" curve.
   * - ``Air/Databases.json``
     - Phelps, Trinity, Biagi and Lisbon cross-section sets, plus the
       no-detachment variant.
   * - ``Air/IonSensitivity.json``
     - Baseline; no :math:`\mathrm{O}^- \to \mathrm{O}_2^-` conversion; no
       :math:`\mathrm{O}^- \to \mathrm{O}_3^-` conversion; no detachment.
   * - ``Air/SEE.json``
     - :math:`\gamma_0 \in \{10^{-4}, 10^{-3}, 10^{-2}\}` and
       :math:`\theta_\mathrm{cone} \in \{0^\circ, 45^\circ, 90^\circ\}`.
   * - ``Air/Paschen.json``
     - The textbook limit: no detachment, no ion conversion, no photon
       feedback.  Reproduces :eq:`eq_standard_paschen`.
   * - ``Air/example_config.json``
     - Annotated example for the two-body mechanism ``Air_2body.py``.

How configurations are applied
------------------------------

The mechanism directory's ``Config.py`` defines a ``Config`` dataclass that
implements a three-step protocol used by :func:`Inception.load_mechanism`:

.. code-block:: python

   config.pre_exec_vars()      # -> dict injected into the module namespace
                               #    BEFORE the mechanism file is executed
   config.post_exec_init(mod)  # called AFTER execution (e.g. photoionization setup)
   config.mechanism_params()   # -> kwargs forwarded to Mechanism(...)
   config.label                # display label

``pre_exec_vars`` is how ``cross_sections`` and the ``gamma*`` parameters
reach the mechanism: the file reads them with ``globals().get("BOLSIG_FILE",
default)`` etc. at import time (:ref:`Chap:NewMechanisms`).
``post_exec_init`` calls ``init_photoionization(ngroups, cone_angle)``.
``mechanism_params`` supplies the run-time scalings (multipliers,
``xi_photo``, ``xi_emit``, polarity overrides) that
:class:`Inception.Mechanism` applies inside its accessor methods, so no
solver code ever sees a raw, un-configured module.

.. literalinclude:: ../../../Air/Config.py
   :language: python
   :pyobject: Config.from_dict

A new mechanism family gets its own ``Config.py`` implementing the same
protocol; copy ``Air/Config.py`` and adapt the key list.
