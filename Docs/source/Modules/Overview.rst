.. _Chap:ModulesOverview:

Overview of the scripts
=======================

.. contents:: On this page
   :local:
   :depth: 1

``Incept1D`` has no package or build system.  Every module lives at the
repository root and imports its siblings directly; every script is run from
the repository root with ``python3 <script>.py``.  The modules fall into
three layers:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Role
   * - :mod:`Inception`
     - **Core solver and Paschen-curve CLI.**  Loads a mechanism, integrates
       the augmented ODE across the gap, evaluates
       :math:`\det\bm{Q}(\lambda)`, finds and tracks the :math:`E/N` roots
       over a :math:`pd` sweep.  See :ref:`Chap:Inception`.
   * - :mod:`Eigenvalues`
     - Diagnostic: eigenvalues of the *local* transport matrix
       :math:`\bm{R}\bm{V}^{-1}` vs. :math:`E/N`, no gap integration.  See
       :ref:`Chap:Eigenvalues`.
   * - :mod:`IonizationIntegral`
     - The classical ionization integral :math:`I_\alpha` and the apparent
       ionization integral :math:`I_\lambda` vs. applied voltage, for
       comparison with the full criterion.  See
       :ref:`Chap:IonizationIntegral`.
   * - :mod:`Lambda`
     - Temporal growth rate :math:`\lambda > 0` vs. voltage above the
       inception voltage.  See :ref:`Chap:Lambda`.
   * - :mod:`CreateChomboDischargeData`
     - Exports transport/rate-coefficient tables from a mechanism for the 3-D
       solver ``chombo-discharge``.  Not part of the inception solve.  See
       :ref:`Chap:CreateChomboDischargeData`.
   * - :mod:`FieldDistributions`
     - Gap geometry: uniform / sphere-plane / sphere-sphere / tabulated field
       line profiles :math:`f(\xi)` and the shared ``--field`` CLI parsing.
       See :ref:`Chap:FieldDistributions`.
   * - :mod:`Reactions`
     - Declarative reaction-string parser that assembles :math:`\bm{R}`.
       Used by mechanism files.  See :ref:`Chap:Reactions`.
   * - :mod:`Constants`
     - Physical constants from ``scipy.constants``.  See :ref:`Chap:Constants`.

Call graph
----------

.. code-block:: text

   mechanism.py (+ Config.py, *.json) ──► Inception.load_mechanism ──► Mechanism
                                                                          │
                     ┌──────────────────────┬──────────────────┬──────────┴────────────┐
                     ▼                      ▼                  ▼                       ▼
             Inception.main         Eigenvalues.main   IonizationIntegral.main    Lambda.main
           (Paschen curves)       (local eigenvalues)   (ionization integrals)   (growth rate)

:mod:`Reactions` and :mod:`FieldDistributions` sit underneath everything;
:mod:`Constants` sits under all of them.  :mod:`CreateChomboDischargeData`
has its own lightweight mechanism loader since it needs the raw rate
functions rather than the :math:`\bm{R}` matrix.

Common command-line conventions
-------------------------------

All four solver scripts share the same first arguments and most options:

.. code-block:: console

   python3 <script>.py MECHANISM.py [CONFIG.json ...] [options]

``MECHANISM.py``
   Path to a mechanism file, e.g. ``Air/Air_Pancheshnyi.py``
   (:ref:`Chap:NewMechanisms`).

``CONFIG.json ...``
   Zero or more JSON configuration files.  Each may contain a single
   configuration object or a list under a ``"configurations"`` key; all
   configurations across all files are run in sequence and plotted in the
   same figure.  With no file, a single baseline configuration with the
   mechanism's defaults is used.  See :ref:`Chap:Configuration`.

``--p P``, ``--T T``
   Gas pressure in bar and temperature in K (defaults 1.0 bar, 293 K).

``--field SPEC``
   Gap geometry: ``uniform`` (default), ``sphere-plane R_mm``,
   ``sphere-sphere R_mm``, or ``fieldline FILE [UNIT]``.  See
   :ref:`Chap:FieldDistributions` and :ref:`Chap:FieldLines`.

``--dx [N_min [N_max [tol]]]``, ``--method {midpoint,magnus2}``
   Integration grid and propagator for non-uniform fields; see
   :ref:`Chap:Numerics`.

``--write-to-file FILE``
   Write results as a tab-separated table with a metadata header (date, git
   commit, full command line, configuration labels, one column group per
   curve).

``--no-plot``
   Skip the matplotlib figure (``Inception.py``, ``IonizationIntegral.py``,
   ``Lambda.py``).  On headless machines also set ``MPLBACKEND=Agg``.

Units on the command line are bar, mm, kV, K and Td; the output files use
the same units and label every column.

Typical workflows
-----------------

**Paschen curve for a gas** — ``Inception.py`` with a mechanism and,
optionally, configuration files for sensitivity variants:

.. code-block:: bash

   python3 Inception.py Air/Air_Pancheshnyi.py Air/Databases.json --p 1 --pd-min 1e-2 --pd-max 1e3

**Understand why the curve looks the way it does** — inspect the local
eigenvalues and the ionization integrals:

.. code-block:: bash

   python3 Eigenvalues.py Air/Air_Pancheshnyi.py --p 1 --EN-lo 20 --EN-hi 300
   python3 IonizationIntegral.py Air/Air_Pancheshnyi.py --p 10 --d 10 --voltage-lo 20 --voltage-hi 60

**How fast does the discharge grow above threshold?**

.. code-block:: bash

   python3 Lambda.py Air/Air_Pancheshnyi.py --pd 10 --p 1 --v-max-factor 1.5

**Hand the chemistry to a 3-D simulation:**

.. code-block:: bash

   python3 CreateChomboDischargeData.py Air/Air_Pancheshnyi.py --write-to-file air_transport.dat

The remaining pages of this section document each script, the configuration format
(:ref:`Chap:Configuration`), and how to modify or write mechanisms
(:ref:`Chap:ModifyingReactions`, :ref:`Chap:NewMechanisms`).
