.. _Chap:ModulesOverview:

Overview of the package
=======================

.. contents:: On this page
   :local:
   :depth: 1

``Incept1D`` is a Python package, ``incept1d`` (sources in
``src/incept1d/``), installed with ``pip install -e .`` and driven by one
console command with subcommands:

.. code-block:: console

   incept1d COMMAND [options]      # COMMAND: pdiv | eigenvalues | ionization
                                   #          | growth | field | chombo

Every subcommand is a thin front end in ``incept1d.cli.<command>`` that
parses arguments, calls the library and plots or writes tables; the physics
lives in the library modules, which can equally be imported from your own
scripts (``from incept1d.solver import inception_det``).  The modules fall
into three layers:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Role
   * - :mod:`incept1d.mechanism`
     - Loads a mechanism file and its JSON configurations into a
       :class:`~incept1d.mechanism.Mechanism`.  See :ref:`Chap:NewMechanisms`
       and :ref:`Chap:Configuration`.
   * - :mod:`incept1d.solver`
     - **Core solver.**  Builds the augmented ODE, integrates it across the
       gap with the midpoint / Magnus propagators and evaluates
       :math:`\det\bm{Q}(\lambda)`.  See :ref:`Chap:Inception` and
       :ref:`Chap:Numerics`.
   * - :mod:`incept1d.inception`
     - Finds the :math:`E/N` roots of :math:`\det\bm{Q} = 0` and tracks them
       over a :math:`pd` sweep (inception-curve branches).  Front end:
       ``incept1d pdiv``.  See :ref:`Chap:Inception`.
   * - :mod:`incept1d.eigenvalues`
     - Diagnostic: eigenvalues of the *local* transport matrix
       :math:`\bm{R}\bm{V}^{-1}` vs. :math:`E/N`, no gap integration.  See
       :ref:`Chap:Eigenvalues`.
   * - :mod:`incept1d.ionization`
     - The classical ionization integral :math:`I_\alpha` and the apparent
       ionization integral :math:`I_\lambda` vs. applied voltage, for
       comparison with the full criterion.  See
       :ref:`Chap:IonizationIntegral`.
   * - :mod:`incept1d.growth`
     - Temporal growth rate :math:`\lambda > 0` vs. voltage above the
       inception voltage.  See :ref:`Chap:Lambda`.
   * - :mod:`incept1d.chombo`
     - Exports transport/rate-coefficient tables from a mechanism for the 3-D
       solver ``chombo-discharge``.  Not part of the inception solve.  See
       :ref:`Chap:ThirdParty`.
   * - :mod:`incept1d.fields`
     - Gap geometry: uniform / sphere-plane / sphere-sphere / coaxial / tabulated field
       line profiles :math:`f(\xi)` and the shared ``--field`` CLI parsing.
       See :ref:`Chap:FieldDistributions`.
   * - :mod:`incept1d.reactions`
     - Declarative reaction-string parser that assembles :math:`\bm{R}`.
       Used by mechanism files.  See :ref:`Chap:Reactions`.
   * - :mod:`incept1d.constants`
     - Physical constants from ``scipy.constants``.  See :ref:`Chap:Constants`.
   * - :mod:`incept1d.output`
     - The metadata header (date, git revision, command line) shared by all
       ``--write-to-file`` outputs.
   * - :mod:`incept1d.cli`
     - The ``incept1d`` console command and one module per subcommand.

Call graph
----------

.. code-block:: text

   mechanism file (+ config.py, *.json) ──► incept1d.mechanism.load_mechanism ──► Mechanism
                                                                                    │
                        ┌───────────────────┬───────────────────┬───────────────────┤
                        ▼                   ▼                   ▼                   ▼
                incept1d.solver     incept1d.eigenvalues  incept1d.ionization  incept1d.growth
                incept1d.inception     (local eigenvalues)  (ionization ints.)   (growth rate)
                (inception curves)
                        ▲                   ▲                   ▲                   ▲
                incept1d.cli.pdiv  cli.eigenvalues       cli.ionization       cli.growth

:mod:`incept1d.reactions` and :mod:`incept1d.fields` sit underneath everything;
:mod:`incept1d.constants` sits under all of them.  :mod:`incept1d.chombo`
has its own lightweight mechanism loader
(:func:`~incept1d.chombo.load_raw_mechanism`) since it needs the raw rate
functions rather than the :math:`\bm{R}` matrix.

Common command-line conventions
-------------------------------

The four solver commands share the same first arguments and most options:

.. code-block:: console

   incept1d COMMAND MECHANISM.py [CONFIG.json ...] [options]

``MECHANISM.py``
   Path to a mechanism file, e.g. ``mechanisms/air/air_pancheshnyi.py``
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
   ``sphere-sphere R_mm``, ``coaxial A_mm B_mm``, or ``fieldline FILE [UNIT]``.  See
   :ref:`Chap:FieldDistributions` and :ref:`Chap:FieldLines`.

``--dx [N_min [N_max [tol]]]``, ``--method {midpoint,magnus2}``
   Integration grid and propagator for non-uniform fields; see
   :ref:`Chap:Numerics`.

``--write-to-file FILE``
   Write results as a tab-separated table with a metadata header (date, git
   commit, full command line, configuration labels, one column group per
   curve).

``--no-plot``
   Skip the matplotlib figure (``incept1d pdiv``, ``incept1d ionization``,
   ``incept1d growth``).  On headless machines also set ``MPLBACKEND=Agg``.

Units on the command line are bar, mm, kV, K and Td; the output files use
the same units and label every column.

Typical workflows
-----------------

Each command answers a different question about the same mechanism:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Question
     - Command
   * - At what voltage does the gap break down?
     - ``incept1d pdiv``
   * - Why does the curve look the way it does?
     - ``incept1d eigenvalues``, ``incept1d ionization``
   * - How fast does the discharge grow above threshold?
     - ``incept1d growth``
   * - What does the imported field profile look like?
     - ``incept1d field``
   * - How do I hand this chemistry to a 3-D solver?
     - ``incept1d chombo``

:ref:`Chap:QuickStart` works two of these through end to end.  The remaining
pages of this chapter document one command each; the mechanism and
configuration files they take as input are described in
:ref:`Chap:ConfigurationOverview`.
