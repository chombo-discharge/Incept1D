.. _Chap:Inception:

Inception curves — ``incept1d pdiv``
=======================================

.. contents:: On this page
   :local:
   :depth: 1

``incept1d pdiv`` is the main command.  It answers: **at what voltage does
this gap break down?**  For a given gas, geometry and pressure it sweeps
:math:`pd` and reports the reduced field :math:`(E/N)^*` and inception
voltage :math:`U^*` at which the criterion
:math:`\det\bm{Q}(\lambda) = 0` of :ref:`Chap:InceptionCriterion` is met —
a generalized Paschen curve.

How it works
------------

For every point of a logarithmic :math:`pd` sweep the command takes four
steps:

1. Evaluate :math:`\det\bm{Q}(0; E/N)` on a coarse logarithmic scan of
   :math:`E/N` to find every sign change.
2. Refine each bracket with Brent's method to obtain the root(s).
3. Match the roots to those found at the previous :math:`pd` point, so that
   each solution *branch* is tracked continuously.
4. Convert :math:`(E/N)^*` to the breakdown voltage
   :math:`U^* = (E/N)^*\,N\,d`.

The criterion is derived in :ref:`Chap:InceptionCriterion` and how it is
evaluated is described in :ref:`Chap:Numerics`.

The sweep runs in either of two modes: **fixed pressure** (``--p``, vary
:math:`d`) or **fixed gap distance** (``--d``, vary :math:`p`).  Both may be
combined, and several pressures or distances may be given.  For a
non-uniform field that is not symmetric — sphere-plane, or a field line —
both polarities are computed.

Inputs
------

.. code-block:: console

   incept1d pdiv MECHANISM.py [CONFIG.json ...]
       [--p P [P ...]] [--d D [D ...]] [--T T]
       [--pd-min PD] [--pd-max PD] [--pd-num N]
       [--field SPEC] [--dx [N_min [N_max [tol]]]] [--method {midpoint,magnus2}]
       [--lam LAM] [--all-branches] [--plot-separate-branches]
       [--plot-ionization-integral] [--streamer-criterion C]
       [--write-to-file FILE] [--save-subplots] [--no-plot]

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Meaning
   * - ``--p P [P ...]``
     - Fixed-pressure mode: pressure(s) in bar.  Default 1.0 when ``--d``
       is not given.
   * - ``--d D [D ...]``
     - Fixed-distance mode: gap distance(s) in mm.
   * - ``--pd-min``, ``--pd-max``, ``--pd-num``
     - The :math:`pd` grid in bar·mm (default :math:`10^{-2}`–:math:`10^{3}`,
       50 log-spaced points).
   * - ``--T T``
     - Gas temperature in K (default 293).
   * - ``--field SPEC``
     - Geometry: ``uniform`` | ``sphere-plane R_mm`` | ``sphere-sphere R_mm``
       | ``fieldline FILE [UNIT]``.  See :ref:`Chap:FieldDistributions`.
   * - ``--dx N_min [N_max [tol]]``
     - Adaptive integration grid (defaults 5, 200, 0.03).  See
       :ref:`Chap:Numerics:Propagator`.
   * - ``--method``
     - Propagator: ``midpoint`` (default, adaptive) or ``magnus2``.
   * - ``--lam LAM``
     - Temporal growth rate :math:`\lambda` in s\ :sup:`-1` for the
       generalized criterion :math:`\det\bm{Q}(\lambda) = 0` (default 0).
       :math:`\lambda > 0` lowers the curve, :math:`\lambda < 0` raises it.
   * - ``--all-branches``
     - Find and plot *every* root in :math:`E/N`, not only the lowest one.
   * - ``--plot-separate-branches``
     - Give each branch its own line style and legend entry.
   * - ``--plot-ionization-integral``
     - Overlay :math:`I_\alpha` along the solution on a second axis.
       Requires ``alpha``/``eta`` in the mechanism.
   * - ``--streamer-criterion C``
     - Also solve :math:`I_\alpha = C` (streamer criterion) and plot/write
       that curve.
   * - ``--write-to-file FILE``
     - Tab-separated output with a metadata header.
   * - ``--save-subplots``
     - Save each panel as a separate PDF.
   * - ``--no-plot``
     - No figure.

Examples
--------

Uniform-field inception curve at three pressures, all branches:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py --p 0.1 1 10 --all-branches --plot-separate-branches

Fixed 10 mm gap, varying pressure, comparing the four cross-section
databases:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py mechanisms/air/databases.json --d 10 --pd-min 0.1 --pd-max 200

Sphere-plane gap (50 mm sphere), both polarities, second-order Magnus
propagator on a fixed 40-step grid:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py --field sphere-plane 50 --method magnus2 --dx 40 40

Inception curve with the growth-rate contour :math:`\lambda = 10^{8}` s\ :sup:`-1`:

.. code-block:: bash

   incept1d pdiv mechanisms/air/air_pancheshnyi.py --lam 1e8

Outputs
-------

Unless ``--no-plot`` is given, a two-panel figure is shown: the inception
voltage :math:`U^*(pd)` and the reduced field :math:`(E/N)^*(pd)`, one curve
per configuration and polarity.  The same numbers are printed as a table.

``--write-to-file`` writes one row per :math:`pd` point.  The header records
the date, git commit, full command line, mechanism, temperature, pressures,
configuration labels, field type, :math:`\lambda`, propagator and grid
settings, and then names every column.  For each configuration and polarity
there is a group of five columns:

.. code-block:: text

   # Column 1: pd_bar_mm
   # Column 2: p_bar[Baseline, p=1.0 bar (sphere=positive)]
   # Column 3: d_mm[Baseline, p=1.0 bar (sphere=positive)]
   # Column 4: EN_Td[Baseline, p=1.0 bar (sphere=positive)]
   # Column 5: U_kV[Baseline, p=1.0 bar (sphere=positive)]
   # Column 6: E_Vm[Baseline, p=1.0 bar (sphere=positive)]
   ...

followed, when requested, by ionization-integral columns, the
:math:`\alpha = \eta` reference curve, and streamer-criterion columns.
Points where no root was found are written as ``nan``.

API reference
-------------

Loading the mechanism and its configurations is
:mod:`incept1d.mechanism` (:ref:`Chap:ConfigurationOverview`); the modules
below are the solve itself.

``incept1d.solver``
~~~~~~~~~~~~~~~~~~~

.. automodule:: incept1d.solver
   :members:
   :undoc-members:
   :private-members: _build_A_aug, _assemble_det_Q, _adaptive_midpoint_segment, _expm_shifted

``incept1d.inception``
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: incept1d.inception
   :members:
   :undoc-members:
   :private-members: _accept_root

``incept1d.cli.pdiv``
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: incept1d.cli.pdiv
   :members: add_arguments, run
