.. _Chap:Inception:

Inception.py — Paschen curves
=============================

.. contents:: On this page
   :local:
   :depth: 1

``Inception.py`` is the core solver.  It implements the boundary-value
inception criterion :math:`\det\bm{Q}(\lambda) = 0` of
:ref:`Chap:InceptionCriterion` and provides the command-line interface that
computes generalized Paschen curves — the reduced field :math:`(E/N)^*` and
voltage :math:`U^*` at which the criterion is met, as a function of
:math:`pd`.

What it does
------------

For every point of a logarithmic :math:`pd` sweep the script

1. evaluates :math:`\det\bm{Q}(0; E/N)` on a coarse logarithmic scan of
   :math:`E/N` to find every sign change,
2. refines each bracket with Brent's method to obtain the root(s),
3. matches the roots to those found at the previous :math:`pd` point so that
   each solution *branch* is tracked continuously,
4. converts :math:`(E/N)^*` to the breakdown voltage :math:`U^* = (E/N)^*
   N d`.

Steps 1–2 are :func:`Inception.find_all_breakdown_EN`; steps 3–4 are
:func:`Inception.compute_paschen_curve`.  The sweep can be run in two modes:
**fixed pressure** (``--p``, vary :math:`d`) or **fixed gap distance**
(``--d``, vary :math:`p`); both may be combined, and several pressures or
distances may be given.  For a non-uniform field that is not symmetric
(sphere-plane, field line) both polarities are computed.

Command-line interface
----------------------

.. code-block:: console

   python3 Inception.py MECHANISM.py [CONFIG.json ...]
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

Uniform-field Paschen curve at three pressures, all branches:

.. code-block:: bash

   python3 Inception.py Air/Air_Pancheshnyi.py --p 0.1 1 10 --all-branches --plot-separate-branches

Fixed 10 mm gap, varying pressure, comparing the four cross-section
databases:

.. code-block:: bash

   python3 Inception.py Air/Air_Pancheshnyi.py Air/Databases.json --d 10 --pd-min 0.1 --pd-max 200

Sphere-plane gap (50 mm sphere), both polarities, second-order Magnus
propagator on a fixed 40-step grid:

.. code-block:: bash

   python3 Inception.py Air/Air_Pancheshnyi.py --field sphere-plane 50 --method magnus2 --dx 40 40

Paschen curve with the growth-rate contour :math:`\lambda = 10^{8}` s\ :sup:`-1`:

.. code-block:: bash

   python3 Inception.py Air/Air_Pancheshnyi.py --lam 1e8

Output file format
------------------

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

Loading a mechanism
-------------------

:func:`load_mechanism` reads a mechanism ``.py`` file with ``importlib`` —
mechanism files are *not* imported as packages, they are executed as
standalone modules, with optional override variables injected into their
namespace beforehand by a companion ``Config.py`` living next to the
mechanism file (:ref:`Chap:Configuration`).  The loaded module is validated
against the required interface (:ref:`Chap:NewMechanisms`) and wrapped in a
:class:`Mechanism` object that bakes in every configuration parameter
(photoionization/photoemission scale factors, reaction-rate multipliers,
per-polarity SEE overrides), so every other function can treat
``mod.get_R(EN, p, T)`` etc. as the final, fully-configured answer.

.. literalinclude:: ../../../Inception.py
   :language: python
   :pyobject: load_mechanism

Assembling the augmented ODE
----------------------------

:func:`_build_A_aug` builds :math:`\bm{\mathcal{A}}` of
:eq:`eq_augmented_ode` at a given reduced field, pressure, temperature and
:math:`\lambda`.

.. literalinclude:: ../../../Inception.py
   :language: python
   :pyobject: _build_A_aug

Boundary conditions and the determinant
---------------------------------------

:func:`_assemble_det_Q` builds :math:`\bm{Q} = [\bm{Q}_0; \bm{Q}_d]`
(:eq:`eq_Q0`, :eq:`eq_Qd`) from the propagator :math:`\bm{M}(d)` and the
mechanism's row selectors and SEE yields, and returns :math:`\det\bm{Q}`.

.. literalinclude:: ../../../Inception.py
   :language: python
   :pyobject: _assemble_det_Q

:func:`inception_det` ties the two together: it integrates
:math:`\bm{\mathcal{A}}(x)` across the gap to build :math:`\bm{M}(d)`, then
calls ``_assemble_det_Q``.  Its root in :math:`E/N` at fixed :math:`pd` *is*
the inception field.

.. literalinclude:: ../../../Inception.py
   :language: python
   :pyobject: inception_det

API reference
-------------

.. automodule:: Inception
   :members:
   :undoc-members:
   :private-members: _build_A_aug, _assemble_det_Q, _adaptive_midpoint_segment, _expm_shifted, _accept_root
