.. _Chap:Numerics:Cost:

Cost and accuracy
=================

.. contents:: On this page
   :local:
   :depth: 1

What a calculation costs
------------------------

The cost of a calculation is the number of propagator evaluations times the
cost of one.  For the six-species dry-air scheme with three photon groups
(:math:`\dim\bm{\mathcal{A}} = 12`), measured on a laptop:

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Evaluation
     - Time
     - Notes
   * - :math:`\det\bm{Q}`, uniform field
     - ≈ 5 ms
     - one ``expm``, dominated by Python overhead in ``get_R``
   * - :math:`\det\bm{Q}`, sphere-sphere, default ``--dx``
     - ≈ 50 ms
     - adaptive midpoint, :math:`N_{\min} = 5`, :math:`N_{\max} = 200`
   * - one :math:`pd` point of an inception curve
     - 20–100 evaluations
     - warm start ≈ 20, full scan ≈ 200 + refinement
   * - 100-point inception curve, sphere gap, two configurations
     - minutes
     - the :ref:`Chap:Examples:IEC60052` runs

Non-uniform fields are 10–100× more expensive than uniform ones because the
adaptive grid places many small steps near the high-field electrode.  The
following controls trade accuracy for time:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Control
     - Effect
   * - ``--pd-num N``
     - Number of :math:`pd` points; cost is linear in it.
   * - ``--dx N_min N_max tol``
     - Fewer initial segments, a smaller refinement budget or a looser
       tolerance reduce the number of ``expm`` calls per evaluation.
       ``--dx 5 25 0.05`` (used for the ELECTRA example) is 5–10× cheaper
       than the default ``5 200 0.03`` at small :math:`pd`.
   * - ``--method magnus2 --dx N N``
     - Fourth-order propagator on a fixed grid; no adaptive refinement.
       Cheap and accurate for smooth profiles when :math:`N` is chosen
       adequately (see below).
   * - ``--all-branches``
     - Disables warm starts, so every :math:`pd` point pays for a full
       coarse scan.
   * - ``ngroups`` (configuration)
     - Each explicitly propagated photon group adds two rows and columns to
       :math:`\bm{\mathcal{A}}`; groups with :math:`\kappa_jd > 12` are
       folded into :math:`\bm{A}` at no cost (:ref:`Chap:Numerics:Propagator`).

Checking accuracy
-----------------

The solver has no built-in error estimate for the final :math:`E/N`; the
adaptive tolerance controls the propagator, not the root.  To check a
result:

* **Grid convergence.**  Tighten ``--dx`` (e.g. ``5 400 0.01``) or double
  a fixed ``magnus2`` grid and compare the curves; differences below a
  percent in :math:`E/N` are typical for the default settings.
* **Propagator cross-check.**  ``--method midpoint`` and ``--method
  magnus2`` converge to the same answer; a persistent difference points to
  an under-resolved profile.
* **Closed-form limit.**  ``mechanisms/Air/Paschen.json`` reduces the chemistry to the
  case with an analytic answer (:eq:`eq_standard_paschen`).
* **Root residual.**  The ``[det Q check]`` messages printed during a run
  flag roots whose residual is large relative to the bracket
  (:ref:`Chap:Numerics:RootFinding`); a clean run prints none.

Limits of the geometry model
----------------------------

A ``max d/R`` warning is printed when the gap exceeds four sphere radii.
Beyond that the on-axis bispherical field is still exact, but the gap is
so non-uniform that the single-field-line model itself (no transverse
structure, photons along the axis only) becomes a poor description of the
discharge; see :ref:`Chap:TheoryOverview`.  Sphere gaps in the examples
stay well inside :math:`d/R < 1`.
