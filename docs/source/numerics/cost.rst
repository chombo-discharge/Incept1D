.. _Chap:Numerics:Cost:

Cost and accuracy
=================

.. contents:: On this page
   :local:
   :depth: 1

What a calculation costs
------------------------

The cost of a calculation is the number of criterion evaluations times the
cost of one.  Measured on a desktop, for a chemistry with six species and
three photon groups, all explicit (:math:`\dim\bm{\mathcal{A}} = 12`):

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Evaluation
     - Time
     - Notes
   * - reflection criterion, uniform field
     - ≈ 1 ms
     - One exponent, split into well-conditioned substeps combined by
       adding–doubling; ``--dx`` and ``--method`` do not apply
   * - reflection criterion, sphere gap, default ``--dx``
     - ≈ 10 ms
     - Adaptive midpoint from a field-following start, :math:`N_{\min} = 5`,
       :math:`N_{\max} = 200`, tolerance 0.03.  The same grid costs about as
       much for the direct :math:`\det\bm{Q}`
   * - :math:`\det\bm{Q}` with an optically thick photon group
     - 0.1 s to minutes
     - Takes the compound-matrix route (:ref:`Chap:Numerics:Determinant`);
       why ``--criterion detq`` and ``--verify`` are not the default
   * - one :math:`pd` point of an inception curve
     - 30–300 evaluations
     - Warm start plus the guard scan below it, or a bottom-up scan that
       stops at the first root, plus refinement
   * - 40-point inception curve, strongly non-uniform gap, both polarities
     - ≈ 10 s on 16 cores
     - ≈ 30 s with ``--jobs 1``

Non-uniform fields are roughly an order of magnitude more expensive than
uniform ones: the uniform propagator is a single exact ``expm``, whereas a
non-uniform gap needs one :math:`\bm{\mathcal{A}}` assembly and one ``expm``
per grid point, and the adaptive grid places extra points near the
high-field electrode.  The
following controls trade accuracy for time:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Control
     - Effect
   * - ``--pd-num N``
     - Number of :math:`pd` points; cost is linear in it.
   * - ``--jobs N``
     - Worker processes (default: the number of physical cores).  The
       :math:`pd` points of ``pdiv``, and the voltages and cases of
       ``growth``, are solved concurrently; the answer does not change.
   * - ``--criterion detq``
     - Evaluates :math:`\det\bm{Q}` instead of the reflection criterion:
       the same roots, but far slower wherever a photon group is optically
       thick.
   * - ``--verify``
     - Adds two :math:`\det\bm{Q}` evaluations per root; with a thick
       photon group these can cost far more than the root search itself.
   * - ``--dx N_min N_max tol``
     - Fewer initial segments, a smaller refinement budget or a looser
       tolerance reduce the number of ``expm`` calls per evaluation.
       ``--dx 5 25 0.05`` (used for the ELECTRA example) is 5–10× cheaper
       than the default ``5 200 0.03`` at small :math:`pd`, and less accurate.
   * - ``--method magnus2 --dx N N``
     - Fourth-order propagator on a fixed grid; no adaptive refinement.
       Cheap and accurate for smooth profiles when :math:`N` is chosen
       adequately (see below).
   * - ``--all-branches``
     - Disables warm starts, so every :math:`pd` point pays for a full
       coarse scan.
   * - ``ngroups`` (configuration)
     - Each photon group adds two rows and columns to
       :math:`\bm{\mathcal{A}}`.  The reflection criterion's cost grows
       mildly with them; the compound route of :math:`\det\bm{Q}` grows
       combinatorially.

Checking accuracy
-----------------

The solver has no built-in error estimate for the final :math:`E/N`; the
adaptive tolerance controls the propagator, not the root.  To check a
result:

* **Grid convergence.**  Tighten ``--dx`` (e.g. ``5 1000 1e-3``, about
  0.01 %) or double a fixed ``magnus2`` grid and compare the curves.  The
  defaults keep :math:`E/N` within about 0.1–0.2 % of the converged value
  on sphere-plane and thin-wire coaxial gaps.
* **Propagator cross-check.**  ``--method midpoint`` and ``--method
  magnus2`` converge to the same answer; a persistent difference points to
  an under-resolved profile.
* **Closed-form limit.**  Reducing the chemistry to a single ionizing
  reaction with no attachment gives a case with an analytic answer
  (:eq:`eq_standard_paschen`) that the solver must reproduce.
* **Independent check.**  ``--verify`` evaluates :math:`\det\bm{Q}` on the
  same field and grid just below and just above every root, and marks the
  root ✓ if it changes sign (:ref:`Chap:Numerics:Riccati`).
* **Root messages.**  ``[root check]`` and ``[det Q check]`` messages
  printed during a run flag rejected or suspect roots
  (:ref:`Chap:Numerics:RootFinding`); a clean run prints none.
