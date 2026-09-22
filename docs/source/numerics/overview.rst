.. _Chap:Numerics:

Overview
========

The theory pages end with a determinant condition, :math:`\det\bm{Q}(\lambda)
= 0`, whose ingredients are a matrix exponential (or a product of them) and
a handful of small linear-algebra operations.  This chapter describes how
that condition is actually evaluated and solved:

* :ref:`Chap:Numerics:Propagator` — how the propagator :math:`\bm{M}(d)` of
  :eq:`eq_theta_soln` is built for uniform and non-uniform fields, including
  the adaptive grid, the second-order Magnus alternative, eigenvalue
  shifting against overflow, and the local treatment of strongly absorbed
  photon groups.
* :ref:`Chap:Numerics:Determinant` — how :math:`\det\bm{Q}` is evaluated so
  that its *sign* is reliable even when :math:`\bm{Q}` is nearly singular,
  and the ``NaN`` convention that the root finders rely on.
* :ref:`Chap:Numerics:RootFinding` — locating the roots in :math:`E/N` at
  fixed :math:`pd` (coarse scan, bracketing, Brent refinement, acceptance
  tests, warm starts), tracking solution branches across a :math:`pd`
  sweep, and solving for the growth rate :math:`\lambda` above threshold.
* :ref:`Chap:Numerics:Cost` — what a calculation costs and which knobs
  trade accuracy for time.

All of it lives in :mod:`incept1d.solver` (and :mod:`incept1d.growth` for the growth
rate); the pages cite the implementing function for every step, and the
command-line options that control it are summarised in
:ref:`Chap:Inception`.

Conventions
-----------

Throughout this chapter :math:`h` is the width of an integration segment,
:math:`\|\cdot\|_F` the Frobenius norm, and ``expm`` the matrix
exponential as implemented by :func:`scipy.linalg.expm` (Padé
approximation with scaling and squaring).  The augmented matrix
:math:`\bm{\mathcal{A}}` is small — :math:`N_s + 2N_\gamma^\mathrm{eff}`,
i.e. 6 to 18 for the air schemes — so every linear-algebra operation is
cheap; the cost of a calculation is dominated by the *number* of
propagator evaluations.
