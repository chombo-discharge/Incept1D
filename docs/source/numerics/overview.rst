.. _Chap:Numerics:

Overview
========

The theory pages end with a determinant condition, :math:`\det\bm{Q}(\lambda)
= 0`, whose ingredients are a matrix exponential (or a product of them) and
a handful of small linear-algebra operations.  This chapter describes how
that condition is actually evaluated and solved:

* :ref:`Chap:Numerics:Propagator` — how the propagator :math:`\bm{M}(d)` of
  :eq:`eq_theta_soln` is built for uniform and non-uniform fields: the
  adaptive grid, the second-order Magnus alternative and eigenvalue
  shifting against overflow.
* :ref:`Chap:Numerics:Riccati` — the default way of evaluating the
  condition: the anode condition carried back to the cathode as a
  reflection operator, which never forms the propagator and stays
  well-conditioned however optically thick the photon groups.
* :ref:`Chap:Numerics:Determinant` — evaluating :math:`\det\bm{Q}` itself,
  the reference the default is checked against, including the
  compound-matrix route for a numerically singular :math:`\bm{Q}`.
* :ref:`Chap:Numerics:RootFinding` — locating the roots in :math:`E/N` at
  fixed :math:`pd` (scan, bracketing, Brent refinement, acceptance tests,
  warm starts), tracking solution branches across a :math:`pd` sweep, and
  solving for the growth rate :math:`\lambda` above threshold.
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
:math:`\bm{\mathcal{A}}` is small — :math:`N_s + 2N_\gamma` is
of order ten for a realistic chemistry — so every linear-algebra operation is
cheap; the cost of a calculation is dominated by the *number* of propagator
evaluations.
