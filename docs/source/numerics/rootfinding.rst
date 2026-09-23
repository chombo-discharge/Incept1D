.. _Chap:Numerics:RootFinding:

Root finding and branch tracking
================================

.. contents:: On this page
   :local:
   :depth: 1

Root finding in :math:`E/N`
---------------------------

:func:`incept1d.inception.find_all_breakdown_EN` locates the :math:`E/N` roots of
:math:`\det\bm{Q}(0; E/N)` at a single :math:`pd`:

1. **Coarse scan.**  :math:`\det\bm{Q}` is evaluated on 200 log-spaced
   points between 10 Td and :math:`3\times10^5` Td and every sign change is
   bracketed.  For non-uniform fields the scan uses a cheap proxy — a
   single-step uniform-field determinant, with a 5-step version as fallback
   — and the bracket is re-validated (and widened if needed) with the full
   determinant before refinement.
2. **Refinement.**  Each bracket is refined with
   :func:`scipy.optimize.brentq` (``xtol=1e-8``, ``rtol=1e-14``).  ``NaN``
   values are mapped to a tiny negative sentinel (:math:`-10^{-300}`) so
   that Brent's method converges *through* the singular region rather than
   failing.
3. **Acceptance.**  :func:`incept1d.inception._accept_root` distinguishes a genuine
   root (finite bracket endpoints of opposite sign; ``NaN`` at the root
   itself is then the expected consequence of :math:`\bm{Q}` being exactly
   singular) from an artefact created by the sentinel at a bracket
   endpoint (discarded with a warning).  A finite but large residual is
   accepted with a warning.
4. **Hidden roots.**  A sign change lying just below the ``NaN`` boundary
   would be invisible to the coarse scan; when the scan goes from positive
   directly to ``NaN`` the interval is re-scanned finely with the full
   determinant.

By default only the lowest root is kept (``first_only``); with
``--all-branches`` every root is returned.

Warm starts
...........

Along a :math:`pd` sweep the root moves smoothly, so in single-branch mode
the previous root is used as a hint: the bracket
:math:`[E/N_\mathrm{prev}/2, 2E/N_\mathrm{prev}]` is tried first.  If the
upper end of the hint bracket is in the ``NaN`` region it is narrowed
geometrically until a finite negative value is found.  Warm starts are
disabled in ``--all-branches`` mode because new branches can appear at any
:math:`pd`.

A hint only records where the root was at the *previous* :math:`pd` point,
so refining it is not by itself enough: a new, lower root may have appeared
since.  This is not hypothetical — near the left-branch asymptote the lowest
root can drop by two orders of magnitude between two adjacent :math:`pd`
points, and in an electronegative gas it does.  The
warm-start root is therefore accepted, and the coarse scan skipped, only
when a second scan over :math:`[E/N_\mathrm{lo}, E/N_\mathrm{warm}]` finds
no sign change beneath it.  That scan samples at the same points per decade
as the full scan, so accepting a warm-start root is exactly as reliable as
running the scan it replaces — only cheaper, because the interval is
shorter.  Without this check a single spurious root at one end of the sweep
propagates through every subsequent :math:`pd` point.

Branch tracking
---------------

:func:`incept1d.inception.compute_inception_curve` repeats the root search over the
:math:`pd` grid and assigns each new root to an existing branch by greedy
nearest-neighbour matching in :math:`\log(E/N)`; unmatched roots start new
branches.  This handles saddle-node bifurcations, where a pair of new
low-:math:`E/N` roots appears mid-sweep, without corrupting the
pre-existing branch.  Branch 0 is the branch that first appeared at the
lowest :math:`pd`.

Root finding in :math:`\lambda`
-------------------------------

:func:`incept1d.growth.find_lambda_for_voltage` solves
:math:`\det\bm{Q}(\lambda; E/N) = 0` at fixed :math:`E/N` above threshold,
in three steps:

1. Evaluate :math:`\det\bm{Q}(0)` to establish the reference sign
   (``NaN`` is treated as negative, per the convention above).
2. Expand an upper bracket geometrically, :math:`\lambda_\mathrm{hi} = 1,
   10, 100, \ldots` s\ :sup:`-1`, until the sign flips — a large positive
   :math:`\lambda` always damps the solution, since :math:`\bm{A} =
   (\bm{R} - \lambda\bm{I})\bm{V}^{-1}`.
3. Refine :math:`[0, \lambda_\mathrm{hi}]` with Brent's method.
