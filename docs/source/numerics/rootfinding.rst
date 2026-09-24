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
   single-step uniform-field determinant — and each bracket is confirmed
   with the full determinant before refinement.
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

Confirming a proxy bracket
..........................

The proxy is a different geometry, so it is reliable for *detecting* that a
root is nearby but not for locating it: for a strongly non-uniform field it
places the lowest root tens of per cent away from where the full
determinant has it, which is several scan steps.  A bracket is therefore
confirmed by re-scanning :math:`[E/N_a/1.6,\; E/N_b \times 1.6]` with the
full determinant, a window sized to that error rather than to the scan
resolution.

If the sign change is still not found, the bracket is *unconfirmed*: either
the proxy invented a root, or the real one lies outside even that window,
and only the full determinant can tell which.  The whole scan is then
repeated with it.  Skipping the bracket instead would be worse than slow —
it would silently promote the next root up to "lowest", which in an
electronegative gas can be two orders of magnitude higher in :math:`E/N`
and shows up as isolated points far above an otherwise smooth curve.

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
:math:`\det\bm{Q}(\lambda; E/N) = 0` at fixed :math:`E/N` above threshold
for the growth rate of the dominant mode, which is the **largest** real
root (:ref:`Sec:Lambda:Range` explains why).  Every evaluation uses the
compound-matrix route where :math:`\bm{Q}` is singular
(:ref:`Chap:Numerics:Determinant`), so no ``NaN`` convention is involved:

1. Evaluate :math:`\det\bm{Q}(0)`.  If it has the sub-threshold sign, the
   voltage is below inception and :math:`\lambda^* = 0`.
2. Start from :math:`\lambda_\mathrm{hi} = 10\,\nu_\mathrm{ion}`, ten times
   the fastest local ionization rate in the gap — well above any rate at
   which the gap as a whole can grow — and raise it by decades until
   :math:`\det\bm{Q}(\lambda_\mathrm{hi})` has the sub-threshold sign.
3. Scan **down** by a factor of three until the sign changes.  Since no
   root lies above the dominant one, the first sign change met from above
   brackets it; scanning up from zero could stop at a slower mode instead.
   Two roots closer than the scan factor can be missed, but roots from
   different feedback loops are usually orders of magnitude apart.
4. Refine the bracket with Brent's method.
5. Check that :math:`\det\bm{Q}` changes sign across
   :math:`\lambda^*(1 \pm 10^{-4})`.  Brent's method also converges on a
   discontinuity — :math:`\bm{Q}` changes size where a photon group crosses
   :math:`\kappa_j d = 12` — which is reported as ``suspect``.
