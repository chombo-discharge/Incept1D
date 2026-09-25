.. _Chap:Numerics:RootFinding:

Root finding and branch tracking
================================

.. contents:: On this page
   :local:
   :depth: 1

Root finding in :math:`E/N`
---------------------------

:func:`incept1d.inception.find_all_breakdown_EN` locates the :math:`E/N` roots of
the inception criterion at a single :math:`pd` — the reflection criterion
:math:`g` by default (:ref:`Chap:Numerics:Riccati`), or
:math:`\det\bm{Q}(0; E/N)` with ``--criterion detq``:

1. **Below the scan range.**  :math:`g` is positive below inception and
   negative above for every mechanism, so one evaluation at the bottom of
   the scan range, 10 Td, tells whether inception lies below it.  If it
   does — a strongly non-uniform gap can put its gap-averaged inception
   field there — :math:`E/N` is halved until :math:`g` turns positive
   (down to 0.1 Td) and that bracket is refined; no scan above it could
   find the root.  (:math:`\det\bm{Q}` has no such invariant sign and
   skips this step.)
2. **Coarse scan.**  The criterion is evaluated on 200 log-spaced points
   between 10 Td and :math:`3\times10^5` Td, from the bottom up, and every
   sign change is bracketed.  When only the lowest root is wanted the scan
   stops at the first confirmed root instead of paying for the points
   above it.  For non-uniform fields the scan uses a cheap proxy — a
   single-step uniform-field evaluation — and each bracket is confirmed
   with the full criterion before refinement.
3. **Refinement.**  Each bracket is refined with
   :func:`scipy.optimize.brentq` (``xtol=1e-8``, ``rtol=1e-14``).  For
   :math:`\det\bm{Q}`, ``NaN`` values are mapped to a tiny negative
   sentinel (:math:`-10^{-300}`) so that Brent's method converges
   *through* a singular region rather than failing.
4. **Acceptance.**  :func:`incept1d.inception._accept_root` decides whether a
   refined root is genuine.  For :math:`g` a root must be a crossing from
   :math:`+` to :math:`-` in increasing :math:`E/N`: any such crossing is
   accepted — a zero, or the pole at which :math:`g` jumps from about
   :math:`+1` to :math:`-1` when the cathode contributes little — and a
   crossing from :math:`-` to :math:`+` is rejected, since it cannot be
   inception.  For :math:`\det\bm{Q}` a root bracketed by the ``NaN``
   sentinel is discarded, and a finite but large residual is accepted
   with a warning.
5. **Hidden roots.**  For :math:`\det\bm{Q}`, a sign change lying just
   below a ``NaN`` boundary would be invisible to the coarse scan; when the
   scan goes from positive directly to ``NaN`` the interval is re-scanned
   finely with the full determinant.

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

The warm start matters most in non-uniform gaps, where the uniform-field
proxy misleads a cold search into unconfirmed brackets and fallback scans;
there, starting every point cold cost more than twice the work.  A
parallel sweep (``--jobs``, :func:`incept1d.parallel.parallel_map`)
therefore splits the :math:`pd` points into contiguous blocks, one per
worker, each solved in order with warm starts, so only the first point of
each block starts cold; the roots are assigned to branches afterwards in
:math:`pd` order, and the result is the same as a sequential sweep.

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

:func:`incept1d.growth.find_lambda_for_voltage` solves the inception
criterion at fixed :math:`E/N` above threshold, as a function of
:math:`\lambda`, for the growth rate of the dominant mode, which is the
**largest** real root (:ref:`Sec:Lambda:Range` explains why).  With the
default reflection criterion :math:`g` (:ref:`Chap:Numerics:Riccati`) that
root is also the only one: every loop gain falls as :math:`\lambda` grows,
so :math:`g` rises through zero once.  The search does not rely on that,
and works for :math:`\det\bm{Q}` as well:

1. Evaluate the criterion at :math:`\lambda = 0`.  If it has the
   sub-threshold sign, the voltage is below inception and
   :math:`\lambda^* = 0`.
2. Start from :math:`\lambda_\mathrm{hi} = 10\,\nu_\mathrm{ion}`, ten
   times the fastest local ionization rate in the gap — well above any
   rate at which the gap as a whole can grow — and raise it by decades
   until the criterion has the sub-threshold sign there.
3. Scan **down** by a factor of three until the sign changes.  Since no
   root lies above the dominant one, the first sign change met from above
   brackets it; scanning up from zero could stop at a slower mode instead.
   Two roots closer than the scan factor can be missed, but roots from
   different feedback loops are usually orders of magnitude apart.
4. Refine the bracket with Brent's method.
5. Check that the criterion changes sign across
   :math:`\lambda^*(1 \pm 10^{-4})`; a bracket that closes on a jump
   rather than a zero is reported as ``suspect``.

The voltages of a sweep are independent, and are solved concurrently on
``--jobs`` worker processes (:func:`incept1d.growth.solve_voltage`).
