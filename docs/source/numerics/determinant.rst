.. _Chap:Numerics:Determinant:

Evaluating the determinant
==========================

This page describes ``--criterion detq``
(:func:`incept1d.solver.inception_det`): :math:`\bm{Q}` of
:eq:`eq_Q_system` is assembled from :math:`\bm{M}(d)` and its determinant
evaluated.  It is the formulation of the Theory chapter and the reference
the default criterion (:ref:`Chap:Numerics:Riccati`) is checked against,
by ``--verify`` and by the test suite.  Only the *sign* of
:math:`\det\bm{Q}` matters for root finding, and the evaluation is
designed around keeping that sign reliable.

:func:`incept1d.solver._det_Q_norm` does not return ``numpy.linalg.det(Q)``
directly.  Four things happen instead:

1. Any non-finite entry of :math:`\bm{Q}` returns ``NaN`` immediately — that
   is an overflow in the propagator.
2. Every row of :math:`\bm{Q}` is normalised to unit norm, which removes the
   arbitrary scale of :math:`\bm{M}(d)` relative to :math:`\bm{Q}_0`.
3. A condition number above :math:`10^{14}` for the normalised matrix
   returns ``NaN``.
4. The result is
   :math:`\mathrm{sign}\cdot\min(e^{\log|\det|}, 10^{300})`, taken from
   ``numpy.linalg.slogdet``.

This direct evaluation fails whenever :math:`\bm{Q}` is numerically
singular: far above inception, at :math:`\lambda > 0`, and — since every
photon group is explicit — wherever a group is optically thick, below
inception as well as above.  It then returns ``NaN``, and by default
(``resolve=True``) the determinant is computed by the compound-matrix
route below instead.  With ``resolve=False`` the ``NaN`` is returned; the
root finders then treat it as "above threshold", which is only safe when
no thick photon group is present.

When Q is singular: the compound-matrix route
---------------------------------------------

The growth rate needs the value, not only the sign, far above threshold and
at :math:`\lambda > 0`, where the direct evaluation fails for every
:math:`\lambda`.  The reason is the rows :math:`\bm{\Pi}\bm{M}` at the anode.
The modes they have to keep apart differ in growth across the gap by more
than the :math:`\sim 35` e-folds that double precision resolves — with
:math:`\lambda > 0` at once, because :math:`\lambda\bm{V}^{-1}` separates the
slow ion modes by thousands of e-folds — so every column of
:math:`\bm{M}` is parallel to the dominant mode to machine precision.

With ``resolve=True``, :func:`incept1d.solver.inception_det` then computes
the same row-normalised determinant without forming :math:`\bm{M}`
(:func:`incept1d.solver._det_Q_compound`).  Write
:math:`\bm{Q} = [\bm{Q}_0;\ \bm{S}\bm{M}]`, with :math:`\bm{Q}_0` the
cathode rows and :math:`\bm{S}` the anode row selection, let the columns of
:math:`\bm{N}_0` span the null space of :math:`\bm{Q}_0`, and
:math:`\bm{T} = [\bm{Q}_0^\intercal\,|\,\bm{N}_0]`.  Then :math:`\bm{Q}\bm{T}`
is block triangular, and

.. math::

   \det\bm{Q}\,\det\bm{T} = \det(\bm{Q}_0\bm{Q}_0^\intercal)\,
   \det(\bm{S}\bm{M}\bm{N}_0).

Only the last factor involves the propagator, and it is one of the
:math:`k\times k` minors — the Plücker coordinates — of
:math:`\bm{M}\bm{N}_0`, with :math:`k` the number of anode conditions.
Those minors obey a *linear* equation whose generator is the :math:`k`-th
additive compound of :math:`\bm{\mathcal{A}}`, so they are carried from the
cathode to the anode by :math:`\exp(\bm{\Omega}_i^{[k]})` step by step, with
a running scale factor.  The compound generator keeps every structural zero
of :math:`\bm{\mathcal{A}}`, so each minor keeps its own relative precision
however far it falls below the largest.

That last property is the reason for the compound matrices rather than the
more common re-orthonormalisation of the propagated subspace
(Godunov–Conte).  Orthonormalising mixes all components.  A component that
is transiently tiny — electrons attaching in a low-field region before an
avalanche near the other electrode — is then buried under the rounding
error of the large ones, and the avalanche amplifies that error to order
one.  Each step is split into substeps whose propagator has a condition
number below :math:`10^{8/k}`; the bound is on the condition number rather
than on the spread of the eigenvalues, because :math:`\bm{\mathcal{A}}` is
strongly non-normal, electron and ion speeds differing by orders of
magnitude.

The compound route costs :math:`10^2`–:math:`10^4` times the direct one
and grows combinatorially with the number of photon groups, so it runs
only where the direct evaluation returns ``NaN``.  This is why the default
criterion is the reflection operator (:ref:`Chap:Numerics:Riccati`), which
avoids the singular :math:`\bm{Q}` altogether at the cost of the direct
route.
