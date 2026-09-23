.. _Chap:Numerics:Determinant:

Evaluating the determinant
==========================

With :math:`\bm{M}(d)` in hand, :math:`\bm{Q}` of :eq:`eq_Q_system` is
assembled and its determinant evaluated.  Only the *sign* of
:math:`\det\bm{Q}` matters for root finding, and the evaluation is designed
around keeping that sign reliable.

:func:`incept1d.solver._assemble_det_Q` does not return ``numpy.linalg.det(Q)``
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

The ``NaN`` convention is central to the root finders: **above** the
inception field the propagator grows so fast that :math:`\bm{Q}` becomes
numerically singular, so ``NaN`` reliably means "above threshold", where
the true sign of :math:`\det\bm{Q}` is negative.
