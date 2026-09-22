.. _Chap:Numerics:Determinant:

Evaluating the determinant
==========================

With :math:`\bm{M}(d)` in hand, :math:`\bm{Q}` of :eq:`eq_Q_system` is
assembled and its determinant evaluated.  Only the *sign* of
:math:`\det\bm{Q}` matters for root finding, and the evaluation is designed
around keeping that sign reliable.

:func:`Inception._assemble_det_Q` does not return ``numpy.linalg.det(Q)``
directly.  It

1. returns ``NaN`` if any entry of :math:`\bm{Q}` is non-finite (overflow
   in the propagator);
2. normalises every row of :math:`\bm{Q}` to unit norm, which removes the
   arbitrary scale of :math:`\bm{M}(d)` relative to :math:`\bm{Q}_0`;
3. returns ``NaN`` if the condition number of the normalised matrix
   exceeds :math:`10^{14}`;
4. returns :math:`\mathrm{sign}\cdot\min(e^{\log|\det|}, 10^{300})` from
   ``numpy.linalg.slogdet``.

The ``NaN`` convention is central to the root finders: **above** the
inception field the propagator grows so fast that :math:`\bm{Q}` becomes
numerically singular, so ``NaN`` reliably means "above threshold", where
the true sign of :math:`\det\bm{Q}` is negative.
