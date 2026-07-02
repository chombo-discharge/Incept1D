Field distributions
====================

``FieldDistributions.py`` abstracts the gap electrode geometry away from the
solvers: everything downstream (:mod:`Inception`, :mod:`IonizationIntegral`,
:mod:`Lambda`) only ever needs a normalized field profile
:math:`f(\xi)`, :math:`\xi \in [0, 1]`, satisfying :math:`\int_0^1 f(\xi)\,
d\xi = 1`, where :math:`\xi = 0` is the high-field electrode (sphere surface,
or either plate for a uniform field) and :math:`\xi = 1` is the low-field
side (plane, or the far sphere).

:class:`FieldDistribution` is a small dataclass holding the geometry type
(``'uniform'``, ``'sphere-plane'``, ``'sphere-sphere'``) and, where
applicable, the sphere radius; :meth:`FieldDistribution.build` turns that
static description into the callable :math:`f(\xi)` for a specific gap
length :math:`d`. The non-uniform field profiles are computed *exactly* via
a bispherical image-charge series (:func:`_sphere_sphere_axial_field`), not
approximated.

``add_field_argument`` / ``parse_field_spec`` implement the shared
``--field SPEC`` command-line syntax used by every CLI script in the
repository (``uniform`` | ``sphere-plane R_mm`` | ``sphere-sphere R_mm``), so
the geometry parsing only has to be written and tested once.

.. literalinclude:: ../../FieldDistributions.py
   :language: python
   :pyobject: FieldDistribution.build

API reference
--------------

.. automodule:: FieldDistributions
   :members:
   :undoc-members:
