Ionization integral
=====================

``IonizationIntegral.py`` compares the full boundary-value inception
criterion against two classical, purely-local surrogates, evaluated along
the (possibly non-uniform) field profile and plotted vs. applied voltage:

- the classical ionization integral
  :math:`\int_0^d \max(\alpha(x) - \eta(x), 0)\,dx`, i.e. the textbook
  Townsend criterion generalized to a non-uniform field, ignoring
  detachment, ion transit, and photoionization entirely;
- :math:`\int_0^d \max(\mathrm{Re}\,\lambda_\max(\bm{R}\bm{V}^{-1}), 0)\,dx`,
  the same integral but using the leading eigenvalue of the full local
  transport matrix (see :mod:`Eigenvalues`) in place of :math:`\alpha -
  \eta` — only available when the mechanism exposes ``get_R``/``get_V``.

Comparing these two curves against the full :func:`Inception.inception_det`
result quantifies how much the *local* approximation (either form) misses
once negative-ion transit and detachment matter — exactly the "sensitivity
to cross sections" / "comparison with canonical data sources" analysis in
the manuscript.

.. literalinclude:: ../../IonizationIntegral.py
   :language: python
   :pyobject: aed_integral

.. literalinclude:: ../../IonizationIntegral.py
   :language: python
   :pyobject: eig_integral

Command-line usage
--------------------

.. code-block:: console

   python IonizationIntegral.py <mechanism.py> [CONFIG.json ...] \
       --p P [P ...] --d D [D ...] \
       [--voltage-lo V] [--voltage-hi V] [--voltage-num N] \
       [--single-voltage V] [--data-file FILE] \
       [--field SPEC] [--write-to-file FILE]

``--data-file`` reads experimental (pressure, voltage) pairs from an ASCII
table and overlays the corresponding integral values, instead of sweeping a
synthetic voltage range — this is how the model-vs-experiment comparison
figures in the manuscript were generated.

API reference
--------------

.. automodule:: IonizationIntegral
   :members:
   :undoc-members:
