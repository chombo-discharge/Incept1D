Lambda (temporal growth rate)
================================

.. note::

   This page documents ``Lambda.py``. The file name is a Python builtin
   keyword clash risk (``lambda``) resolved by the module simply being
   called ``Lambda`` (capitalized); the physical quantity it solves for is
   the temporal growth rate :math:`\lambda` from :doc:`../model`.

Once the applied voltage exceeds the inception voltage :math:`V^*` found by
:mod:`Inception`, the discharge does not merely "exist" — it grows at some
rate. ``Lambda.py`` answers *how fast*: for a fixed geometry
(:math:`p\cdot d`, pressure, temperature, field profile), it first solves
:math:`\det\bm{Q}(\lambda=0, E/N^*) = 0` for the inception field
:math:`E/N^*` (same criterion as :mod:`Inception`), then for each voltage
:math:`V \in [V^*, F\cdot V^*]` solves :math:`\det\bm{Q}(\lambda^*, E/N) =
0` for the growth rate :math:`\lambda^* > 0`.

The root in :math:`\lambda` is found by exploiting the same ``NaN`` = "above
threshold" convention used in :mod:`Inception`: :math:`\det\bm{Q}(\lambda{=}0)`
is ``NaN`` above inception (near-singular :math:`\bm{Q}`), so the code
expands an upper bracket geometrically until it finds a finite, positive
value of :math:`\det\bm{Q}`, then brackets :math:`[0, \lambda_\mathrm{hi}]`
and refines with :func:`scipy.optimize.brentq`.

.. literalinclude:: ../../Lambda.py
   :language: python
   :pyobject: find_lambda_for_voltage

.. literalinclude:: ../../Lambda.py
   :language: python
   :pyobject: compute_lambda_curve

Command-line usage
--------------------

.. code-block:: console

   python Lambda.py <mechanism.py> [CONFIG.json ...] \
       --pd PD [--p P] [--T T] \
       [--n-voltages N] [--v-max-factor F] [--field SPEC] \
       [--dx [N_min [N_max [tol]]]] [--method midpoint|magnus2] \
       [--write-to-file FILE]

API reference
--------------

.. automodule:: Lambda
   :members:
   :undoc-members:
