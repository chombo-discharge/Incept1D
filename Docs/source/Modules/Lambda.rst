.. _Chap:Lambda:

Lambda.py — temporal growth rate
================================

.. contents:: On this page
   :local:
   :depth: 1

Once the applied voltage exceeds the inception voltage :math:`U^*`, the
discharge does not merely "exist" — it grows at some rate.  ``Lambda.py``
answers *how fast*.  For a fixed geometry (:math:`pd`, :math:`p`, :math:`T`,
field profile) it first solves :math:`\det\bm{Q}(0; E/N) = 0` for the
inception field :math:`(E/N)^*` (the same criterion as ``Inception.py``),
then for each voltage :math:`U \in [U^*, F U^*]` solves
:math:`\det\bm{Q}(\lambda; E/N) = 0` for the growth rate
:math:`\lambda^* > 0`.

The root in :math:`\lambda` is found by exploiting the ``NaN`` convention of
:func:`Inception.inception_det` (see :ref:`Chap:Numerics:Determinant`): above threshold
:math:`\det\bm{Q}(0)` is typically ``NaN`` because :math:`\bm{Q}` is nearly
singular, so the code expands an upper bracket geometrically (×10 per step
from 1 s\ :sup:`-1`) until a finite value with the opposite sign appears,
then refines :math:`[0, \lambda_\mathrm{hi}]` with Brent's method.

.. literalinclude:: ../../../Lambda.py
   :language: python
   :pyobject: find_lambda_for_voltage

.. literalinclude:: ../../../Lambda.py
   :language: python
   :pyobject: compute_lambda_curve

Command-line interface
----------------------

.. code-block:: console

   python3 Lambda.py MECHANISM.py [CONFIG.json ...]
       --pd PD [--p P] [--T T] [--n-voltages N] [--v-max-factor F]
       [--field SPEC] [--dx [N_min [N_max [tol]]]] [--method {midpoint,magnus2}]
       [--no-plot] [--write-to-file FILE]

``--pd`` is the :math:`pd` product in bar·mm (required); ``--n-voltages``
(default 20) and ``--v-max-factor`` (default 2.0) define the log-spaced
voltage sweep from :math:`U^*` to :math:`F U^*`.  The remaining options are
shared with ``Inception.py``.  The output table lists, for each voltage, the
over-voltage ratio :math:`U/U^*`, :math:`E/N`, :math:`\lambda` and the
corresponding e-folding time :math:`1/\lambda`.

Example
-------

.. code-block:: bash

   python3 Lambda.py Air/Air_Pancheshnyi.py --pd 10 --p 1 --n-voltages 30 --v-max-factor 1.5

API reference
-------------

.. automodule:: Lambda
   :members:
   :undoc-members:
