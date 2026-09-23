.. _Chap:Lambda:

Temporal growth rate — ``incept1d growth``
==========================================

.. contents:: On this page
   :local:
   :depth: 1

Once the applied voltage exceeds the inception voltage :math:`U^*`, the
discharge does not merely "exist" — it grows at some rate.  ``incept1d growth``
answers *how fast*.  For a fixed geometry (:math:`pd`, :math:`p`, :math:`T`,
field profile) it first solves :math:`\det\bm{Q}(0; E/N) = 0` for the
inception field :math:`(E/N)^*` (the same criterion as ``incept1d pdiv``),
then for each voltage :math:`U \in [U^*, F U^*]` solves
:math:`\det\bm{Q}(\lambda; E/N) = 0` for the growth rate
:math:`\lambda^* > 0`.

The root in :math:`\lambda` is found by exploiting the ``NaN`` convention of
:func:`incept1d.solver.inception_det` (see :ref:`Chap:Numerics:Determinant`): above threshold
:math:`\det\bm{Q}(0)` is typically ``NaN`` because :math:`\bm{Q}` is nearly
singular, so the code expands an upper bracket geometrically (×10 per step
from 1 s\ :sup:`-1`) until a finite value with the opposite sign appears,
then refines :math:`[0, \lambda_\mathrm{hi}]` with Brent's method.  The
resulting root is accepted only if :math:`\det\bm{Q}(\lambda^*)` is finite —
see :ref:`Sec:Lambda:Range` for why that check matters.

.. _Sec:Lambda:Range:

Range of validity
-----------------

.. warning::

   ``incept1d growth`` only answers over a limited band of overvoltage.  For
   dry air at :math:`pd = 10` bar·mm that band is roughly
   :math:`U \lesssim 1.2\,U^*`; beyond it the command reports
   ``det_Q_unresolved`` and returns no growth rate.

The limit is not in the root finder but in the determinant it is given.  The
propagator is formed as :math:`\bm{M} = \exp(\bm{\mathcal{A}}d)` after
shifting by the largest real eigenvalue
(:func:`incept1d.solver._expm_shifted`), so the subdominant modes appear as
:math:`e^{(\mu_i - \mu_\mathrm{max})d}`.  Once the spectral spread
:math:`(\mu_\mathrm{max} - \mu_\mathrm{min})d` exceeds the
double-precision underflow limit of about 709, those modes become *exactly*
zero: at :math:`2U^*` in the case above the spread is 830, and
:math:`\bm{M}` has numerical rank 2 out of 11.  The anode rows
:math:`\bm{\Pi}_+\bm{M}` then carry no independent information,
:math:`\bm{Q}` is rank deficient, and
:func:`incept1d.solver._assemble_det_Q` returns ``NaN`` from its
conditioning guard — for *every* :math:`\lambda`, not only near a root.

This matters because ``NaN`` is mapped to a negative sentinel inside the
bracketing function, so a run of ``NaN`` below a finite positive value is
indistinguishable from a sign change.  Brent's method then converges on the
edge of the ``NaN`` region.  That edge is the :math:`\lambda` at which a
photon group crosses :math:`\kappa d = 12` and collapses into
:math:`\bm{A}` (:ref:`Chap:Numerics:Determinant`), which does not depend on
:math:`E/N` — so before the acceptance check was added, the same
:math:`\lambda \approx 1.9\times10^{11}` s\ :sup:`-1` was reported at every
overvoltage.

Lifting the restriction means not forming :math:`\bm{M}` at all.  The
standard remedy for a linear boundary-value problem with this dynamic range
is to propagate the *subspace* that satisfies the cathode conditions with
periodic re-orthonormalisation (Godunov–Conte shooting), or to propagate the
corresponding exterior product directly (the compound-matrix or Evans-function
method), so that no individual mode is ever represented.  Until then, treat a
``det_Q_unresolved`` row as "not answered" rather than "no growth".

Inputs
------

.. code-block:: console

   incept1d growth MECHANISM.py [CONFIG.json ...]
       --pd PD [--p P] [--T T] [--n-voltages N] [--v-max-factor F]
       [--field SPEC] [--dx [N_min [N_max [tol]]]] [--method {midpoint,magnus2}]
       [--no-plot] [--write-to-file FILE]

``--pd`` is the :math:`pd` product in bar·mm (required); ``--n-voltages``
(default 20) and ``--v-max-factor`` (default 2.0) define the log-spaced
voltage sweep from :math:`U^*` to :math:`F U^*`.  The remaining options are
shared with ``incept1d pdiv``.  The output table lists, for each voltage, the
over-voltage ratio :math:`U/U^*`, :math:`E/N`, :math:`\lambda` and the
corresponding e-folding time :math:`1/\lambda`.

Example
-------

.. code-block:: bash

   incept1d growth mechanisms/air/air_pancheshnyi.py --pd 10 --p 1 --n-voltages 30 --v-max-factor 1.5

Outputs
-------

The inception point is reported first (:math:`U^*` and :math:`(E/N)^*`),
then a table with one row per voltage:

.. code-block:: text

   V (kV)    V/V*    EN (Td)    λ (s⁻¹)    τ (ns)    λ/ν_ion

:math:`\tau = 1/\lambda` is the e-folding time, and the last column measures
the growth rate against the ionization frequency, which says whether the
discharge grows on the avalanche timescale or far more slowly.  The figure
has two panels, :math:`\lambda` and :math:`\tau` against :math:`U/U^*`.

A row that cannot be resolved is reported as ``NaN`` and flagged with a
status — ``det_Q_unresolved`` means the determinant could not be evaluated
at any :math:`\lambda`, not that the gap does not grow; see *Range of
validity* above.

``--write-to-file`` writes one row per voltage with the columns ``V_kV``
and ``V_ratio``, then ``lambda_s-1[LABEL]``, ``tau_ns[LABEL]`` and
``nu_ion_s-1[LABEL]`` for every configuration side by side.

API reference
-------------

.. automodule:: incept1d.growth
   :members:
   :undoc-members:
