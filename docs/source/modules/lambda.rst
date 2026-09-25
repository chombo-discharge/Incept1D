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

:math:`\det\bm{Q}(\lambda) = 0` has many roots, one for every mode of the
gap; :math:`\lambda^*` is the one that sets the growth, as explained in
:ref:`Sec:Lambda:Range`.  How it is bracketed and refined is described in
:ref:`Chap:Numerics:RootFinding`.

.. _Sec:Lambda:Range:

Which root, and range of validity
---------------------------------

A mode :math:`\vec{n}(x)\,e^{\lambda t}` of the linearised gap exists for
every root of :math:`\det\bm{Q}(\lambda) = 0`, and above :math:`U^*` several
of them can grow at once.  The feedback loops (ions returning to the
cathode, photons) act like delays, so besides real roots there are complex
ones, oscillating modes whose real part can also be positive just above
threshold.

Every coupling in the linear model is a non-negative source — ionization,
attachment into negative ions, detachment, photoionization, secondary
emission — so the gain of each feedback loop, the number of cathode
electrons one generation produces with every path discounted by
:math:`e^{-\lambda t}`, falls as :math:`\lambda` increases.  The dominant
mode is then the unique real :math:`\lambda^*` at which the strongest loop
gain is one; its density profile is positive everywhere, and no root, real
or complex, lies above it.  :math:`\lambda^*` is therefore the **largest real
root** of :math:`\det\bm{Q}`, and that is what ``incept1d growth`` reports:
it scans down from well above any possible growth rate to the first sign
change, rather than up from zero, which could stop at a slower mode.

Above threshold the propagator :math:`\bm{M}(d)` cannot be formed reliably:
the modes it has to keep separate differ in growth by hundreds to thousands
of e-folds, far beyond double precision, and the anode rows of
:math:`\bm{Q}` become parallel.  ``incept1d growth`` evaluates
:math:`\det\bm{Q}` by the compound-matrix method there
(:ref:`Chap:Numerics:Determinant`), so it resolves the growth rate at any
overvoltage.

The limit that remains is physical.  The model is linear: it has no space
charge.  At :math:`1.5\,U^*` an avalanche across a centimetre-scale gap can
multiply by :math:`e^{400}` and more, far beyond the :math:`e^{18}`–:math:`e^{20}`
at which space charge takes over, so :math:`\lambda^*` is a meaningful growth
rate only close to :math:`U^*`, and the further above it the more it is a
property of the equations rather than of the discharge.

Inputs
------

.. code-block:: console

   incept1d growth MECHANISM.py [CONFIG.json ...]
       [--pressure P] [--distance D] [--T T] [--n-voltages N] [--v-max-factor F]
       [--field SPEC] [--dx [N_min [N_max [tol]]]] [--method {midpoint,magnus2}]
       [--criterion {riccati,detq}] [--jobs N] [--no-plot] [--write-to-file FILE]

``--pressure`` is the gas pressure in bar (default 1) and ``--distance`` the
gap length in mm.  ``--distance`` is required unless the geometry fixes the
gap itself: a coaxial gap is :math:`b - a` long and a field line its arc
length, and a ``--distance`` given with those is ignored, with a message
saying so.  ``--n-voltages`` (default 20) and ``--v-max-factor`` (default
2.0) define the log-spaced voltage sweep from :math:`U^*` to :math:`F U^*`.
The voltages are independent and are solved concurrently on ``--jobs``
worker processes (default: the number of physical cores); each is printed
as soon as it is solved, and the table at the end is in voltage order.  The
remaining options are shared with ``incept1d pdiv``.  The output table lists, for each voltage, the
over-voltage ratio :math:`U/U^*`, :math:`E/N`, :math:`\lambda` and the
corresponding e-folding time :math:`1/\lambda`.

Example
-------

.. code-block:: bash

   incept1d growth mechanisms/air/pancheshnyi/air_pancheshnyi.py --pressure 1 --distance 10 --n-voltages 30 --v-max-factor 1.5

Outputs
-------

Each configuration is solved once per polarity, with its own inception
point, exactly as ``incept1d pdiv`` does: the two are solved separately
unless the gap is symmetric and the configuration has no per-polarity
overrides, in which case the negative-polarity curve repeats the positive
one.  For each curve the inception point is reported first (:math:`U^*` and
:math:`(E/N)^*`), then a table with one row per voltage:

.. code-block:: text

   V (kV)    V/V*    EN (Td)    λ (s⁻¹)    τ (ns)    λ/ν_ion

:math:`\tau = 1/\lambda` is the e-folding time, and the last column measures
the growth rate against the ionization frequency, which says whether the
discharge grows on the avalanche timescale or far more slowly.  The figure
has two panels, :math:`\lambda` and :math:`\tau` against the voltage in kV,
with a solid line for positive and a dashed line for negative polarity.

With ``--write-to-file`` all curves share the first column, :math:`U/U^*`.
Because :math:`U^*` differs between configurations and polarities, every
curve then has its own four columns: the voltage in kV, :math:`\lambda`,
:math:`\tau` and :math:`\nu_\mathrm{ion}`, named after the curve (for
example ``lambda_s-1[Baseline (negative)]``) in the ``# Column`` lines of the
header.

:math:`\nu_\mathrm{ion}` is the fastest local net ionization rate in the gap,
evaluated at the peak of the field profile; in a non-uniform gap the
gap-averaged field can lie below the ionization threshold.

A row that is not plainly resolved carries a status.  ``suspect`` means
:math:`\det\bm{Q}` does not change sign across :math:`\lambda^*` — a jump,
not a zero, for example where a photon group crosses the optically thick
threshold; ``det_Q_unresolved`` and ``no_bracket_found`` mean no growth rate
could be determined, not that the gap does not grow.

API reference
-------------

.. automodule:: incept1d.growth
   :members:
   :undoc-members:
