.. _Chap:Lambda:

Temporal growth rate — ``incept1d growth``
==========================================

.. contents:: On this page
   :local:
   :depth: 1

Once the applied voltage exceeds the inception voltage :math:`U^*`, the
discharge does not merely "exist" — it grows at some rate.  ``incept1d growth``
answers *how fast*.  For a geometry (field profile, pressure, gap length,
temperature) it first finds the inception field :math:`(E/N)^*` — the same
criterion as ``incept1d pdiv``, at :math:`\lambda = 0` — and then, for each
voltage :math:`U \in [U^*, F U^*]`, solves the criterion at
:math:`\lambda > 0` for the growth rate :math:`\lambda^*`.  Several
pressures and gap lengths can be given at once; every combination is a
case of its own.

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
:math:`\bm{Q}` become parallel.  The default criterion never forms
:math:`\bm{M}` (:ref:`Chap:Numerics:Riccati`), and with it the growth rate
is resolved at any overvoltage; with ``--criterion detq`` the determinant
falls back on the much slower compound-matrix route
(:ref:`Chap:Numerics:Determinant`).  An independent discretisation of the
time-dependent equations, in the test suite, confirms that
:math:`\lambda^*` found this way is the fastest mode, complex modes
included.

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
       [--pressure P [P ...]] [--distance D [D ...]] [--T T]
       [--n-voltages N] [--v-max-factor F]
       [--field SPEC] [--dx [N_min [N_max [tol]]]] [--method {midpoint,magnus2}]
       [--criterion {riccati,detq}] [--jobs N] [--verify] [--silent]
       [--no-plot] [--write-to-file FILE]

``--pressure`` is the gas pressure in bar (default 1) and ``--distance`` the
gap length in mm.  Each takes one or more values, and ``MIN:MAX:N`` expands
to :math:`N` log-spaced values from MIN to MAX; every pressure is combined
with every distance, so ``--pressure 1:10:3 --distance 5 10`` is six cases.
``--distance`` is required unless the geometry fixes the gap itself: a
coaxial gap is :math:`b - a` long and a field line its arc length, and a
``--distance`` given with those is ignored, with a message saying so.  ``--n-voltages`` (default 20) and ``--v-max-factor`` (default
2.0) define the log-spaced voltage sweep from :math:`U^*` to :math:`F U^*`.
The work runs in two parallel phases on ``--jobs`` worker processes
(default: the number of physical cores): first the inception voltage of
every case, configuration and polarity, then all their voltages together.
Each result is printed as soon as it is solved, and the tables at the end
are in voltage order.  ``--verify`` checks every root with
:math:`\det\bm{Q}` — across :math:`(E/N)^*` for each inception voltage and
across :math:`\lambda^*` for each growth rate — and reports ✓ or ✗ with
the progress; it can be expensive.  ``--silent`` prints only the tables,
and so also turns off that check.  The
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
one colour per configuration and case, with a solid line for positive and a
dashed line for negative polarity; a negative polarity that only mirrors
the positive one is drawn once, labelled "both polarities".  Each curve
carries at most 20 markers, the voltage axis becomes logarithmic when the
curves span more than a factor of five, and with more than four curves the
legend sits below the panels.

With ``--write-to-file`` all curves share the first column, :math:`U/U^*`.
Because :math:`U^*` differs between configurations and polarities, every
curve then has its own four columns: the voltage in kV, :math:`\lambda`,
:math:`\tau` and :math:`\nu_\mathrm{ion}`, named after the curve (for
example ``lambda_s-1[Baseline (negative)]``, or with several cases
``lambda_s-1[Baseline, p=2 bar, d=10 mm (negative)]``) in the ``# Column``
lines of the header, which also lists the pressures and distances.

:math:`\nu_\mathrm{ion}` is the fastest local net ionization rate in the gap,
evaluated at the peak of the field profile; in a non-uniform gap the
gap-averaged field can lie below the ionization threshold.

A row that is not plainly resolved carries a status.  ``suspect`` means
the criterion does not change sign across :math:`\lambda^*` — a jump, not
a zero; ``det_Q_unresolved`` and ``no_bracket_found`` mean no growth rate
could be determined, not that the gap does not grow.

API reference
-------------

.. automodule:: incept1d.growth
   :members:
   :undoc-members:
