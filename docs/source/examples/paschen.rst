.. _Chap:Examples:Paschen:

The classical Paschen curve
===========================

.. contents:: On this page
   :local:
   :depth: 1

Every other example in this chapter compares ``Incept1D`` against measurements.
This one compares it against algebra.

The mechanism ``mechanisms/paschen/paschen.py`` is the textbook Townsend model:
one ionizing reaction, no attachment, no detachment, no photoionization, and
electrons released from the cathode only by arriving positive ions.  That is
exactly the reduced model whose inception condition is
:eq:`eq_standard_paschen`, so the answer is known in closed form and the run is
a *verification* — it checks the solver rather than the physics.

The model
---------

The ionization coefficient is the two-parameter Townsend form

.. math::
   :label: eq_townsend_alpha

   \frac{\alpha}{p} = A\,\exp\!\left(-\frac{B p}{E}\right),

and with :math:`\eta = \delta = 0` the criterion
:eq:`eq_standard_paschen` collapses to :math:`\alpha d = \ln(1 + \gamma^{-1})`.
Substituting :eq:`eq_townsend_alpha` and solving for the voltage gives Paschen's
law,

.. math::
   :label: eq_paschen_curve

   U = \frac{B\,(pd)}{\ln\!\big(A\,pd\big) - \ln\!\Big(\ln\big(1 + \gamma^{-1}\big)\Big)},

with the minimum at

.. math::
   :label: eq_paschen_minimum

   (pd)_\mathrm{min} = \frac{e}{A}\ln\!\left(1 + \frac{1}{\gamma}\right),
   \qquad
   U_\mathrm{min} = \frac{e B}{A}\ln\!\left(1 + \frac{1}{\gamma}\right).

The denominator of :eq:`eq_paschen_curve` vanishes at
:math:`pd = \ln(1+\gamma^{-1})/A`; below that the gap cannot break down at any
voltage, which is the vertical asymptote on the left of every Paschen curve.

Coefficients
------------

:math:`A` and :math:`B` are the standard values tabulated by [Raizer1991]_,
converted to SI inside the mechanism.  They are fits, valid over the stated
:math:`E/p` window, not laws:

.. list-table::
   :header-rows: 1
   :widths: 18 16 18 30 18

   * - Gas
     - :math:`A` (cm⁻¹ Torr⁻¹)
     - :math:`B` (V cm⁻¹ Torr⁻¹)
     - :math:`E/p` validity (V cm⁻¹ Torr⁻¹)
     - :math:`U_\mathrm{min}` at :math:`\gamma = 0.01`
   * - Helium
     - 3
     - 34
     - 20 – 150
     - 142 V at 0.056 bar·mm
   * - Argon
     - 12
     - 180
     - 100 – 600
     - 188 V at 0.014 bar·mm
   * - Air
     - 15
     - 365
     - 100 – 800
     - 305 V at 0.011 bar·mm

:math:`\gamma` is not a property of the gas — it depends on the cathode
material and its surface condition — so it is a configuration parameter rather
than a constant of the mechanism.  The value 0.01 used here is representative
of a metal cathode; changing it moves the whole curve, as
:eq:`eq_paschen_minimum` shows.

.. note::

   This is deliberately *not* a quantitative model of a real discharge.  Air is
   electronegative, and with attachment switched off the left branch and the
   high-\ :math:`pd` behaviour are both wrong.  For real air use
   ``mechanisms/air/`` (:ref:`Chap:AirScheme`); the difference between the two
   is the subject of :ref:`Chap:Examples:Electra`.

Running the calculation
-----------------------

The three gases are three configurations of one mechanism, so a single
invocation produces all three curves:

.. literalinclude:: ../../../examples/paschen/run.sh
   :language: bash
   :start-at: set -euo pipefail

``closed_form.py`` then tabulates :eq:`eq_paschen_curve` from the very same
coefficients, so the comparison is self-contained — no external table is
needed, and the reference cannot drift out of step with the mechanism.

Result
------

.. _fig_paschen:

.. figure:: ../figures/paschen.*
   :width: 100%
   :align: center

   Inception voltage of helium, argon and air in a uniform gap at 1 bar.
   Lines are the ``Incept1D`` solution of :math:`\det\bm{Q}(0) = 0`; circles
   are the closed form :eq:`eq_paschen_curve`.  The solver curves stop at the
   left-branch asymptote, where no solution exists at any voltage.

The two agree to solver tolerance — in the test suite the same comparison is
asserted to a relative error below :math:`10^{-9}` at every :math:`pd`.  That
is a stronger statement than it looks: the solver reaches this answer through
the full machinery — assembling :math:`\bm{\mathcal{A}}`, propagating it across
the gap, building the boundary matrix and locating the root of its determinant
— with no shortcut for the analytic case.  Agreement therefore exercises the
whole chain against a result that can be checked by hand.

Three features of the curve are worth naming, because they recur in the
realistic cases:

* **The right branch** rises because a denser gas means more collisions but
  less energy gained between them, so a higher field is needed for each
  electron to ionize.
* **The minimum** is where the two effects balance.  It depends on the cathode
  through :math:`\gamma` only logarithmically, which is why the Paschen
  minimum of a given gas is a robust number.
* **The left branch** rises steeply and then terminates: with too few
  collisions available, no voltage can sustain the avalanche chain.  This is
  the same asymptote that makes the dry-air curve in
  :ref:`Chap:Examples:Electra` turn upward, although there attachment moves it.

Varying the cathode
-------------------

Because :math:`\gamma` is a configuration key, a sensitivity study is one file:

.. code-block:: json

   {
     "configurations": [
       {"label": "γ = 0.001", "gas": "air", "gamma0": 0.001},
       {"label": "γ = 0.01",  "gas": "air", "gamma0": 0.01},
       {"label": "γ = 0.1",   "gas": "air", "gamma0": 0.1}
     ]
   }

:eq:`eq_paschen_minimum` predicts the outcome exactly — both
:math:`U_\mathrm{min}` and :math:`(pd)_\mathrm{min}` scale as
:math:`\ln(1 + \gamma^{-1})` — which makes this a convenient check when
adapting the mechanism to a new gas.
