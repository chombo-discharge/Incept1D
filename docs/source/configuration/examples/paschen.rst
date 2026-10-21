.. _Chap:PaschenMechanism:

Example: The classical Townsend model
=====================================

.. contents:: On this page
   :local:
   :depth: 1

``mechanisms/paschen/`` is the smallest mechanism that does anything useful,
and the easiest place to see a whole mechanism file at once.  It implements
the textbook Townsend model — one ionizing reaction, no attachment, no
detachment, no photons, ion-induced emission at the cathode — and its
configurations give the coefficients of helium, argon and air.

.. code-block:: text

   mechanisms/paschen/
   ├── paschen.json   the mechanism file
   └── gases.json     the three gases as three configurations

There are no data files: the coefficients are formulas, not tables.

The mechanism file
------------------

Two species and one reaction,

.. math::

   \vec{n} = \left(n_\mathrm{e},\; n_{\mathrm{M}^+}\right)^\intercal,
   \qquad
   \mathrm{e} \longrightarrow 2\mathrm{e} + \mathrm{M}^+ ,

with the ionization coefficient in Townsend's two-parameter form

.. math::
   :label: eq_townsend_two_parameter

   \frac{\alpha}{p} = A\,\exp\left(-\frac{B p}{E}\right),

and the whole file is:

.. literalinclude:: ../../../../mechanisms/paschen/paschen.json
   :language: json

Every part of :ref:`Chap:MechanismFiles` that this model needs, and
nothing else:

* **Parameters.**  :math:`A` and :math:`B` in the units they are tabulated
  in, the cathode yield and the two mobilities; these are what a
  configuration can set.
* **Definitions.**  :math:`A` and :math:`B` converted to SI, and the
  electron drift speed.
* **The reaction.**  Its rate is :math:`\alpha v_\mathrm{e}`, in s\ :sup:`-1`,
  since the string names no neutral partner.  Writing it through ``alpha``
  makes the Townsend coefficient and the reaction one and the same, which is
  what makes the closed form below exact.
* **No photoionization block**, so there are no photon groups, and
  ``eta = 0``.

The configurations
------------------

:math:`A` and :math:`B` are tabulated per gas [Raizer1991]_, and each gas
is one configuration:

.. literalinclude:: ../../../../mechanisms/paschen/gases.json
   :language: json

.. list-table::
   :header-rows: 1
   :widths: 20 20 22 38

   * - Gas
     - :math:`A` (cm⁻¹ Torr⁻¹)
     - :math:`B` (V cm⁻¹ Torr⁻¹)
     - :math:`E/p` validity (V cm⁻¹ Torr⁻¹)
   * - Helium
     - 3
     - 34
     - 20 – 150
   * - Argon
     - 12
     - 180
     - 100 – 600
   * - Air
     - 15
     - 365
     - 100 – 800

The cathode yield :math:`\gamma_0` is not a property of the gas but of the
cathode, so each configuration sets it explicitly; 0.01 is representative of
a metal.

Why it is worth having
----------------------

Because the answer is known.  With :math:`\eta = \delta = 0` the criterion
:eq:`eq_standard_paschen` collapses to
:math:`\alpha d = \ln(1 + \gamma^{-1})`, and substituting
:eq:`eq_townsend_two_parameter` gives Paschen's law in closed form together
with its minimum.

That is what :ref:`Chap:Examples:Paschen` compares against, and what the
test suite asserts to a relative error below :math:`10^{-9}`.  When you
write a mechanism for a new gas, reproducing a limit you can compute by
hand is the fastest way to find out whether it says what you meant — and
this file is a working template for doing so.
