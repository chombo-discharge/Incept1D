.. _Chap:PaschenMechanism:

Example: the classical Townsend model
=====================================

.. contents:: On this page
   :local:
   :depth: 1

``mechanisms/paschen/`` is the smallest mechanism that does anything useful,
and the easiest place to see the whole interface at once.  It implements the
textbook Townsend model — one ionizing reaction, no attachment, no
detachment, no photons, ion-induced emission at the cathode — for helium,
argon and air.

Everything it needs fits in three files:

.. code-block:: text

   mechanisms/paschen/
   ├── paschen.py     the mechanism module
   ├── config.py      selects the gas and the cathode yield
   └── gases.json     the three gases as three configurations

There are no data files at all: the coefficients are functions, not tables.

What the module defines
-----------------------

Two species and one reaction:

.. math::

   \vec{n} = \left(n_\mathrm{e},\; n_{\mathrm{M}^+}\right)^\intercal,
   \qquad
   \mathrm{e} \longrightarrow 2\mathrm{e} + \mathrm{M}^+ ,

with the ionization coefficient in Townsend's two-parameter form

.. math::
   :label: eq_townsend_two_parameter

   \frac{\alpha}{p} = A\,\exp\left(-\frac{B p}{E}\right).

:math:`A` and :math:`B` are tabulated per gas [Raizer1991]_ and converted to
SI inside the module:

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

The attachment coefficient returns zero, and the photon accessors return
zero-width arrays — which is how a mechanism declares that it has no
photoionization:

.. code-block:: python

   def get_B(EN, p=1.0, T=293.0):
       """Photon absorption coupling, shape (2, 0)."""
       return np.zeros((2, 0))

Choosing between the two override paths
---------------------------------------

This mechanism is a clean illustration of the rule in :ref:`Chap:ConfigPy`.

The **gas** must be known before the module body runs, because selecting it
fixes :math:`A` and :math:`B` at module level.  It therefore goes through
``pre_exec_vars``:

.. literalinclude:: ../../../mechanisms/paschen/config.py
   :language: python
   :pyobject: Config.pre_exec_vars

The cathode yield :math:`\gamma_0` travels the same way here, because the
module reads it at import; a mechanism that evaluated it per call could
equally have passed it through ``mechanism_params``.

The configurations themselves are three lines each:

.. literalinclude:: ../../../mechanisms/paschen/gases.json
   :language: json
   :lines: 1-10

Why it is worth having
----------------------

Because the answer is known.  With :math:`\eta = \delta = 0` the criterion
:eq:`eq_standard_paschen` collapses to
:math:`\alpha d = \ln(1 + \gamma^{-1})`, and substituting
:eq:`eq_townsend_two_parameter` gives Paschen's law in closed form together
with its minimum.  The module exposes both as ordinary functions,
``paschen_voltage(pd)`` and ``paschen_minimum()``, so a calculation can be
compared against them directly.

That is what :ref:`Chap:Examples:Paschen` does, and what the test suite
asserts to a relative error below :math:`10^{-9}`.  When you write a
mechanism for a new gas, reproducing a limit you can compute by hand is the
fastest way to find out whether the interface has been implemented
correctly — and this mechanism is a working template for doing so.
