.. _Chap:ModifyingReactions:

Modifying reactions
===================

.. contents:: On this page
   :local:
   :depth: 1

There are three levels at which the chemistry of an existing mechanism can
be changed, from least to most invasive.

1. Scaling or disabling reactions (no code changes)
---------------------------------------------------

Use ``reaction_multipliers`` in a JSON configuration
(:ref:`Chap:Configuration`).  The keys are the reaction strings exactly as
they appear in the mechanism's ``REACTIONS`` list (spacing is ignored):

.. code-block:: json

   {
     "label": "Half detachment",
     "reaction_multipliers": {
       "O2- + O2 -> e + O2 + O2": 0.5,
       "O- + N2 -> e + N2O":      0.5
     }
   }

The multiplier scales the rate callable of that reaction wherever it is
used — in :math:`\bm{R}`, and hence in the eigenvalues, the ionization
integrals, and the inception determinant alike.  It does **not** affect
:math:`\alpha` and :math:`\eta` reported by the mechanism's ``alpha``/``eta``
helper functions, which are computed directly from the BOLSIG+ tables.

The reaction strings of ``Air_Pancheshnyi.py`` are:

.. code-block:: text

   e + N2 -> 2e + N2+
   e + O2 -> 2e + O2+
   e + O2 -> O- + O
   e + 2O2 -> O2- + O2
   O2- + O2 -> e + O2 + O2
   O- + N2 -> e + N2O
   O- + O2 -> O + O2-
   O- + O2 + M -> O3- + M

Run ``incept1d pdiv`` once with ``--write-to-file`` and inspect the header
if you are unsure which mechanism and configuration are active.

2. Changing a rate coefficient
------------------------------

Each reaction in the mechanism file is a pair of a string and a rate
callable, and the callables delegate to small named functions
``k1(EN)`` … ``k8(EN)``:

.. literalinclude:: ../../../mechanisms/Air/Air_Pancheshnyi.py
   :language: python
   :pyobject: k5

To change a rate, edit the corresponding ``k`` function.  The callable in
``REACTIONS`` multiplies it by the neutral densities it consumes, so a
two-body rate coefficient in m\ :sup:`3`/s becomes a first-order rate in
s\ :sup:`-1`:

.. literalinclude:: ../../../mechanisms/Air/Air_Pancheshnyi.py
   :language: python
   :start-at: REACTIONS = [
   :end-at: ]

Keep the docstring's units and reference in sync with the new expression.
If the rate depends on the electron temperature, use the mechanism's
``ElectronTemperature(EN)`` helper (from the BOLSIG+ mean energy) as ``k4``
does.

3. Adding or removing a reaction
--------------------------------

Add a ``(string, callable)`` pair to ``REACTIONS``.  The rules of
:ref:`Chap:Reactions` apply: exactly one *tracked* species on the left,
neutrals anywhere, ``" + "``-separated tokens, optional integer
coefficients.  Because :math:`\bm{R}` is assembled automatically from the
string, no other code changes are required as long as every species you
mention is already in ``SPECIES``.

If the reaction produces or consumes a **new** species, you are extending
the species set — see :ref:`Chap:NewMechanisms`, since ``get_V``, the
row-selection matrices and (if the species is a positive ion) the SEE
yield vector must all be extended consistently.

Removing a reaction is a matter of deleting its entry (or, for a temporary
study, setting its multiplier to 0).

Sanity checks after a change
----------------------------

* ``python3 mechanisms/Air/Air_Pancheshnyi.py`` runs the mechanism standalone and
  plots the entries of :math:`\bm{R}` vs. :math:`E/N`.
* ``incept1d eigenvalues mechanisms/Air/Air_Pancheshnyi.py`` shows whether the
  leading eigenvalue still crosses zero where you expect it.
* ``incept1d chombo mechanisms/Air/Air_Pancheshnyi.py`` prints
  every rate coefficient over the :math:`E/N` grid and plots them.
* ``incept1d pdiv mechanisms/Air/Air_Pancheshnyi.py mechanisms/Air/Paschen.json`` must
  still reproduce the classical Paschen curve — a good regression check
  when only detachment/conversion chemistry was touched.
