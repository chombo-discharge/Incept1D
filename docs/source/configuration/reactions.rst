.. _Chap:ModifyingReactions:

Modifying reactions
===================

.. contents:: On this page
   :local:
   :depth: 1

There are three levels at which the chemistry of an existing mechanism can
be changed, from least to most invasive.

1. Scaling or disabling reactions
---------------------------------

Use ``reaction_multipliers`` in a configuration (:ref:`Chap:Configuration`).
The keys are the reaction strings exactly as the mechanism writes them
(spacing is ignored):

.. code-block:: json

   {
     "label": "Half detachment",
     "reaction_multipliers": {
       "B- + M -> e + B + M": 0.5,
       "B- + A -> B + A-":    0.0
     }
   }

The multiplier scales the rate of that reaction in :math:`\bm{R}`, and hence
in the eigenvalues and the inception criterion alike.  A multiplier of ``0``
removes the reaction without deleting it, which is how a "with / without"
comparison is set up.  A key that names no reaction of the mechanism is an
error, so a typo cannot pass for a reaction that was switched off.

Two things a multiplier does **not** change:

* **The Townsend coefficients.**  A mechanism gives ``alpha`` and ``eta``
  itself rather than deriving them from the reactions
  (:ref:`Chap:MechanismFiles`), so the ionization integral and the streamer
  criterion, which use them, are unaffected.
* **The photon source.**  Photon emission is proportional to the rates of the
  source reactions before multipliers, so switching an ionization off does
  not switch off the photons it would have emitted.  Scale the source with a
  parameter of the photoionization block instead.

2. Changing a rate coefficient
------------------------------

Edit the ``rate`` formula of the reaction in the mechanism file.  The rate is
the coefficient :math:`k` in SI units, and the densities of the neutral
partners the reaction string names are multiplied in by the loader:

.. code-block:: json

   "reactions": [
       {"reaction": "e + A -> 2e + A+",   "rate": "k_ion"},
       {"reaction": "e + 2B -> B- + B",   "rate": "1.4e-41*(300/Te)"}
   ]

so the first is in m\ :sup:`3`/s and is multiplied by :math:`x_A N`, the
second in m\ :sup:`6`/s and multiplied by :math:`(x_B N)^2`.  A coefficient
used in more than one place belongs in ``definitions``, and one a study will
vary belongs in ``parameters``, so that a configuration can set it.  Write
its units and literature reference next to it as a ``_`` comment.

3. Adding or removing a reaction
--------------------------------

Add an entry to ``reactions``.  Exactly one *tracked* species on the left,
neutral partners that are gas components or ``M``, ``" + "``-separated
tokens, optional integer coefficients; products that are not followed are
listed in ``untracked``.  :math:`\bm{R}` is assembled from the string, so
nothing else has to change as long as every species it mentions is already
declared.  A new species is one more entry in ``species``, and, if it is
positive, one more yield in ``emission``; the loader says so if either is
missing.

Removing a reaction is a matter of deleting its entry, or, for a temporary
study, setting its multiplier to 0.  Configurations that name the removed
reaction in ``reaction_multipliers`` are then refused, which finds them.

Sanity checks after a change
----------------------------

In rough order of how quickly they fail, for a mechanism at ``MECH``:

* ``incept1d eigenvalues MECH`` — Loads and checks the file first, so a
  reaction that does not balance, an unknown name or a missing yield is
  reported with its place in the file; then shows whether the leading
  eigenvalue still crosses zero where you expect it.
* ``incept1d pdiv MECH CONFIG.json`` — Run with a configuration that reduces
  the chemistry to a limit with a known answer.  If only detachment or
  conversion chemistry was touched, that limit must be unchanged, which makes
  this a regression check rather than just a sanity check.

For a Python mechanism (:ref:`Chap:NewMechanisms`), the same applies to its
``REACTIONS`` list, where each entry is a reaction string and a callable
returning the rate in s\ :sup:`-1` with the neutral densities already folded
in.
