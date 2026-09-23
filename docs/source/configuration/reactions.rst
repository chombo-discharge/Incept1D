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
       "A- + M -> e + A + M": 0.5,
       "A- + B -> A + B-":    0.0
     }
   }

The multiplier scales the rate callable of that reaction wherever it is
used — in :math:`\bm{R}`, and hence in the eigenvalues, the ionization
integrals, and the inception determinant alike.  A multiplier of ``0``
removes the reaction without deleting it, which is how a "with / without"
comparison is set up.

One thing it does **not** do: a mechanism's optional ``alpha`` and ``eta``
helpers are computed directly from its swarm data, not from
:math:`\bm{R}`, so a multiplier does not move them.  A quantity derived from
those two — the ionization integral, the streamer criterion — is therefore
unaffected by a multiplier, while the inception criterion is.

A key that matches no reaction is reported on stderr rather than ignored
silently, so a typo is visible.  To see which strings a mechanism actually
defines, read its ``REACTIONS`` list, or run ``incept1d pdiv`` once with
``--write-to-file`` and inspect the header.

2. Changing a rate coefficient
------------------------------

Each entry of ``REACTIONS`` is a pair: the reaction string, and a callable
``rate(EN, p, T)`` returning a **first-order** rate in s\ :sup:`-1`.  The
convention that makes this work is that the callable folds in the densities
of the neutral partners it consumes:

.. code-block:: python

   # A rate coefficient in m^3/s becomes a first-order rate in s^-1 by
   # multiplying in the neutral density it acts on.
   ("e + M -> 2e + M+", lambda EN, p, T: k_ion(EN) * x_M * N(p, T)),

   # A three-body process folds in the density twice.
   ("e + 2M -> M- + M",  lambda EN, p, T: k_att(EN, T) * (x_M * N(p, T)) ** 2),

Keeping the coefficient itself in a small named function, and the density
bookkeeping in the ``REACTIONS`` entry, is what makes a rate easy to change
later: edit one function, and every use of it follows.  Keep its docstring's
units and literature reference in step with the expression.

If a rate depends on the electron temperature rather than on :math:`E/N`
directly, derive it from the mechanism's own mean-energy data rather than
re-deriving it at each call site.

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

In rough order of how quickly they fail, for a mechanism at ``MECH``:

* ``python3 MECH`` — If the mechanism has a standalone entry point, it plots
  the entries of :math:`\bm{R}` against :math:`E/N`, which catches a sign or
  magnitude error immediately.
* ``incept1d chombo MECH`` — Prints and plots every rate coefficient over an
  :math:`E/N` grid, so a rate that has gone to zero or to infinity is
  obvious.
* ``incept1d eigenvalues MECH`` — Shows whether the leading eigenvalue still
  crosses zero where you expect it.
* ``incept1d pdiv MECH CONFIG.json`` — Run with a configuration that reduces
  the chemistry to a limit with a known answer.  If only detachment or
  conversion chemistry was touched, that limit must be unchanged, which makes
  this a regression check rather than just a sanity check.
