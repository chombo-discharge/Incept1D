.. _Chap:Reactions:

Declarative reactions — ``incept1d.reactions``
==============================================

.. contents:: On this page
   :local:
   :depth: 1

:mod:`incept1d.reactions` assembles the reaction-rate matrix :math:`\bm{R}` of
:eq:`eq_drift_reaction` from a *declarative* list of reaction strings, so a
mechanism file states its chemistry as human-readable equations,

.. code-block:: python

   REACTIONS = [
       ("e + N2 -> 2e + N2+",        lambda EN, p, T: k1(EN) * xN2 * N(p, T)),
       ("e + 2O2 -> O2- + O2",       lambda EN, p, T: k4(EN, T) * (xO2 * N(p, T))**2),
       ("O- + N2 -> e + N2O",        lambda EN, p, T: k6(EN) * xN2 * N(p, T)),
       ...
   ]

instead of hand-writing stoichiometry into a matrix.  Each string is paired
with a rate callable ``rate(EN, p, T)`` that returns the **first-order rate
in s**\ :sup:`-1` for the driving species — i.e. the rate coefficient
already multiplied by the neutral densities it consumes.

The reaction matrix
-------------------

:math:`R_{ij}` is the rate at which one particle of species :math:`j`
produces (positive) or destroys (negative) one particle of species
:math:`i`.  For a reaction with driver :math:`j`, reactant counts
:math:`r_i` and product counts :math:`p_i`, the contribution is

.. math::

   R_{ij} \mathrel{+}= \left(p_i - r_i\right)\,k_\mathrm{eff},

so ``"e + N2 -> 2e + N2+"`` adds :math:`+k` to :math:`R_{\mathrm{e},\mathrm{e}}`
(net gain of one electron) and :math:`+k` to
:math:`R_{\mathrm{N}_2^+,\mathrm{e}}`, while ``"e + O2 -> O- + O"`` adds
:math:`-k` to :math:`R_{\mathrm{e},\mathrm{e}}` and :math:`+k` to
:math:`R_{\mathrm{O}^-,\mathrm{e}}`.

Parsing rules
-------------

* Species tokens are separated by ``" + "`` (spaces required) so that
  cation names like ``N2+`` are never split.
* An optional leading integer is a stoichiometric coefficient: ``2e``,
  ``2O2``.
* The arrow is ``->`` (spaces optional).
* Background/neutral species (``N2``, ``O2``, ``M``, ``O``, ``N2O``, …) may
  appear freely; anything not in ``SPECIES`` is ignored for matrix
  assembly.  Their densities must instead be folded into the rate callable.
* **Exactly one** tracked species must appear on the reactant side — the
  *driver*, which selects the matrix column.  A reaction between two
  tracked species (e.g. ion-ion recombination) is non-linear in
  :math:`\vec{n}` and cannot be represented; the parser raises an error.

Multiplier keys
---------------

JSON configuration files can scale individual reactions with a
``reaction_multipliers`` dictionary keyed by the reaction string
(:ref:`Chap:ModifyingReactions`).  Keys are normalised — spaces stripped,
``→`` replaced by ``->`` — before comparison, so ``"e+N2->2e+N2+"`` matches
the ``REACTIONS`` entry ``"e + N2 -> 2e + N2+"``.

.. warning::

   A key that matches no reaction is *ignored*.  :func:`build_R` prints a
   warning when this happens, but the compiled fast path used by the shipped
   mechanisms does not — a typo in a JSON multiplier key silently leaves the
   reaction at its default rate.  Check the ``Configs:`` line and the
   reaction list in the output-file header when in doubt.

Fast path
---------

:func:`build_R` re-parses every string on every call and is the simplest
entry point.  :func:`compile_reactions` pre-parses the stoichiometry once at
module load into ``(driver_column, delta_vector, rate_fn, key)`` tuples;
:func:`build_R_from_compiled` then assembles :math:`\bm{R}` with no string
handling on the hot path.  The shipped mechanisms use the compiled pair,
since ``get_R`` is called at every :math:`E/N` sample of every integration
step of every determinant evaluation.

.. literalinclude:: ../../../src/incept1d/reactions.py
   :language: python
   :pyobject: build_R

API reference
-------------

.. automodule:: incept1d.reactions
   :members:
   :undoc-members:
