Reactions
=========

``Reactions.py`` assembles the reaction-rate matrix :math:`\bm{R}` (see
:doc:`../model`) from a *declarative* list of reaction strings, so a
mechanism file states its chemistry as human-readable equations like

.. code-block:: text

   "e + N2 -> 2e + N2+"
   "e + O2 + O2 -> O2- + O2"

instead of hand-writing stoichiometry into a matrix. Each reaction is paired
with a rate callable ``rate(EN, p, T) -> float`` giving the rate coefficient
in SI units.

Parsing rules (see the module docstring for the full specification):

- species tokens are separated by ``" + "`` (spaces required, so that cation
  names like ``N2+`` are never split);
- an optional leading integer is a stoichiometric coefficient (``2e``);
- exactly one tracked species (one listed in the mechanism's ``SPECIES``)
  must appear on the reactant side — this is the matrix column ("driver")
  that the reaction is keyed on;
- background/neutral species (``N2``, ``O2``, ``M``, ...) may appear freely
  and are ignored for matrix assembly.

:func:`build_R` re-parses every reaction string on every call and is the
simplest entry point. :func:`compile_reactions` pre-parses the stoichiometry
once at module load time into a list of ``(driver_column, delta_vector,
rate_fn, normalized_key)`` tuples; :func:`build_R_from_compiled` then
assembles :math:`\bm{R}` from that pre-compiled form with no regex parsing
on the hot path — this is the pair used by the example mechanisms
(``Air/Air_Pancheshnyi.py``) since ``get_R`` is called at every ``E/N``
sample point of every solve.

.. literalinclude:: ../../Reactions.py
   :language: python
   :pyobject: build_R

API reference
--------------

.. automodule:: Reactions
   :members:
   :undoc-members:
