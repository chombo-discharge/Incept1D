# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Reactions.py — declarative R-matrix assembly for mechanism files.

Mechanism files define reactions as a list of (reaction_string, rate_callable)
pairs.  The R-matrix is assembled automatically by parsing each string to
determine which tracked species are produced and consumed.  Reaction multipliers
(for Config/JSON overrides) are keyed by the same reaction strings.

Reaction string format
----------------------
    "e + N2 -> 2e + N2+"

Rules:

- Species tokens must be separated by ` + ` (spaces around the plus sign) so
  that cation names like `N2+` are never broken apart by the tokeniser.
- Optional leading integer stoichiometric coefficient: `2e`, `3N2`.
- Arrow ` -> ` (spaces optional but the `->` is required).
- Background/neutral species (N2, O2, M, O, N, …) may appear freely; they are
  ignored during R-matrix assembly since they are not in SPECIES.
- Exactly one SPECIES member must appear in the reactants — this is the
  "driver" column j.

Multiplier key normalisation
-----------------------------
Keys in a `reaction_multipliers` dict are normalised (spaces stripped,
`→` replaced with `->`) before comparison, so JSON config files may use
`"e+N2->2e+N2+"` to match the REACTIONS entry `"e + N2 -> 2e + N2+"`.
"""

import re
import sys

import numpy as np


def compile_reactions(reactions, species):
    """Pre-parse stoichiometric topology once; return list used by build_R_from_compiled.

    Each entry is (j, delta, rate_fn, norm_key):
      j        — driver column index (the single tracked species in the reactants)
      delta    — shape-(n,) net stoichiometry vector (products − reactants per species)
      rate_fn  — the original callable rate(EN, p, T) -> float
      norm_key — normalised reaction string for multiplier lookup

    Call once at module load; pass the result to build_R_from_compiled on every
    get_R invocation to avoid repeated regex parsing.
    """
    sp_idx = {s: i for i, s in enumerate(species)}
    n = len(species)
    compiled = []
    for rxn_str, rate_fn in reactions:
        reactants, products = parse_reaction(rxn_str)
        drivers = [s for s in reactants if s in sp_idx]
        if len(drivers) != 1:
            raise ValueError(
                f"Reaction '{rxn_str}': expected exactly one tracked species "
                f"in reactants, got {drivers!r}."
            )
        j = sp_idx[drivers[0]]
        delta = np.zeros(n)
        for sp, idx in sp_idx.items():
            delta[idx] = products.get(sp, 0) - reactants.get(sp, 0)
        compiled.append((j, delta, rate_fn, _normalize(rxn_str)))
    return compiled


def build_R_from_compiled(compiled, n, EN, p, T, multipliers=None):
    """Assemble R from pre-compiled topology — no regex, no string parsing.

    Drop-in replacement for build_R when the caller holds a compiled list
    produced by compile_reactions.  multipliers is handled identically to
    build_R: keys are normalised before comparison, and a key matching no
    reaction is reported to stderr once.
    """
    norm_mult = (
        {_normalize(k): v for k, v in multipliers.items()} if multipliers else {}
    )
    if norm_mult:
        _warn_unknown_multipliers(norm_mult, {c[3] for c in compiled})
    R = np.zeros((n, n))
    for j, delta, rate_fn, norm_key in compiled:
        K = rate_fn(EN, p, T) * norm_mult.get(norm_key, 1.0)
        R[:, j] += delta * K
    return R


def _normalize(s: str) -> str:
    """Canonical form of a reaction string for dict-key comparison."""
    return s.replace(" ", "").replace("→", "->")


# Multiplier keys already reported, so that a typo is announced once rather
# than at every E/N sample of every integration step.
_WARNED_KEYS = set()


def _warn_unknown_multipliers(norm_mult, known_keys):
    """Report multiplier keys that match no reaction in this mechanism.

    Parameters
    ----------
    norm_mult : dict
        Normalised ``{reaction_string: multiplier}`` mapping.
    known_keys : set of str
        Normalised reaction strings the mechanism actually defines.

    Notes
    -----
    A key is reported to stderr the first time it is seen; a typo in a JSON
    configuration is therefore visible without flooding the output.
    """
    for k in norm_mult:
        if k not in known_keys and k not in _WARNED_KEYS:
            _WARNED_KEYS.add(k)
            print(
                f"WARNING: reaction multiplier key '{k}' does not match "
                f"any reaction in this mechanism — ignored.",
                file=sys.stderr,
                flush=True,
            )


def _parse_side(text: str) -> dict:
    """Parse one side of a reaction equation into {species_name: count}.

    Splits on whitespace-surrounded `+` so that trailing `+` in cation names
    (e.g. `N2+`) is never treated as a separator.
    """
    result = {}
    for token in re.split(r"\s+\+\s+", text.strip()):
        token = token.strip()
        if not token:
            continue
        m = re.match(r"^(\d+)?(.*)", token)
        coeff = int(m.group(1)) if m.group(1) else 1
        name = m.group(2).strip()
        if name:
            result[name] = result.get(name, 0) + coeff
    return result


def parse_reaction(rxn_str: str):
    """Return (reactants_dict, products_dict) for a reaction string."""
    lhs, rhs = re.split(r"\s*->\s*", rxn_str, maxsplit=1)
    return _parse_side(lhs), _parse_side(rhs)


def build_R(reactions, species, EN, p, T, multipliers=None):
    """
    Assemble the R-matrix from a declarative reaction list.

    Parameters
    ----------
    reactions : list of (rxn_str, rate_callable)
        Each rate_callable has signature rate(EN, p, T) -> float and should
        return the volumetric rate coefficient [s^-1] for that reaction.
    species : list of str
        Ordered species names; index i is row/column i of R.
    EN : float
        Reduced electric field in Townsend.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.
    multipliers : dict, optional
        {reaction_string: float}.  Keys are normalised before comparison so
        spaces are optional.  A warning is printed to stderr for any key that
        does not match a reaction in this mechanism.

    Returns
    -------
    numpy.ndarray, shape (n, n)
        Reaction matrix R; R[i,j] is the rate [s^-1] at which one carrier of
        species j produces (positive) or destroys (negative, diagonal) one
        particle of species i.
    """
    if multipliers is None:
        multipliers = {}

    norm_mult = {_normalize(k): v for k, v in multipliers.items()}
    _warn_unknown_multipliers(
        norm_mult, {_normalize(rxn_str) for rxn_str, _ in reactions}
    )

    n = len(species)
    sp_idx = {s: i for i, s in enumerate(species)}
    R = np.zeros((n, n))

    for rxn_str, rate_fn in reactions:
        scale = norm_mult.get(_normalize(rxn_str), 1.0)
        K = rate_fn(EN, p, T) * scale

        reactants, products = parse_reaction(rxn_str)

        drivers = [s for s in reactants if s in sp_idx]
        if len(drivers) != 1:
            raise ValueError(
                f"Reaction '{rxn_str}': expected exactly one tracked species "
                f"in reactants, got {drivers!r}.  Check SPECIES and the "
                f"reaction string."
            )
        j = sp_idx[drivers[0]]

        for sp, i in sp_idx.items():
            net = products.get(sp, 0) - reactants.get(sp, 0)
            if net:
                R[i, j] += net * K

    return R
