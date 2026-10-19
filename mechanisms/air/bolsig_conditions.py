# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
bolsig_conditions.py — gas conditions recorded in a BOLSIG+ output file.

The electron transport tables are computed for one gas mixture and one gas
temperature, and the neutral densities of a mechanism must use the same mixture.
This helper reads both from the table and checks the mechanism against them:

    resolve_composition(path, x_n2, x_o2, default)  → (xN2, xO2)
    read_conditions(path)                           → ({species: fraction}, T_gas)
    check_temperature(T, T_table, path)             → warns once on a large mismatch

BOLSIG+ output written by the stand-alone solver records the mixture in its
"Conditions" block ("Mole fraction N2   0.7900") and the gas temperature as
"Gas temperature (K)". Tables generated on LXCat record only "Tgas = 300 K",
without the mixture.
"""

import os
import re
import warnings

#: Largest allowed difference between configured and tabulated mole fractions.
COMPOSITION_TOL = 1e-3

#: Gas temperature difference (K) above which check_temperature() warns.
TEMPERATURE_TOL = 20.0

_CACHE = {}
_WARNED = set()


def read_conditions(path):
    """
    Return the mole fractions and gas temperature recorded in a BOLSIG+ file.

    Returns
    -------
    fractions : dict
        {species: mole fraction}, empty if the file does not record the mixture.
    T_gas : float or None
        Gas temperature in K, or None if the file does not record it.
    """
    path = os.path.realpath(path)
    if path not in _CACHE:
        with open(path, "rb") as f:
            text = f.read().decode("ascii", errors="replace").replace("\r", "\n")
        fractions = {
            m.group(1): float(m.group(2))
            for m in re.finditer(r"Mole fraction\s+(\S+)\s+([-+0-9.Ee]+)", text)
        }
        m = re.search(r"Gas temperature \(K\)\s+([-+0-9.Ee]+)", text) or re.search(
            r"Tgas\s*=\s*([-+0-9.Ee]+)\s*K", text
        )
        _CACHE[path] = (fractions, float(m.group(1)) if m else None)
    return _CACHE[path]


def resolve_composition(path, x_n2, x_o2, default):
    """
    N2 and O2 mole fractions for the neutral densities, consistent with the table.

    - The table records its mixture: without x_n2/x_o2 the table's values are
      used; configured values that differ by more than COMPOSITION_TOL raise.
    - The table does not record its mixture: configured values are used as
      given; without them the mechanism default is used, with a warning, since
      the consistency cannot be checked.

    Parameters
    ----------
    path : str
        BOLSIG+ output file of the electron transport data.
    x_n2, x_o2 : float or None
        Configured mole fractions (config keys x_n2 and x_o2), or None.
    default : tuple of float
        The mechanism's own (xN2, xO2).

    Returns
    -------
    tuple of float
        (xN2, xO2).

    Raises
    ------
    ValueError
        If the configured mixture differs from the mixture of the table.
    """
    fractions, _ = read_conditions(path)
    name = os.path.basename(path)
    if "N2" in fractions and "O2" in fractions:
        table = (fractions["N2"], fractions["O2"])
        if x_n2 is None:
            return table
        if abs(x_n2 - table[0]) > COMPOSITION_TOL or abs(x_o2 - table[1]) > COMPOSITION_TOL:
            raise ValueError(
                f"x_n2 = {x_n2:g}, x_o2 = {x_o2:g} differ from the mixture of the "
                f"transport data in {name} (N2 {table[0]:g}, O2 {table[1]:g}). "
                f"Use transport data computed for the configured mixture."
            )
        return (x_n2, x_o2)
    if x_n2 is None:
        key = ("composition", os.path.realpath(path))
        if key not in _WARNED:
            _WARNED.add(key)
            warnings.warn(
                f"{name} does not record the gas mixture of its transport data; "
                f"using the mechanism default N2 {default[0]:g}, O2 {default[1]:g}. "
                f"Set x_n2 and x_o2 to the mixture the table was computed for.",
                stacklevel=2,
            )
        return tuple(default)
    return (x_n2, x_o2)


def check_temperature(T, T_table, path):
    """Warn once if the gas temperature T differs from the table's by more than TEMPERATURE_TOL."""
    if T_table is None or abs(T - T_table) <= TEMPERATURE_TOL:
        return
    key = ("temperature", os.path.realpath(path), float(T))
    if key not in _WARNED:
        _WARNED.add(key)
        warnings.warn(
            f"Gas temperature T = {T:g} K differs from the {T_table:g} K of the "
            f"transport data in {os.path.basename(path)}.",
            stacklevel=3,
        )
