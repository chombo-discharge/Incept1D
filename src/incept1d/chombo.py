"""
Transport-coefficient tables for the external 3-D plasma solver
`chombo-discharge`, generated from a mechanism file.

The command-line front end is :mod:`incept1d.cli.chombo`.

Output columns
--------------
    1      E/N [Townsend]
    2      alpha/N [m^2]
    3      eta/N [m^2]
    4      Electron mean energy [eV]
    5-6    mu*N, D*N for e
    5-6    mu*N, D*N for N2+
    ...    (two columns per species in SPECIES order)
    next   Raw rate coefficient for each reaction in REACTIONS order

Rate coefficients are in m^3/s (two-body) or m^6/s (three-body).
Three-body reactions are flagged in the column header.
"""

import datetime
import json
import os
import sys
import types

import numpy as np

from incept1d.constants import kB, Q
from incept1d.reactions import parse_reaction

# ---------------------------------------------------------------------------
# Modifier / JSON config loading
# ---------------------------------------------------------------------------


def _load_modifier(modifier_path, config_label=None):
    """Return (pre_exec_vars, reaction_multipliers) from a modifier JSON file.

    The JSON may use the `configurations` array format.  If config_label is
    given, the matching entry is used; otherwise the first entry is used.
    If the JSON is flat (no `configurations` key), it is used directly.
    """
    mod_dir = os.path.dirname(os.path.abspath(modifier_path))
    with open(modifier_path) as f:
        data = json.load(f)

    if "configurations" in data:
        configs = data["configurations"]
        if config_label is not None:
            matches = [c for c in configs if c.get("label") == config_label]
            if not matches:
                raise ValueError(
                    f"No configuration with label '{config_label}' in {modifier_path}. "
                    f"Available: {[c.get('label') for c in configs]}"
                )
            cfg = matches[0]
        else:
            cfg = configs[0]
    else:
        cfg = data

    pre_vars = {}
    cs = cfg.get("cross_sections")
    if cs is not None:
        if not os.path.isabs(cs):
            candidate = os.path.join(mod_dir, cs)
            cs = candidate if os.path.isfile(candidate) else cs
        pre_vars["BOLSIG_FILE"] = cs

    for attr, key in [
        ("_GAMMA0", "gamma0"),
        ("_GAMMA1", "gamma1"),
        ("_EREF", "eref"),
        ("_BETA", "beta"),
    ]:
        if cfg.get(key) is not None:
            pre_vars[attr] = cfg[key]

    multipliers = cfg.get("reaction_multipliers", {})
    return pre_vars, multipliers


# ---------------------------------------------------------------------------
# Mechanism loading
# ---------------------------------------------------------------------------


def load_raw_mechanism(mech_path, pre_exec_vars=None):
    """
    Execute a mechanism file and return its namespace as a SimpleNamespace.

    Unlike :func:`incept1d.mechanism.load_mechanism` this applies no
    configuration protocol and does not wrap the result in a ``Mechanism``;
    the raw ``REACTIONS`` list and rate functions are needed to tabulate
    individual rate coefficients.
    """
    with open(mech_path) as f:
        code = f.read()

    namespace = {"__file__": os.path.abspath(mech_path)}
    namespace.update(pre_exec_vars or {})
    exec(compile(code, mech_path, "exec"), namespace)
    mod = types.SimpleNamespace(**namespace)

    for attr in (
        "SPECIES",
        "ELECTRON_INDEX",
        "REACTIONS",
        "ElectronMeanEnergy",
        "get_V",
    ):
        if not hasattr(mod, attr):
            raise AttributeError(
                f"Mechanism '{mech_path}' does not define '{attr}'. "
                "Ensure the file exposes the standard interface."
            )
    return mod


# ---------------------------------------------------------------------------
# Neutral-density factor helpers
# ---------------------------------------------------------------------------


def _neutral_factor(rxn_str, species_set, mod, N):
    """Return the neutral density factor from the LHS of a reaction string.

    Product of (xX * N)^stoich for each neutral reactant, where 'M' means
    the total gas density N.
    """
    reactants, _ = parse_reaction(rxn_str)
    factor = 1.0
    for name, count in reactants.items():
        if name in species_set:
            continue
        if name == "M":
            factor *= N**count
        else:
            mf_attr = f"x{name}"
            xX = getattr(mod, mf_attr, None)
            if xX is None:
                raise RuntimeError(
                    f"Reaction '{rxn_str}': cannot find mole fraction attribute "
                    f"'{mf_attr}' on the mechanism module."
                )
            factor *= (xX * N) ** count
    return factor


def _is_three_body(rxn_str, species_set):
    """Return True if any neutral reactant has stoichiometry > 1 or 'M' is present."""
    reactants, _ = parse_reaction(rxn_str)
    for name, count in reactants.items():
        if name in species_set:
            continue
        if name == "M" or count > 1:
            return True
    return False


def _neutral_reactant_description(rxn_str, species_set, mod, N):
    """Human-readable description of the neutral density factor."""
    reactants, _ = parse_reaction(rxn_str)
    parts = []
    for name, count in reactants.items():
        if name in species_set:
            continue
        if name == "M":
            exp = f"^{count}" if count > 1 else ""
            parts.append(f"N(p,T){exp}")
        else:
            mf_attr = f"x{name}"
            xX = getattr(mod, mf_attr, None)
            xstr = f"x{name}" if xX is None else f"x{name}={xX}"
            exp = f"^{count}" if count > 1 else ""
            parts.append(f"({xstr}*N(p,T)){exp}")
    return " * ".join(parts) if parts else "1"


# ---------------------------------------------------------------------------
# Data generation
# ---------------------------------------------------------------------------


def generate(mod, EN_arr, p, T, reaction_multipliers):
    """Return (columns, headers) for all EN points.

    Column order:
        0          E/N [Townsend]
        1          alpha/N  [m^2]
        2          eta/N    [m^2]
        3          Electron mean energy [eV]
        4..3+2n    mu*N, D*N per species
        4+2n..     Raw rate coefficients per reaction

    columns : list of 1-D numpy arrays, one per output column
    headers : list of str, one per column (no leading '#')
    """
    N = p * 1e5 / (kB * T)
    n_species = len(mod.SPECIES)
    species_set = set(mod.SPECIES)
    electron = mod.SPECIES[mod.ELECTRON_INDEX]
    n_EN = len(EN_arr)
    DN_factor = kB * T / Q

    # --- Compute everything first ---

    energy_col = np.array([mod.ElectronMeanEnergy(en) for en in EN_arr])

    muN_all = np.zeros((n_EN, n_species))
    for idx, en in enumerate(EN_arr):
        V = mod.get_V(en, p, T)
        muN_all[idx] = np.abs(np.diag(V)) / (en * 1e-21)

    norm_mult = {
        rxn_str.replace(" ", "").replace("→", "->"): v
        for rxn_str, v in reaction_multipliers.items()
    }

    ki_cols = []
    ki_hdrs = []
    for rxn_str, rate_fn in mod.REACTIONS:
        three_body = _is_three_body(rxn_str, species_set)
        units = "m^6/s" if three_body else "m^3/s"
        nf = _neutral_factor(rxn_str, species_set, mod, N)
        mult = norm_mult.get(rxn_str.replace(" ", ""), 1.0)
        ki_col = np.zeros(n_EN)
        for idx, en in enumerate(EN_arr):
            assembled = rate_fn(en, p, T)
            raw_ki = assembled / nf if nf != 0.0 else 0.0
            ki_col[idx] = raw_ki * mult
        ki_cols.append(ki_col)
        neutral_desc = _neutral_reactant_description(rxn_str, species_set, mod, N)
        if three_body:
            ki_hdrs.append(
                f"Rate coefficient for '{rxn_str}' [{units}]  [THREE-BODY]  "
                f"Multiply by {neutral_desc} to get volumetric rate [s^-1]. "
                f"Volumetric rate scales as p^2."
            )
        else:
            ki_hdrs.append(
                f"Rate coefficient for '{rxn_str}' [{units}]  "
                f"Multiply by {neutral_desc} to get volumetric rate [s^-1]."
            )

    v_e = muN_all[:, mod.ELECTRON_INDEX] * EN_arr * 1e-21

    alpha_col = np.zeros(n_EN)
    eta_col = np.zeros(n_EN)
    alpha_rxns = []
    eta_rxns = []
    for rxn_idx, (rxn_str, _) in enumerate(mod.REACTIONS):
        reactants, products = parse_reaction(rxn_str)
        drivers = [s for s in reactants if s in species_set]
        if len(drivers) != 1 or drivers[0] != electron:
            continue
        nf = _neutral_factor(rxn_str, species_set, mod, N)
        assembled = ki_cols[rxn_idx] * nf
        if products.get(electron, 0) > 0:
            alpha_col += assembled
            alpha_rxns.append(rxn_str)
        else:
            eta_col += assembled
            eta_rxns.append(rxn_str)
    alpha_col /= v_e
    eta_col /= v_e
    alpha_col /= N
    eta_col /= N

    # --- Assemble in the specified column order ---

    columns = []
    headers = []

    columns.append(EN_arr.copy())
    headers.append("E/N [Townsend]")

    columns.append(alpha_col)
    headers.append(
        f"Reduced Townsend ionization coefficient alpha/N [m^2] = "
        f"sum(ionization volumetric rates) / (v_e * N); "
        f"reactions: {'; '.join(alpha_rxns)}"
    )

    columns.append(eta_col)
    headers.append(
        f"Reduced Townsend attachment coefficient eta/N [m^2] = "
        f"sum(attachment volumetric rates) / (v_e * N); "
        f"reactions: {'; '.join(eta_rxns)}"
    )

    columns.append(energy_col)
    headers.append("Electron mean energy [eV]")

    for i, sp in enumerate(mod.SPECIES):
        columns.append(muN_all[:, i])
        headers.append(f"mu*N for {sp} [m^-1 V^-1 s^-1]")
        columns.append(muN_all[:, i] * DN_factor)
        headers.append(
            f"D*N for {sp} [m^-1 s^-1] (Einstein: D*N = mu*N * kB*T/Q, T={T} K)"
        )

    for ki_col, hdr in zip(ki_cols, ki_hdrs):
        columns.append(ki_col)
        headers.append(hdr)

    return columns, headers


def build_header(
    mech_path,
    modifier_path,
    config_label,
    p,
    T,
    headers,
    pre_exec_vars=None,
    EN_min=None,
    EN_max=None,
    num_EN=None,
    species=None,
):
    """Assemble the full file header string."""
    N = p * 1e5 / (kB * T)
    bolsig = (pre_exec_vars or {}).get("BOLSIG_FILE", "(module default)")
    lines = [
        "chombo-discharge transport coefficient table",
        f"Generated  : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Command    : {' '.join(sys.argv)}",
        f"Python     : {sys.version.split()[0]}",
        "",
        f"Mechanism  : {os.path.abspath(mech_path)}",
        f"Database   : {bolsig}",
    ]
    if modifier_path:
        lab = f" (label: {config_label})" if config_label else ""
        lines.append(f"Modifier   : {os.path.abspath(modifier_path)}{lab}")
    lines += [
        f"Pressure   : {p} bar",
        f"Temperature: {T} K",
        f"N(p,T)     : {N:.6e} m^-3",
    ]
    if species is not None:
        lines.append(f"Species    : {', '.join(species)}")
    if EN_min is not None:
        lines += [
            f"E/N range  : {EN_min} -- {EN_max} Td  ({num_EN} points, log-spaced)",
        ]
    lines += [
        "",
        "Column layout (0-indexed):",
    ]
    for col_idx, hdr in enumerate(headers):
        lines.append(f"  Column {col_idx:2d}: {hdr}")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
