# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
``incept1d chombo`` — export transport-coefficient / rate-coefficient tables
for the external 3-D solver `chombo-discharge`, and plot them.

Output columns
--------------
    1      E/N [Townsend]
    2      alpha/N [m^2]
    3      eta/N [m^2]
    4      Electron mean energy [eV]
    5-6    mu*N, D*N for e
    7-8    mu*N, D*N for the next species
    ...    (two columns per species in SPECIES order)
    next   Raw rate coefficient for each reaction in REACTIONS order

Rate coefficients are in m^3/s (two-body) or m^6/s (three-body).
Three-body reactions are flagged in the column header.
"""

import math
import os

import numpy as np
import matplotlib.pyplot as plt

from incept1d.chombo import (
    _load_modifier,
    _neutral_factor,
    build_header,
    generate,
    load_raw_mechanism,
)
from incept1d.constants import kB


def plot_results(mod, EN_arr, columns, p, T, mech_path):
    """Generate a 5-subplot transport coefficient figure."""
    N = p * 1e5 / (kB * T)
    n_species = len(mod.SPECIES)
    species_set = set(mod.SPECIES)

    fig = plt.figure(figsize=(14, 13))
    gs = fig.add_gridspec(3, 2, hspace=0.38, wspace=0.28)
    ax_rates = fig.add_subplot(gs[0, :])  # full top row
    ax_mu = fig.add_subplot(gs[1, 0])
    ax_diff = fig.add_subplot(gs[1, 1])
    ax_energy = fig.add_subplot(gs[2, 0])
    ax_townsend = fig.add_subplot(gs[2, 1])

    # ------------------------------------------------------------------ #
    # Subplot 1: Reaction rates  k * [neutrals] / N  vs E/N              #
    # Plotted quantity:  rate_fn(EN,p,T) / N  [m^3/s]                    #
    # Two-body:   k * xA               (pressure-independent)            #
    # Three-body: k * xA^n * N^(n-1)   (scales with p^(n-1))            #
    # ------------------------------------------------------------------ #
    # col 0: EN, 1: alpha, 2: eta, 3: energy, 4..3+2n: mu/D, 4+2n..: ki
    rate_col_start = 4 + 2 * n_species
    for rxn_idx, (rxn_str, rate_fn) in enumerate(mod.REACTIONS):
        raw_ki = columns[rate_col_start + rxn_idx]
        nf = _neutral_factor(rxn_str, species_set, mod, N)
        rate_to_plot = raw_ki * nf / N  # m^3/s

        mask = rate_to_plot > 0
        if mask.any():
            ax_rates.loglog(EN_arr[mask], rate_to_plot[mask], label=rxn_str)

    all_rates = np.concatenate(
        [
            columns[rate_col_start + i]
            * _neutral_factor(rxn_str, species_set, mod, N)
            / N
            for i, (rxn_str, _) in enumerate(mod.REACTIONS)
        ]
    )
    y_max = all_rates[all_rates > 0].max() * 2.0
    y_min = y_max * 1e-10

    ax_rates.set_xlim(EN_arr[0], EN_arr[-1])
    ax_rates.set_ylim(y_min, y_max)
    ax_rates.set_xlabel("E/N  [Td]")
    ax_rates.set_ylabel(
        r"$k \cdot [\mathrm{neutrals}] \,/\, N \;\;[\mathrm{m}^3\,\mathrm{s}^{-1}]$"
    )
    ax_rates.set_title("Reaction rates")
    ax_rates.legend(fontsize=7, loc="best")
    ax_rates.grid(True, which="both", alpha=0.3)

    # ------------------------------------------------------------------ #
    # Subplot 2: Reduced mobilities  mu*N  [m^-1 V^-1 s^-1]              #
    # ------------------------------------------------------------------ #
    mu_col_start = 4
    for i, sp in enumerate(mod.SPECIES):
        ax_mu.loglog(EN_arr, columns[mu_col_start + 2 * i], label=sp)

    ax_mu.set_xlim(EN_arr[0], EN_arr[-1])
    ax_mu.set_xlabel("E/N  [Td]")
    ax_mu.set_ylabel(r"$\mu N \;\;[\mathrm{m}^{-1}\,\mathrm{V}^{-1}\,\mathrm{s}^{-1}]$")
    ax_mu.set_title("Reduced mobilities")
    ax_mu.legend(loc="best")
    ax_mu.grid(True, which="both", alpha=0.3)

    # ------------------------------------------------------------------ #
    # Subplot 3: Reduced diffusion coefficients  D*N  [m^-1 s^-1]        #
    # ------------------------------------------------------------------ #
    dn_col_start = 5
    for i, sp in enumerate(mod.SPECIES):
        ax_diff.loglog(EN_arr, columns[dn_col_start + 2 * i], label=sp)

    ax_diff.set_xlim(EN_arr[0], EN_arr[-1])
    ax_diff.set_xlabel("E/N  [Td]")
    ax_diff.set_ylabel(r"$D N \;\;[\mathrm{m}^{-1}\,\mathrm{s}^{-1}]$")
    ax_diff.set_title("Reduced diffusion coefficients")
    ax_diff.legend(loc="best")
    ax_diff.grid(True, which="both", alpha=0.3)

    # ------------------------------------------------------------------ #
    # Subplot 4: Electron mean energy  [eV]                               #
    # ------------------------------------------------------------------ #
    ax_energy.loglog(EN_arr, columns[3])
    ax_energy.set_xlim(EN_arr[0], EN_arr[-1])
    ax_energy.set_xlabel("E/N  [Td]")
    ax_energy.set_ylabel("Mean electron energy  [eV]")
    ax_energy.set_title("Electron mean energy")
    ax_energy.grid(True, which="both", alpha=0.3)

    # ------------------------------------------------------------------ #
    # Subplot 5: Reduced Townsend coefficients alpha/N and eta/N  [m^2]  #
    # ------------------------------------------------------------------ #
    alpha_over_N = columns[1]
    eta_over_N = columns[2]

    mask_a = alpha_over_N > 0
    mask_e = eta_over_N > 0
    if mask_a.any():
        ax_townsend.loglog(
            EN_arr[mask_a], alpha_over_N[mask_a], label=r"$\alpha/N$  (ionization)"
        )
    if mask_e.any():
        ax_townsend.loglog(
            EN_arr[mask_e], eta_over_N[mask_e], label=r"$\eta/N$  (attachment)"
        )

    ax_townsend.set_xlim(EN_arr[0], EN_arr[-1])
    ax_townsend.set_xlabel("E/N  [Td]")
    ax_townsend.set_ylabel(r"$\alpha/N,\,\eta/N \;\;[\mathrm{m}^{2}]$")
    ax_townsend.set_title("Reduced Townsend ionization and attachment coefficients")
    ax_townsend.legend(loc="best")
    ax_townsend.grid(True, which="both", alpha=0.3)

    mech_name = os.path.basename(mech_path)
    fig.suptitle(
        f"Transport coefficients — {mech_name},  p = {p} bar,  T = {T} K",
        fontsize=12,
    )

    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


HELP = "Transport/rate-coefficient tables for chombo-discharge."
DESCRIPTION = "Generate chombo-discharge transport data from a mechanism file."


def add_arguments(parser):
    """Register the command-line arguments on *parser*."""
    parser.add_argument(
        "mechanism",
        help="Path to mechanism Python file (e.g. mechanisms/air/air_pancheshnyi.py)",
    )
    parser.add_argument("--modifier", default=None, help="Path to modifier JSON file")
    parser.add_argument(
        "--config-label",
        default=None,
        help="Label of configuration to use in modifier JSON",
    )
    parser.add_argument(
        "--min-EN",
        type=float,
        default=1.0,
        metavar="TD",
        help="Minimum E/N in Townsend (default: 1)",
    )
    parser.add_argument(
        "--max-EN",
        type=float,
        default=1e6,
        metavar="TD",
        help="Maximum E/N in Townsend (default: 1e6)",
    )
    parser.add_argument(
        "--num-EN", type=int, default=1000, help="Number of E/N points (default: 1000)"
    )
    parser.add_argument(
        "--pressure", type=float, default=1.0, metavar="BAR", help="Gas pressure in bar"
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        metavar="K",
        help="Gas temperature in K",
    )
    parser.add_argument(
        "--write-to-file",
        default=None,
        metavar="FILE",
        help="Write data to FILE (default: do not write)",
    )


def run(args, parser):
    """Run the command with parsed *args*; *parser* is used for ``parser.error``."""
    write = args.write_to_file is not None
    args.output = args.write_to_file

    pre_vars = {}
    reaction_multipliers = {}
    if args.modifier:
        pre_vars, reaction_multipliers = _load_modifier(
            args.modifier, args.config_label
        )

    mod = load_raw_mechanism(args.mechanism, pre_vars)

    EN_arr = np.logspace(
        math.log10(args.min_EN),
        math.log10(args.max_EN),
        args.num_EN,
    )

    columns, headers = generate(
        mod, EN_arr, args.pressure, args.temperature, reaction_multipliers
    )

    if write:
        header_str = build_header(
            args.mechanism,
            args.modifier,
            args.config_label,
            args.pressure,
            args.temperature,
            headers,
            pre_exec_vars=pre_vars,
            EN_min=args.min_EN,
            EN_max=args.max_EN,
            num_EN=args.num_EN,
            species=mod.SPECIES,
        )
        data = np.column_stack(columns)
        np.savetxt(args.output, data, delimiter="\t", header=header_str)
        print(f"Written: {args.output}  ({len(EN_arr)} rows, {data.shape[1]} columns)")

    plot_results(mod, EN_arr, columns, args.pressure, args.temperature, args.mechanism)
