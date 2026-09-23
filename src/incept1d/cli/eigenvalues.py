# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
``incept1d eigenvalues`` — plot eigenvalues of the transport matrix A = R V⁻¹
vs. E/N.

Default mode plots all N eigenvalue tracks (Re λ_j / N) vs E/N for a fixed
pressure, one subplot per configuration.  With ``--pressure-scan`` the
leading (or a chosen) eigenvalue is plotted vs E/N for several log-spaced
pressures.  See :mod:`incept1d.eigenvalues` for the computation.
"""

import datetime
import os

import numpy as np
import matplotlib.pyplot as plt

from incept1d.constants import kB as _kB
from incept1d.eigenvalues import compute_eigenvalues, compute_pressure_scan
from incept1d.mechanism import load_mechanism, read_json_configs


def _write_default(
    outfile, EN_range, eigvals_red, mech_name, p, T, N_density, label, n
):
    """Write all eigenvalue tracks to a tab-separated file (default mode)."""
    SEP = "\t"
    columns = [("EN_Td", lambda i: EN_range[i])]
    for j in range(n):
        columns.append(
            (
                f"Re_lam{j}_N_m2",
                lambda i, _ev=eigvals_red, _j=j: np.real(_ev[i, _j]),
            )
        )
    W = max(18, max(len(c) for c, _ in columns) + 2)

    with open(outfile, "w") as fh:
        fh.write(f"# Mechanism:   {mech_name}\n")
        fh.write(f"# Pressure:    {p} bar\n")
        fh.write(f"# Temperature: {T} K\n")
        fh.write(f"# N:           {N_density:.6e} m^-3\n")
        fh.write(f"# Config:      {label}\n")
        fh.write(f"# Generated:   {datetime.date.today().isoformat()}\n")
        fh.write("#\n")
        fh.write("# " + SEP.join(f"{c:<{W}}" for c, _ in columns) + "\n")
        for i in range(len(EN_range)):
            row = [f"{fn(i):<{W}.6e}" for _, fn in columns]
            fh.write(SEP.join(row) + "\n")

    print(f"Written: {outfile}")


def _write_pressure_scan(
    outfile, EN_range, scan_data, pressures, mech_name, T, eig_index, label
):
    """Write pressure-scan eigenvalue data to a tab-separated file."""
    SEP = "\t"
    p_labels = [f"Re_lam{eig_index}_N_p{p:.3g}bar" for p in pressures]
    columns = ["EN_Td"] + p_labels
    W = max(22, max(len(c) for c in columns) + 2)

    with open(outfile, "w") as fh:
        fh.write(f"# Mechanism:        {mech_name}\n")
        fh.write(f"# Temperature:      {T} K\n")
        fh.write(f"# Eigenvalue index: {eig_index}\n")
        fh.write(
            f"# Pressure range:   {pressures[0]:.3g} – {pressures[-1]:.3g} bar"
            f"  ({len(pressures)} log-spaced)\n"
        )
        fh.write(f"# Config:           {label}\n")
        fh.write(f"# Generated:        {datetime.date.today().isoformat()}\n")
        fh.write("#\n")
        fh.write("# " + SEP.join(f"{c:<{W}}" for c in columns) + "\n")
        for i in range(len(EN_range)):
            row = [f"{EN_range[i]:<{W}.6e}"]
            row += [f"{scan_data[i, pi]:<{W}.6e}" for pi in range(len(pressures))]
            fh.write(SEP.join(row) + "\n")

    print(f"Written: {outfile}")


HELP = "Eigenvalues of the local transport matrix A = R V⁻¹ vs. E/N."
DESCRIPTION = (
    "Plot eigenvalues of A = R V^{-1} as a function of E/N. "
    "Positive real eigenvalues indicate net spatial growth of the "
    "discharge (ionisation exceeds attachment).  "
    "One subplot is produced per configuration."
)


def add_arguments(parser):
    """Register the command-line arguments on *parser*."""
    parser.add_argument(
        "mechanism",
        help="Path to mechanism Python file (e.g. Air/Air_Hosl.py).",
    )
    parser.add_argument(
        "configs",
        nargs="*",
        metavar="CONFIG.json",
        help=(
            "One or more JSON configuration files.  Each file may contain a "
            'single configuration object or a list under a "configurations" '
            "key.  If omitted, a single baseline configuration is used.  "
            "Each configuration produces one subplot."
        ),
    )
    parser.add_argument(
        "--p",
        type=float,
        default=1.0,
        help="Gas pressure in bar — used in default mode (default: 1.0).",
    )
    parser.add_argument(
        "--T",
        type=float,
        default=293.0,
        help="Gas temperature in Kelvin (default: 293.0).",
    )
    parser.add_argument(
        "--EN-lo",
        type=float,
        default=10.0,
        help="Lower E/N bound in Td (default: 10).",
    )
    parser.add_argument(
        "--EN-hi",
        type=float,
        default=500.0,
        help="Upper E/N bound in Td (default: 500).",
    )
    parser.add_argument(
        "--EN-num",
        type=int,
        default=500,
        help="Number of log-spaced E/N points (default: 500).",
    )
    parser.add_argument(
        "--pressure-scan",
        action="store_true",
        help=(
            "Activate pressure-scan mode: plot the selected eigenvalue track "
            "vs E/N for several log-spaced pressures."
        ),
    )
    parser.add_argument(
        "--p-min",
        type=float,
        default=1e-3,
        help="Minimum pressure in bar for pressure-scan mode (default: 1e-3).",
    )
    parser.add_argument(
        "--p-max",
        type=float,
        default=10.0,
        help="Maximum pressure in bar for pressure-scan mode (default: 10).",
    )
    parser.add_argument(
        "--p-num",
        type=int,
        default=5,
        help="Number of log-spaced pressures in pressure-scan mode (default: 5).",
    )
    parser.add_argument(
        "--eig-index",
        type=int,
        default=0,
        help=(
            "Eigenvalue track index to show in pressure-scan mode "
            "(default: 0 = leading/maximum-real mode)."
        ),
    )
    parser.add_argument(
        "--write-to-file",
        type=str,
        default=None,
        metavar="FILE",
        help="Write eigenvalue data to tab-separated file(s).",
    )


def run(args, parser):
    """Run the command with parsed *args*; *parser* is used for ``parser.error``."""
    mech_name = os.path.basename(args.mechanism)
    raw_dicts = read_json_configs(args.configs) if args.configs else [{}]
    mods = [load_mechanism(args.mechanism, d) for d in raw_dicts]
    n_mods = len(mods)

    EN_range = np.logspace(np.log10(args.EN_lo), np.log10(args.EN_hi), args.EN_num)

    base_out, ext_out = (
        os.path.splitext(args.write_to_file) if args.write_to_file else (None, ".txt")
    )

    # ------------------------------------------------------------------
    # Pressure-scan mode
    # ------------------------------------------------------------------
    if args.pressure_scan:
        pressures = np.logspace(np.log10(args.p_min), np.log10(args.p_max), args.p_num)

        n = len(mods[0].SPECIES)
        if args.eig_index >= n:
            raise ValueError(
                f"--eig-index {args.eig_index} out of range for {n} species."
            )

        cmap = plt.cm.viridis
        p_colors = cmap(np.linspace(0.1, 0.9, len(pressures)))
        p_mid = np.sqrt(args.p_min * args.p_max)
        N_mid = p_mid * 1e5 / (_kB * args.T)
        linthresh = 1e-21 / N_mid

        fig, axs = plt.subplots(
            n_mods,
            1,
            figsize=(11, 5 * n_mods),
            sharex=True,
            squeeze=False,
        )
        axs = axs[:, 0]

        for ki, mod in enumerate(mods):
            label = mod.label or "Baseline"
            scan_red = compute_pressure_scan(
                mod, EN_range, pressures, args.T, args.eig_index
            )
            ax = axs[ki]
            for pi, p in enumerate(pressures):
                ax.plot(
                    EN_range,
                    scan_red[:, pi],
                    color=p_colors[pi],
                    label=f"p = {p:.3g} bar",
                )
            ax.set_xscale("log")
            ax.set_yscale("symlog", linthresh=linthresh)
            ax.axhline(0, color="k", lw=0.8, ls="--")
            ax.set_ylabel(rf"Re($\lambda_{{{args.eig_index}}}/N$)  (m$^2$)")
            ax.set_title(
                rf"$\lambda_{{{args.eig_index}}}$ of $A = RV^{{-1}}$  —  "
                f"{label},  {mech_name},  T = {args.T} K"
            )
            ax.legend(loc="best", fontsize=9)
            ax.grid(True, which="both", ls="--", alpha=0.4)

            if args.write_to_file:
                outfile = f"{base_out}_{label.lower()}{ext_out}"
                _write_pressure_scan(
                    outfile,
                    EN_range,
                    scan_red,
                    pressures,
                    mech_name,
                    args.T,
                    args.eig_index,
                    label,
                )

        axs[-1].set_xlabel("E/N (Td)")
        plt.tight_layout()
        plt.show()
        return

    # ------------------------------------------------------------------
    # Default mode: all eigenvalue tracks vs E/N at fixed pressure
    # ------------------------------------------------------------------
    n = len(mods[0].SPECIES)
    N_density = args.p * 1e5 / (_kB * args.T)
    _colors = plt.cm.tab10(np.linspace(0, 0.9, n))

    all_eigvals_red = []
    for mod in mods:
        ev = compute_eigenvalues(mod, EN_range, args.p, args.T)
        all_eigvals_red.append(ev / N_density)

    fig, axs = plt.subplots(
        n_mods,
        1,
        figsize=(11, 5 * n_mods),
        sharex=True,
        squeeze=False,
    )
    axs = axs[:, 0]

    for ki, mod in enumerate(mods):
        label = mod.label or "Baseline"
        eigvals_red = all_eigvals_red[ki]
        ax = axs[ki]
        for j in range(n):
            ax.plot(
                EN_range,
                np.real(eigvals_red[:, j]),
                color=_colors[j],
                label=f"$\\lambda_{{{j}}}/N$",
            )
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1.0 / N_density)
        ax.axhline(0, color="k", lw=0.8, ls="--")
        ax.set_ylabel(r"Re($\lambda/N$)  (m$^{2}$)")
        ax.set_title(
            rf"Eigenvalues of $A = RV^{{-1}}$  —  "
            f"{label},  {mech_name},  p = {args.p} bar,  T = {args.T} K"
        )
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, which="both", ls="--", alpha=0.4)

        if args.write_to_file:
            outfile = f"{base_out}_{label.lower()}{ext_out}"
            _write_default(
                outfile,
                EN_range,
                eigvals_red,
                mech_name,
                args.p,
                args.T,
                N_density,
                label,
                n,
            )

    axs[-1].set_xlabel("E/N (Td)")
    plt.tight_layout()
    plt.show()
