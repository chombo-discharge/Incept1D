"""
Eigenvalues.py — compute and plot eigenvalues of the transport matrix A = R V^{-1}.

For a given reaction mechanism the matrix A = R V^{-1} is formed at each value
of E/N and all eigenvalues are computed.  A positive real eigenvalue at a given
E/N means the discharge has a spatially growing mode (net ionisation exceeds
attachment) and the electron flux amplifies as it travels from cathode to anode.

Two operating modes are available:

Default mode
------------
Plot all N eigenvalue tracks (Re λ_j / N) vs E/N for a fixed pressure.
One subplot per configuration; subplots share the x-axis.

Pressure-scan mode  (--pressure-scan)
--------------------------------------
Plot the leading eigenvalue (or a user-selected track) vs E/N for several
log-spaced pressures.  One subplot per configuration, one line per pressure.

Usage
-----
    python Eigenvalues.py <mechanism> [CONFIG.json ...]
                          [--p PRESSURE] [--T TEMPERATURE]
                          [--EN-lo LO] [--EN-hi HI] [--EN-num N]
                          [--pressure-scan] [--p-min MIN] [--p-max MAX]
                          [--p-num N] [--eig-index IDX]
                          [--write-to-file FILE]

Arguments
---------
mechanism
    Path to a Python file implementing the standard mechanism interface.
    Example: Air/Air_Hosl.py
CONFIG.json
    One or more JSON configuration files.  Each may contain a single object or a
    list under a "configurations" key.  If omitted, a single baseline
    configuration is used.  Each configuration produces one subplot.
--p
    Gas pressure in bar used in default mode (default: 1.0).
--T
    Gas temperature in Kelvin (default: 293.0).
--EN-lo
    Lower E/N bound in Td (default: 10).
--EN-hi
    Upper E/N bound in Td (default: 500).
--EN-num
    Number of log-spaced E/N points (default: 500).
--pressure-scan
    Activate pressure-scan mode.
--p-min
    Minimum pressure in bar for pressure-scan mode (default: 1e-3).
--p-max
    Maximum pressure in bar for pressure-scan mode (default: 10).
--p-num
    Number of log-spaced pressures (default: 5).
--eig-index
    Which eigenvalue track to show in pressure-scan mode (default: 0,
    the leading/maximum-real mode).
--write-to-file FILE
    Write reduced eigenvalues (Re(λ/N)) to tab-separated text file(s).
    Default mode: one file per configuration, all N tracks.
    Pressure-scan mode: one file per configuration, one column per pressure.
"""

import os
import sys
import argparse
import datetime

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Inception import _kB, load_mechanism, _read_json_configs


def _track_step(prev, curr):
    """
    Reorder *curr* eigenvalues to best continue the tracks in *prev*.

    Uses the Hungarian algorithm to find the permutation of *curr* that
    minimises total squared distance in the complex plane from *prev*.
    """
    cost = np.abs(prev[:, np.newaxis] - curr[np.newaxis, :]) ** 2
    _, col_ind = scipy.optimize.linear_sum_assignment(cost)
    return curr[col_ind]


def compute_eigenvalues(mod, EN_range, p, T):
    """
    Compute tracked eigenvalues of A = R V^{-1} for each value of E/N.

    At each E/N the transport matrix A = R(EN) @ inv(V(EN)) is formed and all
    N eigenvalues are computed.  Eigenvalues are continuously tracked across
    E/N using the Hungarian algorithm to prevent spurious reordering jumps.
    The initial ordering is by descending real part (most-growing mode in
    slot 0).

    Parameters
    ----------
    mod : Mechanism
        Loaded mechanism (from Inception.load_mechanism).
    EN_range : array-like, shape (M,)
        E/N values in Townsend.
    p : float
        Gas pressure in bar.
    T : float
        Gas temperature in Kelvin.

    Returns
    -------
    numpy.ndarray, shape (M, N), dtype complex
        Tracked eigenvalues.
    """
    n = len(mod.SPECIES)
    n_EN = len(EN_range)
    eigvals = np.zeros((n_EN, n), dtype=complex)

    # Scan high→low so the ionisation mode (unambiguously the largest real
    # eigenvalue at high E/N) seeds track 0, avoiding misassignment to the two
    # structural zero eigenvalues that exist because columns 1 and 2 of R are
    # identically zero (N2+, O2+ have no off-diagonal source reactions).
    for i, EN in enumerate(EN_range[::-1]):
        R = mod.get_R(EN, p, T)
        V = mod.get_V(EN, p, T)
        A = R @ np.linalg.inv(V)
        eigs = np.linalg.eigvals(A)

        if i == 0:
            eigs = eigs[np.argsort(-np.real(eigs))]
        else:
            eigs = _track_step(eigvals[i - 1], eigs)

        eigvals[i] = eigs

    return eigvals[::-1]


def compute_pressure_scan(mod, EN_range, pressures, T, eig_index):
    """
    Compute Re(λ_{eig_index} / N) vs E/N for each pressure.

    Parameters
    ----------
    mod : Mechanism
    EN_range : array-like, shape (M,)
    pressures : array-like, shape (P,)
        Pressures in bar.
    T : float
        Gas temperature in Kelvin.
    eig_index : int
        Which eigenvalue track to extract.

    Returns
    -------
    numpy.ndarray, shape (M, P)
        Re(λ_{eig_index} / N) for each (E/N, pressure) pair, in m².
    """
    n_EN = len(EN_range)
    n_p = len(pressures)
    result = np.zeros((n_EN, n_p))

    for pi, p in enumerate(pressures):
        N_density = p * 1e5 / (_kB * T)
        ev = compute_eigenvalues(mod, EN_range, p, T)
        result[:, pi] = np.real(ev[:, eig_index]) / N_density

    return result


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


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot eigenvalues of A = R V^{-1} as a function of E/N. "
            "Positive real eigenvalues indicate net spatial growth of the "
            "discharge (ionisation exceeds attachment).  "
            "One subplot is produced per configuration."
        )
    )
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
    args = parser.parse_args()

    mech_name = os.path.basename(args.mechanism)
    raw_dicts = _read_json_configs(args.configs) if args.configs else [{}]
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


if __name__ == "__main__":
    main()
