# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
``incept1d growth`` — temporal growth rate λ vs. voltage above inception.

Finds the inception voltage V* from det Q(λ=0) = 0, then solves
det Q(λ, E/N) = 0 for λ at each voltage V ∈ [V*, F·V*].  See
:mod:`incept1d.growth` for the solver.
"""

import functools
import os

import numpy as np
import matplotlib.pyplot as plt

from incept1d.constants import kB as _kB
from incept1d.fields import add_field_argument, parse_field_spec
from incept1d.growth import compute_lambda_curve
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import load_mechanism, read_json_configs
from incept1d.output import write_metadata_header
from incept1d.solver import (
    inception_det,
    midpoint_propagator,
    magnus2_propagator,
    parse_dx_spec,
    polarities_equivalent,
)

_POLARITIES = (("positive", True), ("negative", False))


def _write_results(out_path, curve_records, args, raw_dicts, mech_name, field_dist):
    """Write λ vs V results to a tab-separated file with metadata header."""
    SEP = "\t"
    # Build column headers from all curve labels
    col_headers = ["V_ratio"]
    for label, *_ in curve_records:
        col_headers += [
            f"V_kV[{label}]",
            f"lambda_s-1[{label}]",
            f"tau_ns[{label}]",
            f"nu_ion_s-1[{label}]",
        ]
    W = max(14, max(len(h) for h in col_headers) + 2)

    with open(out_path, "w") as fh:
        write_metadata_header(fh, extra_paths=[args.mechanism])
        fh.write(f"# Mechanism:   {mech_name}\n")
        fh.write(f"# pd:          {args.pd} bar·mm\n")
        fh.write(f"# Pressure:    {args.p} bar\n")
        fh.write(f"# Temperature: {args.T} K\n")
        fh.write(f"# Field type:  {field_dist.label}\n")
        fh.write(f"# V_max_factor:{args.v_max_factor}\n")
        fh.write(f"# n_voltages:  {args.n_voltages}\n")
        fh.write(
            f'# Configs:     {", ".join(d.get("label", "Baseline") for d in raw_dicts)}\n'
        )
        fh.write("#\n")
        for i, h in enumerate(col_headers, start=1):
            fh.write(f"# Column {i}: {h}\n")
        fh.write("#\n")

        # Every curve is sampled on the same V/V* grid, but V* differs between
        # configurations and polarities, so each curve carries its own V_kV.
        first_rows = curve_records[0][-1]
        for j in range(len(first_rows)):
            row_vals = [first_rows[j][1]]  # V_ratio
            for *_, rows in curve_records:
                V_kV, _, _, lam, tau, nu, _ = rows[j]
                row_vals += [V_kV, lam, tau, nu]
            fh.write(SEP.join(f"{v:<{W}.6e}" for v in row_vals) + "\n")

    print(f"\nResults written to: {out_path}")


# ── Entry point ───────────────────────────────────────────────────────────────


HELP = "Temporal growth rate λ vs. voltage above inception."
DESCRIPTION = (
    "Compute temporal growth rate λ vs. voltage above inception.  "
    "Finds V* from det Q(λ=0) = 0, then solves det Q(λ, EN) = 0 for "
    "λ at each voltage V ∈ [V*, F·V*]."
)


def add_arguments(parser):
    """Register the command-line arguments on *parser*."""
    parser.add_argument(
        "mechanism",
        help="Path to mechanism Python file (e.g. mechanisms/air/pancheshnyi/air_pancheshnyi.py).",
    )
    parser.add_argument(
        "configs",
        nargs="*",
        metavar="CONFIG.json",
        help=(
            "One or more JSON configuration files.  Each may contain a single "
            'object or a list under "configurations".  If omitted, a single '
            "baseline configuration is used."
        ),
    )
    parser.add_argument(
        "--pd",
        type=float,
        required=True,
        metavar="PD",
        help="Pressure × gap product in bar·mm.",
    )
    parser.add_argument(
        "--p",
        type=float,
        default=1.0,
        metavar="P",
        help="Gas pressure in bar (default: 1.0).",
    )
    parser.add_argument(
        "--T",
        type=float,
        default=293.0,
        help="Gas temperature in Kelvin (default: 293.0).",
    )
    parser.add_argument(
        "--n-voltages",
        type=int,
        default=20,
        metavar="N",
        help="Number of logarithmically spaced voltage points (default: 20).",
    )
    parser.add_argument(
        "--v-max-factor",
        type=float,
        default=2.0,
        metavar="F",
        help="Upper voltage as multiple of V* (default: 2.0).",
    )
    add_field_argument(parser)
    parser.add_argument(
        "--dx",
        nargs="*",
        default=None,
        metavar="SPEC",
        help=(
            "Adaptive integration stepping: N_min [N_max [tol]].  "
            "Defaults: N_min=5, N_max=200, tol=0.03."
        ),
    )
    parser.add_argument(
        "--method",
        choices=["midpoint", "magnus2"],
        default="midpoint",
        help="Propagator algorithm (default: 'midpoint').",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        default=False,
        help="Suppress the matplotlib figure.",
    )
    parser.add_argument(
        "--write-to-file",
        type=str,
        default=None,
        metavar="FILE",
        help="Write tab-separated results to FILE.",
    )


def run(args, parser):
    """Run the command with parsed *args*; *parser* is used for ``parser.error``."""
    if args.v_max_factor <= 1.0:
        parser.error("--v-max-factor must be > 1.0")
    if args.n_voltages < 2:
        parser.error("--n-voltages must be >= 2")

    _field_dist = parse_field_spec(
        args.field, parser, applied_voltage_kv=args.fieldline_voltage
    )
    _N_min, _N_max, _tol = parse_dx_spec(args.dx, parser)
    _propagator = (
        midpoint_propagator if args.method == "midpoint" else magnus2_propagator
    )

    pd_m = args.pd * 1e-3  # bar·mm → bar·m
    mech_name = os.path.basename(args.mechanism)

    raw_dicts = read_json_configs(args.configs) if args.configs else [{}]

    # ── Find inception and compute λ curve for each config and polarity ────
    # list of (label, config_label, polarity, rows)
    # where label names the curve and rows = compute_lambda_curve(...)
    curve_records = []

    for cfg_dict in raw_dicts:
        mod = load_mechanism(args.mechanism, cfg_dict)
        label = mod.label or "Baseline"
        same = polarities_equivalent(mod, _field_dist)
        positive_rows = None

        for polarity, positive in _POLARITIES:
            pol_label = f"{label} ({_field_dist.polarity_label(polarity)})"
            if not positive and same:
                # Symmetric gap, identical electrodes: negative mirrors positive.
                if positive_rows is not None:
                    curve_records.append((pol_label, label, polarity, positive_rows))
                continue

            # Build det_fn for the inception solve (λ=0).
            det_fn = functools.partial(
                inception_det,
                field_dist=_field_dist,
                N_min=_N_min,
                N_max=_N_max,
                tol=_tol,
                lam=0.0,
                positive_polarity=positive,
                propagator=_propagator,
            )

            print(
                f"\nFinding inception voltage: {pol_label}  "
                f"(pd = {args.pd} bar·mm, p = {args.p} bar)"
            )
            roots = find_all_breakdown_EN(
                pd_m, mod, args.p, args.T, first_only=True, det_fn=det_fn
            )
            if not roots:
                print(f"  [{pol_label}] No inception found — skipping.")
                continue

            EN_star = roots[0]
            V_star = EN_star * pd_m * 1e-16 / (_kB * args.T)
            print(f"  V* = {V_star * 1e-3:.4f} kV   EN* = {EN_star:.2f} Td")
            print(f"\nComputing λ(V) from V* to {args.v_max_factor:.2g}×V*:")

            rows = compute_lambda_curve(
                EN_star,
                pd_m,
                mod,
                args.p,
                args.T,
                _field_dist,
                _N_min,
                _N_max,
                _tol,
                _propagator,
                args.n_voltages,
                args.v_max_factor,
                positive_polarity=positive,
            )
            curve_records.append((pol_label, label, polarity, rows))
            if positive:
                positive_rows = rows

    if not curve_records:
        print("\nNo curves computed.")
        return

    # ── Console summary table ───────────────────────────────────────────────
    print()
    for label, _, _, rows in curve_records:
        print(f"  {label}")
        print("  " + "-" * 78)
        print(
            f"  {'V (kV)':>10}  {'V/V*':>7}  {'EN (Td)':>10}  "
            f"{'λ (s⁻¹)':>14}  {'τ (ns)':>10}  {'λ/ν_ion':>10}"
        )
        for V_kV, V_ratio, EN_ref, lam, tau, nu, status in rows:
            lam_str = f"{lam:.4e}" if np.isfinite(lam) else "     NaN"
            tau_str = f"{tau:.4f}" if np.isfinite(tau) else "     NaN"
            ratio_str = f"{lam/nu:.4f}" if (np.isfinite(lam) and nu > 0) else "     NaN"
            flag = f"  [{status}]" if status != "ok" else ""
            print(
                f"  {V_kV:>10.4f}  {V_ratio:>7.4f}  {EN_ref:>10.2f}  "
                f"{lam_str:>14}  {tau_str:>10}  {ratio_str:>10}{flag}"
            )
        print()

    # ── File output ─────────────────────────────────────────────────────────
    if args.write_to_file:
        _write_results(
            args.write_to_file, curve_records, args, raw_dicts, mech_name, _field_dist
        )

    # ── Plot ─────────────────────────────────────────────────────────────────
    if args.no_plot:
        return

    plt.rcParams["font.size"] += 2
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"Growth rate vs. voltage  —  {mech_name},  "
        f"pd = {args.pd} bar·mm,  p = {args.p} bar,  T = {args.T} K,  "
        f"{_field_dist.label}",
        fontsize=12,
    )

    _markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    # One colour and marker per configuration; linestyle marks the polarity.
    _config_idx = {}
    for label, config, polarity, rows in curve_records:
        idx = _config_idx.setdefault(config, len(_config_idx))
        V_ratio_arr = np.array([r[1] for r in rows])
        lam_arr = np.array([r[3] for r in rows])
        tau_arr = np.array([r[4] for r in rows])
        kw = dict(
            color=f"C{idx}",
            marker=_markers[idx % len(_markers)],
            markersize=5,
            ls="-" if polarity == "positive" else "--",
            label=label,
        )

        mask = np.isfinite(lam_arr)
        ax1.semilogy(V_ratio_arr[mask], lam_arr[mask], **kw)
        mask2 = np.isfinite(tau_arr)
        ax2.semilogy(V_ratio_arr[mask2], tau_arr[mask2], **kw)

    ax1.set_xlabel(r"$V / V^*$")
    ax1.set_ylabel(r"$\lambda$  (s$^{-1}$)")
    ax1.set_title("Temporal growth rate")
    ax1.grid(True, which="both", ls="--", alpha=0.4)

    ax2.set_xlabel(r"$V / V^*$")
    ax2.set_ylabel(r"$\tau = 1/\lambda$  (ns)")
    ax2.set_title("e-folding time")
    ax2.grid(True, which="both", ls="--", alpha=0.4)

    if len(curve_records) > 1:
        ax1.legend(loc="best", framealpha=1.0)
        ax2.legend(loc="best", framealpha=1.0)

    plt.tight_layout()
    plt.show()
