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
import math
import time
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from incept1d.constants import kB as _kB
from incept1d.fields import add_field_argument, parse_field_spec
from incept1d.growth import solve_voltage, voltage_sweep
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import load_mechanism, read_json_configs
from incept1d.output import write_metadata_header
from incept1d.parallel import parallel_map, physical_cores
from incept1d.solver import (
    inception_det,
    DX_N_MAX_DEFAULT,
    DX_N_MIN_DEFAULT,
    DX_TOL_DEFAULT,
    CRITERIA,
    add_criterion_argument,
    midpoint_propagator,
    magnus2_propagator,
    parse_dx_spec,
    polarities_equivalent,
)

_POLARITIES = (("positive", True), ("negative", False))


#: Relative distance from a root at which det Q is checked for a sign change;
#: not smaller, since the adaptive grid is chosen per evaluation.
_VERIFY_EPS = 1e-3


def _listing(values):
    """'2', or '1, 2.5, 10' for a short list, or '1–10 (5 values)'."""
    if len(values) == 1:
        return f"{values[0]:.4g}"
    if len(values) <= 4:
        return ", ".join(f"{v:.4g}" for v in values)
    return f"{min(values):.4g}–{max(values):.4g} ({len(values)} values)"


def _parse_values(tokens, option, parser):
    """
    Values of --pressure / --distance: numbers, or MIN:MAX:N ranges.

    MIN:MAX:N expands to N log-spaced values from MIN to MAX inclusive.
    Returns the values in the order given, each positive.
    """
    values = []
    for token in tokens:
        parts = token.split(":")
        try:
            if len(parts) == 1:
                values.append(float(parts[0]))
            elif len(parts) == 3:
                lo, hi, n = float(parts[0]), float(parts[1]), int(parts[2])
                if n < 1 or (n == 1 and lo != hi):
                    raise ValueError("N must be >= 2 unless MIN = MAX")
                if lo <= 0.0 or hi <= 0.0:
                    raise ValueError("MIN and MAX must be > 0")
                values.extend(float(v) for v in np.geomspace(lo, hi, n))
            else:
                raise ValueError("expected a number or MIN:MAX:N")
        except ValueError as exc:
            parser.error(f"{option} {token!r}: {exc}")
    if any(v <= 0.0 for v in values):
        parser.error(f"{option} values must be > 0")
    return values


#: At most this many markers on a plotted curve; the line carries the rest.
_MAX_MARKERS = 20


def _markevery(mask):
    """Every k-th point, k chosen so a curve gets at most _MAX_MARKERS markers."""
    return max(1, math.ceil(int(np.count_nonzero(mask)) / _MAX_MARKERS))


def _verdict(check):
    """' det Q(...) = a / b ✓' for a (below, above) pair, '' for None."""
    if check is None:
        return ""
    what, below, above = check
    ok = (
        np.isfinite(below)
        and np.isfinite(above)
        and np.sign(below) * np.sign(above) < 0.0  # the product can underflow
    )
    return f"   det Q({what}·(1∓{_VERIFY_EPS:g})) = {below:+.2e} / {above:+.2e} " + (
        "✓" if ok else "✗"
    )


def _progress(label, i, n, row, seconds, check=None):
    """Print one voltage of one curve as soon as it is solved."""
    V_kV, V_ratio, EN_ref, lam, tau_ns, nu_ion, status = row
    head = (
        f"  [{label}] [{i + 1:{len(str(n))}d}/{n}]  "
        f"V = {V_kV:10.4f} kV ({V_ratio:.3f}×V*)"
    )
    tag = f"  [{status}]" if status != "ok" else ""
    if i == 0:
        body = "λ = 0  (inception)"
    elif np.isfinite(lam):
        body = f"λ = {lam:.4e} s⁻¹   τ = {tau_ns:10.3f} ns"
    else:
        body = "λ = NaN"
    print(f"{head}   {body}{tag}   ({seconds:.1f} s){_verdict(check)}", flush=True)


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
        fh.write(f"# Pressure:    {_listing(args.pressures)} bar\n")
        fh.write(f"# Distance:    {_listing(args.distances)} mm\n")
        if len(args.pressures) == len(args.distances) == 1:
            fh.write(
                f"# pd:          {args.pressures[0] * args.distances[0]:.6g} bar·mm\n"
            )
        else:
            fh.write("# Cases:       every pressure with every distance (see labels)\n")
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
        "--pressure",
        nargs="+",
        default=["1"],
        metavar="P",
        help=(
            "Gas pressure(s) in bar (default: 1).  One or more values; "
            "MIN:MAX:N expands to N log-spaced values, e.g. 1:10:5."
        ),
    )
    parser.add_argument(
        "--distance",
        nargs="+",
        default=None,
        metavar="D",
        help=(
            "Gap length(s) in mm, as for --pressure; every pressure is combined "
            "with every distance.  Required unless the geometry fixes the gap: "
            "a coaxial gap is b - a long and a field line its arc length, and "
            "--distance is then ignored (with a message saying so)."
        ),
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
    add_criterion_argument(parser)
    parser.add_argument(
        "--silent",
        action="store_true",
        default=False,
        help="Print only the result tables, no progress.",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        default=False,
        help=(
            "Check every root with the independent det Q criterion on the same "
            "field and grid: det Q must change sign across E/N* for each "
            "inception voltage and across lambda* for each growth rate.  "
            "Printed with the progress (so not with --silent).  Expensive: det "
            "Q with every photon group explicit takes the compound route."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Worker processes; the voltages are solved concurrently (default: "
            f"the number of physical cores, {physical_cores()} here, or 1 if it "
            "cannot be determined).  --jobs 1 solves them one after another."
        ),
    )
    parser.add_argument(
        "--dx",
        nargs="*",
        default=None,
        metavar="SPEC",
        help=(
            "Adaptive integration stepping: N_min [N_max [tol]].  "
            f"Defaults: N_min={DX_N_MIN_DEFAULT}, N_max={DX_N_MAX_DEFAULT}, "
            f"tol={DX_TOL_DEFAULT:g}."
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

    if args.jobs is None:
        args.jobs = physical_cores()
    if args.jobs < 1:
        parser.error("--jobs must be >= 1")
    pressures = _parse_values(args.pressure, "--pressure", parser)

    # The gap length(s): given, or fixed by the geometry.
    fixed = _field_dist.fixed_gap_length
    if fixed is not None:
        why = "b - a" if _field_dist.field_type == "coaxial" else "its arc length"
        given = (
            _parse_values(args.distance, "--distance", parser) if args.distance else []
        )
        if any(abs(d * 1e-3 / fixed - 1.0) > 1e-9 for d in given):
            shown = " ".join(args.distance)
            print(
                f"Note: --distance {shown} mm is ignored: the "
                f"{_field_dist.field_type} geometry fixes the gap at "
                f"{fixed * 1e3:.6g} mm ({why})."
            )
        distances = [fixed * 1e3]
    elif args.distance is None:
        parser.error(f"--distance is required for --field {_field_dist.field_type}")
    else:
        distances = _parse_values(args.distance, "--distance", parser)
    # Every pressure with every distance.
    cases = [(p, d) for p in pressures for d in distances]
    args.pressures, args.distances = pressures, distances
    mech_name = os.path.basename(args.mechanism)

    raw_dicts = read_json_configs(args.configs) if args.configs else [{}]

    # ── The curves: one per configuration and polarity ─────────────────────
    # A symmetric gap with identical electrodes solves positive polarity only
    # and mirrors it for negative.
    # With several (pressure, distance) cases every label names its case;
    # "config" then groups a configuration at one case (colour in the plot).
    curves = []
    for cfg_dict in raw_dicts:
        mod = load_mechanism(args.mechanism, cfg_dict)
        name = mod.label or "Baseline"
        same = polarities_equivalent(mod, _field_dist)
        for p_bar, d_mm in cases:
            label = (
                name
                if len(cases) == 1
                else f"{name}, p={p_bar:.4g} bar, d={d_mm:.4g} mm"
            )
            for polarity, positive in _POLARITIES:
                curves.append(
                    dict(
                        label=f"{label} ({_field_dist.polarity_label(polarity)})",
                        config=label,
                        polarity=polarity,
                        positive=positive,
                        mod=mod,
                        mirror=(not positive and same),
                        p=p_bar,
                        d=d_mm,
                        pd_m=p_bar * d_mm * 1e-3,  # bar·m
                    )
                )
    solved = [c for c in curves if not c["mirror"]]

    def criterion_for(c):
        return functools.partial(
            CRITERIA[args.criterion],
            field_dist=_field_dist,
            N_min=_N_min,
            N_max=_N_max,
            tol=_tol,
            positive_polarity=c["positive"],
            propagator=_propagator,
        )

    def detq_for(c):
        """det Q on the same field, grid and polarity, for --verify."""
        return functools.partial(
            inception_det,
            field_dist=_field_dist,
            N_min=_N_min,
            N_max=_N_max,
            tol=_tol,
            positive_polarity=c["positive"],
            propagator=_propagator,
        )

    verify = args.verify and not args.silent and args.criterion != "detq"

    def say(*a, **k):
        if not args.silent:
            print(*a, **k)

    # ── Phase 1: the inception voltage of every curve, in parallel ─────────
    if len(cases) == 1:
        ((p_bar, d_mm),) = cases
        say(
            f"\nFinding the inception voltages (p = {p_bar:g} bar, "
            f"d = {d_mm:g} mm, pd = {p_bar * d_mm:.6g} bar·mm):"
        )
    else:
        say(
            f"\nFinding the inception voltages ({len(cases)} cases: "
            f"{len(pressures)} pressure(s) × {len(distances)} distance(s)):"
        )

    def find_star(k, previous=None):
        t_start = time.perf_counter()
        c = solved[k]
        roots = find_all_breakdown_EN(
            c["pd_m"],
            c["mod"],
            c["p"],
            args.T,
            first_only=True,
            det_fn=criterion_for(c),
        )
        seconds = time.perf_counter() - t_start
        EN_star = roots[0] if roots else None
        check = None
        if verify and EN_star is not None:
            detq = detq_for(c)
            m, p, pd = c["mod"], c["p"], c["pd_m"]
            check = (
                "E/N*",
                detq(EN_star * (1.0 - _VERIFY_EPS), pd, m, p, args.T),
                detq(EN_star * (1.0 + _VERIFY_EPS), pd, m, p, args.T),
            )
        return EN_star, seconds, check

    def report_star(k, result):
        if args.silent:
            return
        EN_star, seconds, check = result
        if EN_star is None:
            print(f"  [{solved[k]['label']}]  no inception found  ({seconds:.1f} s)")
            return
        V_star = EN_star * solved[k]["pd_m"] * 1e-16 / (_kB * args.T)
        print(
            f"  [{solved[k]['label']}]  V* = {V_star * 1e-3:.4f} kV   "
            f"E/N* = {EN_star:.2f} Td   ({seconds:.1f} s){_verdict(check)}",
            flush=True,
        )

    stars = parallel_map(find_star, len(solved), args.jobs, on_result=report_star)

    # ── Phase 2: every (curve, voltage) in one pool ─────────────────────────
    sweeps = {}
    for k, (EN_star, *_) in enumerate(stars):
        if EN_star is not None:
            sweeps[k] = voltage_sweep(
                EN_star, solved[k]["pd_m"], args.T, args.n_voltages, args.v_max_factor
            )
    tasks = [(k, i) for k in sweeps for i in range(args.n_voltages)]
    say(
        f"\nComputing λ(V) from V* to {args.v_max_factor:.2g}×V* "
        f"({len(tasks)} voltages on {min(args.jobs, max(len(tasks), 1))} workers):"
    )

    def solve_task(t, previous=None):
        k, i = tasks[t]
        V_star, voltages = sweeps[k]
        c = solved[k]
        row, seconds = solve_voltage(
            voltages[i],
            V_star,
            c["pd_m"],
            c["mod"],
            c["p"],
            args.T,
            _field_dist,
            _N_min,
            _N_max,
            _tol,
            _propagator,
            positive_polarity=c["positive"],
            criterion=CRITERIA[args.criterion],
            at_inception=(i == 0),
        )
        lam, EN_ref = row[3], row[2]
        check = None
        if verify and i > 0 and np.isfinite(lam) and lam > 0.0:
            detq = detq_for(c)
            m, p, pd = c["mod"], c["p"], c["pd_m"]
            check = (
                "λ*",
                detq(EN_ref, pd, m, p, args.T, lam=lam * (1.0 - _VERIFY_EPS)),
                detq(EN_ref, pd, m, p, args.T, lam=lam * (1.0 + _VERIFY_EPS)),
            )
        return row, seconds, check

    def report_task(t, result):
        if args.silent:
            return
        k, i = tasks[t]
        _progress(solved[k]["label"], i, args.n_voltages, *result)

    results = parallel_map(solve_task, len(tasks), args.jobs, on_result=report_task)
    rows_of = {k: [None] * args.n_voltages for k in sweeps}
    for (k, i), (row, *_) in zip(tasks, results):
        rows_of[k][i] = row

    # list of (label, config_label, polarity, rows), in curve order; a mirrored
    # negative polarity takes the rows of its configuration's positive one.
    curve_records = []
    positive_rows = {}
    k = 0
    for c in curves:
        if c["mirror"]:
            rows = positive_rows.get(c["config"])
        else:
            rows = rows_of.get(k)
            if c["positive"]:
                positive_rows[c["config"]] = rows
            k += 1
        if rows is not None:
            curve_records.append((c["label"], c["config"], c["polarity"], rows))

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
    # Room below the panels for a shared legend when there are many curves.
    n_drawn = len({id(rec[3]) for rec in curve_records})
    extra = 0.2 * math.ceil(n_drawn / (2 if n_drawn > 6 else 1)) if n_drawn > 4 else 0.0
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5 + extra))
    fig.suptitle(
        f"Growth rate vs. voltage  —  {mech_name},  "
        f"p = {_listing(args.pressures)} bar,  d = {_listing(args.distances)} mm,  "
        f"T = {args.T} K,  "
        f"{_field_dist.label}",
        fontsize=12,
    )

    _markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    # One colour and marker per configuration (and case); linestyle marks the
    # polarity.  A mirrored negative polarity (symmetric gap, identical
    # electrodes) has the very same rows as its positive one: draw it once.
    plotted = []  # (label, config, polarity, rows) actually drawn
    for label, config, polarity, rows in curve_records:
        twin = next((k for k, rec in enumerate(plotted) if rec[3] is rows), None)
        if twin is not None:
            name = plotted[twin][0].rsplit(" (", 1)[0]
            plotted[twin] = (f"{name} (both polarities)",) + plotted[twin][1:]
            continue
        plotted.append((label, config, polarity, rows))

    _config_idx = {}
    all_V = []
    for label, config, polarity, rows in plotted:
        idx = _config_idx.setdefault(config, len(_config_idx))
        V_arr = np.array([r[0] for r in rows])  # kV
        lam_arr = np.array([r[3] for r in rows])
        tau_arr = np.array([r[4] for r in rows])
        all_V.extend(V_arr)
        kw = dict(
            color=f"C{idx % 10}",
            marker=_markers[idx % len(_markers)],
            markersize=5,
            ls="-" if polarity == "positive" else "--",
            label=label,
        )

        mask = np.isfinite(lam_arr) & (lam_arr > 0.0)  # λ = 0 at V* has no log
        ax1.semilogy(V_arr[mask], lam_arr[mask], markevery=_markevery(mask), **kw)
        mask2 = np.isfinite(tau_arr)
        ax2.semilogy(V_arr[mask2], tau_arr[mask2], markevery=_markevery(mask2), **kw)

    # Cases at very different voltages: a log axis gives each of them room.
    if all_V and max(all_V) / min(all_V) > 5.0:
        for ax in (ax1, ax2):
            ax.set_xscale("log")
            # Plain numbers at 1, 2 and 5 per decade, not 2×10^1 3×10^1 ...
            ax.xaxis.set_major_locator(
                mticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0))
            )
            ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
            ax.xaxis.set_minor_formatter(mticker.NullFormatter())

    ax1.set_xlabel("Voltage (kV)")
    ax1.set_ylabel(r"$\lambda$  (s$^{-1}$)")
    ax1.set_title("Temporal growth rate")
    ax1.grid(True, which="both", ls="--", alpha=0.4)

    ax2.set_xlabel("Voltage (kV)")
    ax2.set_ylabel(r"$\tau = 1/\lambda$  (ns)")
    ax2.set_title("e-folding time")
    ax2.grid(True, which="both", ls="--", alpha=0.4)

    if len(plotted) > 4:
        # Many curves: one legend below the panels, not over the data.
        handles, labels = ax1.get_legend_handles_labels()
        ncol = 2 if len(plotted) > 6 else 1
        n_rows = math.ceil(len(plotted) / ncol)
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=ncol,
            framealpha=1.0,
            fontsize="small",
        )
        bottom = min(0.5, (0.25 + 0.19 * n_rows) / fig.get_size_inches()[1])
        plt.tight_layout(rect=(0.0, bottom, 1.0, 1.0))
    else:
        if len(plotted) > 1:
            ax1.legend(loc="best", framealpha=1.0)
            ax2.legend(loc="best", framealpha=1.0)
        plt.tight_layout()
    plt.show()
