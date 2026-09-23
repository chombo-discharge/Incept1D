# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
``incept1d pdiv`` — compute and plot the inception curve det Q(λ=0) = 0,
i.e. the (partial) discharge inception voltage PDIV = V*(p·d), over a p·d
sweep.

Usage
-----
    incept1d pdiv <mechanism> [CONFIG.json ...] [--p P [P ...]] [--d D [D ...]]
                     [--T T] [--field SPEC] [--dx SPEC] [--write-to-file FILE] ...

One curve is produced for every (pressure or distance, configuration)
combination; both polarities are solved for asymmetric field geometries.
"""

import functools
import os

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as _mticker

from incept1d.constants import kB as _kB
from incept1d.fields import (
    FieldDistribution,
    add_field_argument,
    parse_field_spec,
)
from incept1d.mechanism import load_mechanism, read_json_configs
from incept1d.solver import (
    inception_det,
    midpoint_propagator,
    magnus2_propagator,
    parse_dx_spec,
)
from incept1d.inception import compute_inception_curve
from incept1d.output import write_metadata_header

HELP = "Inception curve PDIV(p·d): roots of det Q(E/N, p·d) = 0 over a p·d sweep."
DESCRIPTION = (
    "Compute and plot the inception curve det Q(E/N, pd) = 0 (PDIV vs. pd).  "
    "Sweeps pd over [--pd-min, --pd-max] in two modes: fixed pressure "
    "(--p) or fixed gap distance (--d).  Both may be combined."
)


def add_arguments(parser):
    """Register the ``pdiv`` command-line arguments on *parser*."""
    parser.add_argument(
        "mechanism",
        help="Path to mechanism Python file (e.g. mechanisms/air/air_pancheshnyi.py).",
    )
    parser.add_argument(
        "configs",
        nargs="*",
        metavar="CONFIG.json",
        help=(
            "One or more JSON configuration files for the mechanism.  "
            "Each file may contain a single configuration object or a list "
            'of objects under a "configurations" key.  All configurations '
            "across all files are run in sequence.  If omitted, a single "
            "baseline configuration with default parameters is used."
        ),
    )
    parser.add_argument(
        "--p",
        type=float,
        nargs="+",
        default=[],
        metavar="P",
        help="Fixed-p mode: pressure(s) in bar (default: 1.0 when --d is not given).",
    )
    parser.add_argument(
        "--d",
        type=float,
        nargs="+",
        default=[],
        metavar="D",
        help="Fixed-d mode: gap distance(s) in mm.",
    )
    parser.add_argument(
        "--pd-min",
        type=float,
        default=1e-2,
        metavar="PD_MIN",
        help="Minimum p*d in bar·mm (default: 1e-2).",
    )
    parser.add_argument(
        "--pd-max",
        type=float,
        default=1e3,
        metavar="PD_MAX",
        help="Maximum p*d in bar·mm (default: 1e3).",
    )
    parser.add_argument(
        "--pd-num",
        type=int,
        default=50,
        metavar="PD_NUM",
        help="Number of logarithmically spaced p*d grid points (default: 50).",
    )
    parser.add_argument(
        "--T",
        type=float,
        default=293.0,
        help="Gas temperature in Kelvin (default: 293.0).",
    )
    parser.add_argument(
        "--write-to-file",
        type=str,
        default=None,
        metavar="FILE",
        help="Write all curves to a tab-separated file.",
    )
    parser.add_argument(
        "--save-subplots",
        action="store_true",
        default=False,
        help="Save each subplot as a separate PDF (no titles).",
    )
    parser.add_argument(
        "--all-branches",
        action="store_true",
        default=False,
        help=(
            "Compute and plot all branches of the breakdown curve.  "
            "By default only branch 1 (lowest E/N root) is computed."
        ),
    )
    add_field_argument(parser)
    parser.add_argument(
        "--dx",
        nargs="*",
        default=None,
        metavar="SPEC",
        help=(
            "Adaptive integration stepping: N_min [N_max [tol]].  "
            "N_min (default 5): minimum number of integration segments.  "
            "N_max (default 200): maximum total fine steps; N_max = N_min gives "
            "a constant uniform grid with no adaptive refinement.  "
            "tol (default 0.03): relative Frobenius error threshold for the "
            "midpoint step-halving check (fraction, not percent).  "
            "Example: --dx 10 400 0.01  sets N_min=10, N_max=400, tol=1%%."
        ),
    )
    parser.add_argument(
        "--lam",
        type=float,
        default=0.0,
        metavar="LAM",
        help=(
            "Temporal growth rate λ in s⁻¹ for the generalised inception criterion "
            "det Q(λ) = 0 (default: 0.0 = standard inception threshold).  "
            "λ > 0 → growing discharge (lower breakdown voltage); "
            "λ < 0 → decaying discharge (higher breakdown voltage)."
        ),
    )
    parser.add_argument(
        "--plot-separate-branches",
        action="store_true",
        default=False,
        help=(
            "Give each branch its own line style (solid / dashed / dotted / "
            "dash-dot) and legend entry ('branch 1', 'branch 2', …).  "
            "By default all branches share the same style and a single legend entry."
        ),
    )
    parser.add_argument(
        "--plot-ionization-integral",
        action="store_true",
        default=False,
        help=(
            "Overlay the ionization integral ∫max(α−η,0)dx on a second y-axis "
            "in the voltage subplot.  Only segments where α > η contribute.  "
            "Requires the mechanism to expose both alpha and eta."
        ),
    )
    parser.add_argument(
        "--streamer-criterion",
        dest="streamer_criterion",
        type=float,
        default=None,
        metavar="C",
        help=(
            "Solve for the streamer criterion ∫max(α−η,0)dx = C.  "
            "C must be > 0.  Plots the streamer curve on both panels and "
            "adds it to --write-to-file output.  "
            "Requires the mechanism to expose alpha and eta."
        ),
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        default=False,
        help="Skip the matplotlib figure entirely (useful for batch/scripted runs).",
    )
    parser.add_argument(
        "--method",
        choices=["midpoint", "magnus2"],
        default="midpoint",
        help=(
            "Propagator algorithm for the path-ordered matrix exponential.  "
            "'midpoint' (default): zeroth-order Magnus / midpoint rule.  "
            "'magnus2': second-order Magnus with 2-point Gauss-Legendre quadrature; "
            "reduces to midpoint for uniform fields."
        ),
    )


def run(args, parser):
    """
    Solve for breakdown E/N across a p*d sweep and display a two-panel
    inception curve figure.

    The p*d sweep range and resolution are controlled by --pd-min, --pd-max,
    and --pd-num (defaults: 1e-2 to 1e3 bar·mm, 50 points).

    Panel 1 — V* vs p*d (log-log):
        Shows the inception voltage.  The classical Paschen minimum appears
        as the lowest point on each curve.

    Panel 2 — E/N* vs p*d (log-log):
        Shows the critical reduced electric field at breakdown for each curve.

    One curve is produced for every (pressure, configuration) combination.

    A summary table is printed to stdout for each (p, tag, pd, E/N*, V*) triple
    where a solution was found.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments (see :func:`add_arguments`).
    parser : argparse.ArgumentParser
        The (sub)parser that produced *args*; used for ``parser.error``.
    """
    if args.streamer_criterion is not None and args.streamer_criterion <= 0:
        parser.error("--streamer-criterion value must be > 0")

    _field_dist = parse_field_spec(
        args.field, parser, applied_voltage_kv=args.fieldline_voltage
    )

    if not args.p and not args.d:
        if _field_dist.field_type == "fieldline":
            # A tabulated field line has an intrinsic length: sweep p at d = L.
            args.d = [float(f"{_field_dist.fieldline_length * 1e3:.6g}")]
            print(
                f"--field fieldline: no --p/--d given, using d = L = "
                f"{args.d[0]:.4g} mm (arc length of the tabulated line)."
            )
        else:
            args.p = [1.0]
    elif _field_dist.field_type == "fieldline" and args.p:
        print(
            "Note: --field fieldline with fixed --p sweeps the gap length d; "
            "d ≠ L corresponds to the same electrode arrangement scaled "
            "geometrically by d/L (f(xi) unchanged)."
        )
    _N_min, _N_max, _dx_tol = parse_dx_spec(args.dx, parser)

    _propagator = (
        midpoint_propagator if args.method == "midpoint" else magnus2_propagator
    )

    mech_name = os.path.basename(args.mechanism)

    # Read JSON config files as raw dicts; no gas-specific class needed here.
    raw_dicts = read_json_configs(args.configs) if args.configs else [{}]

    # Determine n (number of species) from the first config.
    _mod0 = load_mechanism(args.mechanism, raw_dicts[0])
    n = len(_mod0.SPECIES)
    del _mod0

    if args.streamer_criterion is not None:
        _mod_check = load_mechanism(args.mechanism, raw_dicts[0])
        if not (hasattr(_mod_check, "alpha") and hasattr(_mod_check, "eta")):
            parser.error(
                "--streamer-criterion requires the mechanism to expose alpha and eta"
            )
        del _mod_check

    # Fixed pd sweep
    pd_arr = np.logspace(
        np.log10(args.pd_min * 1e-3), np.log10(args.pd_max * 1e-3), args.pd_num
    )

    # Build curve specs: (label, p_arr, d_arr, mod)
    # Each config dict gets its own freshly loaded Mechanism instance.
    curve_specs = []
    for cfg_dict in raw_dicts:
        mod = load_mechanism(args.mechanism, cfg_dict)
        cfg_label = mod.label
        if len(mod.SPECIES) != n:
            raise ValueError(
                f"Config '{cfg_label}' has {len(mod.SPECIES)} species; "
                f"expected {n} (from first config)."
            )
        for p in args.p:
            curve_specs.append(
                (
                    f"{cfg_label}, p={p} bar",
                    np.full_like(pd_arr, p),
                    pd_arr / p,
                    mod,
                )
            )
        for d_mm in args.d:
            d_m = d_mm * 1e-3
            curve_specs.append(
                (
                    f"{cfg_label}, d={d_mm} mm",
                    pd_arr / d_m,
                    np.full_like(pd_arr, d_m),
                    mod,
                )
            )

    # Compute alpha=eta crossover E/N from the first config's module (for annotation).
    _en_cross = {}
    if curve_specs:
        _ref_mod = curve_specs[0][3]
        if hasattr(_ref_mod, "alpha") and hasattr(_ref_mod, "eta"):
            EN_scan = np.logspace(np.log10(10.0), np.log10(1e5), 200)
            for p in args.p:

                def f_cross(EN, _p=p):
                    return _ref_mod.alpha(EN, _p, args.T) - _ref_mod.eta(EN, _p, args.T)

                f_vals = np.array([f_cross(en) for en in EN_scan])
                sign_changes = np.where(f_vals[:-1] * f_vals[1:] < 0)[0]
                if len(sign_changes):
                    i = sign_changes[0]
                    _en_cross[p] = scipy.optimize.brentq(
                        f_cross, EN_scan[i], EN_scan[i + 1], xtol=1e-3, rtol=1e-5
                    )
                else:
                    _en_cross[p] = None

    _first_mod = curve_specs[0][3] if curve_specs else None
    _plot_aed = (
        args.plot_ionization_integral
        and _first_mod is not None
        and hasattr(_first_mod, "alpha")
        and hasattr(_first_mod, "eta")
    )

    def _aed_integral(EN_ref, p_val, d_val, _mod):
        """Compute ∫max(α−η,0)dx respecting the actual field profile."""
        if _field_dist.field_type == "uniform":
            return (
                max(
                    0.0,
                    _mod.alpha(EN_ref, p_val, args.T) - _mod.eta(EN_ref, p_val, args.T),
                )
                * d_val
            )
        f = _field_dist.build(d_val)
        xis = (np.arange(_N_min) + 0.5) / _N_min
        EN_arr = EN_ref * np.array([f(xi) for xi in xis])
        ds = d_val / _N_min
        diff = np.array(
            [
                _mod.alpha(en, p_val, args.T) - _mod.eta(en, p_val, args.T)
                for en in EN_arr
            ]
        )
        return float(np.sum(np.maximum(0.0, diff)) * ds)

    def _polarity_desc(polarity):
        if _field_dist.field_type == "uniform":
            return polarity  # "positive" / "negative"
        if _field_dist.field_type == "fieldline":
            return f"start={polarity}"  # xi = 0 (first data row) is anode/cathode
        return f"sphere={polarity}"  # "sphere=positive" / "sphere=negative"

    def _branch0_grids(branches):
        V_g = np.full_like(pd_arr, np.nan)
        EN_g = np.full_like(pd_arr, np.nan)
        p_g = np.full_like(pd_arr, np.nan)
        d_g = np.full_like(pd_arr, np.nan)
        if branches:
            br0 = branches[0]
            V_g[br0["idx"]] = br0["V"]
            EN_g[br0["idx"]] = br0["EN"]
            p_g[br0["idx"]] = br0["p"]
            d_g[br0["idx"]] = br0["d"]
        return V_g, EN_g, p_g, d_g

    # Collects (label, p_arr, d_arr, EN_star, V_star) for every solved curve.
    _file_records = []
    _streamer_records = []  # (label, p_arr, d_arr, EN_arr, V_arr) per streamer setting
    _all_V_pos = {}  # label -> full 200-point V array, positive polarity
    _all_V_neg = {}  # label -> full 200-point V array, negative polarity
    _curve_color = (
        {}
    )  # label -> matplotlib line colour (for consistent colouring on ax3)

    if not args.no_plot:
        plt.rcParams["font.size"] += 2

        # Set up figure — panels: V*, E/N*, optionally modifier ratio, optionally field-error %
        _show_ratio = len(raw_dicts) > 1
        _n_panels = 2 + int(_show_ratio)
        fig, _axes = plt.subplots(1, _n_panels, figsize=(7 * _n_panels, 6))
        ax1, ax2 = _axes[0], _axes[1]
        _next_ax = 2
        if _show_ratio:
            ax3 = _axes[_next_ax]
            _next_ax += 1
        else:
            ax3 = None
        _field_str = _field_dist.label
        _method_str = "" if args.method == "midpoint" else f",  {args.method}"
        _suptitle = fig.suptitle(
            f"Inception curve  —  {mech_name},  T = {args.T} K,  {_field_str}{_method_str}",
            fontsize=13,
        )

        _AED_ALPHA = 0.4
        ax1b = ax1.twinx() if _plot_aed else None
        if ax1b is not None:
            ax1b.set_yscale("symlog", linthresh=1)
            ax1b.yaxis.set_major_locator(
                _mticker.SymmetricalLogLocator(
                    linthresh=1, base=10, subs=[1.0, 2.0, 5.0]
                )
            )
            ax1b.yaxis.set_major_formatter(
                _mticker.LogFormatter(minor_thresholds=(np.inf, np.inf))
            )
            ax1b.set_ylabel(
                r"$\int_0^d \max(\alpha-\eta,\,0)\,\mathrm{d}x$", alpha=_AED_ALPHA
            )
            ax1b.spines["right"].set_alpha(_AED_ALPHA)
            ax1b.tick_params(axis="y", colors=(0, 0, 0, _AED_ALPHA))

            def _aed_format_coord(x, y, _a1=ax1, _a1b=ax1b):
                _, y1 = _a1.transData.inverted().transform(
                    _a1b.transData.transform((x, y))
                )
                return f"pd = {x:.3g} bar·mm    U = {y1:.4g} kV    (α−η)d = {y:.4g}"

            ax1b.format_coord = _aed_format_coord

        _markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
        _BRANCH_LS = ["-", "--", ":", "-."]
    else:
        ax1 = ax2 = ax3 = ax1b = None
        _markers = _BRANCH_LS = None

    for curve_idx, (label, p_arr, _d_arr, mod) in enumerate(curve_specs):

        if _field_dist.field_type == "sphere-sphere":
            max_dR = np.max(_d_arr) / _field_dist.sphere_R
            if max_dR > 4.0:
                print(
                    f"  Warning [{label}]: max d/R = {max_dR:.3g} > 4 — "
                    f"sphere-sphere field approximation may be inaccurate."
                )

        # Determine det functions for both polarities.
        # fast_det_fn  — uniform field, always N=1 (coarse sign-change scan)
        # med_det_fn   — N=min(5,N), modest resolution (fallback scan)
        # det_fn       — full field_dist accuracy (Brentq)
        if _field_dist.field_type == "uniform":
            det_pos = det_neg = functools.partial(
                inception_det,
                field_dist=_field_dist,
                N_min=_N_min,
                N_max=_N_max,
                tol=_dx_tol,
                lam=args.lam,
                positive_polarity=True,
                propagator=_propagator,
            )
            fast_det_pos = fast_det_neg = None
            med_det_pos = med_det_neg = None
        else:
            _fast_fd = FieldDistribution("uniform")
            _med_fd = _field_dist  # same geometry; only the step count differs
            _N_med = min(
                5, _N_min
            )  # cheap constant scan; N_med = N_med disables adaptation
            det_pos = functools.partial(
                inception_det,
                field_dist=_field_dist,
                N_min=_N_min,
                N_max=_N_max,
                tol=_dx_tol,
                lam=args.lam,
                positive_polarity=True,
                propagator=_propagator,
            )
            det_neg = (
                det_pos
                if _field_dist.is_symmetric
                else functools.partial(
                    inception_det,
                    field_dist=_field_dist,
                    N_min=_N_min,
                    N_max=_N_max,
                    tol=_dx_tol,
                    lam=args.lam,
                    positive_polarity=False,
                    propagator=_propagator,
                )
            )
            fast_det_pos = fast_det_neg = functools.partial(
                inception_det,
                field_dist=_fast_fd,
                N_min=1,
                N_max=1,
                tol=_dx_tol,
                lam=args.lam,
                propagator=_propagator,
            )
            med_det_pos = functools.partial(
                inception_det,
                field_dist=_med_fd,
                N_min=_N_med,
                N_max=_N_med,
                tol=_dx_tol,
                lam=args.lam,
                positive_polarity=True,
                propagator=_propagator,
            )
            med_det_neg = (
                med_det_pos
                if _field_dist.is_symmetric
                else functools.partial(
                    inception_det,
                    field_dist=_med_fd,
                    N_min=_N_med,
                    N_max=_N_med,
                    tol=_dx_tol,
                    lam=args.lam,
                    positive_polarity=False,
                    propagator=_propagator,
                )
            )

        print(f"\nSolving inception curve: {label} ({_polarity_desc('positive')})")
        branches_pos = compute_inception_curve(
            pd_arr,
            mod,
            p_arr,
            args.T,
            all_branches=args.all_branches,
            det_fn=det_pos,
            fast_det_fn=fast_det_pos,
            med_det_fn=med_det_pos,
        )
        if _field_dist.is_symmetric:
            branches_neg = branches_pos  # symmetric; reuse same object
        else:
            print(f"\nSolving inception curve: {label} ({_polarity_desc('negative')})")
            branches_neg = compute_inception_curve(
                pd_arr,
                mod,
                p_arr,
                args.T,
                all_branches=args.all_branches,
                det_fn=det_neg,
                fast_det_fn=fast_det_neg,
                med_det_fn=med_det_neg,
            )

        if not branches_pos and not branches_neg:
            print(f"  [{label}] No breakdown found for any pd value.")
            _all_V_pos[label] = np.full_like(pd_arr, np.nan)
            _all_V_neg[label] = np.full_like(pd_arr, np.nan)
            continue

        polarity_pairs = [
            ("positive", branches_pos, "-"),
            ("negative", branches_neg, "--"),
        ]
        for polarity, branches, ls_pol in polarity_pairs:
            if not branches:
                continue
            pol_label = f"{label} ({_polarity_desc(polarity)})"
            for b_idx, br in enumerate(branches):
                order = np.argsort(br["pd"])
                br_pd = br["pd"][order]
                br_EN = br["EN"][order]
                br_V = br["V"][order]
                br_p = br["p"][order]
                br_d = br["d"][order]

                if not args.no_plot:
                    markevery = [0, len(br_pd) - 1]
                    if args.plot_separate_branches:
                        ls = _BRANCH_LS[b_idx % len(_BRANCH_LS)]
                        b_label = f"{pol_label} (branch {b_idx + 1})"
                    else:
                        ls = ls_pol
                        b_label = pol_label if b_idx == 0 else "_nolegend_"
                    is_first = polarity == "positive" and b_idx == 0
                    color_kw = {} if is_first else {"color": _curve_color[label]}
                    plot_kw = dict(
                        label=b_label,
                        marker=_markers[curve_idx % len(_markers)],
                        markevery=markevery,
                        markersize=6,
                        ls=ls,
                        **color_kw,
                    )
                    ax1.loglog(br_pd * 1e3, br_V / 1000, **plot_kw)
                    if is_first:
                        _curve_color[label] = ax1.get_lines()[-1].get_color()
                    col = _curve_color[label]
                    if ax1b is not None:
                        aed = np.array(
                            [
                                (
                                    _aed_integral(en, p, d, mod)
                                    if np.isfinite(en)
                                    else np.nan
                                )
                                for en, p, d in zip(br_EN, br_p, br_d)
                            ]
                        )
                        ax1b.plot(
                            br_pd * 1e3,
                            aed,
                            color=col,
                            alpha=_AED_ALPHA,
                            ls=ls_pol,
                            lw=1.5,
                            zorder=1,
                        )
                    ax2.loglog(
                        br_pd * 1e3,
                        br_EN,
                        **{**plot_kw, "color": col, "label": b_label},
                    )

                # Print summary table for this branch.  Only the shape of a
                # tabulated field line enters the solve, so the useful extra
                # quantity is how far the declared excitation is from
                # inception: U* / U_applied.  Without --fieldline-voltage
                # there is nothing to divide by -- the field units of the
                # file are unknown -- and the column is omitted.
                _vfile = _field_dist.fieldline_applied_voltage
                _scale_hdr = f"  {'U*/U_applied':>14}" if _vfile else ""
                print(
                    f"\n  {'pd (bar·mm)':>14}  {'p (bar)':>10}  {'d (mm)':>8}  "
                    f"{'E/N (Td)':>12}  {'U (kV)':>12}  {'E (V/m)':>14}{_scale_hdr}  "
                    f"[{pol_label} (branch {b_idx + 1})]"
                )
                print("  " + "-" * (90 + (16 if _vfile else 0)))
                step = max(1, len(br_pd) // 20)
                for j in range(0, len(br_pd), step):
                    _scale = f"  {br_V[j]/_vfile:>14.4f}" if _vfile else ""
                    print(
                        f"  {br_pd[j]*1e3:>14.4e}  "
                        f"{br_p[j]:>10.4g}  "
                        f"{br_d[j]*1e3:>8.4g}  "
                        f"{br_EN[j]:>12.4f}  "
                        f"{br_V[j]/1000:>12.4f}  "
                        f"{br_V[j]/br_d[j]:>14.4e}{_scale}"
                    )

        # Map branch 0 onto the full pd grid for file output and ratio plot
        V_pos, EN_pos, p_grid, d_grid = _branch0_grids(branches_pos)
        V_neg, EN_neg, _, _ = _branch0_grids(branches_neg)

        _all_V_pos[label] = V_pos
        _all_V_neg[label] = V_neg

        _file_records.append(
            (
                f"{label} ({_polarity_desc('positive')})",
                p_grid,
                d_grid,
                EN_pos,
                V_pos,
                mod,
            )
        )
        _file_records.append(
            (
                f"{label} ({_polarity_desc('negative')})",
                p_grid,
                d_grid,
                EN_neg,
                V_neg,
                mod,
            )
        )

    if args.streamer_criterion is not None:
        _C = args.streamer_criterion
        _EN_sc = np.logspace(1.0, np.log10(3e5), 200)

        for _cs_label, _p_arr_cs, _d_arr_cs, _mod_cs in curve_specs:
            _s_label = f"Streamer (C={_C}), {_cs_label}"
            _EN_s = np.full_like(pd_arr, np.nan)
            _V_s = np.full_like(pd_arr, np.nan)

            print(f"\nSolving streamer criterion: {_s_label}")
            for _i, (_pd_i, _p_i, _d_i) in enumerate(zip(pd_arr, _p_arr_cs, _d_arr_cs)):
                _fvals = np.array(
                    [_aed_integral(_en, _p_i, _d_i, _mod_cs) - _C for _en in _EN_sc]
                )
                _idx = np.where(_fvals[:-1] * _fvals[1:] < 0)[0]
                if _idx.size == 0:
                    continue
                _k = _idx[0]
                try:
                    _root = scipy.optimize.brentq(
                        lambda _en, __p=_p_i, __d=_d_i, __m=_mod_cs: _aed_integral(
                            _en, __p, __d, __m
                        )
                        - _C,
                        _EN_sc[_k],
                        _EN_sc[_k + 1],
                        xtol=1e-6,
                        rtol=1e-10,
                    )
                    _EN_s[_i] = _root
                    _V_s[_i] = _root * _pd_i * 1e-21 / (_kB * args.T) * 1e5
                except ValueError:
                    pass

            _streamer_records.append((_s_label, _p_arr_cs, _d_arr_cs, _EN_s, _V_s))

            _smask = np.isfinite(_EN_s)
            if not np.any(_smask):
                print(f"  [{_s_label}] No solution found for any pd value.")
                continue

            if not args.no_plot:
                _sm_kw = dict(
                    color="k",
                    ls="-.",
                    lw=1.5,
                    marker="x",
                    markersize=5,
                    markevery=max(1, int(_smask.sum()) // 10),
                    label=_s_label,
                )
                ax1.loglog(pd_arr[_smask] * 1e3, _V_s[_smask] / 1000, **_sm_kw)
                ax2.loglog(pd_arr[_smask] * 1e3, _EN_s[_smask], **_sm_kw)

            print(
                f"\n  {'pd (bar·mm)':>14}  {'p (bar)':>10}  {'d (mm)':>8}  "
                f"{'E/N (Td)':>12}  {'U (kV)':>12}  {'E (V/m)':>14}  [{_s_label}]"
            )
            print("  " + "-" * 90)
            _sstep = max(1, int(_smask.sum()) // 20)
            _sindices = np.where(_smask)[0][::_sstep]
            for _j in _sindices:
                print(
                    f"  {pd_arr[_j]*1e3:>14.4e}  "
                    f"{_p_arr_cs[_j]:>10.4g}  "
                    f"{_d_arr_cs[_j]*1e3:>8.4g}  "
                    f"{_EN_s[_j]:>12.4f}  "
                    f"{_V_s[_j]/1000:>12.4f}  "
                    f"{_V_s[_j]/_d_arr_cs[_j]:>14.4e}"
                )

    if args.write_to_file and _file_records:

        SEP = "\t"

        # Build flat column list: (name, value_fn) where value_fn(j) -> float
        columns = [("pd_bar_mm", lambda j: pd_arr[j] * 1e3)]
        for lbl, p_arr_c, d_arr_c, EN_c, V_c, _rec_mod in _file_records:
            columns += [
                (f"p_bar[{lbl}]", lambda j, a=p_arr_c: a[j]),
                (f"d_mm[{lbl}]", lambda j, a=d_arr_c: a[j] * 1e3),
                (f"EN_Td[{lbl}]", lambda j, a=EN_c: a[j]),
                (f"U_kV[{lbl}]", lambda j, a=V_c: a[j] / 1e3),
                (
                    f"E_Vm[{lbl}]",
                    lambda j, a=V_c, b=d_arr_c: (
                        a[j] / b[j] if np.isfinite(a[j]) else np.nan
                    ),
                ),
            ]

        if _plot_aed:
            for lbl, p_arr_c, d_arr_c, EN_c, V_c, _rec_mod in _file_records:
                columns.append(
                    (
                        f"ionization_integral[{lbl}]",
                        lambda j, _EN=EN_c, _p=p_arr_c, _d=d_arr_c, _m=_rec_mod: (
                            _aed_integral(_EN[j], _p[j], _d[j], _m)
                            if np.isfinite(_EN[j])
                            else np.nan
                        ),
                    )
                )

        # Alpha=eta crossover columns (one group per pressure where a root exists)
        for p_val, en_c in _en_cross.items():
            if en_c is None:
                continue
            lbl = f"alpha=eta, p={p_val} bar"
            columns += [
                (f"p_bar[{lbl}]", lambda j, _p=p_val: _p),
                (f"d_mm[{lbl}]", lambda j, _p=p_val: pd_arr[j] / _p * 1e3),
                (f"EN_Td[{lbl}]", lambda j, _en=en_c: _en),
                (
                    f"U_kV[{lbl}]",
                    lambda j, _en=en_c, _T=args.T: _en
                    * pd_arr[j]
                    * 1e-16
                    / (_kB * _T)
                    / 1e3,
                ),
                (
                    f"E_Vm[{lbl}]",
                    lambda j, _en=en_c, _p=p_val, _T=args.T: _en
                    * _p
                    * 1e5
                    / (_kB * _T)
                    * 1e-21,
                ),
            ]

        # Streamer criterion columns
        for _slbl, _sp_arr, _sd_arr, _sEN, _sV in _streamer_records:
            columns += [
                (f"EN_Td[{_slbl}]", lambda j, a=_sEN: a[j]),
                (f"U_kV[{_slbl}]", lambda j, a=_sV: a[j] / 1e3),
                (
                    f"E_Vm[{_slbl}]",
                    lambda j, a=_sV, b=_sd_arr: (
                        a[j] / b[j] if np.isfinite(a[j]) else np.nan
                    ),
                ),
            ]

        # Column width: wide enough for every header name plus a 2-char margin,
        # and at least 14 to fit a 12-char scientific-notation value.
        W = max(14, max(len(name) for name, _ in columns) + 2)

        with open(args.write_to_file, "w") as fh:
            write_metadata_header(fh, extra_paths=[args.mechanism])
            fh.write(f"# Mechanism:   {mech_name}\n")
            fh.write(f"# Temperature: {args.T} K\n")
            # In fixed-d mode --p is empty and the pressure varies along the
            # sweep, so report the range actually solved rather than the
            # (unused) argument list.
            if args.p:
                p_str = ", ".join(f"{p} bar" for p in args.p)
            else:
                _p_all = (
                    np.concatenate(
                        [
                            rec[1][np.isfinite(rec[1])]
                            for rec in _file_records
                            if len(rec[1])
                        ]
                    )
                    if _file_records
                    else np.array([])
                )
                p_str = (
                    f"{_p_all.min():.4g} – {_p_all.max():.4g} bar "
                    f"(varying; fixed-d sweep)"
                    if _p_all.size
                    else "n/a"
                )
            fh.write(f"# Pressures:   {p_str}\n")
            if args.d:
                d_str = ", ".join(f"{d} mm" for d in args.d)
                fh.write(f"# Distances:   {d_str}\n")
            cfg_labels = ", ".join(d.get("label", "Baseline") for d in raw_dicts)
            fh.write(f"# Configs:     {cfg_labels}\n")
            fh.write(f"# Field type:  {_field_dist.field_type}\n")
            fh.write(f"# Lambda:      {args.lam} s⁻¹\n")
            fh.write(f"# Method:      {args.method}\n")
            fh.write(
                f"# Stepping:    N_min={_N_min}, N_max={_N_max}, tol={_dx_tol:.3g}\n"
            )
            if _field_dist.sphere_R is not None:
                fh.write(f"# Sphere R:    {_field_dist.sphere_R*1e3:.4g} mm\n")
            if _field_dist.field_type == "sphere-plane":
                fh.write(
                    "# Polarity:    sphere=positive → sphere is anode (+),  "
                    "sphere=negative → sphere is cathode (−)\n"
                )
            if _field_dist.field_type == "fieldline":
                fh.write(f"# Field line:  {_field_dist.fieldline_path}\n")
                fh.write(f"# Arc length:  {_field_dist.fieldline_length*1e3:.6g} mm\n")
                fh.write(
                    f"# ∫|E| ds:     {_field_dist.fieldline_integral:.6g} "
                    f"(field units of the file × m; a voltage only if the "
                    f"file tabulates |E| in V/m)\n"
                )
                if _field_dist.fieldline_applied_voltage is not None:
                    fh.write(
                        f"# U_applied:   "
                        f"{_field_dist.fieldline_applied_voltage/1e3:.6g} kV "
                        f"(--fieldline-voltage; only the shape of the profile "
                        f"enters the solve, so U*/U_applied is the factor the "
                        f"excitation must be scaled by to reach inception)\n"
                    )
                fh.write(
                    "# Polarity:    start=positive → first tabulated point is "
                    "anode (+),  start=negative → cathode (−)\n"
                )
            if _streamer_records:
                fh.write(f"# Streamer C:  {args.streamer_criterion}\n")
            fh.write("#\n")

            # Per-column descriptions
            for col_idx, (name, _) in enumerate(columns, start=1):
                fh.write(f"# Column {col_idx}: {name}\n")
            fh.write("#\n")

            # Data rows — one per pd point
            for j in range(len(pd_arr)):
                vals = [fn(j) for _, fn in columns]
                fh.write(SEP.join(f"{v:<{W}.6e}" for v in vals) + "\n")

        print(f"\nResults written to: {args.write_to_file}")

    if not args.no_plot:
        if _en_cross:
            multi_p = len(args.p) > 1
            for p, en_c in _en_cross.items():
                if en_c is None:
                    continue
                lbl = (
                    rf"$\alpha=\eta$, p={p} bar ({en_c:.1f} Td)"
                    if multi_p
                    else rf"$\alpha=\eta$  ({en_c:.1f} Td)"
                )
                ax2.axhline(en_c, color="k", ls=":", lw=1.5, label=lbl)

        if ax3 is not None and _all_V_pos:
            _ref_cfg_label = raw_dicts[0].get("label", "Baseline")
            n_settings = len(args.p) + len(args.d)
            for curve_idx, (label, _p_arr, _d_arr, _cm) in enumerate(curve_specs):
                # Extract config label and p/d setting from the curve label.
                # Label format: "<config.label>, p=X bar" or "<config.label>, d=X mm"
                for _sep in [", p=", ", d="]:
                    if _sep in label:
                        cfg_lbl, setting = label.split(_sep, 1)
                        setting = _sep.lstrip(", ") + setting  # e.g. "p=1.0 bar"
                        break
                else:
                    continue
                if cfg_lbl == _ref_cfg_label:
                    continue
                # Find the matching reference (first config) curve.
                ref_label = f"{_ref_cfg_label}, {setting}"
                if ref_label not in _all_V_pos or label not in _all_V_pos:
                    continue
                V_k = _all_V_pos[label]
                V_b = _all_V_pos[ref_label]
                with np.errstate(invalid="ignore", divide="ignore"):
                    ratio = V_k / V_b
                mask_r = np.isfinite(ratio)
                if not np.any(mask_r):
                    continue
                ratio_label = f"{cfg_lbl} ({setting})" if n_settings > 1 else cfg_lbl
                marker = _markers[curve_idx % len(_markers)]
                markevery = max(1, mask_r.sum() // 15)
                ax3.semilogx(
                    pd_arr[mask_r] * 1e3,
                    ratio[mask_r],
                    label=ratio_label,
                    color=_curve_color.get(label),
                    marker=marker,
                    markevery=markevery,
                    markersize=6,
                )
            ax3.axhline(1.0, color="k", ls="--", lw=1.0)
            ax3.set_xlabel("p·d  (bar·mm)")
            ax3.set_ylabel(rf"$U / U_{{\mathrm{{{_ref_cfg_label}}}}}$")
            ax3.set_title(f"Voltage ratio vs. {_ref_cfg_label!r}")
            ax3.legend(loc="best", framealpha=1.0)
            ax3.grid(True, which="both", ls="--", alpha=0.4)

        ax1.set_xlabel("p·d  (bar·mm)")
        ax1.set_ylabel("U  (kV)")
        ax1.set_title("Breakdown voltage")
        ax1.legend(loc="upper left", framealpha=1.0)
        ax1.grid(True, which="both", ls="--", alpha=0.4)

        ax2.set_xlabel("p·d  (bar·mm)")
        ax2.set_ylabel("E/N  (Td)")
        ax2.set_title("Critical reduced field (average)")
        ax2.legend(loc="best", framealpha=1.0)
        ax2.grid(True, which="both", ls="--", alpha=0.4)

        plt.tight_layout()

        if args.save_subplots:
            import matplotlib.transforms as _mtrans

            mech_stem = os.path.splitext(mech_name)[0]

            # Panels: (primary_ax, twin_ax_or_None, output_filename)
            panels = [
                (ax1, None, f"{mech_stem}_U.pdf"),
                (ax2, None, f"{mech_stem}_EN.pdf"),
            ]
            if ax3 is not None:
                panels.append((ax3, None, f"{mech_stem}_ratio.pdf"))

            # Hide all titles before saving
            _suptitle.set_visible(False)
            _saved_titles = {ax: ax.get_title() for ax, *_ in panels}
            for ax, *_ in panels:
                ax.set_title("")

            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()

            for ax_p, ax_t, fname in panels:
                axes_in_panel = [a for a in (ax_p, ax_t) if a is not None]
                bboxes = [a.get_tightbbox(renderer) for a in axes_in_panel]
                bbox_px = _mtrans.Bbox.union(bboxes)
                bbox_in = bbox_px.transformed(fig.dpi_scale_trans.inverted())
                fig.savefig(fname, bbox_inches=bbox_in)
                print(f"Saved: {fname}")

            # Restore titles for the interactive window
            _suptitle.set_visible(True)
            for ax, *_ in panels:
                ax.set_title(_saved_titles[ax])

        plt.show()
