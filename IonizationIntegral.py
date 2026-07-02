"""
IonizationIntegral.py — plot ∫ max(α−η, 0) dx and ∫ max(Re(λ_max), 0) dx vs voltage.

Sweeps a voltage range and plots, for every (config, p, d) combination:

- ∫ max(α−η, 0) dx — classical ionisation integral
- ∫ max(Re(λ_max(RV⁻¹)), 0) dx — integral of the leading eigenvalue of the
  transport matrix (shown only when the mechanism exposes get_R and get_V)

The field distribution is specified with the same --field syntax as Inception.py.
Use --single-voltage to evaluate at a single voltage instead, or --data-file to
read (pressure, voltage) pairs from an experimental data file.

Usage
-----
::

    python IonizationIntegral.py <mechanism> [CONFIG.json ...]
        --p P [P ...]              pressure(s) in bar (required unless --data-file)
        --d D [D ...]              gap distance(s) in mm (required)
        --voltage-lo V             lower voltage bound in kV (default: 0.0)
        --voltage-hi V             upper voltage bound in kV (default: 1.0)
        --voltage-num N            number of points in sweep (default: 100)
        --single-voltage V         evaluate at one voltage [kV]; print result, no plot
        --data-file FILE           read (p, V) pairs from ASCII data file
        --pressure-column COL      column name or 0-based index for pressure
        --voltage-column COL       column name or 0-based index for voltage
        --T T                      temperature in K (default 293.0)
        --field SPEC               'uniform' | 'sphere-plane R_mm [tol [steps]]'
                                   | 'sphere-sphere R_mm [tol [steps]]'
        --no-plot                  suppress the matplotlib window
        --write-to-file FILE       write tab-separated results to FILE
"""

import os
import sys
import argparse
import subprocess
import datetime

import numpy as np

# ------------------------------------------------------------------
# Import shared infrastructure from Inception.py (same directory).
# ------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Inception import _kB, load_mechanism, _read_json_configs
from FieldDistributions import FieldDistribution, add_field_argument, parse_field_spec


# ------------------------------------------------------------------
# Core integrators
# ------------------------------------------------------------------


def _max_real_eigenvalue(mod, EN, p, T):
    """Return Re(λ_max) of A = R V^{-1} at a single E/N value [m^{-1}]."""
    R = mod.get_R(EN, p, T)
    V = mod.get_V(EN, p, T)
    return float(np.max(np.real(np.linalg.eigvals(R @ np.linalg.inv(V)))))


def aed_integral(EN_ref, p_val, d_val, mod, T, field_dist: FieldDistribution, N: int):
    """Return ∫ max(α−η, 0) dx [m⁻¹], or NaN on overflow.

    Parameters
    ----------
    EN_ref     : float            — E/N in Td at the reference point
    p_val      : float            — pressure in bar
    d_val      : float            — gap distance in metres
    mod        :                  — loaded mechanism module (must expose alpha, eta)
    T          : float            — temperature in K
    field_dist : FieldDistribution — field geometry specification
    N          : int              — number of midpoint-rule integration steps
    """
    try:
        if field_dist.field_type == "uniform":
            val = (mod.alpha(EN_ref, p_val, T) - mod.eta(EN_ref, p_val, T)) * d_val
            result = max(0.0, float(val))
            return result if np.isfinite(result) else float("nan")

        f = field_dist.build(d_val)

        def _net(xi):
            en = EN_ref * f(xi)
            return mod.alpha(en, p_val, T) - mod.eta(en, p_val, T)

        # If there is no net ionisation even at the sphere surface, return 0 immediately.
        if _net(0.0) <= 0.0:
            return 0.0

        # Locate xi_cross ∈ (0, 1] where α(E) = η(E).  The field is monotonically
        # decreasing from sphere (xi=0) to plane/far-sphere (xi=1), so _net is also
        # decreasing and the active region is exactly [0, xi_cross].
        # Restricting integration to this interval avoids the partial-cell error that
        # arises when a coarse quadrature cell straddles the threshold.
        if _net(1.0) > 0.0:
            xi_cross = 1.0  # active all the way to the far electrode
        else:
            lo, hi = 0.0, 1.0
            for _ in range(52):  # 52 bisections → ~1e-15 precision in xi
                mid = 0.5 * (lo + hi)
                if _net(mid) > 0.0:
                    lo = mid
                else:
                    hi = mid
            xi_cross = 0.5 * (lo + hi)

        # Midpoint-rule quadrature over [0, xi_cross·d].  Every cell is within the
        # active region so all contributions are positive — no max(0,·) needed.
        xis = (np.arange(N) + 0.5) / N * xi_cross
        EN_arr = EN_ref * np.array([f(xi) for xi in xis])
        ds = xi_cross * d_val / N
        diff = np.array(
            [mod.alpha(en, p_val, T) - mod.eta(en, p_val, T) for en in EN_arr]
        )
        result = float(np.sum(diff) * ds)

    except (OverflowError, ValueError):
        return float("nan")

    return result if np.isfinite(result) else float("nan")


def eig_integral(EN_ref, p_val, d_val, mod, T, field_dist: FieldDistribution, N: int):
    """Return ∫ max(Re(λ_max(RV^{-1})), 0) dx [m⁻¹], or NaN on overflow.

    Parameters
    ----------
    EN_ref     : float            — E/N in Td at the reference point
    p_val      : float            — pressure in bar
    d_val      : float            — gap distance in metres
    mod        :                  — loaded mechanism module (must expose get_R, get_V)
    T          : float            — temperature in K
    field_dist : FieldDistribution
    N          : int              — number of midpoint-rule integration steps
    """
    try:
        if field_dist.field_type == "uniform":
            val = _max_real_eigenvalue(mod, EN_ref, p_val, T) * d_val
            result = max(0.0, float(val))
            return result if np.isfinite(result) else float("nan")

        f = field_dist.build(d_val)

        def _lmax(xi):
            return _max_real_eigenvalue(mod, EN_ref * f(xi), p_val, T)

        if _lmax(0.0) <= 0.0:
            return 0.0

        if _lmax(1.0) > 0.0:
            xi_cross = 1.0
        else:
            lo, hi = 0.0, 1.0
            for _ in range(52):
                mid = 0.5 * (lo + hi)
                if _lmax(mid) > 0.0:
                    lo = mid
                else:
                    hi = mid
            xi_cross = 0.5 * (lo + hi)

        xis = (np.arange(N) + 0.5) / N * xi_cross
        EN_arr = EN_ref * np.array([f(xi) for xi in xis])
        ds = xi_cross * d_val / N
        eigs = np.array([_max_real_eigenvalue(mod, en, p_val, T) for en in EN_arr])
        result = float(np.sum(eigs) * ds)

    except (OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan")

    return result if np.isfinite(result) else float("nan")


# ------------------------------------------------------------------
# Data-file helpers
# ------------------------------------------------------------------


def _read_data_file(path, pressure_col, voltage_col):
    """Read (pressure [bar], voltage [kV]) columns from a whitespace/tab ASCII file.

    pressure_col, voltage_col: column name string or 0-based integer index.
    Returns (p_arr, V_arr, p_col_name, V_col_name).
    """
    with open(path) as fh:
        header_line = fh.readline()
    col_names = header_line.split()

    def _resolve(spec):
        try:
            return int(spec)
        except (ValueError, TypeError):
            try:
                return col_names.index(str(spec))
            except ValueError:
                raise ValueError(
                    f"Column '{spec}' not found in {os.path.basename(path)}.  "
                    f"Available: {col_names}"
                )

    p_idx = _resolve(pressure_col)
    v_idx = _resolve(voltage_col)
    data = np.loadtxt(path, skiprows=1)
    return data[:, p_idx], data[:, v_idx], col_names[p_idx], col_names[v_idx]


def _print_data_file_table(records, has_matrix):
    from itertools import groupby

    for lbl, rows in groupby(records, key=lambda r: r[0]):
        rows = list(rows)
        print(f"\n  {lbl}")
        print("  " + "-" * (80 + (26 if has_matrix else 0)))
        hdr = (
            f"  {'p (bar)':>10}    {'V (kV)':>12}    {'E/N (Td)':>12}"
            f"    {'∫max(α−η,0)dx (m⁻¹)':>22}"
        )
        if has_matrix:
            hdr += f"    {'∫max(λ_max,0)dx (m⁻¹)':>22}"
        print(hdr)
        for _, p_bar, d_mm, V_kV, EN_ref, aed_val, eig_val in rows:
            aed_str = f"{aed_val:>22.6e}" if np.isfinite(aed_val) else f"{'NaN':>22}"
            row = f"  {p_bar:>10.4f}    {V_kV:>12.4f}    {EN_ref:>12.2f}    {aed_str}"
            if has_matrix:
                eig_str = (
                    f"{eig_val:>22.6e}" if np.isfinite(eig_val) else f"{'NaN':>22}"
                )
                row += f"    {eig_str}"
            print(row)
    print()


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot ∫ max(α−η, 0) dx vs. applied voltage for a fixed (p, d) geometry.  "
            "Multiple pressures, distances, and configurations each produce a separate curve."
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
            "key.  If omitted, a single baseline configuration is used."
        ),
    )
    parser.add_argument(
        "--p",
        type=float,
        nargs="*",
        required=False,
        default=[],
        metavar="P",
        help="Pressure(s) in bar (required unless --data-file is given).",
    )
    parser.add_argument(
        "--data-file",
        type=str,
        default=None,
        metavar="FILE",
        help="ASCII data file containing (pressure, voltage) columns.",
    )
    parser.add_argument(
        "--pressure-column",
        type=str,
        default=None,
        metavar="COL",
        help="Column name or 0-based index for pressure in --data-file.",
    )
    parser.add_argument(
        "--voltage-column",
        type=str,
        default=None,
        metavar="COL",
        help="Column name or 0-based index for voltage in --data-file.",
    )
    parser.add_argument(
        "--d",
        type=float,
        nargs="+",
        required=True,
        metavar="D",
        help="Gap distance(s) in mm (required, one or more).",
    )
    parser.add_argument(
        "--voltage-lo",
        type=float,
        default=0.0,
        metavar="V",
        help="Lower bound of voltage sweep in kV (default: 0.0).",
    )
    parser.add_argument(
        "--voltage-hi",
        type=float,
        default=1.0,
        metavar="V",
        help="Upper bound of voltage sweep in kV (default: 1.0).",
    )
    parser.add_argument(
        "--voltage-num",
        type=int,
        default=100,
        metavar="N",
        help="Number of linearly spaced voltage points in sweep (default: 100).",
    )
    parser.add_argument(
        "--single-voltage",
        type=float,
        default=None,
        metavar="V",
        help="Evaluate at a single voltage V [kV] and print result; no plot.",
    )
    parser.add_argument(
        "--T",
        type=float,
        default=293.0,
        help="Gas temperature in Kelvin (default: 293.0).",
    )
    add_field_argument(parser)
    parser.add_argument(
        "--no-plot",
        action="store_true",
        default=False,
        help="Suppress the matplotlib window (useful for batch/scripted runs).",
    )
    parser.add_argument(
        "--write-to-file",
        type=str,
        default=None,
        metavar="FILE",
        help="Write tab-separated results to FILE.",
    )
    args = parser.parse_args()

    # ---- Validation -----------------------------------------------------
    _data_file_mode = args.data_file is not None

    if _data_file_mode:
        if args.pressure_column is None or args.voltage_column is None:
            parser.error(
                "--data-file requires both --pressure-column and --voltage-column."
            )
    else:
        if not args.p:
            parser.error("--p is required unless --data-file is given.")
        if args.single_voltage is None and args.voltage_hi <= args.voltage_lo:
            parser.error("--voltage-hi must be greater than --voltage-lo.")

    # ---- Field specification --------------------------------------------
    _field_dist, _N = parse_field_spec(args.field, parser)
    _field_str = _field_dist.label

    # ---- Load configurations --------------------------------------------
    mech_dir = os.path.dirname(os.path.abspath(args.mechanism))
    mech_name = os.path.basename(args.mechanism)

    raw_dicts = _read_json_configs(args.configs) if args.configs else [{}]

    # ---- Load modules ---------------------------------------------------
    curve_specs = []  # (label, mod)
    for cfg_dict in raw_dicts:
        mod = load_mechanism(args.mechanism, cfg_dict)
        if not (hasattr(mod, "alpha") and hasattr(mod, "eta")):
            parser.error(
                f"Mechanism '{args.mechanism}' does not expose 'alpha' and 'eta'; "
                "these are required by IonizationIntegral.py."
            )
        curve_specs.append((mod.label, mod))

    has_matrix = all(
        hasattr(m, "get_R") and hasattr(m, "get_V") for _, m in curve_specs
    )

    # When there is exactly one pressure and one distance they go in the title;
    # otherwise they stay in the legend to disambiguate curves.
    single_pd = (not _data_file_mode) and (len(args.p) == 1 and len(args.d) == 1)

    # ---- Single-voltage mode --------------------------------------------
    if args.single_voltage is not None:
        records = []
        for cfg_label, mod in curve_specs:
            for p_bar in args.p:
                for d_mm in args.d:
                    d_m = d_mm * 1e-3
                    label = f"{cfg_label}, p={p_bar} bar, d={d_mm} mm"
                    V_kV = args.single_voltage
                    EN_ref = V_kV * 1e3 * _kB * args.T / (p_bar * d_m * 1e-16)
                    aed_val = aed_integral(
                        EN_ref, p_bar, d_m, mod, args.T, _field_dist, _N
                    )
                    eig_val = (
                        eig_integral(EN_ref, p_bar, d_m, mod, args.T, _field_dist, _N)
                        if has_matrix
                        else float("nan")
                    )
                    records.append((label, p_bar, d_mm, V_kV, EN_ref, aed_val, eig_val))

        # Print table
        from itertools import groupby

        eig_col = "  ∫max(λ_max,0)dx (m⁻¹)" if has_matrix else ""
        for lbl, rows in groupby(records, key=lambda r: r[0]):
            rows = list(rows)
            print(f"\n  {lbl}")
            print("  " + "-" * (64 + (26 if has_matrix else 0)))
            hdr = f"  {'V (kV)':>12}    {'E/N (Td)':>12}    {'∫max(α−η,0)dx (m⁻¹)':>22}"
            if has_matrix:
                hdr += f"    {'∫max(λ_max,0)dx (m⁻¹)':>22}"
            print(hdr)
            for _, p_bar, d_mm, V_kV, EN_ref, aed_val, eig_val in rows:
                aed_str = (
                    f"{aed_val:>22.6e}"
                    if np.isfinite(aed_val)
                    else f"{'NaN (overflow)':>22}"
                )
                row = f"  {V_kV:>12.4f}    {EN_ref:>12.2f}    {aed_str}"
                if has_matrix:
                    eig_str = (
                        f"{eig_val:>22.6e}"
                        if np.isfinite(eig_val)
                        else f"{'NaN (overflow)':>22}"
                    )
                    row += f"    {eig_str}"
                print(row)
        print()

        if args.write_to_file:
            sv_curve_data = [
                (
                    lbl,
                    np.array([V_kV]),
                    np.array([aed_val]),
                    np.array([eig_val]) if has_matrix else None,
                )
                for lbl, p_bar, d_mm, V_kV, EN_ref, aed_val, eig_val in records
            ]
            _write_results(
                args, sv_curve_data, raw_dicts, mech_name, _field_dist, _N, has_matrix
            )
        return

    # ---- Data-file mode -------------------------------------------------
    if _data_file_mode:
        try:
            p_arr, V_arr_file, p_col_name, V_col_name = _read_data_file(
                args.data_file, args.pressure_column, args.voltage_column
            )
        except (ValueError, OSError) as exc:
            parser.error(str(exc))
        single_d = len(args.d) == 1

        curve_data = []
        file_curve_data = []
        all_records = []

        for cfg_label, mod in curve_specs:
            for d_mm in args.d:
                d_m = d_mm * 1e-3
                label = f"{cfg_label}, d={d_mm} mm"
                legend_label = cfg_label if single_d else label

                aed_arr = np.array(
                    [
                        aed_integral(
                            V_kV * 1e3 * _kB * args.T / (p_bar * d_m * 1e-16),
                            p_bar,
                            d_m,
                            mod,
                            args.T,
                            _field_dist,
                            _N,
                        )
                        for p_bar, V_kV in zip(p_arr, V_arr_file)
                    ]
                )
                eig_arr = (
                    np.array(
                        [
                            eig_integral(
                                V_kV * 1e3 * _kB * args.T / (p_bar * d_m * 1e-16),
                                p_bar,
                                d_m,
                                mod,
                                args.T,
                                _field_dist,
                                _N,
                            )
                            for p_bar, V_kV in zip(p_arr, V_arr_file)
                        ]
                    )
                    if has_matrix
                    else None
                )
                curve_data.append((legend_label, V_arr_file, aed_arr, eig_arr))
                file_curve_data.append((label, V_arr_file, aed_arr, eig_arr))

                for p_bar, V_kV, aed_val, eig_val in zip(
                    p_arr,
                    V_arr_file,
                    aed_arr,
                    (
                        eig_arr
                        if eig_arr is not None
                        else np.full_like(aed_arr, float("nan"))
                    ),
                ):
                    EN_ref = V_kV * 1e3 * _kB * args.T / (p_bar * d_m * 1e-16)
                    all_records.append(
                        (label, p_bar, d_mm, V_kV, EN_ref, aed_val, eig_val)
                    )

        _print_data_file_table(all_records, has_matrix)

        if args.write_to_file:
            _write_results(
                args,
                file_curve_data,
                raw_dicts,
                mech_name,
                _field_dist,
                _N,
                has_matrix,
                data_file=args.data_file,
            )

        if not args.no_plot:
            import matplotlib.pyplot as plt

            plt.rcParams["font.size"] += 2
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.set_yscale("symlog", linthresh=1.0)
            ax.set_xlabel(f"{V_col_name}  (kV)")
            ax.set_ylabel(r"Ionisation integral  (m$^{-1}$)")
            title = f"{mech_name},  T = {args.T} K,  {_field_str}"
            if single_d:
                title += f",  d = {args.d[0]} mm"
            ax.set_title(title)
            ax.grid(True, which="both", ls="--", alpha=0.4)

            _markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
            for idx, (lbl, V_file, aed_arr, eig_arr) in enumerate(curve_data):
                color = f"C{idx}"
                mask = np.isfinite(aed_arr)
                ax.plot(
                    V_file[mask],
                    aed_arr[mask],
                    color=color,
                    ls="-",
                    marker=_markers[idx % len(_markers)],
                    markersize=5,
                    label=rf"{lbl}  ($\alpha-\eta$)",
                )
                if eig_arr is not None:
                    mask2 = np.isfinite(eig_arr)
                    ax.plot(
                        V_file[mask2],
                        eig_arr[mask2],
                        color=color,
                        ls="--",
                        marker=_markers[idx % len(_markers)],
                        markersize=5,
                        label=rf"{lbl}  ($\lambda_\mathrm{{max}}$)",
                    )
            ax.legend(loc="best", framealpha=1.0)
            plt.tight_layout()
            plt.show()
        return

    # ---- Sweep mode -----------------------------------------------------
    V_arr = np.linspace(args.voltage_lo, args.voltage_hi, args.voltage_num)

    # curve_data: list of (label, V_arr, aed_arr, eig_arr|None) — one per (config, p, d)
    curve_data = []
    file_curve_data = []
    all_records = []

    for cfg_label, mod in curve_specs:
        for p_bar in args.p:
            for d_mm in args.d:
                d_m = d_mm * 1e-3
                label = f"{cfg_label}, p={p_bar} bar, d={d_mm} mm"
                legend_label = cfg_label if single_pd else label
                EN_arr = V_arr * 1e3 * _kB * args.T / (p_bar * d_m * 1e-16)

                aed_arr = np.array(
                    [
                        aed_integral(EN_ref, p_bar, d_m, mod, args.T, _field_dist, _N)
                        for EN_ref in EN_arr
                    ]
                )
                eig_arr = (
                    np.array(
                        [
                            eig_integral(
                                EN_ref, p_bar, d_m, mod, args.T, _field_dist, _N
                            )
                            for EN_ref in EN_arr
                        ]
                    )
                    if has_matrix
                    else None
                )
                curve_data.append((legend_label, V_arr, aed_arr, eig_arr))
                file_curve_data.append((label, V_arr, aed_arr, eig_arr))

                # Summary line
                finite_vals = aed_arr[np.isfinite(aed_arr)]
                nan_mask = ~np.isfinite(aed_arr)
                nan_info = ""
                if np.any(nan_mask):
                    nan_info = f"  [first NaN at U = {V_arr[nan_mask][0]:.4g} kV]"
                if len(finite_vals):
                    print(
                        f"  {label} (α−η) :  "
                        f"min={finite_vals.min():.3e}  "
                        f"max={finite_vals.max():.3e}"
                        f"{nan_info}"
                    )
                else:
                    print(f"  {label} (α−η) :  all NaN (overflow at all voltages)")

                if eig_arr is not None:
                    finite_eig = eig_arr[np.isfinite(eig_arr)]
                    nan_eig = ~np.isfinite(eig_arr)
                    nan_info2 = ""
                    if np.any(nan_eig):
                        nan_info2 = f"  [first NaN at U = {V_arr[nan_eig][0]:.4g} kV]"
                    if len(finite_eig):
                        print(
                            f"  {label} (λ_max):  "
                            f"min={finite_eig.min():.3e}  "
                            f"max={finite_eig.max():.3e}"
                            f"{nan_info2}"
                        )
                    else:
                        print(f"  {label} (λ_max):  all NaN")

                for V_kV, EN_ref, aed_val, eig_val in zip(
                    V_arr,
                    EN_arr,
                    aed_arr,
                    (
                        eig_arr
                        if eig_arr is not None
                        else np.full_like(aed_arr, float("nan"))
                    ),
                ):
                    all_records.append(
                        (label, p_bar, d_mm, V_kV, EN_ref, aed_val, eig_val)
                    )

    print()

    # ---- Write to file --------------------------------------------------
    if args.write_to_file:
        _write_results(
            args, file_curve_data, raw_dicts, mech_name, _field_dist, _N, has_matrix
        )

    # ---- Plot -----------------------------------------------------------
    if not args.no_plot:
        import matplotlib.pyplot as plt

        plt.rcParams["font.size"] += 2
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_xlabel("U  (kV)")
        ax.set_ylabel(r"Ionisation integral  (m$^{-1}$)")
        if single_pd:
            title = (
                f"{mech_name},  p = {args.p[0]} bar,  d = {args.d[0]} mm,"
                f"  T = {args.T} K,  {_field_str}"
            )
        else:
            title = f"{mech_name},  T = {args.T} K,  {_field_str}"
        ax.set_title(title)
        ax.grid(True, which="both", ls="--", alpha=0.4)

        _markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
        for idx, (label, V_kV_arr, aed_arr, eig_arr) in enumerate(curve_data):
            color = f"C{idx}"
            markevery = max(1, int(np.isfinite(aed_arr).sum()) // 15)

            # α−η curve (solid)
            mask_aed = np.isfinite(aed_arr)
            (line,) = ax.plot(
                V_kV_arr[mask_aed],
                aed_arr[mask_aed],
                color=color,
                ls="-",
                marker=_markers[idx % len(_markers)],
                markevery=markevery,
                markersize=5,
                label=rf"{label}  ($\alpha-\eta$)",
            )
            if np.any(~mask_aed):
                ax.axvline(
                    V_kV_arr[~mask_aed][0], color=color, ls=":", lw=1.0, alpha=0.7
                )

            # λ_max curve (dashed), same color
            if eig_arr is not None:
                mask_eig = np.isfinite(eig_arr)
                ax.plot(
                    V_kV_arr[mask_eig],
                    eig_arr[mask_eig],
                    color=color,
                    ls="--",
                    marker=_markers[idx % len(_markers)],
                    markevery=markevery,
                    markersize=5,
                    label=rf"{label}  ($\lambda_\mathrm{{max}}$)",
                )
                if np.any(~mask_eig):
                    ax.axvline(
                        V_kV_arr[~mask_eig][0], color=color, ls=":", lw=1.0, alpha=0.4
                    )

        ax.legend(loc="best", framealpha=1.0)
        plt.tight_layout()
        plt.show()


# ------------------------------------------------------------------
# File output
# ------------------------------------------------------------------


def _write_results(
    args,
    curve_data,
    raw_dicts,
    mech_name,
    field_dist: FieldDistribution,
    N: int,
    has_matrix: bool = False,
    data_file=None,
):
    """Write results in wide format: one row per voltage point, one column pair per curve."""
    SEP = "\t"

    col_headers = ["U_kV"]
    for lbl, _, _, eig_arr in curve_data:
        col_headers.append(f"streamer_integral[{lbl}]")
        if eig_arr is not None:
            col_headers.append(f"eigenvalue_integral[{lbl}]")

    W = max(14, max(len(h) for h in col_headers) + 2)

    try:
        _repo_dir = os.path.dirname(os.path.abspath(__file__))
        _git_hash = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=_repo_dir,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        _dirty = (
            subprocess.check_output(
                [
                    "git",
                    "status",
                    "--porcelain",
                    os.path.abspath(__file__),
                    os.path.abspath(args.mechanism),
                ],
                cwd=_repo_dir,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        git_str = _git_hash + (" (dirty)" if _dirty else "")
    except Exception:
        git_str = "unavailable"

    with open(args.write_to_file, "w") as fh:
        fh.write("# --- METADATA ---\n")
        fh.write(
            f"# Date:    {datetime.datetime.now().isoformat(timespec='seconds')}\n"
        )
        fh.write(f"# Git:     {git_str}\n")
        fh.write(f"# Command: {' '.join(sys.argv)}\n")
        fh.write("# ---\n")
        fh.write(f"# Mechanism:   {mech_name}\n")
        fh.write(f"# Temperature: {args.T} K\n")
        if data_file:
            fh.write(f"# Data file:   {os.path.basename(data_file)}\n")
            fh.write(f"# p column:    {args.pressure_column}\n")
            fh.write(f"# V column:    {args.voltage_column}\n")
        else:
            fh.write(f"# Pressures:   {', '.join(f'{p} bar' for p in args.p)}\n")
        fh.write(f"# Distances:   {', '.join(f'{d} mm' for d in args.d)}\n")
        fh.write(
            f"# Configs:     {', '.join(d.get('label', 'Baseline') for d in raw_dicts)}\n"
        )
        fh.write(f"# Field type:  {field_dist.field_type}\n")
        if field_dist.field_type != "uniform":
            fh.write(f"# Sphere R:    {field_dist.sphere_R*1e3:.4g} mm,  N = {N}\n")
        fh.write("#\n")
        for i, h in enumerate(col_headers, start=1):
            fh.write(f"# Column {i}: {h}\n")
        fh.write("#\n")

        fh.write(SEP.join(f"{h:<{W}}" for h in col_headers) + "\n")

        V_ref = curve_data[0][1]
        for j in range(len(V_ref)):
            row = [V_ref[j]]
            for _, _, aed_arr, eig_arr in curve_data:
                row.append(aed_arr[j])
                if eig_arr is not None:
                    row.append(eig_arr[j])
            fh.write(SEP.join(f"{v:<{W}.6e}" for v in row) + "\n")

    print(f"\nResults written to: {args.write_to_file}")


if __name__ == "__main__":
    main()
