"""
Lambda.py — compute the temporal growth rate λ vs. voltage above inception.

For a fixed geometry (pd, p, T, field), finds the inception voltage V* by solving
det Q(λ=0, EN*) = 0, then for each voltage V ∈ [V*, F·V*] finds λ* > 0 such that
det Q(λ*, EN) = 0.  The result is the discharge growth rate as a function of applied
voltage, with V* as the natural normalisation point.

The root det Q(λ) = 0 is found for fixed EN by:
  1. Evaluating det Q at λ=0 to establish the reference sign.
  2. Expanding an upper bracket geometrically (×10 per step) from 1 s⁻¹ until a
     sign change is detected.
  3. Refining with scipy.optimize.brentq.
NaN returns from inception_det (Q nearly singular near the root) are mapped to the
opposite-sign sentinel so brentq converges through the singularity.

Usage
-----
    python Lambda.py <mechanism> [CONFIG.json ...]
        --pd PD               p·d in bar·mm (required)
        --p  P                gas pressure in bar (default: 1.0)
        --T  T                temperature in K (default: 293.0)
        --n-voltages N        number of voltage points (default: 20)
        --v-max-factor F      upper voltage as multiple of V* (default: 2.0)
        --field SPEC          'uniform' | 'sphere-plane R_mm' | 'sphere-sphere R_mm'
        --dx [N_min [N_max [tol]]]  adaptive stepping parameters
        --method METHOD       propagator: 'midpoint' (default) or 'magnus2'
        --no-plot             suppress plot
        --write-to-file FILE  write tab-separated results to FILE
"""

import os
import sys
import json
import argparse
import datetime
import functools
import subprocess

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Inception import (
    _kB,
    load_mechanism,
    _read_json_configs,
    inception_det,
    find_all_breakdown_EN,
    midpoint_propagator,
    magnus2_propagator,
    parse_dx_spec,
    _DX_N_MIN_DEFAULT,
    _DX_N_MAX_DEFAULT,
    _DX_TOL_DEFAULT,
)
from FieldDistributions import FieldDistribution, add_field_argument, parse_field_spec
from IonizationIntegral import _max_real_eigenvalue


# ── Core solver ───────────────────────────────────────────────────────────────

def find_lambda_for_voltage(EN_ref, pd, mod, p, T, field_dist,
                             N_min, N_max, tol, propagator, lam_scale=None):
    """
    Find λ* > 0 such that det Q(λ*, EN_ref) = 0.

    For EN_ref > EN* (overvoltage), det Q(λ=0) is NaN because Q is near-singular
    (same convention as Inception.py: NaN = "above inception").  The root λ* is
    found by locating a λ_hi where det Q is finite and positive (sub-threshold),
    then bracketing [0, λ_hi] with brentq (NaN → -1e-300 sentinel throughout).

    Parameters
    ----------
    EN_ref     : float   — reduced field in Townsend
    pd         : float   — pressure × gap in bar·m
    mod        : Mechanism
    p          : float   — pressure in bar
    T          : float   — temperature in K
    field_dist : FieldDistribution
    N_min, N_max, tol : int, int, float
    propagator : callable
    lam_scale  : float or None
        Approximate λ scale in s⁻¹ for the upper-bracket search.  If None,
        estimated from the dominant ionisation eigenvalue × electron speed.

    Returns
    -------
    lam_star : float   — growth rate in s⁻¹ (NaN on failure)
    status   : str     — 'ok', 'suspect', or reason for failure
    """
    def _det_raw(lam_val):
        return inception_det(EN_ref, pd, mod, p, T, field_dist=field_dist,
                             N_min=N_min, N_max=N_max, tol=tol,
                             lam=lam_val, positive_polarity=True,
                             propagator=propagator)

    # NaN convention: NaN → -1e-300 (= "above inception", same as Inception.py).
    # This is correct for λ < λ* (where Q is near-singular due to overvoltage),
    # AND is used consistently throughout brentq.
    def _f_brentq(lam_val):
        val = _det_raw(lam_val)
        return val if np.isfinite(val) else -1e-300

    # det Q(0) > 0 means V ≤ V*: already sub-threshold at λ=0.
    f0 = _det_raw(0.0)
    if np.isfinite(f0) and f0 > 0.0:
        return 0.0, 'below_inception'

    # Estimate the scale for the upper-bracket search.
    if lam_scale is None:
        lam_eig   = max(1.0, _max_real_eigenvalue(mod, EN_ref, p, T))
        Vmat      = mod.get_V(EN_ref, p, T)
        v_e       = abs(float(np.diag(Vmat)[mod.ELECTRON_INDEX]))
        lam_scale = max(lam_eig * v_e, 1e4)

    # Probe around lam_scale (and broader) for a λ_hi where det Q is finite & > 0.
    # This point is in the "sub-threshold" region (λ > λ*) where Q is well-conditioned.
    lam_hi = None
    probes = [1.0, 0.3, 3.0, 0.1, 10.0, 30.0, 0.03, 100.0, 300.0, 1000.0]
    for factor in probes:
        test_lam = lam_scale * factor
        v = _det_raw(test_lam)
        if np.isfinite(v) and v > 0.0:
            lam_hi = test_lam
            break
    if lam_hi is None:
        for log_lam in range(3, 14):
            v = _det_raw(10.0 ** log_lam)
            if np.isfinite(v) and v > 0.0:
                lam_hi = 10.0 ** log_lam
                break
    if lam_hi is None:
        return float('nan'), 'no_bracket_found'

    # f_brentq(0) ≤ 0  (NaN → -1e-300, or genuinely negative)
    # f_brentq(lam_hi) > 0
    # → valid bracket for brentq.
    lam_star = scipy.optimize.brentq(
        _f_brentq, 0.0, lam_hi, xtol=1.0, rtol=1e-6, maxiter=60,
    )

    # Very small result means V ≈ V* (Q was singular at λ=0 and root is near zero).
    if lam_star < 1.0:
        return 0.0, 'ok'

    f_check = _det_raw(lam_star)
    f_ref   = _det_raw(lam_hi)
    if np.isfinite(f_check) and np.isfinite(f_ref) and abs(f_check) > 1e-2 * abs(f_ref):
        return lam_star, 'suspect'
    return lam_star, 'ok'


def compute_lambda_curve(EN_star, pd, mod, p, T, field_dist, N_min, N_max, tol,
                         propagator, n_voltages, v_max_factor):
    """
    Sweep voltages from V* to v_max_factor·V* and find λ at each point.

    Parameters
    ----------
    EN_star     : float   — inception E/N in Townsend
    pd          : float   — pressure × gap in bar·m
    ...

    Returns
    -------
    list of (V_kV, V_ratio, EN_ref, lam_star, tau_ns, nu_ion, status)
        V_kV    — voltage in kV
        V_ratio — V / V*
        EN_ref  — reduced field in Townsend
        lam_star — growth rate in s⁻¹
        tau_ns  — e-folding time in ns
        nu_ion  — ionisation frequency scale in s⁻¹ (max eigenvalue × v_e)
        status  — 'ok' or diagnostic string
    """
    V_star     = EN_star * pd * 1e-16 / (_kB * T)      # inception voltage in V
    voltages_V = np.geomspace(V_star, v_max_factor * V_star, n_voltages)

    results = []
    for i, V in enumerate(voltages_V):
        EN_ref  = V * _kB * T / (pd * 1e-16)       # Townsend
        V_kV    = V * 1e-3
        V_ratio = V / V_star

        # Ionisation frequency scale: max real eigenvalue of A = R V⁻¹, times v_e.
        lam_eig   = _max_real_eigenvalue(mod, EN_ref, p, T)   # m⁻¹
        Vmat      = mod.get_V(EN_ref, p, T)
        v_e       = abs(float(np.diag(Vmat)[mod.ELECTRON_INDEX]))  # m/s
        nu_ion    = max(lam_eig, 0.0) * v_e                        # s⁻¹
        lam_scale = max(nu_ion, 1e4)

        if i == 0:
            # V = V* by construction; λ* = 0 by definition.
            lam_star, status = 0.0, 'ok'
        else:
            lam_star, status = find_lambda_for_voltage(
                EN_ref, pd, mod, p, T, field_dist, N_min, N_max, tol, propagator,
                lam_scale=lam_scale,
            )

        tau_ns = 1e9 / lam_star if (np.isfinite(lam_star) and lam_star > 0.0) else float('nan')
        results.append((V_kV, V_ratio, EN_ref, lam_star, tau_ns, nu_ion, status))

        tag = f'  [{status}]' if status != 'ok' else ''
        if i == 0:
            print(f'  V = {V_kV:8.4f} kV  (1.000×V*)  λ = 0  (inception){tag}')
        elif np.isfinite(lam_star):
            print(f'  V = {V_kV:8.4f} kV  ({V_ratio:.3f}×V*)  '
                  f'λ = {lam_star:.4e} s⁻¹  τ = {tau_ns:.3f} ns{tag}')
        else:
            print(f'  V = {V_kV:8.4f} kV  ({V_ratio:.3f}×V*)  λ = NaN{tag}')

    return results


# ── File output ───────────────────────────────────────────────────────────────

def _write_results(out_path, curve_records, args, raw_dicts, mech_name, field_dist):
    """Write λ vs V results to a tab-separated file with metadata header."""
    SEP = '\t'
    try:
        _repo_dir = os.path.dirname(os.path.abspath(__file__))
        _git_hash = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=_repo_dir, stderr=subprocess.DEVNULL,
        ).decode().strip()
        _dirty = subprocess.check_output(
            ['git', 'status', '--porcelain', os.path.abspath(__file__),
             os.path.abspath(args.mechanism)],
            cwd=_repo_dir, stderr=subprocess.DEVNULL,
        ).decode().strip()
        git_str = _git_hash + (' (dirty)' if _dirty else '')
    except Exception:
        git_str = 'unavailable'

    # Build column headers from all curve labels
    col_headers = ['V_kV', 'V_ratio']
    for label, _ in curve_records:
        col_headers += [f'lambda_s-1[{label}]', f'tau_ns[{label}]',
                        f'nu_ion_s-1[{label}]']
    W = max(14, max(len(h) for h in col_headers) + 2)

    with open(out_path, 'w') as fh:
        fh.write('# --- METADATA ---\n')
        fh.write(f'# Date:        {datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'# Git:         {git_str}\n')
        fh.write(f'# Command:     {" ".join(sys.argv)}\n')
        fh.write('# ---\n')
        fh.write(f'# Mechanism:   {mech_name}\n')
        fh.write(f'# pd:          {args.pd} bar·mm\n')
        fh.write(f'# Pressure:    {args.p} bar\n')
        fh.write(f'# Temperature: {args.T} K\n')
        fh.write(f'# Field type:  {field_dist.label}\n')
        fh.write(f'# V_max_factor:{args.v_max_factor}\n')
        fh.write(f'# n_voltages:  {args.n_voltages}\n')
        fh.write(f'# Configs:     {", ".join(d.get("label", "Baseline") for d in raw_dicts)}\n')
        fh.write('#\n')
        for i, h in enumerate(col_headers, start=1):
            fh.write(f'# Column {i}: {h}\n')
        fh.write('#\n')
        fh.write(SEP.join(f'{h:<{W}}' for h in col_headers) + '\n')

        # All curves share the same V array (same pd/p)
        n_pts = len(curve_records[0][1])
        for j in range(n_pts):
            row_vals = [curve_records[0][1][j][0],   # V_kV
                        curve_records[0][1][j][1]]    # V_ratio
            for _, rows in curve_records:
                _, _, _, lam, tau, nu, _ = rows[j]
                row_vals += [lam, tau, nu]
            fh.write(SEP.join(f'{v:<{W}.6e}' for v in row_vals) + '\n')

    print(f'\nResults written to: {out_path}')


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            'Compute temporal growth rate λ vs. voltage above inception.  '
            'Finds V* from det Q(λ=0) = 0, then solves det Q(λ, EN) = 0 for '
            'λ at each voltage V ∈ [V*, F·V*].'
        )
    )
    parser.add_argument(
        'mechanism',
        help='Path to mechanism Python file (e.g. Air/Air_Pancheshnyi.py).',
    )
    parser.add_argument(
        'configs', nargs='*', metavar='CONFIG.json',
        help=(
            'One or more JSON configuration files.  Each may contain a single '
            'object or a list under "configurations".  If omitted, a single '
            'baseline configuration is used.'
        ),
    )
    parser.add_argument(
        '--pd', type=float, required=True, metavar='PD',
        help='Pressure × gap product in bar·mm.',
    )
    parser.add_argument(
        '--p', type=float, default=1.0, metavar='P',
        help='Gas pressure in bar (default: 1.0).',
    )
    parser.add_argument(
        '--T', type=float, default=293.0,
        help='Gas temperature in Kelvin (default: 293.0).',
    )
    parser.add_argument(
        '--n-voltages', type=int, default=20, metavar='N',
        help='Number of logarithmically spaced voltage points (default: 20).',
    )
    parser.add_argument(
        '--v-max-factor', type=float, default=2.0, metavar='F',
        help='Upper voltage as multiple of V* (default: 2.0).',
    )
    add_field_argument(parser)
    parser.add_argument(
        '--dx', nargs='*', default=None, metavar='SPEC',
        help=(
            'Adaptive integration stepping: N_min [N_max [tol]].  '
            'Defaults: N_min=5, N_max=200, tol=0.03.'
        ),
    )
    parser.add_argument(
        '--method', choices=['midpoint', 'magnus2'], default='midpoint',
        help="Propagator algorithm (default: 'midpoint').",
    )
    parser.add_argument(
        '--no-plot', action='store_true', default=False,
        help='Suppress the matplotlib figure.',
    )
    parser.add_argument(
        '--write-to-file', type=str, default=None, metavar='FILE',
        help='Write tab-separated results to FILE.',
    )
    args = parser.parse_args()

    if args.v_max_factor <= 1.0:
        parser.error('--v-max-factor must be > 1.0')
    if args.n_voltages < 2:
        parser.error('--n-voltages must be >= 2')

    _field_dist           = parse_field_spec(args.field, parser)
    _N_min, _N_max, _tol  = parse_dx_spec(args.dx, parser)
    _propagator           = midpoint_propagator if args.method == 'midpoint' else magnus2_propagator

    pd_m      = args.pd * 1e-3              # bar·mm → bar·m
    mech_name = os.path.basename(args.mechanism)

    raw_dicts = _read_json_configs(args.configs) if args.configs else [{}]

    # ── Find inception and compute λ curve for each config ─────────────────
    curve_records = []   # list of (label, rows)  where rows = compute_lambda_curve(...)

    for cfg_dict in raw_dicts:
        mod   = load_mechanism(args.mechanism, cfg_dict)
        label = mod.label or 'Baseline'

        # Build det_fn for the inception solve (λ=0).
        det_fn = functools.partial(
            inception_det, field_dist=_field_dist,
            N_min=_N_min, N_max=_N_max, tol=_tol,
            lam=0.0, positive_polarity=True, propagator=_propagator,
        )

        print(f'\nFinding inception voltage: {label}  '
              f'(pd = {args.pd} bar·mm, p = {args.p} bar)')
        roots = find_all_breakdown_EN(pd_m, mod, args.p, args.T,
                                      first_only=True, det_fn=det_fn)
        if not roots:
            print(f'  [{label}] No inception found — skipping.')
            continue

        EN_star = roots[0]
        V_star  = EN_star * pd_m * 1e-16 / (_kB * args.T)
        print(f'  V* = {V_star * 1e-3:.4f} kV   EN* = {EN_star:.2f} Td')
        print(f'\nComputing λ(V) from V* to {args.v_max_factor:.2g}×V*:')

        rows = compute_lambda_curve(
            EN_star, pd_m, mod, args.p, args.T, _field_dist,
            _N_min, _N_max, _tol, _propagator,
            args.n_voltages, args.v_max_factor,
        )
        curve_records.append((label, rows))

    if not curve_records:
        print('\nNo curves computed.')
        return

    # ── Console summary table ───────────────────────────────────────────────
    print()
    for label, rows in curve_records:
        print(f'  {label}')
        print('  ' + '-' * 78)
        print(f"  {'V (kV)':>10}  {'V/V*':>7}  {'EN (Td)':>10}  "
              f"{'λ (s⁻¹)':>14}  {'τ (ns)':>10}  {'λ/ν_ion':>10}")
        for V_kV, V_ratio, EN_ref, lam, tau, nu, status in rows:
            lam_str    = f'{lam:.4e}'    if np.isfinite(lam)  else '     NaN'
            tau_str    = f'{tau:.4f}'    if np.isfinite(tau)  else '     NaN'
            ratio_str  = f'{lam/nu:.4f}' if (np.isfinite(lam) and nu > 0) else '     NaN'
            flag = f'  [{status}]' if status != 'ok' else ''
            print(f'  {V_kV:>10.4f}  {V_ratio:>7.4f}  {EN_ref:>10.2f}  '
                  f'{lam_str:>14}  {tau_str:>10}  {ratio_str:>10}{flag}')
        print()

    # ── File output ─────────────────────────────────────────────────────────
    if args.write_to_file:
        _write_results(args.write_to_file, curve_records, args, raw_dicts,
                       mech_name, _field_dist)

    # ── Plot ─────────────────────────────────────────────────────────────────
    if args.no_plot:
        return

    plt.rcParams['font.size'] += 2
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f'Growth rate vs. voltage  —  {mech_name},  '
        f'pd = {args.pd} bar·mm,  p = {args.p} bar,  T = {args.T} K,  '
        f'{_field_dist.label}',
        fontsize=12,
    )

    _markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
    for idx, (label, rows) in enumerate(curve_records):
        V_ratio_arr = np.array([r[1] for r in rows])
        lam_arr     = np.array([r[3] for r in rows])
        tau_arr     = np.array([r[4] for r in rows])
        color       = f'C{idx}'
        marker      = _markers[idx % len(_markers)]
        kw          = dict(color=color, marker=marker, markersize=5, label=label)

        mask = np.isfinite(lam_arr)
        ax1.semilogy(V_ratio_arr[mask], lam_arr[mask], **kw)
        mask2 = np.isfinite(tau_arr)
        ax2.semilogy(V_ratio_arr[mask2], tau_arr[mask2], **kw)

    ax1.set_xlabel(r'$V / V^*$')
    ax1.set_ylabel(r'$\lambda$  (s$^{-1}$)')
    ax1.set_title('Temporal growth rate')
    ax1.grid(True, which='both', ls='--', alpha=0.4)

    ax2.set_xlabel(r'$V / V^*$')
    ax2.set_ylabel(r'$\tau = 1/\lambda$  (ns)')
    ax2.set_title('e-folding time')
    ax2.grid(True, which='both', ls='--', alpha=0.4)

    if len(curve_records) > 1:
        ax1.legend(loc='best', framealpha=1.0)
        ax2.legend(loc='best', framealpha=1.0)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
