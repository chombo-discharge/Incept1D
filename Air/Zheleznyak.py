#!/usr/bin/env python3
"""
Fit the Zheleznyak (1982) photoionization absorption curve with a multi-exponential
two-stream approximation.

The Zheleznyak absorption function in reduced units chi = p_O2 * r [Pa m] is:

    K(chi) = (exp(-lambda1 * chi) - exp(-lambda2 * chi)) / chi

with lambda1 and lambda2 being the O2 absorption limits for the UV photoionisation
window in air (Zheleznyak 1982), expressed here in SI pressure-reduced units [m^-1 Pa^-1].

Note that K(chi) = integral_{lambda1}^{lambda2} exp(-kappa*chi) dkappa, so the
two-stream approximation is a quadrature of this integral:

    K_2S(chi) = sum_j g_j * kappa_j * exp(-kappa_j * chi)

The effective weight per group is w_j = g_j * kappa_j, so that

    K_2S(chi) = sum_j w_j * exp(-kappa_j * chi)  ~  K(chi)

Both kappa_j and w_j are optimised jointly by minimising the mean squared relative
error, starting from a non-negative least-squares (NNLS) initial guess.

In the two-stream model (manuscript eq. two_stream / eq. drift_reaction), the
photoionisation rate at x from source at x' is proportional to

    sum_j g_j * kappa_j * exp(-kappa_j * |x - x'|)

so the reported kappa_j and g_j feed directly into the C and B matrices.

Usage
-----
    python Zheleznyak.py [--ngroups N]

    --ngroups N   Number of photon groups (default: 5)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

# numpy >= 2.0 renamed trapz -> trapezoid; support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz
from scipy.optimize import nnls, minimize

# ---------------------------------------------------------------------------
# Zheleznyak (1982) absorption limits for O2 in air UV window
# Original values: 0.035 and 2.0 cm^-1 Torr^-1; converted to SI [m^-1 Pa^-1]
# via  1 cm^-1 Torr^-1 = 100 m^-1 / 133.322 Pa  =  100/133.322 m^-1 Pa^-1
# ---------------------------------------------------------------------------
_CGS_TO_SI = 100.0 / 133.322  # cm^-1 Torr^-1  →  m^-1 Pa^-1
_LAMBDA1 = 0.035 * _CGS_TO_SI  # lower limit (longest photon range) [m^-1 Pa^-1]
_LAMBDA2 = 2.0 * _CGS_TO_SI  # upper limit (shortest photon range) [m^-1 Pa^-1]


def zheleznyak(chi):
    """
    Zheleznyak absorption function K(chi) = (exp(-lam1*chi) - exp(-lam2*chi)) / chi.

    Parameters
    ----------
    chi : array-like
        Reduced path length  p_O2 * r  [Pa m].

    Returns
    -------
    ndarray  [m^-1 Pa^-1]
    """
    chi = np.asarray(chi, dtype=float)
    return np.where(
        chi > 1e-14,
        (np.exp(-_LAMBDA1 * chi) - np.exp(-_LAMBDA2 * chi)) / chi,
        _LAMBDA2 - _LAMBDA1,
    )


def fit_twostream(n_groups, n_pts=2000, threshold=0.0, chi_hi=None):
    """
    Fit kappa_j and g_j to approximate K(chi) with sum_j g_j*kappa_j*exp(-kappa_j*chi).

    Procedure
    ---------
    1. Build an NNLS initial guess with log-spaced kappa on [lambda1, lambda2].
    2. Jointly optimise kappa_j and w_j = g_j*kappa_j by minimising the mean squared
       relative error (L-BFGS-B in log-parameter space for positivity).
    3. Normalise g_j so that sum(g_j) = 1.
    4. Drop any group with g_j < threshold and renormalise the survivors.

    Parameters
    ----------
    n_groups : int
    n_pts : int
        Number of sample points (log-spaced in chi) used for the fit.
    threshold : float
        Groups with g_j < threshold after normalisation are discarded and the
        remaining g_j are renormalised to sum to 1.  Default 0.0 (keep all).
    chi_hi : float or None
        Upper bound of the chi sample range [Pa m].  Points beyond chi_hi are
        excluded from the objective, avoiding the region where K(chi) ≈ 0 and
        relative errors are meaningless.  Default None uses 50/λ₁ (original
        full-range behaviour).

    Returns
    -------
    kappa : ndarray, shape (n_kept,)  [m^-1 Pa^-1]  sorted ascending
    g     : ndarray, shape (n_kept,)  photon group fractions, sum = 1
    rel_rms : float  RMS relative error of the final (pruned) fit
    """
    if chi_hi is None:
        chi_hi = 50.0 / _LAMBDA1
    chi = np.logspace(np.log10(0.02 / _LAMBDA2), np.log10(chi_hi), n_pts)
    target = zheleznyak(chi)

    # --- Step 1: NNLS initial guess ---
    kappa0 = np.logspace(np.log10(_LAMBDA1), np.log10(_LAMBDA2), n_groups)
    A_rel = np.column_stack(
        [np.exp(-kappa0[j] * chi) / target for j in range(n_groups)]
    )
    w0, _ = nnls(A_rel, np.ones_like(chi))
    w0 = np.maximum(w0, 1e-10)

    # --- Step 2: joint optimisation in log space ---
    # Parametrise as x = [log(kappa_1), ..., log(kappa_N), log(w_1), ..., log(w_N)]
    lk_lo = np.log(_LAMBDA1 * 0.5)
    lk_hi = np.log(_LAMBDA2 * 2.0)
    x0 = np.concatenate([np.log(kappa0), np.log(w0)])
    bounds = [(lk_lo, lk_hi)] * n_groups + [(-30.0, 10.0)] * n_groups

    def objective(x):
        lk = np.clip(x[:n_groups], lk_lo, lk_hi)
        lw = np.clip(x[n_groups:], -30.0, 10.0)
        k = np.exp(lk)
        w = np.exp(lw)
        approx = sum(w[j] * np.exp(-k[j] * chi) for j in range(n_groups))
        return float(np.mean(((approx - target) / target) ** 2))

    res = minimize(
        objective,
        x0,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 5000, "ftol": 1e-14, "gtol": 1e-10},
    )

    lk = np.clip(res.x[:n_groups], lk_lo, lk_hi)
    lw = np.clip(res.x[n_groups:], -30.0, 10.0)
    kappa = np.exp(lk)
    w = np.exp(lw)

    # --- Step 3: sort by kappa and normalise g ---
    idx = np.argsort(kappa)
    kappa = kappa[idx]
    w = w[idx]
    g = w / kappa
    g /= g.sum()
    w = g * kappa  # adjusted to match normalised g

    # --- Step 4: discard negligible groups and renormalise ---
    if threshold > 0.0:
        keep = g >= threshold
        kappa = kappa[keep]
        g = g[keep]
        g /= g.sum()
        w = g * kappa

    approx = sum(w[j] * np.exp(-kappa[j] * chi) for j in range(len(kappa)))
    scale = _trapz(target, chi) / max(_trapz(approx, chi), 1e-30)
    rel_rms = float(np.sqrt(np.mean(((approx * scale - target) / target) ** 2)))

    return kappa, g, rel_rms


def _convolution_1d(kernel_fn, source, x):
    """Numerical 1-D convolution  integral K(|x - x'|) S(x') dx'."""
    dx = x[1] - x[0]
    result = np.zeros_like(x)
    for i in range(len(x)):
        R = np.abs(x[i] - x)
        result[i] = np.sum(kernel_fn(R) * source) * dx
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Fit Zheleznyak (1982) absorption curve with two-stream approximation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--ngroups",
        type=int,
        default=5,
        metavar="N",
        help="Number of photon groups",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        metavar="T",
        help="Drop groups with g_j < T after normalisation (default: keep all)",
    )
    parser.add_argument(
        "--pd_max",
        type=float,
        default=None,
        metavar="PD",
        help="Upper fit limit as reduced gap pd [bar cm]; chi_hi = 200*pd_max [Pa m]"
        " (default: 50/lambda1)",
    )
    parser.add_argument(
        "--write-to-file",
        dest="write_to_file",
        type=str,
        default=None,
        metavar="FILE",
        help="Write subplot data to FILE_absorption_curve.dat, FILE_error_terms.dat,"
        " FILE_solutions.dat",
    )
    args = parser.parse_args()

    n = args.ngroups
    chi_hi = 200.0 * args.pd_max if args.pd_max is not None else None
    kappa, g, rel_rms = fit_twostream(n, threshold=args.threshold, chi_hi=chi_hi)
    w = g * kappa
    n_kept = len(kappa)

    # --- Print coefficients ---
    print()
    header = f"{n} groups"
    if n_kept < n:
        header += f", {n_kept} kept after threshold={args.threshold:.1e}"
    print(f"Two-stream fit to Zheleznyak (1982) absorption curve  ({header})")
    print(
        f"  UV window:  lambda1 = {_LAMBDA1:.6e} m^-1 Pa^-1, "
        f"lambda2 = {_LAMBDA2:.6e} m^-1 Pa^-1"
    )
    print()
    col_w = 24
    print(f"  {'j':>4}  {'(kappa/pO2)_j [1/(m Pa)]':>{col_w}}  {'g_j':>14}")
    print("  " + "-" * 50)
    for j in range(n_kept):
        print(f"  {j+1:>4}  {kappa[j]:>{col_w}.6e}  {g[j]:>14.6e}")
    print("  " + "-" * 50)
    print(f"  {'sum':>4}  {'':>{col_w}}  {g.sum():>14.6e}")
    print()
    print(f"  RMS rel. error = {rel_rms*100:.2f} %")

    # --- Absorption function comparison ---
    chi = np.logspace(np.log10(0.01 / _LAMBDA2), np.log10(1e3), 3000)
    K_zh = zheleznyak(chi)
    K_2s = sum(w[j] * np.exp(-kappa[j] * chi) for j in range(n_kept))
    # scale K_2s to same integral as K_zh for visual comparison
    scale = _trapz(K_zh, chi) / max(_trapz(K_2s, chi), 1e-30)
    K_2s_scaled = K_2s * scale

    # --- Convolution with a Gaussian source ---
    x_max = 4.0 / _LAMBDA1
    x = np.linspace(0.0, x_max, 600)
    sigma = 0.04 * x_max
    source = np.exp(-(((x - 0.3 * x_max) / sigma) ** 2))  # off-centre Gaussian

    def K_zh_fn(R):
        return zheleznyak(np.maximum(R, 1e-14))

    def K_2s_fn(R):
        R = np.asarray(R)
        return (
            sum(w[j] * np.exp(-kappa[j] * np.maximum(R, 1e-14)) for j in range(n_kept))
            * scale
        )

    conv_zh = _convolution_1d(K_zh_fn, source, x)
    conv_2s = _convolution_1d(K_2s_fn, source, x)
    norm_c = max(conv_zh.max(), 1e-30)
    conv_zh /= norm_c
    conv_2s /= norm_c

    # Relative error (hoisted so it is available for both file output and plotting)
    err = np.abs((K_2s_scaled - K_zh) / (K_zh + 1e-30)) * 100.0

    # --- Write subplot data to files ---
    if args.write_to_file is not None:
        import datetime, sys as _sys

        stem = args.write_to_file
        generated = datetime.datetime.now().isoformat(timespec="seconds")
        command = " ".join(_sys.argv)

        def _write(fname, header_lines, col_names, columns):
            """Write columns (list of 1-D arrays) to a fixed-width dat file."""
            W = max(22, max(len(c) for c in col_names) + 2)
            with open(fname, "w") as fh:
                for line in header_lines:
                    fh.write(f"# {line}\n")
                fh.write("#\n")
                fh.write("# " + "  ".join(f"{c:>{W}}" for c in col_names) + "\n")
                for row in zip(*columns):
                    fh.write("  " + "  ".join(f"{v:>{W}.8e}" for v in row) + "\n")
            print(f"Written: {fname}")

        # --- FILE 1: absorption curve ---
        group_col_names = [f"group_{j+1}_kappa{kappa[j]:.3e}" for j in range(n_kept)]
        group_cols = [w[j] * scale * np.exp(-kappa[j] * chi) for j in range(n_kept)]
        _write(
            f"{stem}_absorption_curve.dat",
            [
                "Zheleznyak (1982) absorption curve and two-stream fit",
                f"Generated:  {generated}",
                f"Command:    {command}",
                f"ngroups={n},  n_kept={n_kept},  threshold={args.threshold:.1e}",
                f"lambda1 = {_LAMBDA1:.6e} m^-1 Pa^-1,  "
                f"lambda2 = {_LAMBDA2:.6e} m^-1 Pa^-1",
                f"RMS relative error = {rel_rms * 100:.3f} %",
                "",
                "Columns:",
                "  chi            : reduced path length p_O2*r  [Pa m]",
                "  K_zheleznyak   : Zheleznyak (1982) absorption function  [m^-1 Pa^-1]",
                "  K_twostream    : two-stream fit, integral-normalised to K_zheleznyak  [m^-1 Pa^-1]",
            ]
            + [
                f"  {group_col_names[j]:20s}: group {j+1}: "
                f"kappa={kappa[j]:.4e} m^-1 Pa^-1,  g={g[j]:.4e}  [m^-1 Pa^-1]"
                for j in range(n_kept)
            ],
            ["chi", "K_zheleznyak", "K_twostream"] + group_col_names,
            [chi, K_zh, K_2s_scaled] + group_cols,
        )

        # --- FILE 2: error terms ---
        _write(
            f"{stem}_error_terms.dat",
            [
                "Pointwise relative error of two-stream fit to Zheleznyak (1982)",
                f"Generated:  {generated}",
                f"Command:    {command}",
                f"ngroups={n},  n_kept={n_kept},  threshold={args.threshold:.1e}",
                f"RMS relative error = {rel_rms * 100:.3f} %",
                "",
                "Columns:",
                "  chi            : reduced path length p_O2*r  [Pa m]",
                "  rel_error_pct  : |K_twostream - K_zheleznyak| / K_zheleznyak * 100  [%]",
            ],
            ["chi", "rel_error_pct"],
            [chi, err],
        )

        # --- FILE 3: convolution solutions ---
        _write(
            f"{stem}_solutions.dat",
            [
                "Convolution of a Gaussian source with Zheleznyak and two-stream kernels",
                f"Generated:  {generated}",
                f"Command:    {command}",
                f"ngroups={n},  n_kept={n_kept},  threshold={args.threshold:.1e}",
                f"Source: Gaussian centred at 0.3*x_max, sigma = 0.04*x_max",
                f"x_max = {x_max:.6e} Pa m  (= 4 / lambda1)",
                "All convolution values normalised by max(conv_zheleznyak)",
                "",
                "Columns:",
                "  x                  : reduced coordinate p_O2*r  [Pa m]",
                "  source             : normalised Gaussian source (max = 1)  [dimensionless]",
                "  conv_zheleznyak    : photoionisation rate, Zheleznyak kernel, normalised  [dimensionless]",
                "  conv_twostream     : photoionisation rate, two-stream kernel, normalised  [dimensionless]",
            ],
            ["x", "source", "conv_zheleznyak", "conv_twostream"],
            [x, source / source.max(), conv_zh, conv_2s],
        )

    # --- Plot ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: absorption function (log-log)
    ax = axes[0]
    ax.loglog(chi, K_zh, "k-", lw=2, label="Zheleznyak (1982)")
    ax.loglog(chi, K_2s_scaled, "r--", lw=2, label=f"Two-stream ({n_kept} groups)")
    for j in range(n_kept):
        if g[j] > 1e-6:
            ax.loglog(
                chi,
                w[j] * scale * np.exp(-kappa[j] * chi),
                ":",
                color="steelblue",
                alpha=0.55,
                label="Individual groups" if j == 0 else None,
            )
    ax.set_xlabel(r"$\chi = p_{\mathrm{O}_2}\,r$  [Pa m]")
    ax.set_ylabel(r"$K(\chi)$  [Pa$^{-1}$ m$^{-1}$]")
    ax.set_title("Absorption function")
    ax.set_ylim(bottom=1e-14, top=10)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", ls=":", alpha=0.4)

    # Panel 2: relative error
    ax = axes[1]
    ax.semilogx(chi, err, "r-", lw=1.5)
    ax.set_xlabel(r"$\chi = p_{\mathrm{O}_2}\,r$  [Pa m]")
    ax.set_ylabel("Relative error [%]")
    ax.set_title("Pointwise relative error")
    ax.grid(True, ls=":", alpha=0.4)

    # Panel 3: convolution with test source
    ax = axes[2]
    ax.plot(x, conv_zh, "k-", lw=2, label="Zheleznyak (1982)")
    ax.plot(x, conv_2s, "r--", lw=2, label=f"Two-stream ({n} groups)")
    ax.plot(
        x, source / source.max(), color="gray", ls=":", lw=1.5, label="Source (norm.)"
    )
    ax.set_xlabel(r"$\chi = p_{\mathrm{O}_2}\,r$  [Pa m]")
    ax.set_ylabel("Normalised photoionisation rate")
    ax.set_title("Convolution: Gaussian source")
    ax.legend(fontsize=8)
    ax.grid(True, ls=":", alpha=0.4)

    plt.suptitle(
        rf"Zheleznyak vs two-stream ({n_kept} groups, RMS rel. error = {rel_rms*100:.1f}%)",
        fontsize=11,
    )
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
