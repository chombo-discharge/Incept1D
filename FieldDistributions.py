"""
FieldDistributions.py — field geometry abstraction for gap-discharge computations.

Provides:
  FieldDistribution   — geometry spec independent of gap length d; .build(d) → f(xi)
  add_field_argument  — add --field SPEC to an argparse.ArgumentParser
  parse_field_spec    — parse a --field token list → (FieldDistribution, N)

Low-level functions (moved here from Inception.py):
  _sphere_sphere_axial_field   — exact bispherical image-charge series
  _compute_n_steps             — minimum N for ≤ field_tol relative change in first step
  _n_steps_sphere_cached       — LRU-cached version of _compute_n_steps for sphere geometry

Standalone usage
----------------
    python FieldDistributions.py --field sphere-plane 50 [N] --d D_mm
    python FieldDistributions.py --field sphere-sphere 50 25 --d 40
    python FieldDistributions.py --field uniform --d 20

Plots the normalised on-axis field and midpoint-quadrature grid for the specified
geometry.
"""

import argparse
import dataclasses
import functools
import math
from typing import Callable, Optional

import numpy as np
import matplotlib.pyplot as plt


# ── Bispherical image-charge series ──────────────────────────────────────────


def _sphere_sphere_axial_field(xi, alpha, a_focal, d):
    """
    Exact normalised on-axis field between two equal spheres at ±U/2.

    Bispherical image-charge series:
        a  = R · sinh(α),   z_n = a · coth(n·α),   q_n = ½a / sinh(n·α)
        f(ξ) = d · Σ_{n=1}^∞ q_n · [1/(z_n−z)² + 1/(z_n+z)²],   z = d(½−ξ)

    Returns f(ξ) satisfying ∫₀¹ f(ξ) dξ = 1 by construction.  Vectorised
    over n; truncated at n = 30/α (relative tail < 1e-13) or na = 700.
    """
    if alpha < 1e-14:
        return 1.0
    z = d * (0.5 - xi)
    n_max = min(100000, max(50, int(30.0 / alpha) + 10))
    na = np.arange(1, n_max + 1, dtype=float) * alpha
    na = na[na < 700.0]
    sinh_na = np.sinh(na)
    z_n = a_focal * np.cosh(na) / sinh_na
    q_n = 0.5 * a_focal / sinh_na
    return float(np.sum(q_n * (1.0 / (z_n - z) ** 2 + 1.0 / (z_n + z) ** 2))) * d


# ── Step-count helpers (kept for backward compatibility / plot_sphere_plane.py) ──


def _compute_n_steps(field_func, field_tol, max_steps=2000):
    """Return minimum N so that the relative field change across the first step ≤ field_tol."""
    f0 = field_func(0.0)
    if f0 <= 0.0 or field_tol >= 1.0:
        return 1
    lo, hi = 1, max_steps
    if (f0 - field_func(1.0 / hi)) / f0 > field_tol:
        return hi
    while lo < hi:
        mid = (lo + hi) // 2
        if (f0 - field_func(1.0 / mid)) / f0 <= field_tol:
            hi = mid
        else:
            lo = mid + 1
    return lo


@functools.lru_cache(maxsize=1024)
def _n_steps_sphere_cached(alpha, is_plane, field_tol, max_steps):
    """N_steps for sphere geometry — depends only on shape parameter alpha, not scale."""
    sinh_a = math.sinh(float(alpha))
    cosh_a = math.cosh(float(alpha))
    a_norm = sinh_a / (2.0 * (cosh_a - 1.0))
    if is_plane:

        def _f(xi):
            return _sphere_sphere_axial_field(0.5 * xi, float(alpha), a_norm, 1.0)

    else:

        def _f(xi):
            return _sphere_sphere_axial_field(xi, float(alpha), a_norm, 1.0)

    return _compute_n_steps(_f, field_tol, max_steps)


# ── Field distribution class ──────────────────────────────────────────────────


@dataclasses.dataclass
class FieldDistribution:
    """
    Field geometry specification, independent of gap length d and grid resolution.

    Call .build(d) to obtain the normalised field callable f(xi) for a specific
    gap length.  The number of integration steps N is a separate concern and is
    not stored here — it is supplied by the caller (CLI arg, adaptive scheme, etc.).

    Attributes
    ----------
    field_type : str
        One of 'uniform', 'sphere-plane', 'sphere-sphere'.
    sphere_R : float or None
        Sphere radius in metres.  None for uniform field.
    """

    field_type: str
    sphere_R: Optional[float] = None

    @property
    def is_symmetric(self) -> bool:
        """True when both polarities give the same inception result."""
        return self.field_type in ("uniform", "sphere-sphere")

    @property
    def label(self) -> str:
        """Human-readable description for plot titles and file headers."""
        if self.field_type == "uniform":
            return "uniform field"
        R_mm = self.sphere_R * 1e3
        return f"{self.field_type}, R = {R_mm:.4g} mm"

    def build(self, d: float) -> Callable[[float], float]:
        """
        Return the normalised field callable f(xi) for gap length d.

        Parameters
        ----------
        d : float
            Gap length in metres.

        Returns
        -------
        f : callable
            f(xi) → normalised field at fractional position xi ∈ [0,1].
            xi = 0 is the sphere surface (field maximum); xi = 1 is the plane
            / far sphere.  Satisfies ∫₀¹ f(xi) dxi = 1.
            For uniform field, f(xi) = 1.0 everywhere.
        """
        if self.field_type == "uniform":
            return lambda xi: 1.0

        if self.field_type == "sphere-plane":
            alpha = float(np.arccosh(1.0 + d / self.sphere_R))
            a = self.sphere_R * np.sinh(alpha)

            def f_sp(xi, _alpha=alpha, _a=a, _d=d):
                return _sphere_sphere_axial_field(0.5 * xi, _alpha, _a, 2.0 * _d)

            return f_sp

        # sphere-sphere
        alpha = float(np.arccosh(1.0 + d / (2.0 * self.sphere_R)))
        a = self.sphere_R * np.sinh(alpha)

        def f_ss(xi, _alpha=alpha, _a=a, _d=d):
            return _sphere_sphere_axial_field(xi, _alpha, _a, _d)

        return f_ss


# ── CLI utilities ─────────────────────────────────────────────────────────────


def add_field_argument(parser: argparse.ArgumentParser) -> None:
    """Add --field SPEC argument to an argparse.ArgumentParser."""
    parser.add_argument(
        "--field",
        nargs="+",
        default=["uniform"],
        metavar="SPEC",
        help=(
            "Field distribution.  'uniform' (default): spatially uniform field.  "
            "'sphere-plane R_mm': sphere-plane gap with sphere radius R_mm in mm.  "
            "'sphere-sphere R_mm': symmetric sphere-sphere gap.  "
            "Grid resolution is set separately via --dx."
        ),
    )


def parse_field_spec(spec: list, parser=None) -> "FieldDistribution":
    """
    Parse a --field token list into a FieldDistribution.

    Parameters
    ----------
    spec : list of str
        The raw token list from argparse (e.g. ['sphere-plane', '50']).
    parser : argparse.ArgumentParser or None
        If given, error messages are routed through parser.error(); otherwise
        a ValueError is raised.

    Returns
    -------
    FieldDistribution
    """

    def _err(msg):
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    field_type = spec[0]
    if field_type == "uniform":
        return FieldDistribution(field_type="uniform")

    if field_type in ("sphere-plane", "sphere-sphere"):
        if len(spec) < 2:
            _err(f"--field {field_type} requires a sphere radius in mm")
        sphere_R = float(spec[1]) * 1e-3  # mm → m
        return FieldDistribution(field_type=field_type, sphere_R=sphere_R)

    _err(
        f"Unknown field type: {field_type!r}. "
        "Use 'uniform', 'sphere-plane', or 'sphere-sphere'."
    )


# ── Standalone field plotter ──────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the normalised on-axis electric field and midpoint quadrature grid "
            "for a given gap geometry."
        )
    )
    add_field_argument(parser)
    parser.add_argument(
        "--d",
        type=float,
        required=True,
        metavar="D_mm",
        help="Gap distance in mm (required).",
    )
    parser.add_argument(
        "--N",
        type=int,
        default=25,
        metavar="N",
        help="Number of grid points to display (default: 25).",
    )
    args = parser.parse_args()

    fd = parse_field_spec(args.field, parser)
    N = args.N
    d_mm = args.d
    d = d_mm * 1e-3
    f = fd.build(d)
    color = "tab:blue"

    geom_str = f"{fd.label},  d = {d_mm} mm"
    if fd.sphere_R is not None:
        geom_str += f"  (d/R = {d/fd.sphere_R:.3f})"

    print(f"Geometry:  {geom_str}")
    print(f"f(0) = {f(0.0):.6f}")
    print(f"f(1) = {f(1.0):.6f}")
    print(f"N = {N}  (cell width = {d_mm / N:.3f} mm)")
    print()

    xi_plot = np.linspace(0.0, 1.0, 2000)
    f_plot = np.array([f(xi) for xi in xi_plot])
    x_plot = xi_plot * d_mm

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle(f"{geom_str},  N = {N}", fontsize=13)

    ax.plot(x_plot, f_plot, "k-", lw=1.8, zorder=3, label="$f(\\xi)$ (exact)")
    ax.axhline(1.0, color="gray", lw=0.8, ls="--", zorder=1, label="uniform ($f=1$)")

    for k in range(N + 1):
        ax.axvline(k / N * d_mm, color=color, lw=0.7, ls=":", alpha=0.6, zorder=2)

    xi_mids = (np.arange(N) + 0.5) / N
    f_mids = np.array([f(xi) for xi in xi_mids])
    ax.plot(
        xi_mids * d_mm,
        f_mids,
        "o",
        color=color,
        ms=6,
        zorder=5,
        label=f"cell midpoints (N={N})",
    )
    for i in range(N):
        ax.hlines(
            f_mids[i],
            i / N * d_mm,
            (i + 1) / N * d_mm,
            color=color,
            lw=2.5,
            alpha=0.55,
            zorder=4,
        )

    ax.set_xlabel("$x$  (mm)", fontsize=11)
    ax.set_ylabel("Normalised field  $f(\\xi) = E(x)\\,/\\,(V/d)$", fontsize=10)
    ax.set_xlim(-0.3, d_mm + 0.3)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, ls="--", alpha=0.35)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
