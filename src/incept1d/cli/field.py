# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
``incept1d field`` — plot the normalised field profile f(ξ) and the midpoint
quadrature grid for a gap geometry (see :mod:`incept1d.fields`).
"""

import numpy as np
import matplotlib.pyplot as plt

from incept1d.fields import add_field_argument, parse_field_spec

HELP = "Plot the normalised field profile f(ξ) of a gap geometry."
DESCRIPTION = (
    "Plot the normalised on-axis electric field and midpoint quadrature grid "
    "for a given gap geometry."
)


def add_arguments(parser):
    """Register the command-line arguments on *parser*."""
    add_field_argument(parser)
    parser.add_argument(
        "--d",
        type=float,
        default=None,
        metavar="D_mm",
        help=(
            "Gap distance in mm (required unless --field fieldline, where it "
            "defaults to the arc length of the tabulated line)."
        ),
    )
    parser.add_argument(
        "--N",
        type=int,
        default=25,
        metavar="N",
        help="Number of grid points to display (default: 25).",
    )


def run(args, parser):
    """Run the command with parsed *args*; *parser* is used for ``parser.error``."""
    fd = parse_field_spec(args.field, parser)
    N = args.N
    if args.d is None:
        if fd.field_type != "fieldline":
            parser.error("--d is required for this field type")
        args.d = fd.fieldline_length * 1e3
    d_mm = args.d
    d = d_mm * 1e-3
    f = fd.build(d)
    color = "tab:blue"

    geom_str = f"{fd.label},  d = {d_mm} mm"
    if fd.sphere_R is not None:
        geom_str += f"  (d/R = {d/fd.sphere_R:.3f})"

    print(f"Geometry:  {geom_str}")
    if fd.field_type == "fieldline":
        print(
            f"U_file   = {fd.fieldline_voltage/1e3:.6g} kV"
            f"   (= ∫|E| ds in the field units of the file)"
        )
        print(
            "           only the shape of f(ξ) enters the solve, so this is a "
            "units check, not an input"
        )
    print(f"f(0) = {f(0.0):.6f}")
    print(f"f(1) = {f(1.0):.6f}")
    print(f"N = {N}  (cell width = {d_mm / N:.3f} mm)")
    print()

    xi_plot = np.linspace(0.0, 1.0, 2000)
    f_plot = np.array([f(xi) for xi in xi_plot])
    x_plot = xi_plot * d_mm

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle(f"{geom_str},  N = {N}", fontsize=13)

    exact_lbl = "tabulated" if fd.field_type == "fieldline" else "exact"
    ax.plot(x_plot, f_plot, "k-", lw=1.8, zorder=3, label=f"$f(\\xi)$ ({exact_lbl})")
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
