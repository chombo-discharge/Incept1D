# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
An independent reference for the temporal growth rate.

``incept1d growth`` finds λ as a root of the boundary determinant
det Q(λ).  This module finds it a different way: it discretises the
linearised time-dependent equations of the Theory chapter directly,

    ∂_t n = −∂_x(V n) + R n + B (Ψ⁺ + Ψ⁻),
    ∂_x Ψ⁺ = C n − κ Ψ⁺,    ∂_x Ψ⁻ = −C n + κ Ψ⁻,

with the cathode and anode conditions of eq. ``eq_Q_system``, by Chebyshev
collocation, and solves the generalised eigenproblem L u = λ M u densely.
Every mode comes out, complex ones included, so the eigenvalue with the
largest real part *is* the fastest-growing mode — no argument about which
root of det Q to take is needed — and its eigenvector is the mode shape.

Like ``closed_form.py`` it must stay independent of the solver: it uses the
mechanism (the data) and the field profile, never ``incept1d.solver``.

Photon retardation (κ + λ/c) is neglected, which moves λ by a relative
amount of order λ/(cκ), far below the tolerances it is used with.  Every
photon group is kept explicit; the solver's optically thick shortcut
(κd > 12) has no counterpart, so compare only where no group is folded.
"""

import numpy as np
import scipy.linalg


def chebyshev(N, d):
    """Chebyshev–Gauss–Lobatto nodes on [0, d] (x[0] = 0) and d/dx there."""
    k = np.arange(N + 1)
    s = np.cos(np.pi * k / N)  # 1 → −1
    c = np.where((k == 0) | (k == N), 2.0, 1.0) * (-1.0) ** k
    S = np.tile(s, (N + 1, 1)).T
    D = np.outer(c, 1.0 / c) / (S - S.T + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    # s runs 1 → −1 while x = d (1 − s) / 2 runs 0 → d.
    return d * (1.0 - s) / 2.0, -2.0 / d * D


def growth_modes(mod, EN_ref, pd, p, T, field_dist, positive_polarity=True, N=120):
    """
    All modes of the linearised gap, fastest first.

    Parameters
    ----------
    mod : Mechanism
        Loaded mechanism.
    EN_ref : float
        Uniform-equivalent reduced field in Townsend.
    pd : float
        Pressure × gap in bar·m.
    p, T : float
        Pressure in bar and temperature in K.
    field_dist : FieldDistribution
        Gap geometry.
    positive_polarity : bool
        True when the ξ = 0 electrode is the anode.
    N : int
        Polynomial degree; N + 1 collocation nodes.

    Returns
    -------
    lam : ndarray
        Finite eigenvalues in s⁻¹, by decreasing real part.
    n : ndarray
        (len(lam), N + 1, n_species) densities of each mode at the nodes.
    x : ndarray
        The nodes in metres, cathode (0) to anode (d).
    """
    d = pd / p
    x, D = chebyshev(N, d)
    f = field_dist.build(d)
    xi = (d - x) / d if positive_polarity else x / d
    EN = EN_ref * np.array([f(v) for v in xi])
    eff = mod.resolve(positive_polarity)

    ns = len(eff.SPECIES)
    ng = eff.get_B(EN[0], p, T).shape[1]
    kappa = eff.get_kappa(p, T)
    P = N + 1
    size = (ns + 2 * ng) * P

    def sp(s):  # rows/columns of species s at all nodes
        return slice(s * P, (s + 1) * P)

    def ph(j, stream):  # photon group j, stream 0 = forward, 1 = backward
        base = (ns + stream * ng + j) * P
        return slice(base, base + P)

    L = np.zeros((size, size))
    M = np.zeros((size, size))
    R = np.array([eff.get_R(e, p, T) for e in EN])  # (P, ns, ns)
    v = np.array([np.diag(eff.get_V(e, p, T)) for e in EN])  # (P, ns)
    Bm = np.array([eff.get_B(e, p, T) for e in EN])  # (P, ns, ng)
    Cm = np.array([eff.get_C(e, p, T) for e in EN])  # (P, ng, ns)

    # Species: λ n_s = −∂x(v_s n_s) + Σ_r R_sr n_r + Σ_j B_sj (Ψ⁺_j + Ψ⁻_j).
    for s in range(ns):
        L[sp(s), sp(s)] -= D * v[:, s][np.newaxis, :]
        for r in range(ns):
            L[sp(s), sp(r)] += np.diag(R[:, s, r])
        for j in range(ng):
            L[sp(s), ph(j, 0)] += np.diag(Bm[:, s, j])
            L[sp(s), ph(j, 1)] += np.diag(Bm[:, s, j])
        M[sp(s), sp(s)] = np.eye(P)
    # Photons (no time derivative): ∂xΨ⁺ + κΨ⁺ − C n = 0, ∂xΨ⁻ − κΨ⁻ + C n = 0.
    for j in range(ng):
        L[ph(j, 0), ph(j, 0)] = D + kappa[j] * np.eye(P)
        L[ph(j, 1), ph(j, 1)] = D - kappa[j] * np.eye(P)
        for r in range(ns):
            L[ph(j, 0), sp(r)] = -np.diag(Cm[:, j, r])
            L[ph(j, 1), sp(r)] = np.diag(Cm[:, j, r])

    def replace(row, coeffs):
        """Replace equation `row` by the boundary condition Σ c · u = 0."""
        L[row, :] = 0.0
        M[row, :] = 0.0
        for col, value in coeffs:
            L[row, col] += value

    first, last = 0, P - 1
    sel = {
        name: [int(np.argmax(r)) for r in getter()]
        for name, getter in (
            ("e", eff.get_Pi_e),
            ("+", eff.get_Pi_plus),
            ("-", eff.get_Pi_minus),
        )
    }
    g_plus = eff.get_gamma_plus(EN[0], p, T)
    g_psi = eff.get_gamma_Psi(EN[0], p, T) if ng else np.zeros(0)
    # Cathode: electron emission, no negative ions, no forward photons.
    (e,) = sel["e"]
    replace(
        e * P + first,
        [(e * P + first, v[first, e])]
        + [(s * P + first, g * v[first, s]) for s, g in zip(sel["+"], g_plus)]
        + [(ph(j, 1).start + first, -g_psi[j]) for j in range(ng)],
    )
    for s in sel["-"]:
        replace(s * P + first, [(s * P + first, 1.0)])
    for j in range(ng):
        replace(ph(j, 0).start + first, [(ph(j, 0).start + first, 1.0)])
    # Anode: no positive ions, no backward photons.
    for s in sel["+"]:
        replace(s * P + last, [(s * P + last, 1.0)])
    for j in range(ng):
        replace(ph(j, 1).start + last, [(ph(j, 1).start + last, 1.0)])

    lam, vec = scipy.linalg.eig(L, M)
    finite = np.isfinite(lam)
    lam, vec = lam[finite], vec[:, finite]
    order = np.argsort(-lam.real)
    lam, vec = lam[order], vec[:, order]
    n = np.stack(
        [np.stack([vec[sp(s), i] for s in range(ns)], axis=-1) for i in range(len(lam))]
    )
    return lam, n, x


def converged_modes(
    mod,
    EN_ref,
    pd,
    p,
    T,
    field_dist,
    positive_polarity=True,
    Ns=(100, 150, 225),
    rtol=1e-4,
):
    """
    The modes that are resolved: present, to *rtol*, at every resolution.

    Chebyshev collocation of an advection operator also produces spurious
    eigenvalues whose imaginary part grows with N and whose real part does
    not settle; they are excluded by keeping only eigenvalues of the finest
    grid that the coarser grids reproduce.

    Returns
    -------
    lam : ndarray
        Converged eigenvalues, fastest first.
    n : ndarray
        Their densities at the nodes of the finest grid.
    """
    runs = [
        growth_modes(mod, EN_ref, pd, p, T, field_dist, positive_polarity, N)
        for N in Ns
    ]
    lam, n, _ = runs[-1]
    keep = [
        i
        for i, z in enumerate(lam)
        if all(
            np.min(np.abs(coarse - z)) < rtol * max(abs(z), 1.0)
            for coarse, _, _ in runs[:-1]
        )
    ]
    return lam[keep], n[keep]
