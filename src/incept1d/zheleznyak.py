# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The Zheleznyak photoionization absorption function and its two-stream fit.

The absorption function of the UV photoionization window of O2 in air, in
the reduced path length chi = p_O2 r [Pa m] (Zheleznyak 1982), is

    K(chi) = (exp(-lambda1 chi) - exp(-lambda2 chi)) / chi
           = integral_{lambda1}^{lambda2} exp(-kappa chi) dkappa,

and :func:`fit_twostream` approximates it by a sum of exponentials,

    K(chi) ~ sum_j g_j kappa_j exp(-kappa_j chi),

one photon group per term, with the absorption coefficients kappa_j
[m^-1 Pa^-1] and the group fractions g_j (sum 1).  A mechanism file uses it
through ``"absorption": {"model": "zheleznyak", ...}``
(:mod:`incept1d.declarative`); ``mechanisms/air/zheleznyak.py`` plots the
fit.
"""

import functools

import numpy as np
from scipy.optimize import minimize, nnls

# numpy >= 2.0 renamed trapz -> trapezoid; support both.
_trapz = getattr(np, "trapezoid", None) or np.trapz

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


@functools.lru_cache(maxsize=None)
def twostream_groups(n_groups):
    """
    :func:`fit_twostream` for *n_groups*, computed once per process.

    Returns
    -------
    kappa : ndarray  [m^-1 Pa^-1], sorted ascending
    g : ndarray  group fractions, sum 1
    """
    kappa, g, _ = fit_twostream(int(n_groups))
    kappa.setflags(write=False)
    g.setflags(write=False)
    return kappa, g
