# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Closed-form inception condition for the three-species reduced model.

Transcribed directly from ``docs/source/theory/inceptioncriterion.rst``,
equations ``eq_generalized_paschen`` and ``eq_standard_paschen``.  This module
deliberately shares no code with :mod:`incept1d`: it is the independent
reference the solver is checked against, so it must stay a transcription of
the documented algebra rather than a refactoring of the implementation.
"""

import math

import scipy.optimize


def lambda_pm(alpha, eta, delta):
    """
    Eigenvalues of the electron/negative-ion block of A = R V^-1.

        lambda_pm = (alpha - eta - delta +- Delta) / 2,
        Delta     = sqrt((alpha - eta - delta)^2 + 4 alpha delta)

    ``lambda_p`` is the apparent effective ionization coefficient; it reduces
    to ``alpha - eta`` when ``delta = 0``.

    Returns
    -------
    (lambda_p, lambda_m, Delta) : tuple of float
    """
    tr = alpha - eta - delta
    Delta = math.sqrt(tr * tr + 4.0 * alpha * delta)
    return 0.5 * (tr + Delta), 0.5 * (tr - Delta), Delta


def generalized_paschen_lhs(alpha, eta, delta, gamma, d):
    """
    Left-hand side of :eq:`eq_generalized_paschen`; inception where it equals 1.

        gamma alpha [ c_+/lambda_+ (e^{lambda_+ d} - 1)
                    + c_-/lambda_- (e^{lambda_- d} - 1) ] = 1

    with ``c_pm = (Delta +- (alpha - eta + delta)) / (2 Delta)``.
    """
    lam_p, lam_m, Delta = lambda_pm(alpha, eta, delta)
    c_p = (Delta + (alpha - eta + delta)) / (2.0 * Delta)
    c_m = (Delta - (alpha - eta + delta)) / (2.0 * Delta)

    def term(c, lam):
        # (e^{lam d} - 1) / lam, continuous at lam -> 0 where it tends to d.
        if abs(lam * d) < 1e-12:
            return c * d
        return c * math.expm1(lam * d) / lam

    return gamma * alpha * (term(c_p, lam_p) + term(c_m, lam_m))


def standard_paschen_lhs(alpha, eta, gamma, d):
    """
    Left-hand side of :eq:`eq_standard_paschen` (the delta = 0 limit).

        gamma alpha / (alpha - eta) (e^{(alpha - eta) d} - 1) = 1
    """
    ae = alpha - eta
    if abs(ae * d) < 1e-12:
        return gamma * alpha * d
    return gamma * alpha * math.expm1(ae * d) / ae


def solve_EN(alpha_of_EN, eta_of_EN, delta_of_EN, gamma, d, EN_lo=1.0, EN_hi=1e5):
    """
    Solve the closed-form condition for the reduced field E/N in Townsend.

    Brackets the root of ``generalized_paschen_lhs(...) - 1`` by scanning
    ``[EN_lo, EN_hi]`` logarithmically, then refines with Brent's method.

    Returns
    -------
    float
        The lowest E/N at which the closed-form criterion is met.

    Raises
    ------
    RuntimeError
        If no sign change is found in the scan range.
    """

    def f(EN):
        return (
            generalized_paschen_lhs(
                alpha_of_EN(EN), eta_of_EN(EN), delta_of_EN(EN), gamma, d
            )
            - 1.0
        )

    n = 400
    prev_EN = EN_lo
    prev = f(prev_EN)
    for i in range(1, n + 1):
        EN = EN_lo * (EN_hi / EN_lo) ** (i / n)
        cur = f(EN)
        if prev * cur < 0.0:
            return float(
                scipy.optimize.brentq(f, prev_EN, EN, xtol=1e-12, rtol=1e-14)
            )
        prev_EN, prev = EN, cur
    raise RuntimeError(f"no closed-form root in [{EN_lo}, {EN_hi}] Td for d = {d} m")
