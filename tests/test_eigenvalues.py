# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Eigenvalues of the local transport matrix A = R V^-1."""

import numpy as np
import pytest

from closed_form import lambda_pm
from incept1d.constants import kB
from incept1d.eigenvalues import (
    compute_eigenvalues,
    compute_pressure_scan,
    max_real_eigenvalue,
)


class TestTracking:
    def test_shape_and_ordering(self, toy):
        """E3: the initial ordering is by descending real part."""
        EN = np.logspace(1.5, 3.0, 40)
        ev = compute_eigenvalues(toy, EN, 1.0, 293.0)
        assert ev.shape == (len(EN), len(toy.SPECIES))
        first = np.real(ev[0])
        assert np.all(np.diff(first) <= 1e-12)

    def test_leading_track_is_lambda_plus_everywhere(self, toy, toy_raw):
        """
        E2: continuity is what makes track 0 meaningful.

        If the Hungarian assignment ever swapped tracks, track 0 would stop
        agreeing with the analytic lambda_+ partway along the sweep.
        """
        p, T = 1.0, 293.0
        EN = np.logspace(1.5, 3.5, 60)
        ev = compute_eigenvalues(toy, EN, p, T)
        for i, en in enumerate(EN):
            lam_p, _, _ = lambda_pm(
                toy_raw.alpha(en, p, T), toy_raw.eta(en, p, T), toy_raw.delta(en, p, T)
            )
            assert np.real(ev[i, 0]) == pytest.approx(lam_p, rel=1e-8)

    def test_tracks_vary_continuously(self, toy):
        """E2b: no track may jump between adjacent, closely spaced points."""
        EN = np.logspace(2.0, 3.0, 200)
        ev = np.real(compute_eigenvalues(toy, EN, 1.0, 293.0))
        scale = np.ptp(ev) or 1.0
        assert np.max(np.abs(np.diff(ev, axis=0))) < 0.05 * scale


class TestMaxRealEigenvalue:
    def test_agrees_with_the_leading_track(self, toy):
        EN = np.logspace(1.5, 3.0, 25)
        ev = compute_eigenvalues(toy, EN, 1.0, 293.0)
        for i, en in enumerate(EN):
            assert max_real_eigenvalue(toy, en, 1.0, 293.0) == pytest.approx(
                float(np.max(np.real(ev[i]))), rel=1e-10
            )


class TestPressureScan:
    def test_shape_and_reduced_scaling(self, toy):
        """E4/E5: shape (M, P) and values reduced by the number density."""
        p_arr = np.array([0.1, 1.0, 10.0])
        EN = np.logspace(2.0, 3.0, 15)
        T = 293.0
        scan = compute_pressure_scan(toy, EN, p_arr, T, 0)
        assert scan.shape == (len(EN), len(p_arr))
        for pi, p in enumerate(p_arr):
            N = p * 1e5 / (kB * T)
            direct = np.real(compute_eigenvalues(toy, EN, p, T)[:, 0]) / N
            assert scan[:, pi] == pytest.approx(direct, rel=1e-12)

    def test_reduced_eigenvalue_is_pressure_independent(self, toy):
        """
        The toy coefficients scale linearly with p, so lambda/N must not.

        This is the similarity property the reduced plot exists to show.
        """
        EN = np.logspace(2.0, 3.0, 12)
        scan = compute_pressure_scan(toy, EN, np.array([0.5, 5.0]), 293.0, 0)
        assert scan[:, 0] == pytest.approx(scan[:, 1], rel=1e-10)
