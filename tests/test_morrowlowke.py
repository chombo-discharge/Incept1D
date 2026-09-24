# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The Morrow–Lowke air mechanism, ``mechanisms/air/morrowlowke/air_morrowlowke.py``.

The published fits are transcribed here a second time, independently of the
mechanism file, in the paper's own units (E/N in V cm², N in cm⁻³, cm/s), so
that a unit-conversion slip in the mechanism shows up as a mismatch.  The
strongest check is :class:`TestClosedFormLimit`: the scheme has no
detachment, so with photon feedback switched off it *is* the reduced model of
eq_standard_paschen.
"""

import functools
import math
import os

import numpy as np
import pytest
import scipy.optimize

from closed_form import standard_paschen_lhs
from incept1d.constants import kB
from incept1d.fields import FieldDistribution
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import REQUIRED_ATTRS, load_mechanism
from incept1d.solver import riccati_criterion

ML_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "mechanisms",
    "air",
    "morrowlowke",
    "air_morrowlowke.py",
)
P, T = 1.0, 293.0
N_CM3 = P * 1e5 / (kB * T) * 1e-6


# ── Independent transcription of Morrow & Lowke (1997), Appendix ─────────────


def pub_alpha_N(en):
    if en > 1.5e-15:
        return 2e-16 * math.exp(-7.248e-15 / en)
    return 6.619e-17 * math.exp(-5.593e-15 / en)


def pub_eta_N(en, N):
    eta2 = 8.889e-5 * en + 2.567e-19 if en > 1.05e-15 else 6.089e-4 * en - 2.893e-19
    return max(eta2, 0.0) + N * 4.7778e-59 * en**-1.2749


def pub_We(en):
    if en > 2e-15:
        return 7.4e21 * en + 7.1e6
    if en > 1e-16:
        return 1.03e22 * en + 1.3e6
    if en > 2.6e-17:
        return 7.2973e21 * en + 1.63e6
    return 6.87e22 * en + 3.38e4


@pytest.fixture(scope="module")
def ml():
    return load_mechanism(ML_PATH, {})


@pytest.fixture(scope="module")
def ml_raw(ml):
    return ml._mod


class TestInterface:
    def test_required_attributes(self, ml_raw):
        assert all(hasattr(ml_raw, a) for a in REQUIRED_ATTRS)

    def test_three_species(self, ml):
        assert ml.SPECIES == ["e", "M+", "M-"]
        assert ml.SPECIES[ml.ELECTRON_INDEX] == "e"

    def test_velocity_signs_follow_the_convention(self, ml):
        """Electrons and anions towards the anode (+x), cations to the cathode."""
        V = np.diag(ml.get_V(120.0, P, T))
        assert V[0] > 0.0 and V[1] < 0.0 and V[2] > 0.0

    def test_photoionization_makes_electron_ion_pairs(self, ml):
        B = ml.get_B(120.0, P, T)
        assert np.all(B[0] > 0.0)
        assert np.array_equal(B[0], B[1])
        assert np.all(B[2] == 0.0)


class TestPublishedFits:
    @pytest.mark.parametrize("EN_Td", [30.0, 80.0, 120.0, 200.0, 500.0])
    def test_townsend_coefficients(self, ml, EN_Td):
        """α and η in m⁻¹ equal the published per-cm² fits times N."""
        en = EN_Td * 1e-17
        assert ml.alpha(EN_Td, P, T) == pytest.approx(
            100.0 * pub_alpha_N(en) * N_CM3, rel=1e-10
        )
        assert ml.eta(EN_Td, P, T) == pytest.approx(
            100.0 * pub_eta_N(en, N_CM3) * N_CM3, rel=1e-10
        )

    @pytest.mark.parametrize("EN_Td", [1.0, 5.0, 50.0, 500.0])
    def test_electron_drift_velocity(self, ml, EN_Td):
        V = ml.get_V(EN_Td, P, T)
        assert V[0, 0] == pytest.approx(1e-2 * pub_We(EN_Td * 1e-17), rel=1e-12)

    def test_ion_velocities_at_the_reference_density(self, ml):
        """W+ = 2.34 E and W− = 2.7 E (cm/s, E in V/cm) at N = 2.5e19 cm⁻³."""
        EN_Td = 100.0
        E_Vcm = EN_Td * 1e-17 * 2.5e19
        V = ml.get_V(EN_Td, P, T)
        assert -V[1, 1] == pytest.approx(1e-2 * 2.34 * E_Vcm, rel=1e-12)
        assert V[2, 2] == pytest.approx(1e-2 * 2.7 * E_Vcm, rel=1e-12)

    def test_ion_velocity_scales_with_E_over_N(self, ml):
        """μN is the pressure-independent quantity, so W depends on E/N only."""
        assert np.allclose(ml.get_V(150.0, 0.1, T), ml.get_V(150.0, 10.0, T))

    def test_attachment_is_never_negative(self, ml):
        """The low-field η₂ branch crosses zero near 47.5 Td and is clipped."""
        for EN_Td in np.linspace(1.0, 60.0, 60):
            assert ml._mod.eta2(EN_Td, P, T) >= 0.0

    def test_three_body_attachment_scales_with_pressure_squared(self, ml_raw):
        r = ml_raw.eta3(100.0, 2.0, T) / ml_raw.eta3(100.0, 1.0, T)
        assert r == pytest.approx(4.0, rel=1e-12)

    def test_no_detachment(self, ml):
        """Nothing returns an electron from M-: column 2 of R is empty."""
        assert np.all(ml.get_R(120.0, P, T)[:, 2] == 0.0)


class TestContinuity:
    """The piecewise fits at their breakpoints, as published."""

    @pytest.mark.parametrize("EN_Td", [200.0, 2.6])
    def test_electron_drift_is_continuous(self, ml_raw, EN_Td):
        lo = ml_raw.ElectronDriftVelocity(EN_Td * (1 - 1e-12))
        hi = ml_raw.ElectronDriftVelocity(EN_Td * (1 + 1e-12))
        assert hi == pytest.approx(lo, rel=2e-4)

    def test_electron_drift_step_at_10_Td_is_in_the_published_fit(self, ml_raw):
        """|W_e| drops by ~1.3 % across 1e-16 V cm²; inherent in the fit."""
        lo = ml_raw.ElectronDriftVelocity(10.0 * (1 - 1e-12))
        hi = ml_raw.ElectronDriftVelocity(10.0 * (1 + 1e-12))
        assert hi / lo - 1.0 == pytest.approx(-0.0126, abs=1e-3)

    def test_alpha_is_continuous(self, ml):
        lo, hi = ml.alpha(150.0 * (1 - 1e-12), P, T), ml.alpha(
            150.0 * (1 + 1e-12), P, T
        )
        assert hi == pytest.approx(lo, rel=3e-3)

    def test_two_body_attachment_is_continuous(self, ml_raw):
        lo = ml_raw.eta2(105.0 * (1 - 1e-12), P, T)
        hi = ml_raw.eta2(105.0 * (1 + 1e-12), P, T)
        assert hi == pytest.approx(lo, rel=1e-4)

    def test_negative_ion_mobility_step_at_50_Td_is_in_the_published_fit(self, ml):
        lo = ml.get_V(50.0 * (1 - 1e-9), P, T)[2, 2] / (50.0 * (1 - 1e-9))
        hi = ml.get_V(50.0 * (1 + 1e-9), P, T)[2, 2] / (50.0 * (1 + 1e-9))
        assert hi / lo == pytest.approx(2.7 / 1.86, rel=1e-6)


class TestCriticalField:
    def test_alpha_equals_eta_near_108_Td(self, ml):
        root = scipy.optimize.brentq(
            lambda e: ml.alpha(e, P, T) - ml.eta(e, P, T), 50.0, 300.0
        )
        assert 100.0 < root < 120.0

    def test_three_body_attachment_raises_the_critical_field(self, ml):
        """η₃ ∝ N², so the reduced critical field grows with pressure."""

        def crit(p):
            return scipy.optimize.brentq(
                lambda e: ml.alpha(e, p, T) - ml.eta(e, p, T), 50.0, 300.0
            )

        assert crit(0.1) < crit(1.0) < crit(10.0)


@pytest.fixture(scope="module")
def townsend():
    """Photon feedback off: ionization, attachment and ion emission only."""
    return load_mechanism(ML_PATH, {"cone_angle": 0.0, "xi_emit": 0.0})


# Not marked slow: the fits are analytic, so there is no swarm data to load.
class TestClosedFormLimit:
    """
    Without photon feedback there is nothing left but ionization, attachment
    and ion-induced emission, which is the reduced model of
    eq_standard_paschen.  The whole boundary-value solve must reproduce it.
    """

    @pytest.mark.parametrize("pd_mm", [1.0, 10.0, 100.0])
    def test_reproduces_standard_paschen(self, townsend, pd_mm):
        det = functools.partial(
            riccati_criterion, field_dist=FieldDistribution("uniform")
        )
        pd = pd_mm * 1e-3
        roots = find_all_breakdown_EN(pd, townsend, P, T, first_only=True, det_fn=det)
        assert roots, f"no inception found at pd = {pd_mm} bar·mm"
        EN = roots[0]
        a, e = townsend.alpha(EN, P, T), townsend.eta(EN, P, T)
        gamma = float(townsend.get_gamma_plus(EN, P, T)[0])
        assert standard_paschen_lhs(a, e, gamma, pd / P) == pytest.approx(1.0, rel=2e-2)
