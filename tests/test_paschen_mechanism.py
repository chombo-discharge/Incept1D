# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The classical Townsend mechanism, verified against Paschen's law.

``mechanisms/paschen/paschen.py`` has no attachment, no detachment and no
photon feedback, so its inception condition is eq_standard_paschen and the
breakdown voltage has the closed form

    U = B (pd) / [ln(A pd) - ln(ln(1 + 1/gamma))].

The solver reaches the same answer through the full machinery -- A_aug,
propagation, boundary matrix, determinant root -- with no shortcut for the
analytic case, so agreement exercises the whole chain.
"""

import functools
import os

import numpy as np
import pytest

from incept1d.constants import kB
from incept1d.fields import FieldDistribution
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import load_mechanism
from incept1d.solver import riccati_criterion

MECH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "mechanisms",
    "paschen",
    "paschen.py",
)
GASES = ["helium", "argon", "air"]
UNIFORM = FieldDistribution("uniform")


@pytest.fixture(scope="module")
def det():
    return functools.partial(riccati_criterion, field_dist=UNIFORM)


def _load(gas, **cfg):
    return load_mechanism(MECH, {"gas": gas, **cfg})


def _voltage(EN, pd, T=293.0):
    """Breakdown voltage [V] from the reduced field and pd [bar m]."""
    return EN * pd * 1e-21 / (kB * T) * 1e5


class TestClosedForm:
    @pytest.mark.parametrize("gas", GASES)
    @pytest.mark.parametrize("pd_mm", [0.05, 0.5, 5.0, 50.0])
    def test_matches_paschens_law(self, gas, pd_mm, det):
        """The solver reproduces U(pd) exactly, for every gas."""
        mod = _load(gas)
        raw = mod._mod
        p, T, pd = 1.0, 293.0, pd_mm * 1e-3
        roots = find_all_breakdown_EN(pd, mod, p, T, first_only=True, det_fn=det)
        assert roots, f"{gas}: no inception at pd = {pd_mm} bar·mm"
        assert _voltage(roots[0], pd, T) == pytest.approx(
            raw.paschen_voltage(pd), rel=1e-9
        )

    @pytest.mark.parametrize("gas", GASES)
    def test_minimum_matches_the_analytic_minimum(self, gas, det):
        """
        The Paschen minimum, located by sweeping, matches
        (pd)_min = e/A ln(1 + 1/gamma),  U_min = e B/A ln(1 + 1/gamma).
        """
        mod = _load(gas)
        raw = mod._mod
        pd_min, U_min = raw.paschen_minimum()
        p, T = 1.0, 293.0
        pd_arr = np.logspace(np.log10(pd_min * 0.4), np.log10(pd_min * 3.0), 60)
        best = min(
            (
                (_voltage(r[0], pd, T), pd)
                for pd in pd_arr
                for r in [
                    find_all_breakdown_EN(pd, mod, p, T, first_only=True, det_fn=det)
                ]
                if r
            ),
            default=None,
        )
        assert best is not None
        U_found, pd_found = best
        assert U_found == pytest.approx(U_min, rel=1e-3)
        assert pd_found == pytest.approx(pd_min, rel=0.1)

    @pytest.mark.parametrize("gas", GASES)
    def test_no_solution_below_the_left_asymptote(self, gas, det):
        """
        Below pd = ln(1 + 1/gamma)/A the denominator of Paschen's law is not
        positive and the gap cannot break down at any voltage.  The solver must
        agree, rather than inventing a root.
        """
        mod = _load(gas)
        raw = mod._mod
        gamma = float(mod.get_gamma_plus(100.0, 1.0, 293.0)[0])
        pd_asym = np.log1p(1.0 / gamma) / raw.A_SI
        assert np.isnan(raw.paschen_voltage(pd_asym * 0.5))
        roots = find_all_breakdown_EN(
            pd_asym * 0.5, mod, 1.0, 293.0, first_only=True, det_fn=det
        )
        assert not roots


class TestGammaDependence:
    @pytest.mark.parametrize("gamma", [1e-3, 1e-2, 1e-1])
    def test_minimum_scales_with_log_gamma(self, gamma):
        """Both (pd)_min and U_min go as ln(1 + 1/gamma)."""
        raw = _load("air", gamma0=gamma)._mod
        pd_min, U_min = raw.paschen_minimum()
        L = np.log1p(1.0 / gamma)
        assert pd_min == pytest.approx(np.e * L / raw.A_SI, rel=1e-12)
        assert U_min == pytest.approx(np.e * raw.B_SI * L / raw.A_SI, rel=1e-12)

    def test_solver_follows_the_gamma_override(self, det):
        """A configuration key must actually reach the solve."""
        p, T, pd = 1.0, 293.0, 1e-3
        volts = []
        for gamma in (1e-3, 1e-1):
            mod = _load("air", gamma0=gamma)
            r = find_all_breakdown_EN(pd, mod, p, T, first_only=True, det_fn=det)
            assert r
            volts.append(_voltage(r[0], pd, T))
            assert volts[-1] == pytest.approx(mod._mod.paschen_voltage(pd), rel=1e-9)
        assert volts[0] > volts[1], "a better emitter must break down sooner"


class TestMechanismData:
    def test_known_gases(self):
        raw = _load("air")._mod
        assert set(raw.GASES) == set(GASES)

    def test_unknown_gas_is_rejected(self):
        with pytest.raises(ValueError, match="Unknown gas"):
            _load("xenon")

    @pytest.mark.parametrize(
        "gas, A, B",
        [("helium", 3.0, 34.0), ("argon", 12.0, 180.0), ("air", 15.0, 365.0)],
    )
    def test_tabulated_coefficients(self, gas, A, B):
        """Guards the Raizer values against an accidental edit."""
        raw = _load(gas)._mod
        assert raw.GASES[gas][0] == A
        assert raw.GASES[gas][1] == B

    def test_no_attachment_and_no_photons(self):
        """The model's defining simplifications, asserted rather than assumed."""
        mod = _load("air")
        assert mod._mod.eta(300.0, 1.0, 293.0) == 0.0
        assert mod.get_B(300.0, 1.0, 293.0).shape[1] == 0
        assert mod.get_kappa(1.0, 293.0).size == 0
        assert mod.get_Pi_minus().shape[0] == 0
