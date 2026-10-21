# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The classical Townsend mechanism, verified against Paschen's law.

``mechanisms/paschen/paschen.json`` has no attachment, no detachment and no
photon feedback, so its inception condition is eq_standard_paschen and the
breakdown voltage has the closed form

    U = B (pd) / [ln(A pd) - ln(ln(1 + 1/gamma))].

The solver reaches the same answer through the full machinery -- the
compiled mechanism, A_aug, propagation, boundary matrix, criterion root --
with no shortcut for the analytic case, so agreement exercises the whole
chain.  The closed form is written out here from A, B and gamma, not taken
from the mechanism.
"""

import functools
import json
from pathlib import Path

import numpy as np
import pytest

from incept1d.constants import kB
from incept1d.declarative import MechanismFileError
from incept1d.fields import FieldDistribution
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import load_mechanism
from incept1d.solver import riccati_criterion

HERE = Path(__file__).resolve().parent.parent / "mechanisms" / "paschen"
MECH = str(HERE / "paschen.json")
with open(HERE / "gases.json") as _fh:
    GASES = {c["label"]: c for c in json.load(_fh)["configurations"]}
UNIFORM = FieldDistribution("uniform")
#: [cm^-1 Torr^-1] -> [m^-1 bar^-1]
PER_CM_TORR = 100.0 * 750.061682704


@pytest.fixture(scope="module")
def det():
    return functools.partial(riccati_criterion, field_dist=UNIFORM)


def _load(gas, **cfg):
    return load_mechanism(MECH, {**GASES[gas], **cfg})


def _voltage(EN, pd, T=293.0):
    """Breakdown voltage [V] from the reduced field and pd [bar m]."""
    return EN * pd * 1e-21 / (kB * T) * 1e5


def _coefficients(gas, gamma=None):
    """(A [m^-1 bar^-1], B [V m^-1 bar^-1], gamma) of a configuration."""
    c = GASES[gas]
    return c["A"] * PER_CM_TORR, c["B"] * PER_CM_TORR, gamma or c["gamma0"]


def paschen_voltage(pd, A, B, gamma):
    """U = B pd / [ln(A pd) - ln(ln(1 + 1/gamma))], NaN below the asymptote."""
    denom = np.log(A * pd) - np.log(np.log1p(1.0 / gamma))
    return B * pd / denom if denom > 0.0 else float("nan")


def paschen_minimum(A, B, gamma):
    """((pd)_min, U_min) = (e/A, e B/A) ln(1 + 1/gamma)."""
    L = np.log1p(1.0 / gamma)
    return np.e * L / A, np.e * B * L / A


class TestClosedForm:
    @pytest.mark.parametrize("gas", sorted(GASES))
    @pytest.mark.parametrize("pd_mm", [0.05, 0.5, 5.0, 50.0])
    def test_matches_paschens_law(self, gas, pd_mm, det):
        """The solver reproduces U(pd) exactly, for every gas."""
        mod = _load(gas)
        p, T, pd = 1.0, 293.0, pd_mm * 1e-3
        roots = find_all_breakdown_EN(pd, mod, p, T, first_only=True, det_fn=det)
        assert roots, f"{gas}: no inception at pd = {pd_mm} bar·mm"
        assert _voltage(roots[0], pd, T) == pytest.approx(
            paschen_voltage(pd, *_coefficients(gas)), rel=1e-9
        )

    @pytest.mark.parametrize("gas", sorted(GASES))
    def test_minimum_matches_the_analytic_minimum(self, gas, det):
        """
        The Paschen minimum, located by sweeping, matches
        (pd)_min = e/A ln(1 + 1/gamma),  U_min = e B/A ln(1 + 1/gamma).
        """
        mod = _load(gas)
        pd_min, U_min = paschen_minimum(*_coefficients(gas))
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

    @pytest.mark.parametrize("gas", sorted(GASES))
    def test_no_solution_below_the_left_asymptote(self, gas, det):
        """
        Below pd = ln(1 + 1/gamma)/A the denominator of Paschen's law is not
        positive and the gap cannot break down at any voltage.  The solver must
        agree, rather than inventing a root.
        """
        mod = _load(gas)
        A, B, gamma = _coefficients(gas)
        pd_asym = np.log1p(1.0 / gamma) / A
        assert np.isnan(paschen_voltage(pd_asym * 0.5, A, B, gamma))
        roots = find_all_breakdown_EN(
            pd_asym * 0.5, mod, 1.0, 293.0, first_only=True, det_fn=det
        )
        assert not roots


class TestGammaDependence:
    def test_solver_follows_the_gamma_override(self, det):
        """A configuration key must actually reach the solve."""
        p, T, pd = 1.0, 293.0, 1e-3
        volts = []
        for gamma in (1e-3, 1e-1):
            mod = _load("Air", gamma0=gamma)
            r = find_all_breakdown_EN(pd, mod, p, T, first_only=True, det_fn=det)
            assert r
            volts.append(_voltage(r[0], pd, T))
            assert volts[-1] == pytest.approx(
                paschen_voltage(pd, *_coefficients("Air", gamma)), rel=1e-9
            )
        assert volts[0] > volts[1], "a better emitter must break down sooner"


class TestMechanismData:
    @pytest.mark.parametrize(
        "gas, A, B",
        [("Helium", 3.0, 34.0), ("Argon", 12.0, 180.0), ("Air", 15.0, 365.0)],
    )
    def test_tabulated_coefficients(self, gas, A, B):
        """Guards the Raizer values against an accidental edit."""
        assert (GASES[gas]["A"], GASES[gas]["B"]) == (A, B)

    def test_a_configuration_for_another_model_is_refused(self):
        """The old 'gas' key, or any other, is not a parameter of the file."""
        with pytest.raises(MechanismFileError, match="unknown key"):
            load_mechanism(MECH, {"gas": "xenon"})

    def test_no_attachment_and_no_photons(self):
        """The model's defining simplifications, asserted rather than assumed."""
        mod = _load("Air")
        assert mod.eta(300.0, 1.0, 293.0) == 0.0
        assert mod.get_B(300.0, 1.0, 293.0).shape[1] == 0
        assert mod.get_kappa(1.0, 293.0).size == 0
        assert mod.get_Pi_minus().shape[0] == 0
