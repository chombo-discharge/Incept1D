# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Consistency between the gas mixture of the transport data and the neutral
densities of the air mechanisms (mechanisms/air/bolsig_conditions.py).
"""

import os

import pytest

from conftest import AIR_PATH, REPO_ROOT
from incept1d.mechanism import load_helper, load_mechanism

AIR = os.path.join(REPO_ROOT, "mechanisms", "air")
LXCAT = os.path.join(AIR, "lxcat")
LISBON = os.path.join(LXCAT, "lisbon.txt")  # records N2 0.79, O2 0.21, 300 K
BIAGI = os.path.join(LXCAT, "biagi.txt")  # LXCat output: Tgas only, no mixture

conditions = load_helper(os.path.join(AIR, "bolsig_conditions.py"))


class TestReadConditions:
    def test_bolsig_output_records_mixture_and_temperature(self):
        fractions, T = conditions.read_conditions(LISBON)
        assert fractions == {"N2": 0.79, "O2": 0.21}
        assert T == 300.0

    def test_lxcat_output_records_temperature_only(self):
        fractions, T = conditions.read_conditions(BIAGI)
        assert fractions == {}
        assert T == 300.0


class TestResolveComposition:
    def test_table_mixture_is_the_default(self):
        assert conditions.resolve_composition(LISBON, None, None, (0.8, 0.2)) == (0.79, 0.21)

    def test_matching_configuration_is_accepted(self):
        assert conditions.resolve_composition(LISBON, 0.79, 0.21, (0.8, 0.2)) == (0.79, 0.21)

    def test_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ from the mixture"):
            conditions.resolve_composition(LISBON, 0.8, 0.2, (0.8, 0.2))

    def test_unrecorded_mixture_uses_configuration(self):
        assert conditions.resolve_composition(BIAGI, 0.79, 0.21, (0.8, 0.2)) == (0.79, 0.21)

    def test_unrecorded_mixture_without_configuration_warns(self):
        conditions._WARNED.clear()
        with pytest.warns(UserWarning, match="does not record the gas mixture"):
            assert conditions.resolve_composition(BIAGI, None, None, (0.8, 0.2)) == (0.8, 0.2)


class TestCheckTemperature:
    def test_small_difference_is_silent(self, recwarn):
        conditions.check_temperature(293.0, 300.0, LISBON)
        assert not recwarn.list

    def test_large_difference_warns(self):
        conditions._WARNED.clear()
        with pytest.warns(UserWarning, match="differs from the 300 K"):
            conditions.check_temperature(400.0, 300.0, LISBON)


@pytest.mark.slow
class TestMechanism:
    def test_default_mixture_follows_the_table(self):
        mech = load_mechanism(AIR_PATH, {"cross_sections": LISBON})
        assert (mech._mod.xN2, mech._mod.xO2) == (0.79, 0.21)

    def test_conflicting_mixture_is_rejected(self):
        with pytest.raises(ValueError, match="differ from the mixture"):
            load_mechanism(AIR_PATH, {"cross_sections": LISBON, "x_n2": 0.8, "x_o2": 0.2})
