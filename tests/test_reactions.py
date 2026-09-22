# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Reaction-string parsing and R-matrix assembly."""

import numpy as np
import pytest

from incept1d.reactions import (
    build_R,
    build_R_from_compiled,
    compile_reactions,
    parse_reaction,
)

SPECIES = ["e", "N2+", "O-"]


def _const(k):
    """A rate callable with the mechanism signature, returning a constant."""
    return lambda EN, p, T: k


class TestParsing:
    def test_cation_names_survive_the_tokeniser(self):
        """R1: splitting on ' + ' must not break 'N2+' into 'N2' and ''."""
        reactants, products = parse_reaction("e + N2 -> 2e + N2+")
        assert reactants == {"e": 1, "N2": 1}
        assert products == {"e": 2, "N2+": 1}

    def test_leading_stoichiometric_coefficients(self):
        """R2: '2e' and '2O2' carry integer coefficients."""
        reactants, products = parse_reaction("e + 2O2 -> O2- + O2")
        assert reactants == {"e": 1, "O2": 2}
        assert products == {"O2-": 1, "O2": 1}

    def test_anion_names_survive_the_tokeniser(self):
        reactants, products = parse_reaction("O- + N2 -> e + N2O")
        assert reactants == {"O-": 1, "N2": 1}
        assert products == {"e": 1, "N2O": 1}

    def test_repeated_species_accumulate(self):
        reactants, _ = parse_reaction("O- + O2 + M -> O3- + M")
        assert reactants == {"O-": 1, "O2": 1, "M": 1}


class TestRMatrix:
    def test_ionization_signs(self):
        """R3: e + N2 -> 2e + N2+ adds one electron and one ion per electron."""
        R = build_R([("e + N2 -> 2e + N2+", _const(7.0))], SPECIES, 100.0, 1.0, 293.0)
        assert R[0, 0] == pytest.approx(7.0)  # net +1 electron
        assert R[1, 0] == pytest.approx(7.0)  # +1 N2+
        assert R[2, 0] == 0.0

    def test_attachment_signs(self):
        """R4: e + O2 -> O- + O destroys the electron and makes an anion."""
        R = build_R([("e + O2 -> O- + O", _const(3.0))], SPECIES, 100.0, 1.0, 293.0)
        assert R[0, 0] == pytest.approx(-3.0)
        assert R[2, 0] == pytest.approx(3.0)
        assert R[1, 0] == 0.0

    def test_untracked_background_species_are_ignored(self):
        """R5: N2, O2, M are not in SPECIES and contribute no rows."""
        R = build_R([("e + N2 -> 2e + N2+", _const(1.0))], SPECIES, 100.0, 1.0, 293.0)
        assert R.shape == (3, 3)
        assert np.count_nonzero(R) == 2

    def test_detachment_drives_from_the_anion_column(self):
        """The driver column is the tracked reactant, here O-."""
        R = build_R([("O- + N2 -> e + N2O", _const(2.0))], SPECIES, 100.0, 1.0, 293.0)
        assert R[0, 2] == pytest.approx(2.0)  # e produced per O-
        assert R[2, 2] == pytest.approx(-2.0)  # O- consumed

    def test_two_tracked_reactants_are_rejected(self):
        """R6: the one-driver rule is enforced with a helpful message."""
        with pytest.raises(ValueError, match="exactly one tracked species"):
            build_R([("e + O- -> 2e", _const(1.0))], SPECIES, 100.0, 1.0, 293.0)


class TestMultipliers:
    RX = [("e + N2 -> 2e + N2+", _const(5.0))]

    @pytest.mark.parametrize(
        "key", ["e + N2 -> 2e + N2+", "e+N2->2e+N2+", "e + N2 → 2e + N2+"]
    )
    def test_key_normalisation(self, key):
        """R7: spacing and the unicode arrow must not affect matching."""
        R = build_R(self.RX, SPECIES, 100.0, 1.0, 293.0, multipliers={key: 2.0})
        assert R[0, 0] == pytest.approx(10.0)

    def test_multiplier_scales_linearly(self):
        """R9: the multiplier is a plain factor on the rate."""
        for m in (0.5, 2.0, 3.0):
            R = build_R(
                self.RX,
                SPECIES,
                100.0,
                1.0,
                293.0,
                multipliers={"e + N2 -> 2e + N2+": m},
            )
            assert R[0, 0] == pytest.approx(5.0 * m)

    def test_zero_multiplier_switches_a_reaction_off(self):
        """R9b: this is how the configuration files disable reactions."""
        R = build_R(
            self.RX,
            SPECIES,
            100.0,
            1.0,
            293.0,
            multipliers={"e + N2 -> 2e + N2+": 0.0},
        )
        assert np.count_nonzero(R) == 0

    def test_unknown_key_warns_and_is_ignored(self, capsys):
        """R8: a typo in a config file must be visible, not silent."""
        R = build_R(
            self.RX, SPECIES, 100.0, 1.0, 293.0, multipliers={"e + Xe -> 2e + Xe+": 2.0}
        )
        assert "does not match" in capsys.readouterr().out
        assert R[0, 0] == pytest.approx(5.0)


class TestCompiledFastPath:
    RX = [
        ("e + N2 -> 2e + N2+", _const(5.0)),
        ("e + O2 -> O- + O", _const(3.0)),
        ("O- + N2 -> e + N2O", _const(2.0)),
    ]

    def test_matches_build_R(self):
        """R10: the pre-compiled path is a drop-in replacement."""
        compiled = compile_reactions(self.RX, SPECIES)
        a = build_R(self.RX, SPECIES, 120.0, 1.0, 293.0)
        b = build_R_from_compiled(compiled, len(SPECIES), 120.0, 1.0, 293.0)
        assert a == pytest.approx(b)

    def test_matches_build_R_with_multipliers(self):
        """R10b: including multiplier handling."""
        mult = {"e+O2->O-+O": 0.0, "e + N2 -> 2e + N2+": 1.5}
        compiled = compile_reactions(self.RX, SPECIES)
        a = build_R(self.RX, SPECIES, 120.0, 1.0, 293.0, multipliers=mult)
        b = build_R_from_compiled(
            compiled, len(SPECIES), 120.0, 1.0, 293.0, multipliers=mult
        )
        assert a == pytest.approx(b)

    def test_compile_rejects_two_tracked_reactants(self):
        with pytest.raises(ValueError, match="exactly one tracked species"):
            compile_reactions([("e + O- -> 2e", _const(1.0))], SPECIES)
