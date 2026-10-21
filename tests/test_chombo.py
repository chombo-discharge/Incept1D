# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Transport-table export for chombo-discharge: structure, and what each column is."""

from pathlib import Path

import numpy as np
import pytest

from incept1d.chombo import (
    build_header,
    generate,
    generate_declarative,
    load_raw_mechanism,
)
from incept1d.constants import kB
from incept1d.declarative import MechanismFileError
from incept1d.mechanism import load_mechanism, read_json_configs

ROOT = Path(__file__).resolve().parent.parent
AIR = ROOT / "mechanisms" / "air" / "pancheshnyi" / "air_pancheshnyi.json"


@pytest.fixture(scope="module")
def air():
    return load_mechanism(str(AIR))


@pytest.fixture(scope="module")
def table(air):
    EN = np.logspace(1.0, 3.0, 12)
    return generate_declarative(air, EN, 1.0, 300.0)


class TestStructure:
    def test_one_header_per_column(self, table):
        """K1: a mismatch here silently mislabels every downstream column."""
        columns, headers = table
        assert len(columns) == len(headers)

    def test_all_columns_have_the_same_length(self, table):
        columns, _ = table
        assert len({len(c) for c in columns}) == 1

    def test_documented_column_order(self, table, air):
        """K2: E/N, alpha/N, eta/N, energy, then two per species, then rates."""
        columns, headers = table
        n = len(air.SPECIES)
        assert "E/N" in headers[0]
        assert "alpha" in headers[1].lower()
        assert "eta" in headers[2].lower()
        assert len(headers) == 4 + 2 * n + len(air.reaction_strings)

    def test_mobility_and_diffusion_pairs_name_their_species(self, table, air):
        _, headers = table
        block = headers[4 : 4 + 2 * len(air.SPECIES)]
        for i, sp in enumerate(air.SPECIES):
            assert sp in block[2 * i] and sp in block[2 * i + 1]

    def test_reduced_field_column_is_the_input_grid(self, air):
        EN = np.logspace(1.0, 3.0, 7)
        columns, _ = generate_declarative(air, EN, 1.0, 300.0)
        assert columns[0] == pytest.approx(EN)

    def test_three_body_reactions_are_flagged(self, table):
        """K3: units differ (m^6/s), so the distinction must be visible."""
        _, headers = table
        three = [h for h in headers if "THREE-BODY" in h]
        assert three and all("m^6/s" in h for h in three)


class TestValues:
    def test_columns_are_what_the_mechanism_says(self, air):
        """K4: alpha/N, eta/N and every k times its neutral densities give R."""
        EN, p, T = np.array([80.0, 150.0]), 1.0, 300.0
        N = p * 1e5 / (kB * T)
        columns, _ = generate_declarative(air, EN, p, T)
        first = 4 + 2 * len(air.SPECIES)
        for j, en in enumerate(EN):
            assert columns[1][j] == pytest.approx(air.alpha(en, p, T) / N, rel=1e-14)
            assert columns[2][j] == pytest.approx(air.eta(en, p, T) / N, rel=1e-14)
            # The ionization of N2 alone is R[N2+, e].
            K = columns[first][j] * air.neutral_density_factor(0, en, p, T)
            assert K == pytest.approx(air.get_R(en, p, T)[1, 0], rel=1e-14)

    def test_multipliers_reach_the_table(self):
        """K5: the configuration's multiplier scales that column, only."""
        cfg = {"reaction_multipliers": {"O- + N2 -> e + N2O": 0.5}}
        EN = np.logspace(1.5, 2.5, 5)
        base, headers = generate_declarative(load_mechanism(str(AIR)), EN, 1.0, 300.0)
        half, _ = generate_declarative(load_mechanism(str(AIR), cfg), EN, 1.0, 300.0)
        for a, b, h in zip(base, half, headers):
            expected = 0.5 * a if "'O- + N2 -> e + N2O'" in h else a
            np.testing.assert_allclose(b, expected, rtol=1e-15)

    def test_a_file_without_a_mean_energy_is_refused(self):
        """K6: the column is not invented when the file does not define it."""
        paschen = load_mechanism(str(ROOT / "mechanisms" / "paschen" / "paschen.json"))
        with pytest.raises(MechanismFileError, match="mean_energy"):
            generate_declarative(paschen, np.logspace(1, 2, 3), 1.0, 300.0)

    def test_header_names_the_mechanism_and_its_data(self, table, air):
        _, headers = table
        cfg = read_json_configs([str(AIR.parent / "databases.json")])[0]
        mech = load_mechanism(str(AIR), cfg)
        text = build_header(
            str(AIR),
            None,
            None,
            1.0,
            300.0,
            headers,
            EN_min=10.0,
            EN_max=100.0,
            num_EN=4,
            species=air.SPECIES,
            files=mech.files,
        )
        assert "air_pancheshnyi.json" in text and "phelps.txt" in text
        assert "300" in text


_PYTHON_MECHANISM = """
import numpy as np
SPECIES = ["e", "A+"]
ELECTRON_INDEX = 0
xA = 1.0
def k_ion(EN):
    return 1e-15 * np.exp(-400.0 / EN)
REACTIONS = [("e + A -> 2e + A+", lambda EN, p, T: k_ion(EN) * xA * p * 1e5 / (1.380649e-23 * T))]
def ElectronMeanEnergy(EN):
    return 0.05 * EN
def get_V(EN, p, T):
    return np.diag([4e23, -2e21]) * EN * 1e-21
"""


class TestPythonMechanism:
    def test_rejects_a_file_without_reactions(self, tmp_path):
        """K7: the raw loader has its own, different interface requirement."""
        p = tmp_path / "no_reactions.py"
        p.write_text("SPECIES = ['e']\nELECTRON_INDEX = 0\n")
        with pytest.raises(AttributeError, match="REACTIONS"):
            load_raw_mechanism(str(p))

    def test_recovers_the_rate_coefficient_and_derives_alpha(self, tmp_path):
        """K8: k = rate / neutral density; alpha/N = sum of ionization / (v_e N)."""
        path = tmp_path / "tiny.py"
        path.write_text(_PYTHON_MECHANISM)
        mod = load_raw_mechanism(str(path))
        EN = np.array([100.0, 300.0])
        columns, headers = generate(mod, EN, 1.0, 300.0, {})
        k = 1e-15 * np.exp(-400.0 / EN)
        assert columns[-1] == pytest.approx(k, rel=1e-14)
        assert columns[1] == pytest.approx(k / (4e23 * EN * 1e-21), rel=1e-14)
        assert columns[3] == pytest.approx(0.05 * EN)
