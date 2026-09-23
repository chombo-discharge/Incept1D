# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Transport-table export for chombo-discharge: structure, not physics."""

import numpy as np
import pytest

from incept1d.chombo import build_header, generate, load_raw_mechanism

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def raw(air_path_module):
    return load_raw_mechanism(air_path_module)


@pytest.fixture(scope="module")
def air_path_module():
    import os

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(here, "mechanisms", "air", "air_pancheshnyi.py")


@pytest.fixture(scope="module")
def table(raw):
    EN = np.logspace(1.0, 3.0, 12)
    return generate(raw, EN, 1.0, 300.0, {})


class TestStructure:
    def test_one_header_per_column(self, table):
        """K1: a mismatch here silently mislabels every downstream column."""
        columns, headers = table
        assert len(columns) == len(headers)

    def test_all_columns_have_the_same_length(self, table):
        columns, _ = table
        assert len({len(c) for c in columns}) == 1

    def test_documented_column_order(self, table, raw):
        """K2: E/N, alpha/N, eta/N, energy, then two per species, then rates."""
        columns, headers = table
        n = len(raw.SPECIES)
        assert "E/N" in headers[0]
        assert "alpha" in headers[1].lower()
        assert "eta" in headers[2].lower()
        assert len(headers) == 4 + 2 * n + len(raw.REACTIONS)

    def test_mobility_and_diffusion_pairs_name_their_species(self, table, raw):
        _, headers = table
        block = headers[4 : 4 + 2 * len(raw.SPECIES)]
        for i, sp in enumerate(raw.SPECIES):
            assert sp in block[2 * i] and sp in block[2 * i + 1]

    def test_reduced_field_column_is_the_input_grid(self, raw):
        EN = np.logspace(1.0, 3.0, 7)
        columns, _ = generate(raw, EN, 1.0, 300.0, {})
        assert columns[0] == pytest.approx(EN)


class TestHeader:
    def test_header_mentions_the_provenance(self, raw, air_path_module):
        _, headers = generate(raw, np.logspace(1, 2, 4), 1.0, 300.0, {})
        text = build_header(
            air_path_module,
            None,
            None,
            1.0,
            300.0,
            headers,
            pre_exec_vars={},
            EN_min=10.0,
            EN_max=100.0,
            num_EN=4,
            species=raw.SPECIES,
        )
        assert "air_pancheshnyi.py" in text
        assert "300" in text

    def test_three_body_reactions_are_flagged(self, table):
        """K3: units differ (m^6/s), so the distinction must be visible."""
        _, headers = table
        rate_headers = " ".join(headers[4:])
        assert "m^6" in rate_headers or "three" in rate_headers.lower()


class TestRawLoader:
    def test_rejects_a_file_without_reactions(self, tmp_path):
        """K4: the raw loader has its own, different interface requirement."""
        p = tmp_path / "no_reactions.py"
        p.write_text("SPECIES = ['e']\nELECTRON_INDEX = 0\n")
        with pytest.raises(AttributeError, match="REACTIONS"):
            load_raw_mechanism(str(p))

    def test_loads_the_air_mechanism(self, raw):
        assert raw.SPECIES and hasattr(raw, "REACTIONS")
        assert callable(raw.get_V)
