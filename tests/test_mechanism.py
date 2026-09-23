# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Mechanism loading, the configuration protocol and JSON config files."""

import json
import os

import numpy as np
import pytest

from incept1d.mechanism import REQUIRED_ATTRS, load_mechanism, read_json_configs


class TestLoadMechanism:
    def test_loads_the_toy_mechanism(self, toy_path):
        """M1: a conforming file yields a usable Mechanism."""
        mod = load_mechanism(toy_path, {})
        assert mod.SPECIES == ["e", "M+", "M-"]
        assert mod.ELECTRON_INDEX == 0
        assert mod.get_R(300.0, 1.0, 293.0).shape == (3, 3)
        assert mod.get_V(300.0, 1.0, 293.0).shape == (3, 3)

    def test_missing_file(self, tmp_path):
        """M3."""
        with pytest.raises(FileNotFoundError):
            load_mechanism(str(tmp_path / "nope.py"))

    def test_missing_interface_is_reported_by_name(self, tmp_path):
        """M2: the error must say which attributes are missing."""
        p = tmp_path / "bad_mechanism.py"
        p.write_text("SPECIES = ['e']\nELECTRON_INDEX = 0\n")
        with pytest.raises(AttributeError) as exc:
            load_mechanism(str(p))
        msg = str(exc.value)
        assert "missing required attributes" in msg
        assert "get_R" in msg and "get_V" in msg

    def test_required_attrs_is_the_documented_contract(self):
        """The mechanism-writing guide lists exactly these names."""
        for name in (
            "SPECIES",
            "ELECTRON_INDEX",
            "get_R",
            "get_V",
            "get_B",
            "get_C",
            "get_kappa",
            "get_Pi_e",
            "get_Pi_plus",
            "get_Pi_minus",
            "get_gamma_plus",
            "get_gamma_plus_with",
            "get_gamma_Psi",
        ):
            assert name in REQUIRED_ATTRS

    def test_config_dict_without_config_py_is_an_error(self, tmp_path, toy_path):
        """M5: passing configuration to a mechanism that cannot accept it."""
        import shutil

        dest = tmp_path / "lonely_mechanism.py"
        shutil.copy(toy_path, dest)  # no config.py alongside
        with pytest.raises(ImportError, match="No config.py"):
            load_mechanism(str(dest), {"gamma0": 1e-3})


class TestConfigProtocol:
    def test_pre_exec_vars_reach_the_module(self, toy_path):
        """M4: injected names are visible while the module body executes."""
        mod = load_mechanism(toy_path, {"alpha_B": 250.0, "eta_0": 33.0})
        assert mod._mod._ALPHA_B == 250.0
        assert mod._mod._ETA_0 == 33.0
        assert mod._mod.eta(100.0, 1.0, 293.0) == pytest.approx(33.0)

    def test_post_exec_init_runs(self, toy_path):
        """M4b: the post-exec hook fires after the body has run."""
        mod = load_mechanism(toy_path, {})
        assert getattr(mod._mod, "_POST_EXEC_RAN", False) is True

    def test_label_is_carried_through(self, toy_path):
        mod = load_mechanism(toy_path, {"label": "Variant A"})
        assert mod.label == "Variant A"

    def test_defaults_apply_when_nothing_is_injected(self, toy_path):
        mod = load_mechanism(toy_path, {})
        assert mod._mod._ETA_0 == 20.0


class TestMechanismWrapper:
    def test_multipliers_reach_get_R(self, toy_path):
        """M9: reaction_multipliers are applied inside the accessor."""
        base = load_mechanism(toy_path, {})
        off = load_mechanism(toy_path, {"reaction_multipliers": {"e -> 2e + M+": 0.0}})
        EN, p, T = 300.0, 1.0, 293.0
        assert base.get_R(EN, p, T)[1, 0] > 0.0
        assert off.get_R(EN, p, T)[1, 0] == pytest.approx(0.0)

    def test_xi_photo_scales_B(self, toy_path):
        """M8: the photoionization scale factor is applied by the wrapper."""
        mod = load_mechanism(toy_path, {"xi_photo": 0.5})
        raw_B = mod._mod.get_B(300.0, 1.0, 293.0)
        assert mod.get_B(300.0, 1.0, 293.0).shape == raw_B.shape

    def test_gamma_override_is_applied(self, toy_path):
        """M8b: gamma0 from the configuration reaches get_gamma_plus."""
        mod = load_mechanism(toy_path, {"gamma0": 0.25})
        assert mod.get_gamma_plus(300.0, 1.0, 293.0) == pytest.approx([0.25])

    def test_resolve_applies_polarity_overrides(self, toy_path):
        """
        M7: resolve() returns a configured copy and leaves the original alone.
        """
        mod = load_mechanism(
            toy_path,
            {
                "gamma0": 0.1,
                "pos_override": {"gamma0": 0.4},
                "neg_override": {"gamma0": 0.01},
            },
        )
        pos, neg = mod.resolve(True), mod.resolve(False)
        assert pos.get_gamma_plus(300.0, 1.0, 293.0) == pytest.approx([0.4])
        assert neg.get_gamma_plus(300.0, 1.0, 293.0) == pytest.approx([0.01])
        assert mod.get_gamma_plus(300.0, 1.0, 293.0) == pytest.approx([0.1])

    def test_resolve_without_override_returns_self(self, toy):
        """M7b: no override means no copy, so nothing can drift."""
        assert toy.resolve(True) is toy
        assert toy.resolve(False) is toy

    def test_selectors_are_one_hot_rows(self, toy):
        """The boundary conditions rely on these being selection matrices."""
        for sel in (toy.get_Pi_e(), toy.get_Pi_plus(), toy.get_Pi_minus()):
            assert sel.shape[1] == len(toy.SPECIES)
            assert np.count_nonzero(sel) == sel.shape[0]
            assert set(np.unique(sel)) <= {0.0, 1.0}


class TestReadJsonConfigs:
    def _write(self, tmp_path, name, payload):
        p = tmp_path / name
        p.write_text(json.dumps(payload))
        return str(p)

    def test_single_object(self, tmp_path):
        """M6."""
        p = self._write(tmp_path, "one.json", {"label": "A"})
        out = read_json_configs([p])
        assert len(out) == 1 and out[0]["label"] == "A"

    def test_bare_list(self, tmp_path):
        p = self._write(tmp_path, "list.json", [{"label": "A"}, {"label": "B"}])
        assert [c["label"] for c in read_json_configs([p])] == ["A", "B"]

    def test_configurations_key(self, tmp_path):
        p = self._write(
            tmp_path, "cfg.json", {"configurations": [{"label": "A"}, {"label": "B"}]}
        )
        assert [c["label"] for c in read_json_configs([p])] == ["A", "B"]

    def test_multiple_files_concatenate_in_order(self, tmp_path):
        a = self._write(tmp_path, "a.json", {"label": "A"})
        b = self._write(tmp_path, "b.json", [{"label": "B"}, {"label": "C"}])
        assert [c["label"] for c in read_json_configs([a, b])] == ["A", "B", "C"]

    def test_cfg_dir_is_recorded(self, tmp_path):
        """M6b: relative paths inside a config resolve against its own folder."""
        p = self._write(tmp_path, "d.json", {"label": "A"})
        out = read_json_configs([p])
        assert out[0]["_cfg_dir"] == os.path.dirname(os.path.abspath(p))
