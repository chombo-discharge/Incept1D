# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""The shipped mechanism directories: one mechanism each, with its own configs.

A configuration file belongs to the mechanism in its directory.  These tests
check that every shipped configuration loads with that mechanism, and that
the loader refuses one written for another mechanism instead of silently
solving a different problem.  A mechanism is a JSON file (a mechanism file)
or, for a model the format cannot express, a Python module.
"""

import os
import re
from pathlib import Path

import pytest

from incept1d.declarative import is_mechanism_file
from incept1d.mechanism import load_mechanism, read_json_configs

ROOT = Path(__file__).resolve().parent.parent
MECHANISMS = ROOT / "mechanisms"


def _mechanism_files(mech_dir):
    """The mechanism of *mech_dir*, as a module and / or a mechanism file."""
    py = [p for p in mech_dir.glob("*.py") if p.name != "config.py"]
    js = [p for p in mech_dir.glob("*.json") if is_mechanism_file(str(p))]
    return sorted(py + js)


#: Every directory holding a mechanism module or file.
MECH_DIRS = sorted(
    {p.parent for p in MECHANISMS.rglob("config.py")}
    | {p.parent for p in MECHANISMS.rglob("*.json") if is_mechanism_file(str(p))}
)

#: Mechanisms that are slow to load.
_SLOW = set()


def _param(path, mech_dir):
    marks = [pytest.mark.slow] if mech_dir.name in _SLOW else []
    return pytest.param(path, mech_dir, marks=marks, id=str(path.relative_to(ROOT)))


def test_there_are_mechanism_directories():
    """Guard against the collection silently matching nothing."""
    names = {d.name for d in MECH_DIRS}
    assert {"pancheshnyi", "2body", "paschen"} <= names


@pytest.mark.parametrize("mech_dir", MECH_DIRS, ids=lambda d: d.name)
def test_one_mechanism_per_directory(mech_dir):
    """L1: a directory is one mechanism (in one or both forms) and its configs."""
    files = _mechanism_files(mech_dir)
    assert (
        files and len({f.stem for f in files}) == 1
    ), f"{mech_dir} holds {[f.name for f in files]}"


@pytest.mark.parametrize(
    "json_path, mech_dir",
    [
        _param(j, d)
        for d in MECH_DIRS
        for j in sorted(d.glob("*.json"))
        if not is_mechanism_file(str(j))
    ],
)
def test_every_config_loads_with_its_own_mechanism(json_path, mech_dir):
    """L2: each shipped configuration is accepted by the mechanism beside it."""
    for mech in _mechanism_files(mech_dir):
        for cfg in read_json_configs([str(json_path)]):
            mod = load_mechanism(str(mech), cfg)
            assert mod.get_R(120.0, 1.0, 293.0).shape[0] == len(mod.SPECIES)


AIR = MECHANISMS / "air"


@pytest.mark.parametrize(
    "mechanism, config, why",
    [
        pytest.param(
            AIR / "2body" / "air_2body.json",
            AIR / "pancheshnyi" / "nodetachment.json",
            "does not define",
            id="pancheshnyi-reactions-for-2body",
        ),
        pytest.param(
            AIR / "pancheshnyi" / "air_pancheshnyi.json",
            AIR / "2body" / "example_config.json",
            "does not define",
            id="2body-reactions-for-pancheshnyi",
        ),
        pytest.param(
            MECHANISMS / "paschen" / "paschen.json",
            AIR / "pancheshnyi" / "see.json",
            "unknown key",
            id="air-config-for-paschen",
        ),
    ],
)
def test_a_config_for_another_mechanism_is_refused(mechanism, config, why):
    """L3: the pairing is enforced, not merely conventional."""
    with pytest.raises(ValueError, match=why):
        for cfg in read_json_configs([str(config)]):
            load_mechanism(str(mechanism), cfg)


@pytest.mark.parametrize(
    "py", sorted(MECHANISMS.rglob("*.py")), ids=lambda p: str(p.relative_to(ROOT))
)
def test_mechanism_files_do_not_touch_sys_path(py):
    """L4: shared code is loaded by path (load_helper), never via sys.path."""
    assert not re.search(r"sys\.path", py.read_text()), (
        f"{os.path.relpath(py, ROOT)} modifies sys.path; use "
        f"incept1d.mechanism.load_helper"
    )
