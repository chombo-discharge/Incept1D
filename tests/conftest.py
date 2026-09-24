# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Shared fixtures for the Incept1D test suite."""

import functools
import os
import sys

# Every CLI module imports matplotlib at module scope; make sure importing them
# never tries to open a window.  Must happen before the first import.
os.environ.setdefault("MPLBACKEND", "Agg")

import pytest  # noqa: E402

from incept1d.fields import FieldDistribution  # noqa: E402
from incept1d.mechanism import load_mechanism  # noqa: E402
from incept1d.solver import inception_det  # noqa: E402

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
TOY_PATH = os.path.join(TESTS_DIR, "toy", "toy_mechanism.py")
AIR_PATH = os.path.join(
    REPO_ROOT, "mechanisms", "air", "pancheshnyi", "air_pancheshnyi.py"
)

# tests/closed_form.py is imported by name from the test modules.
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)


@pytest.fixture(scope="session")
def toy_path():
    """Path to the three-species toy mechanism file."""
    return TOY_PATH


@pytest.fixture(scope="session")
def air_path():
    """Path to the reference dry-air mechanism."""
    return AIR_PATH


@pytest.fixture(scope="session")
def toy(toy_path):
    """The toy mechanism loaded through the real loader, default configuration."""
    return load_mechanism(toy_path, {})


@pytest.fixture(scope="session")
def toy_raw(toy):
    """The bare toy module, for direct access to alpha/eta/delta and _GAMMA0."""
    return toy._mod


@pytest.fixture(scope="session")
def air():
    """The dry-air mechanism, baseline configuration (slow to load)."""
    return load_mechanism(AIR_PATH, {})


@pytest.fixture(scope="session")
def uniform():
    """A uniform FieldDistribution."""
    return FieldDistribution("uniform")


@pytest.fixture(scope="session")
def det_uniform(uniform):
    """``inception_det`` bound to a uniform field, ready for the root finders."""
    return functools.partial(inception_det, field_dist=uniform)


@pytest.fixture
def counting_det(uniform):
    """
    A det function that records how many times it was evaluated.

    Returns ``(fn, counter)`` where ``counter["n"]`` is the evaluation count.
    Used to assert that the warm start actually saves work.
    """
    counter = {"n": 0}

    def fn(EN, pd, mod, p, T):
        counter["n"] += 1
        return inception_det(EN, pd, mod, p, T, field_dist=uniform)

    return fn, counter
