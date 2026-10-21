# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Declarative mechanism files (``incept1d.declarative``, ``incept1d.expressions``).

A JSON mechanism must compile into exactly what the equivalent Python
mechanism builds by hand, and anything it cannot describe unambiguously must
be refused when it is loaded, not discovered as a wrong curve.
"""

import copy
import functools
from pathlib import Path

import numpy as np
import pytest

from incept1d.declarative import (
    DeclarativeMechanism,
    MechanismFileError,
    is_mechanism_file,
)
from incept1d.expressions import Expression, ExpressionError
from incept1d.fields import FieldDistribution
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import load_mechanism
from incept1d.solver import _build_A_aug, inception_det, riccati_criterion

ROOT = Path(__file__).resolve().parent.parent
PASCHEN_PY = str(ROOT / "mechanisms" / "paschen" / "paschen.py")
PASCHEN_JSON = str(ROOT / "mechanisms" / "paschen" / "paschen.json")
GASES = {"helium": (3.0, 34.0), "argon": (12.0, 180.0), "air": (15.0, 365.0)}


# --- expressions -------------------------------------------------------------


def test_expression_evaluates_a_formula():
    """D1: a rate formula evaluates, and reports the names it uses."""
    f = Expression("1.24e-17*exp(-(259/(14.6 + EN))**2)", ["EN", "p"])
    assert f({"EN": 100.0}) == pytest.approx(1.24e-17 * np.exp(-((259 / 114.6) ** 2)))
    assert f.names == {"EN"}


@pytest.mark.parametrize(
    "text",
    [
        '__import__("os")',
        "x.real",
        "x[0]",
        "exp + 1",
        "x if x else 1",
        "lambda: 1",
        "x < 1",
        '"text"',
        "open(1)",
        "exp(x=1)",
        "exp(*x)",
        "[x]",
        "x and x",
        "True",
        "y",
    ],
)
def test_expression_refuses_anything_but_arithmetic(text):
    """D2: only numbers, known names, arithmetic and the listed functions."""
    with pytest.raises(ExpressionError):
        Expression(text, ["x"])


def test_expression_overflows_instead_of_raising():
    """D3: NumPy arithmetic: 1/0 and huge powers give inf; exp(-inf) = 0."""
    assert Expression("9**9**9", [])({}) == np.inf
    assert Expression("1/x", ["x"])({"x": 0.0}) == np.inf
    assert Expression("exp(-1/x)", ["x"])({"x": 0.0}) == 0.0


# --- a small mechanism to check the rules on -------------------------------

BASE = {
    "format": "incept1d-mechanism/1",
    "parameters": {"k_att": 1e-16, "gamma0": 1e-3},
    "gas": {"N2": "0.8", "O2": "0.2"},
    "untracked": ["O"],
    "species": [
        {"name": "e", "charge": -1, "mobility": "4e23"},
        {"name": "O2+", "charge": 1, "mobility": "2e21"},
        {"name": "O-", "charge": -1, "mobility": "2e21"},
        {"name": "O3-", "charge": -1, "mobility": "2e21"},
    ],
    "reactions": [
        {"reaction": "e + O2 -> 2e + O2+", "rate": "1e-16*exp(-500/EN)"},
        {"reaction": "e + O2 -> O- + O", "rate": "k_att"},
        {"reaction": "O- + O2 + M -> O3- + M", "rate": "1e-42"},
    ],
    "emission": {"ions": {"O2+": "gamma0"}},
}


def _mech(spec=None, config=None):
    return DeclarativeMechanism(copy.deepcopy(spec or BASE), config)


def _with(change):
    spec = copy.deepcopy(BASE)
    change(spec)
    return spec


def test_reaction_matrix_follows_the_reaction_strings():
    """D4: K = k x the density of every neutral reactant; R[:, driver] += delta K."""
    m = _mech()
    EN, p, T = 120.0, 1.0, 300.0
    N = p * 1e5 / (1.380649e-23 * T)
    R = m.get_R(EN, p, T)
    k_ion = 1e-16 * np.exp(-500 / EN) * 0.2 * N
    k_att = 1e-16 * 0.2 * N
    k_3b = 1e-42 * 0.2 * N * N
    assert R[:, 0] == pytest.approx([k_ion - k_att, k_ion, k_att, 0.0])
    assert R[:, 2] == pytest.approx([0.0, 0.0, -k_3b, k_3b])


def test_selectors_and_velocities_follow_the_charges():
    """D5: Pi_e, Pi_+, Pi_- and the drift directions come from the charges."""
    m = _mech()
    assert m.ELECTRON_INDEX == 0
    assert m.get_Pi_e().tolist() == [[1, 0, 0, 0]]
    assert m.get_Pi_plus().tolist() == [[0, 1, 0, 0]]
    assert m.get_Pi_minus().tolist() == [[0, 0, 1, 0], [0, 0, 0, 1]]
    v = np.diag(m.get_V(100.0, 1.0, 300.0))
    assert v[0] > 0 and v[1] < 0 and v[2] > 0 and v[3] > 0


@pytest.mark.parametrize(
    "change, match",
    [
        (lambda s: s.update(colour="red"), "unknown key"),
        (lambda s: s.update(format="something/2"), "format"),
        (lambda s: s["species"][0].update(diffusion="1"), "unknown key"),
        (lambda s: s["species"].pop(0), "electron"),
        (lambda s: s["species"][1].update(charge=0), "non-zero integer"),
        (lambda s: s["species"].append(dict(s["species"][1])), "repeated"),
        (
            lambda s: s["reactions"].append(
                {"reaction": "O- + O2 -> O2+ + O", "rate": "1"}
            ),
            "charge",
        ),
        (
            lambda s: s["reactions"].append({"reaction": "e + O- -> O3-", "rate": "1"}),
            "one tracked",
        ),
        (
            lambda s: s["reactions"].append({"reaction": "e + Ar -> O-", "rate": "1"}),
            "neither",
        ),
        (
            lambda s: s["reactions"].append(
                {"reaction": "O- + O2 -> O3- + Q", "rate": "1"}
            ),
            "untracked",
        ),
        (lambda s: s["reactions"][0].update(rate="k_ion"), "unknown name"),
        (lambda s: s["emission"]["ions"].clear(), "no yield"),
        (lambda s: s.update(gas={"N2": "0.8", "O2": "0.3"}), "sum to 1"),
        (lambda s: s["parameters"].update({"2x": 1.0}), "start with a letter"),
        (lambda s: s["parameters"].update(EN=1.0), "start with a letter"),
        (lambda s: s["parameters"].update(k=True), "number"),
    ],
)
def test_a_file_that_is_not_a_mechanism_is_refused(change, match):
    """D6: every rule of the format is checked when the file is loaded."""
    with pytest.raises((MechanismFileError, ExpressionError), match=match):
        _mech(_with(change))


def test_comments_are_keys_starting_with_underscore():
    """D7: '_' keys are ignored at every level."""
    spec = _with(lambda s: (s.update(_note="x"), s["species"][0].update(_unit="m")))
    _mech(spec)


# --- configurations ----------------------------------------------------------


def test_configuration_sets_parameters_by_name():
    """D8: a configuration overrides the file's parameters, and nothing else."""
    m = _mech(config={"label": "x", "k_att": 2e-16, "_comment": "ok"})
    base = _mech()
    assert m.label == "x"
    assert m.get_R(120.0, 1.0, 300.0)[2, 0] == pytest.approx(
        2 * base.get_R(120.0, 1.0, 300.0)[2, 0]
    )
    with pytest.raises(MechanismFileError, match="unknown key"):
        _mech(config={"gamma": 0.1})


def test_multipliers_name_reactions_and_unknown_ones_are_refused():
    """D9: multipliers are keyed by reaction string; a stray key is an error."""
    m = _mech(config={"reaction_multipliers": {"e + O2 -> O- + O": 0.0}})
    assert m.get_R(120.0, 1.0, 300.0)[2, 0] == 0.0
    with pytest.raises(MechanismFileError, match="does not define"):
        _mech(config={"reaction_multipliers": {"e + N2 -> 2e + N2+": 0.0}})


def test_polarity_overrides_are_parameters_of_one_polarity():
    """D10: 'positive' / 'negative' override parameters for that polarity only."""
    m = _mech(config={"negative": {"gamma0": 0.1}})
    assert m.has_polarity_overrides
    assert m.resolve(True) is m
    assert m.resolve(False).get_gamma_plus(100.0, 1.0, 300.0)[0] == 0.1
    assert m.get_gamma_plus(100.0, 1.0, 300.0)[0] == 1e-3
    with pytest.raises(MechanismFileError, match="unknown key"):
        _mech(config={"positive": {"label": "x"}})


# --- equivalence with the Python Paschen mechanism ---------------------------


def test_loader_dispatches_on_the_file():
    """D11: load_mechanism compiles a .json mechanism; configs are not mechanisms."""
    assert is_mechanism_file(PASCHEN_JSON)
    assert not is_mechanism_file(str(ROOT / "mechanisms" / "paschen" / "gases.json"))
    assert isinstance(load_mechanism(PASCHEN_JSON), DeclarativeMechanism)


@pytest.mark.parametrize("gas", sorted(GASES))
@pytest.mark.parametrize("lam", [0.0, 1e6])
def test_paschen_json_builds_the_same_augmented_matrix(gas, lam):
    """D12: A_aug of the JSON and the Python Paschen mechanism agree to round-off."""
    A, B = GASES[gas]
    py = load_mechanism(PASCHEN_PY, {"gas": gas, "gamma0": 3e-3})
    js = load_mechanism(PASCHEN_JSON, {"A": A, "B": B, "gamma0": 3e-3})
    for EN in (30.0, 150.0, 800.0):
        for p, T in ((1.0, 293.0), (0.2, 350.0)):
            a = _build_A_aug(EN, 1e-3, py, p, T, lam)[0]
            b = _build_A_aug(EN, 1e-3, js, p, T, lam)[0]
            np.testing.assert_allclose(b, a, rtol=1e-14, atol=0.0)
            assert js.get_gamma_plus(EN, p, T) == pytest.approx(
                py.get_gamma_plus(EN, p, T)
            )


@pytest.mark.parametrize("criterion", [riccati_criterion, inception_det])
def test_paschen_json_gives_the_same_inception_field(criterion):
    """D13: the same root, with multipliers and a polarity override."""
    cfg_py = {
        "gas": "argon",
        "reaction_multipliers": {"e -> 2e + M+": 0.9},
        "neg_override": {"gamma0": 0.05},
    }
    cfg_js = {
        "A": 12.0,
        "B": 180.0,
        "reaction_multipliers": {"e -> 2e + M+": 0.9},
        "negative": {"gamma0": 0.05},
    }
    py, js = load_mechanism(PASCHEN_PY, cfg_py), load_mechanism(PASCHEN_JSON, cfg_js)
    for positive in (True, False):
        det = functools.partial(
            criterion,
            field_dist=FieldDistribution("uniform"),
            positive_polarity=positive,
        )
        a = find_all_breakdown_EN(0.01, py, 1.0, 293.0, first_only=True, det_fn=det)
        b = find_all_breakdown_EN(0.01, js, 1.0, 293.0, first_only=True, det_fn=det)
        assert b[0] == pytest.approx(a[0], rel=1e-12)
