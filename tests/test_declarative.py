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
import json
import warnings
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
from incept1d.solver import _build_A_aug, riccati_criterion

ROOT = Path(__file__).resolve().parent.parent
PASCHEN_JSON = str(ROOT / "mechanisms" / "paschen" / "paschen.json")


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
        (lambda s: s["reactions"].append(dict(s["reactions"][1])), "appears twice"),
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


# --- loading -----------------------------------------------------------------


def test_loader_dispatches_on_the_file():
    """D11: load_mechanism compiles a .json mechanism; configs are not mechanisms."""
    assert is_mechanism_file(PASCHEN_JSON)
    assert not is_mechanism_file(str(ROOT / "mechanisms" / "paschen" / "gases.json"))
    assert isinstance(load_mechanism(PASCHEN_JSON), DeclarativeMechanism)


# --- tables, mixture, photoionization ----------------------------------------

AIR = ROOT / "mechanisms" / "air"
PANCHESHNYI_PY = str(AIR / "pancheshnyi" / "air_pancheshnyi.py")
PANCHESHNYI_JSON = str(AIR / "pancheshnyi" / "air_pancheshnyi.json")


def _air_spec():
    with open(PANCHESHNYI_JSON) as fh:
        return json.load(fh)


def _air(change=None, config=None):
    spec = _air_spec()
    if change is not None:
        change(spec)
    return DeclarativeMechanism(spec, config, source=PANCHESHNYI_JSON)


def test_tables_give_columns_and_summed_rates():
    """D14: a BOLSIG+ column by label, a rate summed over its legend entries."""
    from incept1d.tables import read_table

    table = read_table(str(AIR / "lxcat" / "lisbon.txt"), "bolsig+")
    assert table.fractions == {"N2": 0.79, "O2": 0.21} and table.T_gas == 300.0
    EN, k = table.rate("N2 Ionization")
    assert len(EN) == len(k) and np.all(k >= 0) and k[-1] > 0
    with pytest.raises(ValueError, match="no 'N2 Recombination'"):
        table.rate("N2 Recombination")
    with pytest.raises(ValueError, match="no column 'A99'"):
        table.column("A99")


def test_mixture_must_match_the_table_that_records_one():
    """D15: lisbon.txt is 79/21; an 80/20 mixture with it is refused."""
    with pytest.raises(MechanismFileError, match="differs from the one"):
        _air(config={"x_n2": 0.8, "x_o2": 0.2})


def test_a_table_without_a_mixture_warns_unless_one_is_configured():
    """D16: the mixture cannot be checked; say so, unless the user set it."""
    from incept1d import declarative

    declarative._WARNED.clear()
    with pytest.warns(UserWarning, match="does not record the gas mixture"):
        _air(config={"cross_sections": "../lxcat/phelps.txt"})
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _air(config={"cross_sections": "../lxcat/phelps.txt", "x_n2": 0.8, "x_o2": 0.2})


def test_file_parameters_resolve_against_the_configuration_first(tmp_path):
    """D17: a relative file in a configuration is looked for beside it first."""
    import shutil

    shutil.copy(AIR / "lxcat" / "phelps.txt", tmp_path / "mine.txt")
    m = _air(
        config={
            "cross_sections": "mine.txt",
            "x_n2": 0.8,
            "x_o2": 0.2,
            "_cfg_dir": str(tmp_path),
        }
    )
    ref = _air(
        config={"cross_sections": "../lxcat/phelps.txt", "x_n2": 0.8, "x_o2": 0.2}
    )
    np.testing.assert_array_equal(
        m.get_R(120.0, 1.0, 300.0), ref.get_R(120.0, 1.0, 300.0)
    )
    with pytest.raises(MechanismFileError, match="expected a file name"):
        _air(config={"cross_sections": 3.0})


@pytest.mark.parametrize(
    "change, match",
    [
        (
            lambda s: s["photoionization"]["absorption"].update(model="bourdon"),
            "unknown model",
        ),
        (
            lambda s: s["photoionization"]["absorption"].update(absorber="Ar"),
            "not a gas",
        ),
        (lambda s: s["parameters"].update(ngroups=2.5), "not a count"),
        (
            lambda s: s["photoionization"]["source"].update(
                reactions=["e + N2 -> N2+"]
            ),
            "no reaction",
        ),
        (lambda s: s["photoionization"]["products"].update(Xe=1), "unknown key"),
        (lambda s: s["photoionization"].update(direction="up"), "direction"),
        (
            lambda s: s["tables"]["bolsig"]["rates"].update(k_x="N2 Recombination"),
            "no 'N2",
        ),
        (lambda s: s["tables"]["bolsig"].update(format="csv"), "unknown table format"),
    ],
)
def test_photoionization_and_tables_are_checked(change, match):
    """D18: the photoionization block and the tables follow the same rules."""
    with pytest.raises(MechanismFileError, match=match):
        _air(change)


def test_photon_direction_is_declared_per_mechanism():
    """D19: a direction gives PHOTON_DIRECTIONS, one entry per group."""
    m = _air(lambda s: s["photoionization"].update(direction="backward"))
    assert m.PHOTON_DIRECTIONS == ["backward"] * 3
    assert _air().PHOTON_DIRECTIONS is None


def _pancheshnyi_configs():
    from incept1d.mechanism import read_json_configs

    out = []
    for path in sorted((AIR / "pancheshnyi").glob("*.json")):
        if is_mechanism_file(str(path)):
            continue
        for cfg in read_json_configs([str(path)]):
            out.append(pytest.param(cfg, id=f"{path.stem}:{cfg.get('label', '')}"))
    return out


@pytest.mark.slow
@pytest.mark.parametrize("cfg", _pancheshnyi_configs())
def test_pancheshnyi_json_builds_the_same_augmented_matrix(cfg):
    """D20: every shipped configuration, both polarities, to round-off."""
    py = load_mechanism(PANCHESHNYI_PY, dict(cfg))
    js = load_mechanism(PANCHESHNYI_JSON, dict(cfg))
    for positive in (True, False):
        a, b = py.resolve(positive), js.resolve(positive)
        for EN in (40.0, 110.0, 300.0):
            for p, T in ((1.0, 300.0), (3.0, 293.0)):
                for lam in (0.0, 1e7):
                    A = _build_A_aug(EN, 1e-3, a, p, T, lam)[0]
                    B = _build_A_aug(EN, 1e-3, b, p, T, lam)[0]
                    np.testing.assert_allclose(B, A, rtol=1e-13, atol=0.0)
                for f in ("get_gamma_plus", "get_gamma_Psi"):
                    np.testing.assert_allclose(
                        getattr(b, f)(EN, p, T), getattr(a, f)(EN, p, T), rtol=1e-14
                    )
                for f in ("alpha", "eta"):
                    assert getattr(b, f)(EN, p, T) == pytest.approx(
                        getattr(a, f)(EN, p, T), rel=1e-13
                    )


@pytest.mark.slow
@pytest.mark.parametrize(
    "pd_mm, EN_expected", [(1.0, 175.4867), (10.0, 116.2765), (100.0, 100.711)]
)
def test_pancheshnyi_json_meets_the_air_reference_values(pd_mm, EN_expected):
    """D21: the anchors of test_air A3, from the JSON mechanism."""
    det = functools.partial(riccati_criterion, field_dist=FieldDistribution("uniform"))
    roots = find_all_breakdown_EN(
        pd_mm * 1e-3,
        load_mechanism(PANCHESHNYI_JSON),
        1.0,
        293.0,
        first_only=True,
        det_fn=det,
    )
    assert roots[0] == pytest.approx(EN_expected, rel=1e-5)
