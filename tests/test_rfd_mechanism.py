# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The relativistic feedback discharge mechanism, ``mechanisms/air/rfd/``.

Feedback closes entirely in the gas: nothing is emitted at either end, one
photon group radiates forward and one backward, and positrons carry the rest
back.  The threshold is checked against the normal-mode solution of the
same model (:mod:`closed_form_rfd`, independent of the solver), at sea level
and with the authors' scaling off it, against the numbers the paper states,
against the paper's published results when supplied, and against the
invariants of the model.
"""

import functools
from pathlib import Path

import numpy as np
import pytest

import closed_form_rfd as cf
from incept1d.constants import kB
from incept1d.fields import FieldDistribution
from incept1d.inception import find_all_breakdown_EN
from incept1d.mechanism import load_mechanism, read_json_configs
from incept1d.solver import inception_det, riccati_criterion

ROOT = Path(__file__).resolve().parent.parent
RFD = ROOT / "mechanisms" / "air" / "rfd"
MECHANISM = str(RFD / "air_rfd.py")
UNIFORM = FieldDistribution("uniform")
#: Sea level, n~ = 1: the authors' density 2.688e25 m^-3, the mechanism's default.
T0 = 273.15
P0 = 2.688e25 * kB * T0 / 1e5


def _load(cfg=None):
    return load_mechanism(MECHANISM, cfg)


def _radius_configs():
    return {c["label"]: c for c in read_json_configs([str(RFD / "radius.json")])}


def _threshold(mod, L, p=P0, T=T0, criterion=riccati_criterion):
    """Self-sustained field [V/m] of a region of length L [m]."""
    det = functools.partial(criterion, field_dist=UNIFORM)
    roots = find_all_breakdown_EN(p * L, mod, p, T, first_only=True, det_fn=det)
    assert roots, f"no threshold at L = {L} m"
    return mod._mod.E_field(roots[0], p, T)


def _reference(mod, L, n=1.0):
    """
    Normal-mode threshold [V/m] at relative density n, scaled as the authors'
    code does: avalanche length U_R / (E - n E_d) with the branch at n E_t,
    every other coefficient its sea-level function of E, times n.
    """
    m = mod._mod

    def M_of_E(E):
        E_d = 276e3 if E >= 284e3 * n else 320e3
        return cf.mode_matrix(
            (E - n * E_d) / 7.3e6,
            n * m.P_bb(E),
            n * m.P_cp_fwd(E),
            n * m.P_cp_bck(E),
            n * m.P_PP_FWD,
            n * m.inv_positron_attenuation(E),
            n * m.P_BR_FWD,
            n * m.P_BR_BCK,
            n / m._LAMBDA_A_FWD,
            n / m._LAMBDA_A_BCK,
        )

    return cf.threshold(M_of_E, L, 285e3 * n, 6e6, n_scan=2000)


@pytest.mark.parametrize("label", ["R >> L", "R = 100 m"])
@pytest.mark.parametrize("criterion", [riccati_criterion, inception_det])
@pytest.mark.parametrize("L", [25.0, 100.0, 500.0, 5000.0])
def test_threshold_matches_normal_modes(label, criterion, L):
    """RFD1: both criteria give the normal-mode threshold of the same model."""
    mod = _load(_radius_configs()[label])
    assert _threshold(mod, L, criterion=criterion) == pytest.approx(
        _reference(mod, L), rel=1e-8
    )


def test_conventional_breakdown_is_reached_near_23_m():
    """RFD3: the paper's E_s reaches 3 MV/m at L of about 23 m."""
    mod = _load()
    assert _threshold(mod, 22.0) > 3e6 > _threshold(mod, 24.0)


def test_threshold_falls_toward_the_bhabha_onset():
    """
    RFD4: E_s decreases with L, toward the field where positrons start to
    make RREs (305.4 kV/m) and above the RREA threshold E_t.
    """
    mod = _load()
    E = [_threshold(mod, L) for L in (30.0, 100.0, 300.0, 1000.0, 3000.0)]
    assert all(a > b for a, b in zip(E, E[1:]))
    assert 305.4e3 < E[-1] < 310e3


def test_narrow_region_needs_a_slightly_higher_field():
    """RFD5: photons escape a narrow region sooner; E_s rises, but only a little."""
    cfgs = _radius_configs()
    wide, narrow = _load(cfgs["R >> L"]), _load(cfgs["R = 100 m"])
    for L in (100.0, 300.0):
        ratio = _threshold(narrow, L) / _threshold(wide, L)
        assert 1.0 < ratio < 1.05


@pytest.mark.parametrize("L", [100.0, 1000.0])
def test_off_sea_level_follows_the_authors_scaling(L):
    """
    RFD6: at half the sea-level density the solver gives the normal-mode
    threshold of the model scaled as the authors' code scales it.
    """
    mod = _load()
    assert _threshold(mod, L, p=P0 / 2) == pytest.approx(
        _reference(mod, L, n=0.5), rel=1e-8
    )


def test_nothing_enters_and_photons_are_one_way():
    """RFD7: no emission at x = 0; one forward and one backward photon group."""
    mod = _load()
    assert mod.PHOTON_DIRECTIONS == ["forward", "backward"]
    EN = 20.0
    assert not np.any(mod.get_gamma_plus(EN, P0, T0))
    assert not np.any(mod.get_gamma_Psi(EN, P0, T0))
    V = np.diag(mod.get_V(EN, P0, T0))
    assert V[0] > 0.0 > V[1] and V[0] == -V[1]


PAPER = ROOT / "examples" / "rfd" / "paper.dat"


@pytest.mark.skipif(not PAPER.exists(), reason="examples/rfd/paper.dat not supplied")
def test_threshold_matches_the_papers_normal_modes():
    """
    RFD9: the paper's own normal-mode E_s(L), from the authors' archive
    (examples/rfd/README.md), to within its 100 V/m search grid.  Only the
    paper's column is read; ours is solved here.  Beyond 5 km their minimum of
    |det|^2 is too flat to resolve to one grid step.
    """
    mod = _load(_radius_configs()["R >> L"])
    data = np.loadtxt(PAPER)
    for L, theirs in data[data[:, 0] <= 5e3][:, [0, 3]]:
        assert abs(_threshold(mod, L) / 1e3 - theirs) <= 0.1, f"L = {L} m"


def test_unknown_configuration_keys_are_refused():
    """RFD8: a configuration for another mechanism is not solved as this one."""
    with pytest.raises(ValueError, match="gamma0"):
        _load({"gamma0": 0.01})
