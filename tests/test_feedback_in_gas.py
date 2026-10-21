# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Feedback that closes in the gas, with nothing emitted at the cathode.

Two components, f travelling toward the anode and b toward the cathode, with
f(0) = 0 and b(d) = 0 and the flux-form system

    d/dx [f, b] = A [f, b],    A = [[a, -c], [s, l]],

where f grows at a, makes b at s, and b, attenuated at l, makes f back at
c.  There is no cathode emission, so the threshold is reached only when the
gap sustains itself: [exp(A d)]_bb = 0.  With mu = (a + l)/2 and
w^2 = s c - ((a - l)/2)^2,

    exp(A d) = e^{mu d} [cos(w d) I + sin(w d)/w (A - mu I)],

so the threshold is

    cos(w d) + (l - a)/2 * sin(w d)/w = 0,

continued to imaginary w (cosh, sinh) where the modes are real.

The Riccati criterion is identically +1 below this threshold, since Q0 does
not see b at the cathode, and sees inception only as a pole of the reflection
operator.  That pole went undetected when the gap was a single unrefined
step, which put the threshold at 3-4 times the right field.

Two realisations of the same A are checked: a forward and a backward
species (b a positive-ion-like species that regenerates the forward one in
the gas), and electrons with one photon group that radiates backward only
(eq. eq_two_stream_directional).
"""

import functools

import numpy as np
import pytest
import scipy.optimize

from incept1d.fields import FieldDistribution
from incept1d.inception import find_all_breakdown_EN
from incept1d.solver import inception_det, riccati_criterion

UNIFORM = FieldDistribution("uniform")
P, T = 1.0, 293.0
V = 1e5  # m/s, drift speed of both species
S, C, L = 200.0, 200.0, 50.0  # 1/m


def growth(EN):
    """Net growth rate a of the forward component [1/m]."""
    return 1e3 * (EN / 100.0 - 1.0)


def closed_form_threshold(d):
    """Lowest E/N at which cos(w d) + (l - a)/2 sin(w d)/w changes sign."""

    def f(EN):
        a = growth(EN)
        w = np.sqrt(complex(S * C - ((a - L) / 2.0) ** 2))
        return float(
            np.real(np.cos(w * d) + (L - a) / 2.0 * d * np.sinc(w * d / np.pi))
        )

    grid = np.linspace(1.0, 2000.0, 20000)
    values = [f(e) for e in grid]
    for a, b, fa, fb in zip(grid, grid[1:], values, values[1:]):
        if fa * fb < 0.0:
            return scipy.optimize.brentq(f, a, b, xtol=1e-13)
    raise RuntimeError(f"no closed-form threshold for d = {d}")


class _Gap:
    """What both realisations share: no cathode emission, no negative ions."""

    ELECTRON_INDEX = 0

    def resolve(self, positive):
        return self

    def get_Pi_e(self):
        return np.eye(1, len(self.SPECIES))

    def get_gamma_plus(self, EN, p, T):
        return np.zeros(self.get_Pi_plus().shape[0])

    def get_Pi_minus(self):
        return np.zeros((0, len(self.SPECIES)))


class TwoSpecies(_Gap):
    """f is a forward species, b a backward one that regenerates it."""

    SPECIES = ["f", "b"]

    def get_V(self, EN, p, T):
        return np.diag([V, -V])

    def get_R(self, EN, p, T):
        # A = R V^-1, so R = A V; the cross rates are positive.
        return V * np.array([[growth(EN), C], [S, -L]])

    def get_Pi_plus(self):
        return np.array([[0.0, 1.0]])

    def get_B(self, EN, p, T):
        return np.zeros((2, 0))

    def get_C(self, EN, p, T):
        return np.zeros((0, 2))

    def get_kappa(self, p, T):
        return np.zeros(0)

    def get_gamma_Psi(self, EN, p, T):
        return np.zeros(0)


class PhotonGap(_Gap):
    """
    Electrons and photon groups, each ``(direction, b, gamma_Psi)``.

    Every group is emitted at s and absorbed at l, and ionizes at b.  With
    the default single backward group b is its backward stream, and A is
    the same as for :class:`TwoSpecies`.
    """

    SPECIES = ["e"]

    def __init__(self, groups=(("backward", C, 0.0),)):
        self.PHOTON_DIRECTIONS = [g[0] for g in groups]
        self._b = np.array([[g[1] for g in groups]])
        self._gamma_Psi = np.array([g[2] for g in groups])

    def get_V(self, EN, p, T):
        return np.diag([V])

    def get_R(self, EN, p, T):
        return np.array([[V * growth(EN)]])

    def get_Pi_plus(self):
        return np.zeros((0, 1))

    def get_B(self, EN, p, T):
        return self._b

    def get_C(self, EN, p, T):
        # Acts on densities; the solver divides by V.
        return np.full((self._b.shape[1], 1), S * V)

    def get_kappa(self, p, T):
        return np.full(self._b.shape[1], L)

    def get_gamma_Psi(self, EN, p, T):
        return self._gamma_Psi


def _threshold(mod, d, criterion):
    det = functools.partial(criterion, field_dist=UNIFORM)
    roots = find_all_breakdown_EN(P * d, mod, P, T, first_only=True, det_fn=det)
    assert roots, f"no threshold at d = {d} m"
    return roots[0]


@pytest.mark.parametrize("mod", [TwoSpecies(), PhotonGap()], ids=type)
@pytest.mark.parametrize("criterion", [riccati_criterion, inception_det])
@pytest.mark.parametrize("d", [8e-3, 12e-3, 20e-3])
def test_threshold_matches_closed_form(mod, criterion, d):
    """F1: both criteria find [exp(A d)]_bb = 0, with no cathode emission."""
    assert _threshold(mod, d, criterion) == pytest.approx(
        closed_form_threshold(d), rel=1e-8
    )


def test_riccati_is_one_below_and_negative_above():
    """F2: g carries no cathode term here; inception is the pole."""
    d = 12e-3
    EN = closed_form_threshold(d)
    g = functools.partial(riccati_criterion, field_dist=UNIFORM)
    assert g(EN * 0.99, P * d, TwoSpecies(), P, T) == 1.0
    assert g(EN * 1.01, P * d, TwoSpecies(), P, T) < 0.0


def test_forward_stream_lowers_the_threshold():
    """F3: radiating both ways adds the forward stream's ionization."""
    d = 12e-3
    backward = _threshold(PhotonGap(), d, riccati_criterion)
    both = _threshold(PhotonGap((("both", C, 0.0),)), d, riccati_criterion)
    assert both < backward


@pytest.mark.parametrize("criterion", [riccati_criterion, inception_det])
@pytest.mark.parametrize("EN", [60.0, 80.0, 95.0])
def test_both_ways_is_a_forward_and_a_backward_group(criterion, EN):
    """
    F4: a group radiating both ways equals one per direction, same data.

    Below the threshold (near 100 Td), where g is a value and not the pole
    sentinel.
    """
    g = functools.partial(criterion, field_dist=UNIFORM)
    d, b, gamma = 12e-3, 0.4 * C, 0.05
    one = PhotonGap((("both", b, gamma),))
    two = PhotonGap((("forward", b, gamma), ("backward", b, gamma)))
    assert g(EN, P * d, two, P, T) == pytest.approx(g(EN, P * d, one, P, T), rel=1e-10)


@pytest.mark.parametrize("criterion", [riccati_criterion, inception_det])
def test_forward_group_ejects_no_photoelectrons(criterion):
    """F5: no stream of a forward group reaches the cathode."""
    g = functools.partial(criterion, field_dist=UNIFORM)
    d, EN = 12e-3, 80.0  # below the threshold, near 100 Td
    back = ("backward", 0.4 * C, 0.05)
    quiet = PhotonGap((("forward", 0.4 * C, 0.0), back))
    yields = PhotonGap((("forward", 0.4 * C, 0.5), back))
    assert g(EN, P * d, yields, P, T) == g(EN, P * d, quiet, P, T)
