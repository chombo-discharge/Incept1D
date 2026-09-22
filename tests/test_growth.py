# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Temporal growth rate above the inception voltage.

NOTE on scope.  :func:`incept1d.growth.find_lambda_for_voltage` decides whether
a voltage is above inception from the *sign* of det Q at lambda = 0
("det Q(0) > 0 means V <= V*").  That sign is not an invariant of the physics:
it depends on the row ordering and hence the parity of the boundary matrix,
which differs between mechanisms.  The dry-air scheme (6 species, 11x11 Q) has
det Q > 0 below inception; the three-species toy has the opposite sign, so the
lambda solver reports ``below_inception`` for it at every voltage.

The structural behaviour of ``compute_lambda_curve`` is therefore checked here
with the toy, and the physics -- lambda > 0 above threshold, growing with
overvoltage -- is checked against the air mechanism in ``test_air.py``.
"""

import numpy as np
import pytest

from incept1d.fields import FieldDistribution
from incept1d.growth import compute_lambda_curve, find_lambda_for_voltage
from incept1d.inception import find_all_breakdown_EN
from incept1d.solver import inception_det, midpoint_propagator

UNIFORM = FieldDistribution("uniform")
DX = (5, 200, 0.03)


@pytest.fixture(scope="module")
def setup():
    """(mechanism, pd, p, T, EN*) for a gap that has a clean inception point."""
    import functools
    import os

    from incept1d.mechanism import load_mechanism

    path = os.path.join(os.path.dirname(__file__), "toy", "toy_mechanism.py")
    mod = load_mechanism(path, {})
    p, T, pd = 1.0, 293.0, 20e-3
    det = functools.partial(inception_det, field_dist=UNIFORM)
    EN_star = find_all_breakdown_EN(pd, mod, p, T, first_only=True, det_fn=det)[0]
    return mod, pd, p, T, EN_star


class TestSignConvention:
    """The assumption find_lambda_for_voltage rests on, stated explicitly."""

    def test_toy_determinant_sign_is_reversed(self, setup):
        """
        Documents why the toy cannot drive the lambda solver.

        If this test ever fails, the determinant convention has changed and the
        scope note at the top of this module needs revisiting.
        """
        mod, pd, p, T, EN_star = setup
        below = inception_det(EN_star * 0.9, pd, mod, p, T, UNIFORM)
        above = inception_det(EN_star * 1.1, pd, mod, p, T, UNIFORM)
        assert below < 0.0 < above

    def test_below_inception_is_reported_for_the_toy(self, setup):
        """Consequence of the above: the solver short-circuits."""
        mod, pd, p, T, EN_star = setup
        lam, status = find_lambda_for_voltage(
            EN_star * 1.5, pd, mod, p, T, UNIFORM, *DX, midpoint_propagator
        )
        assert status == "below_inception"
        assert lam == 0.0


class TestLambdaCurve:
    def test_curve_starts_at_the_inception_point(self, setup):
        """The first row is V = V*, lambda = 0 by construction."""
        mod, pd, p, T, EN_star = setup
        rows = compute_lambda_curve(
            EN_star, pd, mod, p, T, UNIFORM, *DX, midpoint_propagator, 4, 2.0
        )
        assert len(rows) == 4
        V_kV, V_ratio, EN_ref, lam, tau, nu, status = rows[0]
        assert V_ratio == pytest.approx(1.0)
        assert lam == 0.0
        assert EN_ref == pytest.approx(EN_star, rel=1e-9)

    def test_voltage_ratios_span_the_requested_range(self, setup):
        mod, pd, p, T, EN_star = setup
        rows = compute_lambda_curve(
            EN_star, pd, mod, p, T, UNIFORM, *DX, midpoint_propagator, 5, 3.0
        )
        ratios = [r[1] for r in rows]
        assert ratios[0] == pytest.approx(1.0)
        assert ratios[-1] == pytest.approx(3.0)
        assert all(a < b for a, b in zip(ratios, ratios[1:]))

    def test_reduced_field_tracks_the_voltage(self, setup):
        """E/N is proportional to V for a fixed gap."""
        mod, pd, p, T, EN_star = setup
        rows = compute_lambda_curve(
            EN_star, pd, mod, p, T, UNIFORM, *DX, midpoint_propagator, 4, 2.0
        )
        for V_kV, V_ratio, EN_ref, *_ in rows:
            assert EN_ref == pytest.approx(EN_star * V_ratio, rel=1e-9)
