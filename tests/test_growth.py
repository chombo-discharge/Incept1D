# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Temporal growth rate above the inception voltage.

:func:`incept1d.growth.find_lambda_for_voltage` decides whether a voltage is
above inception from the *sign* of its criterion at lambda = 0.  The default,
:func:`incept1d.solver.riccati_criterion`, returns 1 - (loop gain): positive
below inception and negative above for every mechanism.  det Q has no such
invariant sign -- it depends on the parity of the boundary matrix, and the
three-species toy has the opposite sign to dry air -- which is why the toy
could not drive the lambda solver before.  With the Riccati criterion it can,
so the solver's behaviour is checked here on the toy, cheaply, and the air
physics in ``test_air.py``.
"""

import numpy as np
import pytest

from incept1d.fields import FieldDistribution
from incept1d.growth import compute_lambda_curve, find_lambda_for_voltage
from incept1d.inception import find_all_breakdown_EN
from incept1d.solver import inception_det, midpoint_propagator, riccati_criterion

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
    det = functools.partial(riccati_criterion, field_dist=UNIFORM)
    EN_star = find_all_breakdown_EN(pd, mod, p, T, first_only=True, det_fn=det)[0]
    return mod, pd, p, T, EN_star


class TestSignConvention:
    """The assumption find_lambda_for_voltage rests on, stated explicitly."""

    def test_toy_determinant_sign_is_reversed(self, setup):
        """
        det Q has no invariant sign: the toy's is opposite to dry air's.

        This is why the lambda solver defaults to the Riccati criterion.
        """
        mod, pd, p, T, EN_star = setup
        below = inception_det(EN_star * 0.9, pd, mod, p, T, UNIFORM)
        above = inception_det(EN_star * 1.1, pd, mod, p, T, UNIFORM)
        assert below < 0.0 < above

    def test_riccati_sign_is_physical(self, setup):
        """g = 1 - loop gain: positive below inception, negative above."""
        mod, pd, p, T, EN_star = setup
        assert riccati_criterion(EN_star * 0.9, pd, mod, p, T, UNIFORM) > 0.0
        assert riccati_criterion(EN_star * 1.1, pd, mod, p, T, UNIFORM) < 0.0


class TestLambdaOnTheToy:
    """The lambda solver on a mechanism that can now drive it."""

    def test_growth_rate_is_positive_and_increases(self, setup):
        mod, pd, p, T, EN_star = setup
        lams = []
        for over in (1.05, 1.2, 1.5):
            lam, status = find_lambda_for_voltage(
                EN_star * over, pd, mod, p, T, UNIFORM, *DX, midpoint_propagator
            )
            assert status == "ok", f"{over}x: {status}"
            assert lam > 0.0
            lams.append(lam)
        assert lams[0] < lams[1] < lams[2]

    def test_below_inception_is_reported(self, setup):
        mod, pd, p, T, EN_star = setup
        lam, status = find_lambda_for_voltage(
            EN_star * 0.9, pd, mod, p, T, UNIFORM, *DX, midpoint_propagator
        )
        assert status == "below_inception" and lam == 0.0

    def test_nan_inside_the_bracket_is_reported(self, setup):
        """
        A criterion that is finite at the bracket ends but NaN between them
        (det Q can be) gives "det_Q_unresolved", not a Brent failure.
        """
        mod, pd, p, T, EN_star = setup

        def criterion(*args, lam=0.0, **kwargs):
            if 4e5 < lam < 1e6:
                return float("nan")
            return 1.0 if lam > 5e5 else -1.0

        lam, status = find_lambda_for_voltage(
            EN_star * 1.1,
            pd,
            mod,
            p,
            T,
            UNIFORM,
            *DX,
            midpoint_propagator,
            lam_scale=1e6,
            criterion=criterion,
        )
        assert status == "det_Q_unresolved" and np.isnan(lam)

    def test_det_q_gives_the_same_growth_rate(self, setup):
        """
        The two criteria vanish at the same lambda.

        det Q's sign is reversed for the toy, so its root is located
        independently here rather than through find_lambda_for_voltage.
        """
        import scipy.optimize

        mod, pd, p, T, EN_star = setup
        EN = EN_star * 1.2
        lam, status = find_lambda_for_voltage(
            EN, pd, mod, p, T, UNIFORM, *DX, midpoint_propagator
        )
        assert status == "ok"

        def det(x):
            return inception_det(EN, pd, mod, p, T, UNIFORM, lam=x)

        assert det(0.9 * lam) * det(1.1 * lam) < 0.0
        root = scipy.optimize.brentq(det, 0.9 * lam, 1.1 * lam, rtol=1e-12)
        assert root == pytest.approx(lam, rel=1e-6)


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


class TestPolarity:
    """The polarity reaches every det Q evaluation of the lambda solve."""

    @pytest.mark.parametrize("positive", [True, False])
    def test_polarity_is_forwarded(self, setup, monkeypatch, positive):
        import incept1d.growth

        seen = set()

        def spy(*args, **kwargs):
            seen.add(kwargs["positive_polarity"])
            return riccati_criterion(*args, **kwargs)

        monkeypatch.setattr(incept1d.growth, "riccati_criterion", spy)
        mod, pd, p, T, EN_star = setup
        compute_lambda_curve(
            EN_star,
            pd,
            mod,
            p,
            T,
            UNIFORM,
            *DX,
            midpoint_propagator,
            3,
            2.0,
            positive_polarity=positive,
        )
        assert seen == {positive}
