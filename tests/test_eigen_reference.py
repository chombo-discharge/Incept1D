# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The growth rate is the fastest mode: an independent check.

``incept1d growth`` takes lambda* as the largest real root of its criterion,
on the argument that every coupling is a non-negative source, so the
dominant mode is real, positive everywhere, and grows fastest.  Here the
linearised equations are discretised directly and every mode computed
(``eigen_reference.py``, which does not use the solver); the fastest
resolved mode must be real, equal to lambda*, and non-negative, and no
resolved mode, complex ones included, may grow faster.
"""

import functools

import numpy as np
import pytest

from eigen_reference import converged_modes
from incept1d.fields import FieldDistribution
from incept1d.growth import find_lambda_for_voltage
from incept1d.inception import find_all_breakdown_EN
from incept1d.solver import riccati_criterion

UNIFORM = FieldDistribution("uniform")


@pytest.fixture(scope="module")
def star(toy):
    p, T, pd = 1.0, 293.0, 20e-3
    det = functools.partial(riccati_criterion, field_dist=UNIFORM)
    return (
        p,
        T,
        pd,
        find_all_breakdown_EN(pd, toy, p, T, first_only=True, det_fn=det)[0],
    )


# Up to about 1.2 E*; above it the avalanche is too steep for the
# collocation to resolve the dominant mode at these resolutions.
@pytest.mark.parametrize("over", [1.05, 1.2])
def test_growth_rate_is_the_fastest_mode(toy, star, over):
    p, T, pd, EN_star = star
    EN = EN_star * over
    lam, status = find_lambda_for_voltage(
        EN, pd, toy, p, T, UNIFORM, 5, 1000, 1e-3, None
    )
    assert status == "ok"

    modes, shapes = converged_modes(toy, EN, pd, p, T, UNIFORM)
    assert len(modes) > 0
    fastest = modes[0]
    assert abs(fastest.imag) < 1e-6 * abs(fastest)
    assert fastest.real == pytest.approx(lam, rel=1e-5)
    assert np.all(modes.real <= lam * (1.0 + 1e-5)), "a resolved mode grows faster"

    shape = shapes[0] / shapes[0].flat[np.argmax(np.abs(shapes[0]))]
    assert np.abs(shape.imag).max() < 1e-6
    assert shape.real.min() > -1e-6, "the dominant mode must be non-negative"
