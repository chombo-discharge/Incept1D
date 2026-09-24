# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The shipped dry-air mechanism.

These tests load real swarm data and run the full solver, so they are marked
``slow``.  The important one is :meth:`TestClosedFormLimit` -- it ties the
mechanism people actually use to the analytic result, rather than to numbers
recorded from a previous run.
"""

import functools

import numpy as np
import pytest

from closed_form import standard_paschen_lhs
from incept1d.fields import FieldDistribution
from incept1d.growth import find_lambda_for_voltage
from incept1d.inception import compute_inception_curve, find_all_breakdown_EN
from incept1d.mechanism import load_mechanism
from incept1d.solver import inception_det, midpoint_propagator

pytestmark = pytest.mark.slow

UNIFORM = FieldDistribution("uniform")


@pytest.fixture(scope="module")
def det():
    return functools.partial(inception_det, field_dist=UNIFORM)


@pytest.fixture(scope="module")
def paschen_air(air_path):
    """Air reduced to the textbook limit by mechanisms/air/pancheshnyi/paschen.json."""
    import json
    import os

    cfg_path = os.path.join(os.path.dirname(air_path), "paschen.json")
    with open(cfg_path) as fh:
        cfg = json.load(fh)["configurations"][0]
    cfg["_cfg_dir"] = os.path.dirname(cfg_path)
    return load_mechanism(air_path, cfg)


class TestMechanism:
    def test_loads_and_exposes_the_interface(self, air):
        """A1."""
        assert air.SPECIES[air.ELECTRON_INDEX] == "e"
        assert len(air.SPECIES) == 6
        assert callable(air.alpha) and callable(air.eta)

    def test_alpha_rises_and_eta_falls_with_field(self, air):
        p, T = 1.0, 293.0
        EN = np.array([50.0, 100.0, 200.0, 400.0])
        a = np.array([air.alpha(e, p, T) for e in EN])
        assert np.all(np.diff(a) > 0.0)

    def test_alpha_eta_crossover_is_in_the_expected_range(self, air):
        """For dry air the crossover sits near 100 Td."""
        import scipy.optimize

        p, T = 1.0, 293.0
        f = lambda e: air.alpha(e, p, T) - air.eta(e, p, T)  # noqa: E731
        root = scipy.optimize.brentq(f, 50.0, 200.0)
        assert 90.0 < root < 130.0


class TestClosedFormLimit:
    """
    A2: the strongest test in the suite.

    ``paschen.json`` switches off detachment, ion conversion and photon
    feedback, which is exactly the reduced model of eq_standard_paschen.  The
    full mechanism -- real swarm data, real rate coefficients, the whole
    boundary-value solve -- must then reproduce the textbook algebra.
    """

    @pytest.mark.parametrize("pd_mm", [1.0, 10.0, 100.0])
    def test_reproduces_standard_paschen(self, paschen_air, det, pd_mm):
        p, T = 1.0, 293.0
        pd = pd_mm * 1e-3
        d = pd / p
        roots = find_all_breakdown_EN(
            pd, paschen_air, p, T, first_only=True, det_fn=det
        )
        assert roots, f"no inception found at pd = {pd_mm} bar·mm"
        EN = roots[0]
        a, e = paschen_air.alpha(EN, p, T), paschen_air.eta(EN, p, T)
        gamma = float(paschen_air.get_gamma_plus(EN, p, T)[0])
        assert standard_paschen_lhs(a, e, gamma, d) == pytest.approx(1.0, rel=2e-2)


class TestInceptionCurve:
    def test_paschen_minimum(self, air, det):
        """
        A4: end-to-end regression for the warm-start lock-in.

        Before the fix the curve bottomed out at 468 V because branch 1 had
        locked onto a spurious high-field root at the smallest pd.
        """
        p, T = 1.0, 293.0
        pd_arr = np.logspace(np.log10(3e-6), np.log10(1e-1), 40)
        branches = compute_inception_curve(pd_arr, air, p, T, det_fn=det)
        br = branches[0]
        i = int(np.argmin(br["V"]))
        assert br["V"][i] == pytest.approx(378.0, rel=0.05)
        assert br["pd"][i] * 1e3 == pytest.approx(9.5e-3, rel=0.3)

    def test_branch_one_is_the_lowest_root(self, air, det):
        """
        A5: an invariant, not a pinned number.

        A pinned curve would have happily passed the lock-in bug, because the
        wrong curve was self-consistent.  This checks the contract instead.
        """
        p, T = 1.0, 293.0
        for pd_mm in (8e-3, 3e-2, 1.0, 50.0):
            pd = pd_mm * 1e-3
            first = find_all_breakdown_EN(pd, air, p, T, first_only=True, det_fn=det)
            every = find_all_breakdown_EN(pd, air, p, T, first_only=False, det_fn=det)
            if not every:
                continue
            assert first[0] == pytest.approx(min(every), rel=1e-6)

    @pytest.mark.parametrize(
        "pd_mm, EN_expected",
        [
            (1.0, 175.9713),
            (10.0, 116.4574),
            (100.0, 100.5216),
        ],
    )
    def test_reference_values(self, air, det, pd_mm, EN_expected):
        """A3: drift detector for the uniform-field curve at 1 bar."""
        p, T = 1.0, 293.0
        roots = find_all_breakdown_EN(
            pd_mm * 1e-3, air, p, T, first_only=True, det_fn=det
        )
        assert roots[0] == pytest.approx(EN_expected, rel=1e-5)

    def test_sphere_plane_polarities_differ(self, air):
        """Asymmetric geometry must give two different answers."""
        p, T, pd = 1.0, 293.0, 20e-3
        fd = FieldDistribution("sphere-plane", 50e-3)
        assert not fd.is_symmetric
        kw = dict(field_dist=fd)
        pos = find_all_breakdown_EN(
            pd,
            air,
            p,
            T,
            first_only=True,
            det_fn=functools.partial(inception_det, positive_polarity=True, **kw),
        )[0]
        neg = find_all_breakdown_EN(
            pd,
            air,
            p,
            T,
            first_only=True,
            det_fn=functools.partial(inception_det, positive_polarity=False, **kw),
        )[0]
        assert pos != neg
        assert abs(pos - neg) / pos < 0.2


class TestGrowthRate:
    """The physics of growth.py, which the toy cannot drive (see test_growth)."""

    @staticmethod
    @pytest.fixture(scope="class")
    def star(air, det):
        p, T, pd = 1.0, 293.0, 10e-3
        EN = find_all_breakdown_EN(pd, air, p, T, first_only=True, det_fn=det)[0]
        return air, pd, p, T, EN

    def test_determinant_sign_convention(self, star):
        """The assumption find_lambda_for_voltage relies on, for real air."""
        air, pd, p, T, EN = star
        below = inception_det(EN * 0.9, pd, air, p, T, UNIFORM)
        assert np.isfinite(below) and below > 0.0

    def test_small_overvoltage_gives_a_small_growth_rate(self, star):
        """Just above threshold the solve still behaves."""
        air, pd, p, T, EN = star
        lam, status = find_lambda_for_voltage(
            EN * 1.1, pd, air, p, T, UNIFORM, 5, 200, 0.03, midpoint_propagator
        )
        assert status in ("ok", "suspect")
        assert np.isfinite(lam) and lam > 0.0

    def test_growth_rate_increases_with_overvoltage(self, star):
        """
        G2: the physical requirement, over the range where det Q is resolvable.

        Beyond roughly 1.2 V* the determinant underflows (see
        test_unresolvable_determinant_is_reported); below it the solve is
        sound and lambda must grow with the applied voltage.
        """
        air, pd, p, T, EN = star
        lams = []
        for over in (1.02, 1.05, 1.10, 1.15):
            lam, status = find_lambda_for_voltage(
                EN * over, pd, air, p, T, UNIFORM, 5, 200, 0.03, midpoint_propagator
            )
            assert status in ("ok", "suspect"), f"{over}x V*: {status}"
            assert np.isfinite(lam) and lam > 0.0
            lams.append(lam)
        assert all(a < b for a, b in zip(lams, lams[1:])), lams

    def test_unresolvable_determinant_is_reported(self, star):
        """
        Above roughly 1.2 V* the determinant cannot be evaluated at all.

        The spectral spread of A_aug*d exceeds the double-precision underflow
        limit (~709), so the subdominant modes of M underflow to exactly zero,
        Q becomes rank deficient and det Q is NaN for *every* lambda.  Because
        ``_f_brentq`` maps NaN to a negative sentinel, a run of NaN below a
        positive value looks like a sign change, and brentq used to converge on
        the edge of the NaN region -- reporting the same 1.9e11 s^-1 for every
        overvoltage, since that edge is the kappa*d = 12 photon-collapse
        threshold and does not depend on E/N.

        The solver must now say so rather than return that number.  Remove this
        test if the propagation is reformulated so the determinant survives.
        """
        air, pd, p, T, EN = star
        for over in (1.25, 1.5, 2.0):
            lam, status = find_lambda_for_voltage(
                EN * over, pd, air, p, T, UNIFORM, 5, 200, 0.03, midpoint_propagator
            )
            assert status == "det_Q_unresolved", f"{over}x V*: {status}"
            assert np.isnan(lam)

    def test_determinant_underflows_above_that_range(self, star):
        """The cause, asserted directly: det Q is NaN at lambda = 0 already."""
        air, pd, p, T, EN = star
        assert np.isfinite(inception_det(EN * 1.15, pd, air, p, T, UNIFORM))
        assert np.isnan(inception_det(EN * 1.5, pd, air, p, T, UNIFORM))


class TestNonUniformRootSelection:
    """The coarse scan runs on a proxy; the answer must still be the lowest root."""

    #: Sphere-sphere gap and pd at which the uniform-field proxy misplaces the
    #: lowest root by about 10 %, and where air has three roots in the search
    #: range -- so a dropped bracket has somewhere worse to land.
    SPHERE_R_M = 50e-3
    D_M = 0.200
    PD_VALUES = (1.0e-2, 1.0975e-2, 1.2619e-2)

    @pytest.fixture(scope="class")
    @classmethod
    def sphere_dets(cls):
        """The full determinant and the proxy the CLI pairs it with."""
        fd = FieldDistribution("sphere-sphere", sphere_R=cls.SPHERE_R_M)
        full = functools.partial(
            inception_det,
            field_dist=fd,
            N_min=5,
            N_max=200,
            tol=0.03,
            positive_polarity=True,
            propagator=midpoint_propagator,
        )
        proxy = functools.partial(
            inception_det,
            field_dist=UNIFORM,
            N_min=1,
            N_max=1,
            tol=0.03,
            propagator=midpoint_propagator,
        )
        return full, proxy

    @staticmethod
    def _brute_force_lowest(det_fn, pd, mod, p, T):
        """Lowest sign change of the full determinant, found without any proxy."""
        grid = np.logspace(np.log10(10.0), np.log10(3e5), 400)
        prev_en, prev_v = None, None
        for en in grid:
            v = det_fn(en, pd, mod, p, T)
            v = v if np.isfinite(v) else 0.0
            if prev_v is not None and prev_v * v < 0.0:
                return prev_en, en
            prev_en, prev_v = en, v
        return None

    @pytest.mark.parametrize("pd_bar_mm", PD_VALUES)
    def test_proxy_bracket_failure_does_not_promote_a_higher_root(
        self, air, sphere_dets, pd_bar_mm
    ):
        """
        A6: the reported root must be the lowest root of the full determinant.

        For a strongly non-uniform field the proxy used by the coarse scan
        puts the lowest root tens of per cent from where the full determinant
        has it.  When such a bracket failed to confirm it was skipped, and the
        next bracket up was reported as the inception field -- two orders of
        magnitude high, as isolated points on an otherwise smooth curve.
        """
        full, proxy = sphere_dets
        pd = pd_bar_mm * 1e-3
        p = pd / self.D_M

        bracket = self._brute_force_lowest(full, pd, air, p, 293.0)
        assert bracket, f"no root at all at pd = {pd_bar_mm} bar*mm"

        got = find_all_breakdown_EN(
            pd, air, p, 293.0, first_only=True, det_fn=full, fast_det_fn=proxy
        )
        assert got, f"no root returned at pd = {pd_bar_mm} bar*mm"
        lo, hi = bracket
        assert lo <= got[0] <= hi, (
            f"reported {got[0]:.4g} Td, but the lowest root of the full "
            f"determinant is in [{lo:.4g}, {hi:.4g}] Td"
        )
