# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Verification of the inception criterion against the closed-form reduced model,
and regression tests for the root finder.

The tests in the first class are the core of the suite: they check the whole
stack -- mechanism loading, A_aug assembly, propagation, boundary conditions,
determinant, root finding -- against algebra derived independently in
``tests/closed_form.py``, not against previously recorded output.
"""

import functools

import numpy as np
import pytest

from closed_form import (
    generalized_paschen_lhs,
    lambda_pm,
    solve_EN,
    standard_paschen_lhs,
)
from incept1d.constants import kB
from incept1d.eigenvalues import max_real_eigenvalue
from incept1d.fields import FieldDistribution
from incept1d.inception import (
    _accept_root,
    compute_inception_curve,
    find_all_breakdown_EN,
)
from incept1d.solver import inception_det, riccati_criterion


def _coeffs(raw, EN, p=1.0, T=293.0):
    return raw.alpha(EN, p, T), raw.eta(EN, p, T), raw.delta(EN, p, T)


class TestClosedForm:
    """The solver must reproduce eq_generalized_paschen and its limits."""

    @pytest.mark.parametrize("d_mm", [0.5, 1.0, 5.0, 20.0, 100.0])
    def test_root_matches_generalized_paschen(self, toy, toy_raw, det_uniform, d_mm):
        """I1: det Q = 0 has the same E/N root as eq_generalized_paschen."""
        p, T, d = 1.0, 293.0, d_mm * 1e-3
        roots = find_all_breakdown_EN(
            p * d, toy, p, T, first_only=True, det_fn=det_uniform
        )
        assert roots, f"solver found no root at d = {d_mm} mm"
        expected = solve_EN(
            lambda e: toy_raw.alpha(e, p, T),
            lambda e: toy_raw.eta(e, p, T),
            lambda e: toy_raw.delta(e, p, T),
            toy_raw._GAMMA0,
            d,
        )
        assert roots[0] == pytest.approx(expected, rel=1e-7)

    @pytest.mark.parametrize("d_mm", [1.0, 10.0, 50.0])
    def test_criterion_is_unity_at_the_root(self, toy, toy_raw, det_uniform, d_mm):
        """I1b: the closed-form LHS equals 1 at the root the solver returns."""
        p, T, d = 1.0, 293.0, d_mm * 1e-3
        roots = find_all_breakdown_EN(
            p * d, toy, p, T, first_only=True, det_fn=det_uniform
        )
        lhs = generalized_paschen_lhs(*_coeffs(toy_raw, roots[0]), toy_raw._GAMMA0, d)
        assert lhs == pytest.approx(1.0, rel=1e-6)

    @pytest.mark.parametrize("d_mm", [1.0, 10.0, 50.0])
    def test_no_detachment_reduces_to_standard_paschen(self, toy_path, d_mm):
        """I2: delta = 0 must reproduce eq_standard_paschen."""
        from incept1d.mechanism import load_mechanism

        mod = load_mechanism(toy_path, {"delta_0": 0.0})
        raw = mod._mod
        p, T, d = 1.0, 293.0, d_mm * 1e-3
        det = functools.partial(inception_det, field_dist=FieldDistribution("uniform"))
        roots = find_all_breakdown_EN(p * d, mod, p, T, first_only=True, det_fn=det)
        assert roots
        a, e, dl = _coeffs(raw, roots[0])
        assert dl == 0.0
        assert standard_paschen_lhs(a, e, raw._GAMMA0, d) == pytest.approx(
            1.0, rel=1e-6
        )

    @pytest.mark.parametrize("d_mm", [1.0, 10.0])
    def test_no_attachment_reduces_to_textbook_paschen(self, toy_path, d_mm):
        """I3: with eta = delta = 0, alpha d = ln(1 + 1/gamma)."""
        from incept1d.mechanism import load_mechanism

        mod = load_mechanism(toy_path, {"eta_0": 0.0, "delta_0": 0.0})
        raw = mod._mod
        p, T, d = 1.0, 293.0, d_mm * 1e-3
        det = functools.partial(inception_det, field_dist=FieldDistribution("uniform"))
        roots = find_all_breakdown_EN(p * d, mod, p, T, first_only=True, det_fn=det)
        assert roots
        a = raw.alpha(roots[0], p, T)
        assert a * d == pytest.approx(np.log1p(1.0 / raw._GAMMA0), rel=1e-6)

    @pytest.mark.parametrize("d_mm", [200.0, 1000.0])
    def test_solution_exists_in_the_attachment_regime(self, toy_path, d_mm):
        """
        I4: alpha > eta is not required.

        The documentation notes that eq_standard_paschen also has solutions
        when alpha <= eta: cathode-emitted electrons still make positive ions
        before their avalanches terminate.  A strongly attaching gas in a long
        gap puts the root there.  In that limit e^{(alpha-eta)d} -> 0 and the
        criterion collapses to gamma alpha / (eta - alpha) = 1, i.e.

            alpha -> eta / (1 + gamma),

        which is a sharper check than the inequality alone.
        """
        from incept1d.mechanism import load_mechanism

        mod = load_mechanism(toy_path, {"eta_0": 1000.0, "delta_0": 0.0})
        raw = mod._mod
        p, T, d = 1.0, 293.0, d_mm * 1e-3
        det = functools.partial(inception_det, field_dist=FieldDistribution("uniform"))
        roots = find_all_breakdown_EN(p * d, mod, p, T, first_only=True, det_fn=det)
        assert roots, "no root found in the attachment regime"
        a, e, _ = _coeffs(raw, roots[0])
        assert a <= e, "this configuration was meant to probe alpha <= eta"
        assert standard_paschen_lhs(a, e, raw._GAMMA0, d) == pytest.approx(
            1.0, rel=1e-6
        )
        # Long-gap asymptote; loosest at the shortest gap tested.
        assert a == pytest.approx(e / (1.0 + raw._GAMMA0), rel=3e-3)

    @pytest.mark.parametrize("EN", [50.0, 120.0, 400.0, 2000.0])
    def test_leading_eigenvalue_equals_lambda_plus(self, toy, toy_raw, EN):
        """I5: max Re eig(R V^-1) equals the analytic lambda_+."""
        p, T = 1.0, 293.0
        lam_p, _, _ = lambda_pm(*_coeffs(toy_raw, EN, p, T))
        assert max_real_eigenvalue(toy, EN, p, T) == pytest.approx(lam_p, rel=1e-10)

    def test_lambda_plus_reduces_to_alpha_minus_eta(self, toy_path):
        """I5b: with delta = 0 the apparent coefficient is exactly alpha - eta."""
        from incept1d.mechanism import load_mechanism

        mod = load_mechanism(toy_path, {"delta_0": 0.0})
        raw = mod._mod
        p, T, EN = 1.0, 293.0, 300.0
        a, e, _ = _coeffs(raw, EN, p, T)
        assert max_real_eigenvalue(mod, EN, p, T) == pytest.approx(a - e, rel=1e-10)


class TestRootFinding:
    """Regressions for the warm-start lock-in and the branch contract."""

    def test_lowest_root_is_returned_despite_a_stale_hint(self, toy, det_uniform):
        """
        I6: first_only=True must return the globally lowest root even when the
        caller supplies a hint far above it.

        This is the defect that made the dry-air curve miss its Paschen
        minimum: a hint from the previous pd point was refined and returned
        without checking whether a lower root had appeared.
        """
        p, T, d = 1.0, 293.0, 20e-3
        truth = find_all_breakdown_EN(
            p * d, toy, p, T, first_only=True, det_fn=det_uniform
        )[0]
        for hint in (10.0 * truth, 100.0 * truth, 1000.0 * truth):
            got = find_all_breakdown_EN(
                p * d, toy, p, T, first_only=True, det_fn=det_uniform, EN_hints=[hint]
            )
            assert got, f"no root returned for hint {hint}"
            assert got[0] == pytest.approx(
                truth, rel=1e-8
            ), f"stale hint {hint:.4g} Td changed the answer"

    def test_first_only_agrees_with_min_of_all_roots(self, toy, det_uniform):
        """I6b/I7: branch 1 is exactly min(all roots), and all roots are sorted."""
        p, T = 1.0, 293.0
        for d_mm in (1.0, 10.0, 50.0):
            pd = p * d_mm * 1e-3
            first = find_all_breakdown_EN(
                pd, toy, p, T, first_only=True, det_fn=det_uniform
            )
            every = find_all_breakdown_EN(
                pd, toy, p, T, first_only=False, det_fn=det_uniform
            )
            assert every == sorted(every)
            assert first[0] == pytest.approx(min(every), rel=1e-8)

    def test_every_root_includes_one_below_the_scan_range(self, toy, uniform):
        """
        With a physically signed criterion, a root below EN_lo is found by
        searching down from it — also when every root is wanted, not only
        the lowest.  EN_lo is put above the true root to force that case.
        """
        det = functools.partial(riccati_criterion, field_dist=uniform)
        p, T, pd = 1.0, 293.0, 20e-3
        truth = find_all_breakdown_EN(pd, toy, p, T, first_only=True, det_fn=det)[0]
        for first_only in (True, False):
            got = find_all_breakdown_EN(
                pd, toy, p, T, EN_lo=2.0 * truth, first_only=first_only, det_fn=det
            )
            assert got, f"first_only={first_only}: no root returned"
            assert got == sorted(got)
            assert got[0] == pytest.approx(truth, rel=1e-8)

    def test_lowest_root_does_not_pay_for_the_whole_scan(self, toy, counting_det):
        """
        I8a: the scan runs bottom up and stops at the first confirmed root,
        so a cold search in a uniform gap costs far fewer evaluations than
        the full 200-point scan.
        """
        fn, counter = counting_det
        p, T, pd = 1.0, 293.0, 20e-3
        find_all_breakdown_EN(pd, toy, p, T, first_only=True, det_fn=fn)
        assert counter["n"] < 100, f"{counter['n']} evaluations"

    @pytest.mark.parametrize("pd", [5e-3, 20e-3])
    def test_warm_start_saves_work_in_a_non_uniform_gap(self, toy, pd):
        """
        I8b: the guard must not silently disable the warm start.

        In a non-uniform gap the uniform proxy misleads a cold search into
        unconfirmed brackets and fallback scans; a hint from the neighbouring
        pd point avoids them (on a real sweep, dropping the hints cost 2.3x
        the CPU time), and must find the same root.
        """
        fd, uniform = FieldDistribution("sphere-plane", 1e-3), FieldDistribution(
            "uniform"
        )
        n = {"c": 0}

        def counted(field, N=None):
            def f(EN, pd_, mod, p_, T_):
                n["c"] += 1
                grid = dict(N_min=N, N_max=N) if N else {}
                return riccati_criterion(EN, pd_, mod, p_, T_, field, **grid)

            return f

        full, proxy = counted(fd), counted(uniform, 1)
        p, T = pd / 10e-3, 293.0
        cold = find_all_breakdown_EN(
            pd, toy, p, T, first_only=True, det_fn=full, fast_det_fn=proxy
        )
        n_cold, n["c"] = n["c"], 0
        warm = find_all_breakdown_EN(
            pd,
            toy,
            p,
            T,
            first_only=True,
            det_fn=full,
            fast_det_fn=proxy,
            EN_hints=[cold[0] * 1.03],
        )
        assert warm[0] == pytest.approx(cold[0], rel=1e-6)
        assert n["c"] < n_cold, f"warm {n['c']} vs cold {n_cold} evaluations"

    def test_riccati_root_must_cross_from_plus_to_minus(self, toy, uniform, capsys):
        """
        I10: with a physically signed criterion, + -> - is accepted as is
        (zero or pole), and - -> + is rejected: it cannot be inception.
        """
        ric = functools.partial(riccati_criterion, field_dist=uniform)
        pd, p, T = 20e-3, 1.0, 293.0
        assert _accept_root(80.0, 70.0, 90.0, 1.0, -1.0, pd, ric, toy, p, T) is True
        assert capsys.readouterr().out == ""  # no "suspect root" for a pole
        assert _accept_root(80.0, 70.0, 90.0, -1.0, 1e6, pd, ric, toy, p, T) is False
        assert "not inception" in capsys.readouterr().out

    def test_accept_root_rejects_the_nan_sentinel_artefact(self, toy):
        """I9: a root bracketed by a NaN sentinel is discarded, not reported."""
        p, T, pd = 1.0, 293.0, 20e-3
        nan_det = lambda EN, *a, **k: float("nan")  # noqa: E731
        accepted = _accept_root(
            300.0, 100.0, 900.0, -1e-300, 1.0, pd, nan_det, toy, p, T
        )
        assert accepted is False

    def test_accept_root_keeps_a_genuine_singularity(self, toy):
        """I9b: NaN at the root with finite brackets is the expected Q -> 0."""
        p, T, pd = 1.0, 293.0, 20e-3
        nan_det = lambda EN, *a, **k: float("nan")  # noqa: E731
        assert (
            _accept_root(300.0, 100.0, 900.0, 1.0, -1.0, pd, nan_det, toy, p, T) is True
        )


class TestInceptionCurve:
    """compute_inception_curve: branch bookkeeping and derived quantities."""

    def test_voltage_follows_from_reduced_field(self, toy, det_uniform):
        """I11: V = (E/N) pd 1e-21 / (kB T) * 1e5, and arrays stay consistent."""
        p, T = 1.0, 293.0
        pd_arr = np.logspace(np.log10(1e-3), np.log10(50e-3), 8)
        branches = compute_inception_curve(pd_arr, toy, p, T, det_fn=det_uniform)
        assert branches, "no branch found"
        br = branches[0]
        n = len(br["idx"])
        assert all(len(br[k]) == n for k in ("pd", "EN", "V", "p", "d"))
        expected = br["EN"] * br["pd"] * 1e-21 / (kB * T) * 1e5
        assert br["V"] == pytest.approx(expected, rel=1e-12)
        assert br["d"] == pytest.approx(br["pd"] / br["p"], rel=1e-12)

    def test_branch_one_is_the_lowest_root_at_every_point(self, toy, det_uniform):
        """
        A5 (toy version): an invariant rather than a pinned value.

        At every pd of the sweep, branch 1 must equal the lowest root found by
        an independent all-branches search.
        """
        p, T = 1.0, 293.0
        pd_arr = np.logspace(np.log10(2e-3), np.log10(50e-3), 6)
        branches = compute_inception_curve(pd_arr, toy, p, T, det_fn=det_uniform)
        br = branches[0]
        for i, pd in zip(br["idx"], br["pd"]):
            every = find_all_breakdown_EN(
                pd, toy, p, T, first_only=False, det_fn=det_uniform
            )
            assert br["EN"][list(br["idx"]).index(i)] == pytest.approx(
                min(every), rel=1e-6
            )

    def test_branches_are_tracked_by_continuity(self, toy, det_uniform):
        """I10: a single smooth branch stays one branch across the sweep."""
        p, T = 1.0, 293.0
        pd_arr = np.logspace(np.log10(2e-3), np.log10(100e-3), 12)
        branches = compute_inception_curve(
            pd_arr, toy, p, T, all_branches=True, det_fn=det_uniform
        )
        assert branches
        br = branches[0]
        order = np.argsort(br["pd"])
        EN = br["EN"][order]
        # The toy inception field falls monotonically with pd.
        assert np.all(np.diff(EN) < 0.0)
