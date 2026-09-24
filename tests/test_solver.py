# SPDX-FileCopyrightText: 2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Augmented ODE assembly, propagators and the boundary determinant."""

import math

import numpy as np
import pytest
import scipy.linalg

from incept1d.constants import c_light
from incept1d.fields import FieldDistribution
from incept1d.solver import (
    _additive_compound,
    _boundary_rows,
    _build_A_aug,
    _compound_index,
    _det_Q_compound,
    _det_Q_norm,
    _expm_shifted,
    inception_det,
    magnus2_propagator,
    midpoint_propagator,
    parse_dx_spec,
    polarities_equivalent,
)


class PhotonStub:
    """
    Minimal mechanism exposing exactly what ``_build_A_aug`` consumes.

    Lets the photon blocks be dialled independently of any real chemistry, so
    the block structure and the kappa*d collapse can be checked exactly.
    """

    def __init__(self, n=2, kappa=(1.0,), b=0.3, c=0.7):
        self.SPECIES = [f"s{i}" for i in range(n)]
        self.ELECTRON_INDEX = 0
        self._n = n
        self._kappa = np.asarray(kappa, dtype=float)
        self._b, self._c = b, c

    def get_R(self, EN, p, T):
        return np.arange(self._n * self._n, dtype=float).reshape(self._n, self._n)

    def get_V(self, EN, p, T):
        return np.diag(np.arange(1.0, self._n + 1.0))

    def get_B(self, EN, p, T):
        return np.full((self._n, len(self._kappa)), self._b)

    def get_C(self, EN, p, T):
        return np.full((len(self._kappa), self._n), self._c)

    def get_kappa(self, p, T):
        return self._kappa


class TestBuildAAug:
    def test_without_photons_reduces_to_R_Vinv(self, toy):
        """S1: with no photon groups A_aug is exactly R V^-1."""
        EN, p, T, d = 300.0, 1.0, 293.0, 10e-3
        A_aug, n_gamma, mask, n_aug = _build_A_aug(EN, d, toy, p, T)
        R, V = toy.get_R(EN, p, T), toy.get_V(EN, p, T)
        assert n_gamma == 0 and mask is None and n_aug == len(toy.SPECIES)
        assert A_aug == pytest.approx(R @ np.linalg.inv(V))

    def test_block_structure_with_photons(self):
        """S2: A_aug = [[A, B, B], [C, -D, 0], [-C, 0, D]]."""
        mod = PhotonStub(n=2, kappa=(1.0, 2.0))
        EN, p, T, d = 100.0, 1.0, 293.0, 1e-3  # kappa*d small: nothing collapses
        A_aug, n_gamma, mask, n_aug = _build_A_aug(EN, d, mod, p, T)
        assert n_gamma == 2 and n_aug == 2 + 2 * 2
        n = 2
        A = A_aug[:n, :n]
        B1, B2 = A_aug[:n, n : n + 2], A_aug[:n, n + 2 :]
        C1, C2 = A_aug[n : n + 2, :n], A_aug[n + 2 :, :n]
        D_tl, Z_tr = A_aug[n : n + 2, n : n + 2], A_aug[n : n + 2, n + 2 :]
        Z_bl, D_br = A_aug[n + 2 :, n : n + 2], A_aug[n + 2 :, n + 2 :]

        assert A == pytest.approx(
            mod.get_R(EN, p, T) @ np.linalg.inv(mod.get_V(EN, p, T))
        )
        assert B1 == pytest.approx(B2)
        assert C1 == pytest.approx(-C2)
        assert D_tl == pytest.approx(-np.diag(mod.get_kappa(p, T)))
        assert D_br == pytest.approx(np.diag(mod.get_kappa(p, T)))
        assert Z_tr == pytest.approx(np.zeros((2, 2)))
        assert Z_bl == pytest.approx(np.zeros((2, 2)))

    @pytest.mark.parametrize("kd", [0.01, 11.99, 12.01, 100.0, 1e4])
    def test_every_group_is_propagated(self, kd):
        """
        S3: no photon group is folded into A, however optically thick.

        Folding a thick group as the local source 2 B C / kappa removes the
        upstream seeding by photons, which is the photoionization feedback.
        """
        mod = PhotonStub(n=2, kappa=(1.0, 1e4))
        A_aug, n_gamma, mask, n_aug = _build_A_aug(100.0, kd, mod, 1.0, 293.0)
        assert n_gamma == 2 and n_aug == 2 + 2 * 2
        assert list(mask) == [True, True]
        v_inv = 1.0 / np.diag(mod.get_V(100.0, 1.0, 293.0))
        assert A_aug[:2, :2] == pytest.approx(
            mod.get_R(100.0, 1.0, 293.0) * v_inv[np.newaxis, :]
        )

    def test_lambda_shifts_A_and_kappa(self):
        """S4: A -> (R - lam I) V^-1 and kappa -> kappa + lam / c."""
        mod = PhotonStub(n=2, kappa=(1.0,))
        EN, p, T, d, lam = 100.0, 1.0, 293.0, 1e-3, 5e8
        base, _, _, _ = _build_A_aug(EN, d, mod, p, T, lam=0.0)
        shifted, _, _, _ = _build_A_aug(EN, d, mod, p, T, lam=lam)
        v_inv = 1.0 / np.diag(mod.get_V(EN, p, T))
        assert shifted[:2, :2] == pytest.approx(base[:2, :2] - lam * np.diag(v_inv))
        assert shifted[2, 2] == pytest.approx(-(1.0 + lam / c_light))
        assert shifted[3, 3] == pytest.approx(1.0 + lam / c_light)


class TestPropagators:
    def test_expm_shifted_differs_from_expm_by_a_positive_scalar(self, toy):
        """S5: the dropped factor is exp(lam_max), which cannot change signs."""
        A, _, _, _ = _build_A_aug(300.0, 10e-3, toy, 1.0, 293.0)
        M = A * 1e-3
        lam = max(0.0, float(np.max(np.real(np.linalg.eigvals(M)))))
        assert _expm_shifted(M) * math.exp(lam) == pytest.approx(
            scipy.linalg.expm(M), rel=1e-10
        )

    def test_midpoint_on_constant_A(self, toy):
        """S6: for constant A the midpoint rule is exp(A h), up to the shift."""
        A, _, _, _ = _build_A_aug(300.0, 10e-3, toy, 1.0, 293.0)
        h = 2e-3
        P = midpoint_propagator(lambda x: A, 0.0, h)
        assert P == pytest.approx(_expm_shifted(A * h), rel=1e-12)

    def test_magnus2_reduces_to_midpoint_for_constant_A(self, toy):
        """S7: the commutator term vanishes when A does not vary."""
        A, _, _, _ = _build_A_aug(300.0, 10e-3, toy, 1.0, 293.0)
        h = 2e-3
        assert magnus2_propagator(lambda x: A, 0.0, h) == pytest.approx(
            midpoint_propagator(lambda x: A, 0.0, h), rel=1e-10
        )

    def test_composed_propagators_equal_one_exponential(self, toy):
        """S8: the identity the uniform-field fast path relies on."""
        A, _, _, _ = _build_A_aug(300.0, 10e-3, toy, 1.0, 293.0)
        d, N = 10e-3, 8
        M = np.eye(A.shape[0])
        for i in range(N):
            M = midpoint_propagator(lambda x: A, i * d / N, (i + 1) * d / N) @ M
        assert M == pytest.approx(_expm_shifted(A * d), rel=1e-9)


class TestUniformFastPath:
    """S9/S10: the shortcut must be exact, not merely close."""

    @staticmethod
    def _stepped():
        """A FieldDistribution with f(xi) = 1 that bypasses the fast path."""
        fd = FieldDistribution("uniform")
        fd.field_type = "stepped-uniform"
        fd.build = lambda d: (lambda xi: 1.0)
        return fd

    @pytest.mark.parametrize("pd_mm", [1e-2, 1.0, 100.0])
    @pytest.mark.parametrize("EN", [80.0, 300.0, 3000.0])
    def test_matches_the_stepped_path(self, toy, uniform, pd_mm, EN):
        p, T, pd = 1.0, 293.0, pd_mm * 1e-3
        fast = inception_det(EN, pd, toy, p, T, uniform)
        slow = inception_det(EN, pd, toy, p, T, self._stepped())
        if math.isnan(fast) or math.isnan(slow):
            assert math.isnan(fast) and math.isnan(slow)
        else:
            assert fast == pytest.approx(slow, rel=1e-9)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {},
            {"N_min": 1, "N_max": 1},
            {"N_min": 40, "N_max": 400, "tol": 1e-3},
            {"propagator": magnus2_propagator},
        ],
    )
    def test_grid_and_method_do_not_matter(self, toy, uniform, kwargs):
        """S10: the propagator is exact, so --dx and --method are inert."""
        p, T, pd, EN = 1.0, 293.0, 10e-3, 120.0
        ref = inception_det(EN, pd, toy, p, T, uniform)
        assert inception_det(EN, pd, toy, p, T, uniform, **kwargs) == pytest.approx(
            ref, rel=1e-12
        )

    def test_one_matrix_exponential_per_evaluation(self, toy, uniform, monkeypatch):
        """S9b: guard the saving itself, not only the value."""
        calls = {"n": 0}
        real = scipy.linalg.expm

        def counting(M):
            calls["n"] += 1
            return real(M)

        monkeypatch.setattr("incept1d.solver.scipy.linalg.expm", counting)
        inception_det(120.0, 10e-3, toy, 1.0, 293.0, uniform)
        assert calls["n"] == 1


class TestAdaptiveGrid:
    def test_diagnostics_report_the_segments(self, toy):
        """S11: a non-uniform field exercises the adaptive machinery."""
        fd = FieldDistribution("sphere-plane", 50e-3)
        diag = []
        inception_det(
            120.0, 10e-3, toy, 1.0, 293.0, fd, N_min=5, N_max=200, diag_list=diag
        )
        assert len(diag) == 5
        assert all(entry["halvings"] for entry in diag)

    def test_n_min_equals_n_max_disables_halving(self, toy):
        """S11b: max_depth = 0 means each segment is accepted as one step."""
        fd = FieldDistribution("sphere-plane", 50e-3)
        diag = []
        inception_det(
            120.0, 10e-3, toy, 1.0, 293.0, fd, N_min=4, N_max=4, diag_list=diag
        )
        assert len(diag) == 4
        for entry in diag:
            assert all(h["level"] == 0 for h in entry["halvings"])


class TestParseDxSpec:
    def test_defaults(self):
        assert parse_dx_spec(None) == (5, 1000, 1e-3)
        assert parse_dx_spec([]) == (5, 1000, 1e-3)

    def test_partial_specs(self):
        assert parse_dx_spec(["10"]) == (10, 1000, 1e-3)
        assert parse_dx_spec(["10", "400"]) == (10, 400, 1e-3)
        assert parse_dx_spec(["10", "400", "0.01"]) == (10, 400, 0.01)

    @pytest.mark.parametrize(
        "tokens, match",
        [
            (["0"], "N_min must be"),
            (["10", "5"], "must be ≥ N_min"),
            (["5", "200", "0"], "tol must be in"),
            (["5", "200", "1.5"], "tol must be in"),
            (["x"], "--dx"),
        ],
    )
    def test_validation(self, tokens, match):
        """S12: every malformed --dx is rejected with a specific message."""
        with pytest.raises(ValueError, match=match):
            parse_dx_spec(tokens)


class TestDeterminant:
    def test_sign_changes_across_the_inception_field(self, toy, uniform, det_uniform):
        """S14: det Q has opposite signs either side of the root."""
        from incept1d.inception import find_all_breakdown_EN

        p, T, pd = 1.0, 293.0, 20e-3
        root = find_all_breakdown_EN(
            pd, toy, p, T, first_only=True, det_fn=det_uniform
        )[0]
        below = inception_det(root * 0.9, pd, toy, p, T, uniform)
        above = inception_det(root * 1.1, pd, toy, p, T, uniform)
        assert np.isfinite(below) and np.isfinite(above)
        assert below * above < 0.0

    def test_determinant_is_finite_over_the_working_range(self, toy, uniform):
        """A NaN here means Q was ill-conditioned; it must not be the norm."""
        vals = [
            inception_det(EN, 10e-3, toy, 1.0, 293.0, uniform)
            for EN in np.logspace(1, 3, 40)
        ]
        assert sum(np.isfinite(v) for v in vals) > 30


class TestPolarityEquivalence:
    """When may negative polarity reuse the positive-polarity answer?"""

    @pytest.fixture
    def overridden(self, toy_path):
        from incept1d.mechanism import load_mechanism

        return load_mechanism(
            toy_path,
            {"pos_override": {"gamma0": 0.4}, "neg_override": {"gamma0": 0.01}},
        )

    @pytest.mark.parametrize(
        "field, expected",
        [
            (FieldDistribution("uniform"), True),
            (FieldDistribution("sphere-sphere", 5e-3), True),
            (FieldDistribution("sphere-plane", 5e-3), False),
        ],
    )
    def test_geometry_decides_for_identical_electrodes(self, toy, field, expected):
        assert polarities_equivalent(toy, field) is expected

    def test_overrides_break_the_symmetry_of_a_uniform_gap(self, overridden):
        """Two different cathodes are not a reflection of each other."""
        assert not polarities_equivalent(overridden, FieldDistribution("uniform"))

    def test_the_determinant_agrees(self, toy, overridden):
        """
        The helper's claim, checked against det Q itself.

        In a uniform gap det Q is the same for both polarities exactly when
        the electrodes are identical.
        """
        uniform = FieldDistribution("uniform")
        args = (400.0, 20e-3)

        def det(mod, positive):
            return inception_det(
                *args, mod, 1.0, 293.0, uniform, positive_polarity=positive
            )

        assert det(toy, True) == det(toy, False)
        assert det(overridden, True) != pytest.approx(det(overridden, False))


class TestCompoundRoute:
    """The singular-Q fallback must return the same number as the direct route."""

    @pytest.mark.parametrize("n, k", [(3, 1), (4, 2), (6, 2), (6, 3)])
    def test_additive_compound_generates_the_compound(self, n, k):
        """exp(Ω^[k]) acts on k×k minors exactly as exp(Ω) acts on vectors."""
        rng = np.random.default_rng(n * 10 + k)
        Omega = rng.normal(size=(n, n))
        P = scipy.linalg.expm(Omega)
        subsets = _compound_index(n, k)[0]
        minors = np.array(
            [[np.linalg.det(P[np.ix_(r, c)]) for c in subsets] for r in subsets]
        )
        assert scipy.linalg.expm(_additive_compound(Omega, k)) == pytest.approx(
            minors, rel=1e-9, abs=1e-12
        )

    @staticmethod
    def _routes(mod, EN, pd, lam=0.0, n_steps=1, p=1.0, T=293.0):
        d = pd / p
        A, Ng, mask, n_aug = _build_A_aug(EN, d, mod, p, T, lam=lam)
        Q0, S = _boundary_rows(EN, mod, p, T, Ng, mask, n_aug)
        exponents = [A * d / n_steps] * n_steps
        M = np.linalg.matrix_power(scipy.linalg.expm(A * d / n_steps), n_steps)
        with np.errstate(over="ignore", invalid="ignore"):
            direct = _det_Q_norm(np.vstack([Q0, S @ (M / np.abs(M).max())]))
            compound = _det_Q_compound(exponents, Q0, S)
        return direct, compound

    @pytest.mark.parametrize("EN", [150.0, 300.0, 400.0, 600.0])
    @pytest.mark.parametrize("lam", [0.0, 1e5])
    @pytest.mark.parametrize("n_steps", [1, 7])
    def test_equals_the_direct_route(self, toy, EN, lam, n_steps):
        direct, compound = self._routes(toy, EN, 20e-3, lam, n_steps)
        assert np.isfinite(direct)
        assert compound == pytest.approx(direct, rel=1e-7)

    # A hand-built gap the toy cannot provide: electrons and two positive
    # ions, so the anode block S M has two rows that can become parallel.
    # Ions feed nothing back, as in the air mechanisms.  Rows of Q0/S and
    # the generator in the order (e, ion1, ion2).
    Q0_TWO_IONS = np.array([[1.0, 1e-3, 2e-3]])
    S_TWO_IONS = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    @staticmethod
    def _generator(a, c1=300.0, c2=320.0):
        """Electron growth a, ion production, ion-mode growth c1, c2 (1/m)."""
        return np.array(
            [
                [a, 0.0, 0.0],
                [0.6 * abs(a) + 50.0, c1, 0.0],
                [0.4 * abs(a) + 50.0, 0.0, c2],
            ]
        )

    def _exact(self, exponents):
        """Row-normalised det [Q0; S M] in 300-digit arithmetic."""
        mp = pytest.importorskip("mpmath")
        mp.mp.dps = 300
        M = mp.eye(3)
        for Om in exponents:
            M = mp.expm(mp.matrix(Om.tolist())) * M
        Q = mp.matrix(3, 3)
        rows = [list(self.Q0_TWO_IONS[0])] + [
            [sum(self.S_TWO_IONS[r, j] * M[j, c] for j in range(3)) for c in range(3)]
            for r in range(2)
        ]
        for r, row in enumerate(rows):
            norm = mp.sqrt(sum(mp.mpf(v) ** 2 for v in row))
            for c in range(3):
                Q[r, c] = row[c] / norm
        return float(mp.det(Q))

    @pytest.mark.parametrize("scale", [1.0, 30.0, 300.0])
    def test_uniform_gap_with_widely_separated_modes(self, scale):
        """
        The ion modes run far ahead of the electrons: at scale 300 they are
        ~10^4 e-folds apart and the direct route has long since failed.
        """
        exponents = [self._generator(200.0, 300.0 * scale, 320.0 * scale) * 0.01]
        exact = self._exact(exponents)
        with np.errstate(over="ignore", invalid="ignore"):
            got = _det_Q_compound(exponents, self.Q0_TWO_IONS, self.S_TWO_IONS)
        assert np.sign(got) == np.sign(exact)
        assert got == pytest.approx(exact, rel=1e-6, abs=1e-280)

    @pytest.mark.parametrize("n_steps", [10, 200])
    @pytest.mark.parametrize("c", [6000.0, 9000.0, 10500.0, 11000.0])
    def test_transiently_tiny_electrons(self, n_steps, c):
        """
        Electrons decay by 160 e-folds over the first half of the gap and
        avalanche by 600 over the second, so in mid-gap their component is
        far below the ions' yet decides the sign at the anode.  That defeats
        re-orthonormalisation (Godunov–Conte), which mixes components; and
        for the two smaller ion growth rates c the rows S M are parallel, so
        the direct route fails too.  The root lies between c = 10500 and
        11000, so the sign flips across the parametrisation.
        """
        d = 0.04
        h = d / n_steps
        exponents = [
            self._generator(-8000.0 if i < n_steps // 2 else 30000.0, c, 1.07 * c) * h
            for i in range(n_steps)
        ]
        exact = self._exact(exponents)
        assert exact != 0.0
        with np.errstate(over="ignore", invalid="ignore"):
            got = _det_Q_compound(exponents, self.Q0_TWO_IONS, self.S_TWO_IONS)
        assert np.sign(got) == np.sign(exact)
        assert got == pytest.approx(exact, rel=1e-6)
