import unittest
import numpy as np

from synim.registration.model import GuideStar, DM, WFS, System
from synim.registration.reconstruction import (
    ParameterSpec, jacobian, apply_alpha, local_params_vector,
)
from synim.registration.analysis import (
    normalize_jacobian, condition_number, describe_mode,
    reconstruction_matrix, covariance_from_noise,
    monte_carlo_gauss_newton, analyze,
)


def _build_3wfs_2dm_system(pixel_pitch=0.2):
    system = System(pixel_pitch=pixel_pitch)
    gs_positions = [(0.0, 0.0), (17.5, 0.0), (-12.0, 9.0)]
    for i, pos in enumerate(gs_positions):
        gs = GuideStar(name=f"gs{i}", position=pos,
                        height=90000.0 if pos != (0.0, 0.0) else np.inf)
        system.add_wfs(WFS(name=f"wfs{i}", guide_star=gs))
    system.add_dm(DM(name="dm0", height=0.0))
    system.add_dm(DM(name="dm1", height=6000.0))
    return system


class TestSvdAndConditionNumber(unittest.TestCase):
    def test_condition_number_matches_manual_svd(self):
        Lambda = np.array([[3.0, 0.0], [0.0, 1.0]])
        self.assertAlmostEqual(condition_number(Lambda), 3.0)

    def test_condition_number_is_huge_for_singular_matrix(self):
        Lambda = np.array([[1.0, 2.0], [2.0, 4.0]])  # rank 1
        # The smallest singular value is ~0 up to floating-point round-off
        # (rarely bit-exact 0), so the condition number is huge rather than
        # necessarily inf; both are "hard degeneracy" in practice.
        self.assertGreater(condition_number(Lambda), 1e10)

    def test_condition_number_is_infinite_for_exactly_zero_row(self):
        Lambda = np.array([[1.0, 0.0], [0.0, 0.0]])  # exactly rank 1
        self.assertEqual(condition_number(Lambda), np.inf)


class TestDescribeMode(unittest.TestCase):
    def test_picks_dominant_components(self):
        specs = [ParameterSpec("dm", "dm0", "shift", 0),
                 ParameterSpec("dm", "dm1", "shift", 0),
                 ParameterSpec("wfs", "wfs0", "shift", 0)]
        v = np.array([0.9, 0.1, -0.05])
        described = describe_mode(v, specs, top_n=2)
        self.assertEqual(described[0][0], "dm:dm0:shift_x")
        self.assertAlmostEqual(described[0][1], 0.9)
        self.assertEqual(described[1][0], "dm:dm1:shift_x")


class TestReconstructionMatrix(unittest.TestCase):
    def test_full_rank_matches_pinv(self):
        rng = np.random.default_rng(0)
        Lambda = rng.uniform(-1, 1, size=(6, 3))
        R, n_kept, S = reconstruction_matrix(Lambda)
        self.assertEqual(n_kept, 3)
        np.testing.assert_allclose(R, np.linalg.pinv(Lambda), atol=1e-8)

    def test_truncation_reduces_kept_modes(self):
        Lambda = np.diag([10.0, 1.0, 1e-8])
        R_full, n_full, _ = reconstruction_matrix(Lambda, rcond=1e-10)
        R_trunc, n_trunc, _ = reconstruction_matrix(Lambda, n_modes=2)
        self.assertEqual(n_full, 3)
        self.assertEqual(n_trunc, 2)
        # the truncated reconstruction matrix must not try to invert the
        # near-zero singular value
        self.assertAlmostEqual(R_trunc[2, 2], 0.0)


class TestDegenerateCase(unittest.TestCase):
    """
    Sec. 2/5 of the paper: with a single WFS-DM pair, a GS position error
    and a DM shift produce the same (x-only) local shift and cannot be
    told apart - this must show up as a huge condition number.
    """

    def test_gs_shift_and_dm_shift_are_degenerate_with_one_pair(self):
        gs = GuideStar(name="lgs", position=(0.0, 0.0), height=90000.0)
        system = System(pixel_pitch=0.2)
        system.add_wfs(WFS(name="wfs0", guide_star=gs))
        system.add_dm(DM(name="dm0", height=6000.0))

        specs = [ParameterSpec("dm", "dm0", "shift", 0),
                 ParameterSpec("gs", "lgs", "position_shift", 0)]
        pairs = system.pairs()

        Lambda = jacobian(system, specs, pairs, dof=("shift_x", "shift_y"))
        self.assertGreater(condition_number(Lambda), 1e6)


class TestNoisePropagation(unittest.TestCase):
    """Analytic covariance (through the pseudo-inverse) must match Monte
    Carlo, for a well-conditioned, non-overlapping-shifts setup."""

    def test_analytic_matches_montecarlo(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        specs = [
            ParameterSpec("dm", "dm0", "shift", 0),
            ParameterSpec("dm", "dm1", "shift", 0),
            ParameterSpec("wfs", "wfs0", "shift", 0),
        ]
        dof = ("shift_x", "shift_y")
        sigma = 0.05

        rng = np.random.default_rng(42)
        report = analyze(system, specs, pairs, dof=dof, sigma=sigma,
                          n_trials=8000, rng=rng)

        np.testing.assert_allclose(report["analytic_std"], report["montecarlo_std"],
                                    rtol=0.15)
        # cross-check the full covariance matrices too (not just the diagonal)
        rel_diff = (np.linalg.norm(report["analytic_covariance"] - report["montecarlo_covariance"])
                    / np.linalg.norm(report["analytic_covariance"]))
        self.assertLess(rel_diff, 0.2)

    def test_analytic_covariance_matches_reconstruction_matrix_definition(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        specs = [ParameterSpec("dm", "dm0", "shift", 0),
                 ParameterSpec("dm", "dm1", "shift", 0),
                 ParameterSpec("wfs", "wfs0", "shift", 0)]
        Lambda = jacobian(system, specs, pairs, dof=("shift_x", "shift_y"))
        sigma = 0.03
        cov, std = covariance_from_noise(Lambda, sigma)
        R, _, _ = reconstruction_matrix(Lambda)
        expected_cov = (R * sigma ** 2) @ R.T
        np.testing.assert_allclose(cov, expected_cov)
        np.testing.assert_allclose(std, np.sqrt(np.diag(expected_cov)))

    def test_analyze_exposes_montecarlo_mean(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        specs = [ParameterSpec("dm", "dm0", "shift", 0)]
        report = analyze(system, specs, pairs, dof=("shift_x", "shift_y"), sigma=0.02,
                          n_trials=200, rng=np.random.default_rng(0))
        self.assertIn("montecarlo_mean", report)
        self.assertEqual(report["montecarlo_mean"].shape, (1,))


class TestMonteCarloGaussNewton(unittest.TestCase):
    """
    Unlike `monte_carlo_noise_propagation` (a single linear step
    linearized exactly at the truth, unbiased by construction),
    `monte_carlo_gauss_newton` runs the real iterative estimator starting
    from a given system - here the untouched nominal one, as a real
    estimator would.
    """

    def test_matches_true_alpha_for_a_well_conditioned_case(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        specs = [
            ParameterSpec("dm", "dm0", "shift", 0),
            ParameterSpec("dm", "dm1", "shift", 0),
            ParameterSpec("wfs", "wfs0", "shift", 0),
        ]
        dof = ("shift_x", "shift_y")
        true_alpha = np.array([0.3, -0.4, 0.2])
        true_system = apply_alpha(system, specs, true_alpha)
        D_true = local_params_vector(true_system, pairs, dof=dof)

        rng = np.random.default_rng(1)
        samples, mean, covariance = monte_carlo_gauss_newton(
            system, specs, pairs, D_true, sigma=0.01, dof=dof,
            n_trials=200, n_iter=3, rng=rng)

        self.assertEqual(samples.shape, (200, 3))
        self.assertEqual(covariance.shape, (3, 3))
        # non-overlapping shifts, near-linear: negligible bias expected
        np.testing.assert_allclose(mean, true_alpha, atol=0.01)


class TestNormalizeJacobian(unittest.TestCase):
    def test_scales_columns(self):
        Lambda = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        scaled = normalize_jacobian(Lambda, [None, None, None], [2.0, 0.5, 1.0])
        np.testing.assert_allclose(scaled, [[2.0, 1.0, 3.0], [8.0, 2.5, 6.0]])

    def test_dict_scale_uses_spec_field(self):
        specs = [ParameterSpec("dm", "dm0", "shift", 0),
                 ParameterSpec("dm", "dm0", "rotation")]
        Lambda = np.array([[1.0, 1.0]])
        scaled = normalize_jacobian(Lambda, specs, {"shift": 10.0, "rotation": 0.1})
        np.testing.assert_allclose(scaled, [[10.0, 0.1]])

    def test_missing_field_in_dict_raises(self):
        specs = [ParameterSpec("dm", "dm0", "magnification")]
        with self.assertRaises(KeyError):
            normalize_jacobian(np.array([[1.0]]), specs, {"shift": 1.0})

    def test_does_not_change_rank_or_exact_null_space(self):
        # The "shift every WFS + every DM together" gauge freedom (see
        # 03_sensitivity_and_degeneracy.py in the MORFEO example) is exact
        # only when every WFS's guide star shares the same height, so the
        # DM-side compensation does not depend on which WFS it is paired
        # with - build such a system explicitly here (uniform LGS height).
        system = System(pixel_pitch=0.2)
        for i, angle in enumerate([0.0, 120.0, 240.0]):
            gs = GuideStar(name=f"gs{i}", position=tuple(30.0 * np.array(
                [np.cos(np.radians(angle)), np.sin(np.radians(angle))])), height=90000.0)
            system.add_wfs(WFS(name=f"wfs{i}", guide_star=gs))
        system.add_dm(DM(name="dm0", height=0.0))
        system.add_dm(DM(name="dm1", height=6000.0))
        pairs = system.pairs()

        specs = []
        for w in ("wfs0", "wfs1", "wfs2"):
            specs += [ParameterSpec("wfs", w, "shift", 0), ParameterSpec("wfs", w, "shift", 1)]
        for d in ("dm0", "dm1"):
            specs += [ParameterSpec("dm", d, "shift", 0), ParameterSpec("dm", d, "shift", 1)]
        Lambda = jacobian(system, specs, pairs, dof=("shift_x", "shift_y"))

        # Rescaling columns by an invertible (all-nonzero) diagonal matrix
        # cannot change the rank - a degenerate direction stays degenerate.
        # (Finite-difference noise keeps the smallest singular values from
        # being bit-exact zero, so use the same RELATIVE threshold as
        # elsewhere in this project - e.g. 03_sensitivity_and_degeneracy.py
        # - rather than numpy.linalg.matrix_rank's near-machine-epsilon
        # default, which is too strict to see them as degenerate here.)
        def _rank(Lambda_):
            S = np.linalg.svd(Lambda_, compute_uv=False)
            return int(np.sum(S > 1e-6 * S[0]))

        rank_before = _rank(Lambda)
        scale = {"shift": 0.02}
        Lambda_scaled = normalize_jacobian(Lambda, specs, scale)
        rank_after = _rank(Lambda_scaled)
        self.assertEqual(rank_before, rank_after)
        self.assertLess(rank_before, len(specs))  # this system is exactly degenerate


if __name__ == "__main__":
    unittest.main()
