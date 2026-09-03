import unittest
import numpy as np

from synim.registration.model import GuideStar, DM, WFS, System
from synim.registration.reconstruction import ParameterSpec, jacobian, apply_alpha
from synim.registration.analysis import (
    svd_of_jacobian, condition_number, describe_mode, describe_modes,
    reconstruction_matrix, covariance_from_noise, monte_carlo_noise_propagation,
    analyze,
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


if __name__ == "__main__":
    unittest.main()
