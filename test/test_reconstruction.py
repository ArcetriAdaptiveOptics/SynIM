import unittest
import numpy as np

from synim.registration.model import GuideStar, DM, WFS, System
from synim.registration.reconstruction import (
    ParameterSpec, get_alpha, apply_alpha, local_params_vector,
    jacobian, gauss_newton_invert, LOCAL_DOF,
)


def _build_3wfs_2dm_system(pixel_pitch=0.2):
    """Analogous to the illustrative 3-shift example of Eq. 3 in the SPIE
    paper: 3 WFSs, 2 DMs (one at the pupil, one post-focal)."""
    system = System(pixel_pitch=pixel_pitch)
    gs_positions = [(0.0, 0.0), (17.5, 0.0), (-12.0, 9.0)]
    for i, pos in enumerate(gs_positions):
        gs = GuideStar(name=f"gs{i}", position=pos,
                        height=90000.0 if pos != (0.0, 0.0) else np.inf)
        system.add_wfs(WFS(name=f"wfs{i}", guide_star=gs))
    system.add_dm(DM(name="dm0", height=0.0))
    system.add_dm(DM(name="dm1", height=6000.0))
    return system


class TestJacobianAndVector(unittest.TestCase):
    def test_get_and_apply_alpha_roundtrip(self):
        system = _build_3wfs_2dm_system()
        specs = [
            ParameterSpec("dm", "dm0", "shift", 0),
            ParameterSpec("wfs", "wfs0", "rotation"),
            ParameterSpec("gs", "gs1", "position_shift", 1),
        ]
        alpha = np.array([0.3, 1.5, -2.0])
        new_system = apply_alpha(system, specs, alpha)
        np.testing.assert_allclose(get_alpha(new_system, specs), alpha)
        # original system must be untouched (deep copy)
        np.testing.assert_allclose(get_alpha(system, specs), [0.0, 0.0, 0.0])

    def test_local_params_vector_shape(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        D = local_params_vector(system, pairs, dof=LOCAL_DOF)
        self.assertEqual(D.shape, (len(pairs) * len(LOCAL_DOF),))


class TestLinearShiftReconstruction(unittest.TestCase):
    """
    Minimal case analogous to the paper's Eq. 3 example: 3 pure shifts,
    alpha1 on DM0 (pupil), alpha2 on DM1 (post-focal), alpha3 on WFS0 -
    each affecting a distinct, non-overlapping block of local shifts, so
    the reconstruction is exactly linear (no Gauss-Newton iteration
    needed) and perfectly conditioned.
    """

    def test_exact_recovery_single_iteration(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        specs = [
            ParameterSpec("dm", "dm0", "shift", 0),
            ParameterSpec("dm", "dm1", "shift", 0),
            ParameterSpec("wfs", "wfs0", "shift", 0),
        ]
        true_alpha = np.array([0.35, -0.6, 0.15])
        true_system = apply_alpha(system, specs, true_alpha)
        D_meas = local_params_vector(true_system, pairs, dof=("shift_x", "shift_y"))

        alpha_hat, system_hat, history = gauss_newton_invert(
            system, specs, pairs, D_meas, dof=("shift_x", "shift_y"), n_iter=1
        )

        np.testing.assert_allclose(alpha_hat, true_alpha, atol=1e-6)
        self.assertGreater(history[0], 1e-3)  # residual before the (only) update was nonzero...
        # ... but after applying alpha_hat the residual should now vanish:
        final_residual = D_meas - local_params_vector(system_hat, pairs, dof=("shift_x", "shift_y"))
        np.testing.assert_allclose(final_residual, 0.0, atol=1e-6)

    def test_jacobian_is_well_conditioned_for_this_example(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        specs = [
            ParameterSpec("dm", "dm0", "shift", 0),
            ParameterSpec("dm", "dm1", "shift", 0),
            ParameterSpec("wfs", "wfs0", "shift", 0),
        ]
        Lambda = jacobian(system, specs, pairs, dof=("shift_x", "shift_y"))
        singular_values = np.linalg.svd(Lambda, compute_uv=False)
        condition_number = singular_values[0] / singular_values[-1]
        self.assertLess(condition_number, 100.0,
                         "Non-overlapping shifts should be well conditioned.")


class TestNonlinearReconstructionWithRotation(unittest.TestCase):
    """Rotation + shift together requires a few Gauss-Newton iterations
    (the coupling the paper notes in Sec. 4/5)."""

    def test_converges_with_a_few_iterations(self):
        system = _build_3wfs_2dm_system()
        pairs = system.pairs()
        specs = [
            ParameterSpec("dm", "dm1", "shift", 0),
            ParameterSpec("dm", "dm1", "shift", 1),
            ParameterSpec("dm", "dm1", "rotation"),
            ParameterSpec("wfs", "wfs1", "rotation"),
        ]
        true_alpha = np.array([0.2, -0.15, 1.2, -0.8])
        true_system = apply_alpha(system, specs, true_alpha)
        D_meas = local_params_vector(true_system, pairs, dof=LOCAL_DOF)

        alpha_hat, system_hat, history = gauss_newton_invert(
            system, specs, pairs, D_meas, dof=LOCAL_DOF, n_iter=8
        )

        np.testing.assert_allclose(alpha_hat, true_alpha, atol=1e-5)
        # residual should decrease monotonically (or at least end near zero)
        self.assertLess(history[-1], 1e-6)


if __name__ == "__main__":
    unittest.main()
