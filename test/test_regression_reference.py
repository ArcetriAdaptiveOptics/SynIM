"""
Regression test against test/data/regression_reference.npz.

The reference was generated on commit 12ea3de with
test/data/generate_regression_reference.py (see that script for how the
cases that were wrong on 12ea3de were handled). All the results are computed
here through the public API.
"""
import os
import sys
import unittest

import numpy as np

import synim
import synim.synim as synim_core
import synim.synpm as synpm
import synim.utils as synim_utils
import synim.params_utils as params_utils

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
sys.path.insert(0, DATA_DIR)
import generate_regression_reference as gen  # noqa: E402

# Tolerances relative to max(|reference|). float32 results can differ by a
# few units in the last place between numpy versions and platforms.
RTOL_FLOAT32 = 1e-5
RTOL_FLOAT64 = 1e-8


class TestRegressionReference(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ref = np.load(os.path.join(DATA_DIR, 'regression_reference.npz'))
        cls.pup_mask, cls.dm_mask, cls.dm_array = gen.make_inputs()

    def assert_matches(self, value, key, rtol=RTOL_FLOAT32):
        expected = self.ref[key]
        value = synim.cpuArray(value)
        self.assertEqual(value.shape, expected.shape, msg=key)
        np.testing.assert_array_equal(np.isnan(value), np.isnan(expected), err_msg=key)
        scale = np.nanmax(np.abs(expected))
        self.assertGreater(scale, 0, msg=key)
        max_rel = np.nanmax(np.abs(value.astype(np.float64) - expected)) / scale
        self.assertLessEqual(max_rel, rtol, msg=f'{key}: max relative difference {max_rel:.2e}')

    def test_inputs_unchanged(self):
        np.testing.assert_allclose(
            gen.inputs_checksum(self.pup_mask, self.dm_mask, self.dm_array),
            self.ref['inputs_checksum'], rtol=1e-10)

    def test_interaction_matrix(self):
        for slope_method in ('derivatives', 'telsum'):
            for name, geometry in gen.GEOMETRIES.items():
                with self.subTest(slope_method=slope_method, geometry=name):
                    im = synim_core.interaction_matrix(
                        gen.PUP_DIAM_M, self.pup_mask, self.dm_array, self.dm_mask,
                        wfs_nsubaps=gen.N_SUBAPS, wfs_fov_arcsec=gen.FOV_ARCSEC,
                        slope_method=slope_method, **dict(gen.DEFAULTS, **geometry))
                    self.assert_matches(im, f'im_{slope_method}_{name}')

    def test_interaction_matrix_valid_subapertures(self):
        im = synim_core.interaction_matrix(
            gen.PUP_DIAM_M, self.pup_mask, self.dm_array, self.dm_mask,
            wfs_nsubaps=gen.N_SUBAPS, wfs_fov_arcsec=gen.FOV_ARCSEC,
            idx_valid_sa=self.ref['idx_valid_sa'],
            **dict(gen.DEFAULTS, **gen.GEOMETRIES['dm_and_wfs']))
        self.assert_matches(im, 'im_derivatives_dm_and_wfs_idx_valid')

    def test_interaction_matrices_multi_wfs(self):
        configs = [
            dict(name='a', nsubaps=gen.N_SUBAPS, fov_arcsec=gen.FOV_ARCSEC,
                 gs_pol_coo=(30.0, 0.0), gs_height=90e3),
            dict(name='b', nsubaps=gen.N_SUBAPS, fov_arcsec=gen.FOV_ARCSEC,
                 gs_pol_coo=(30.0, 120.0), gs_height=90e3, rotation=10.0,
                 translation=(0.3, -0.1)),
        ]
        im_dict, _ = synim_core.interaction_matrices_multi_wfs(
            gen.PUP_DIAM_M, self.pup_mask, self.dm_array, self.dm_mask, 8000.0, 1.0, configs)
        for key in ('a', 'b'):
            with self.subTest(wfs=key):
                self.assert_matches(im_dict[key], f'multi_{key}')

    def test_subaperture_illumination(self):
        illumination = synim_core.compute_subaperture_illumination(self.pup_mask, gen.N_SUBAPS)
        self.assert_matches(illumination, 'illumination')

    def test_projection_matrix(self):
        idx = np.where(self.pup_mask > 0.5)
        base = self.dm_array[idx[0], idx[1], :10].astype(np.float64)
        base_inv = np.linalg.pinv(base).astype(np.float32)
        projection = synpm.projection_matrix(
            gen.PUP_DIAM_M, self.pup_mask, self.dm_array, self.dm_mask, base_inv,
            5000.0, 1.0, 0.0, (0, 0), (1, 1), (20, 10), np.inf)
        self.assert_matches(projection, 'projection')

    def test_dm_conversions(self):
        dm_2d = synim_utils.dm3d_to_2d(self.dm_array[:, :, :3].copy(), self.dm_mask,
                                       xp_local=np, float_dtype_local=np.float32)
        self.assert_matches(dm_2d, 'dm3d_to_2d')
        dm_3d = synim.cpuArray(synim_utils.dm2d_to_3d(
            dm_2d.copy(), self.dm_mask, xp_local=np, float_dtype_local=np.float32))
        self.assert_matches(dm_3d[self.dm_mask > 0], 'dm2d_to_3d_on_mask')
        self.assertEqual(float(np.abs(dm_3d[self.dm_mask == 0]).sum()),
                         float(self.ref['dm2d_to_3d_outside_mask_abs_sum']))

    def test_reconstructors(self):
        im = self.ref['im_derivatives_none'].astype(np.float64)
        c_atm = np.diag(1.0 / np.arange(1, gen.N_MODES + 1) ** 1.5)
        c_noise = np.eye(im.shape[0]) * 1e-2
        mmse = params_utils.compute_mmse_reconstructor(
            im, c_atm, C_noise=c_noise, xp=np, dtype=np.float64)
        self.assert_matches(mmse, 'reconstructor_mmse', rtol=RTOL_FLOAT64)
        pinv = params_utils.compute_pseudoinverse_reconstructor(im, xp=np, dtype=np.float64)
        self.assert_matches(pinv, 'reconstructor_pinv', rtol=RTOL_FLOAT64)


if __name__ == '__main__':
    unittest.main()
