"""
interaction_matrix computes the interaction matrix as a sparse linear
operator applied to the DM array (synim.synim._slope_operators,
_transformation_operator, _wfs_operator). The results must be the same, up
to rounding errors, as the explicit computation on the transformed DM array
(apply_dm_transformations_combined followed by
apply_wfs_transformations_combined), for both slope methods, integer and
non-integer numbers of pixels per subaperture, both conventions, and a
selection of the valid subapertures.
"""
import unittest

import numpy as np

import synim
import synim.synim as synim_core
import synim.utils as synim_utils

N_PIX = 96
N_MODES = 10
# Tolerance relative to max(|IM|). On GPU in single precision the explicit
# computation is less accurate: cupyx affine_transform computes the
# interpolation coordinates in single precision.
if synim.float_dtype != np.float32:
    RTOL = 1e-10
elif synim.xp is np:
    RTOL = 1e-5
else:
    RTOL = 1e-4


class TestSparseOperator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        n = N_PIX
        yy, xx = np.mgrid[:n, :n]
        r = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5)
        cls.pup_mask = ((r < n / 2 - 3) & (r > n / 8)).astype(np.float32)
        cls.dm_mask = (r < n / 2 - 1).astype(np.float32)
        cls.dm_array = (np.random.default_rng(0).standard_normal((n, n, N_MODES))
                        * cls.dm_mask[:, :, np.newaxis]).astype(np.float32)

    def _cases(self):
        rng = np.random.default_rng(1)
        for k in range(16):
            n_subaps = (8, 12, 13, 7)[k % 4]
            geometry = dict(
                dm_height=(0.0, 8000.0)[k % 2], dm_rotation=rng.uniform(-20, 20),
                gs_pol_coo=(rng.uniform(0, 40), rng.uniform(0, 360)),
                gs_height=(np.inf, 90e3)[(k // 2) % 2], wfs_rotation=rng.uniform(-30, 30),
                wfs_translation=tuple(rng.normal(0, 1.5, 2)),
                wfs_mag_global=rng.uniform(0.97, 1.03),
                wfs_anamorphosis_90=rng.uniform(0.98, 1.02),
                wfs_anamorphosis_45=rng.uniform(0.98, 1.02))
            idx_valid_sa = None
            if k % 4 == 1:
                idx_valid_sa = np.argwhere(rng.random((n_subaps, n_subaps)) > 0.3)
            yield dict(slope_method=('derivatives', 'telsum')[k % 2], wfs_nsubaps=n_subaps,
                       specula_convention=k % 3 != 0, idx_valid_sa=idx_valid_sa, **geometry)

    def _explicit(self, case):
        case = dict(case)
        idx_valid_sa = case.pop('idx_valid_sa')
        if idx_valid_sa is not None:
            idx_valid_sa = synim.xp.asarray(idx_valid_sa)
        _, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
            synim_core.apply_dm_transformations_combined(
                8.0, self.pup_mask, self.dm_array, self.dm_mask, **case)
        return synim.cpuArray(synim_core.apply_wfs_transformations_combined(
            derivatives_x, derivatives_y, trans_pup_mask, trans_dm_mask,
            case['wfs_nsubaps'], 3.0, 8.0, idx_valid_sa=idx_valid_sa,
            slope_method=case['slope_method'],
            specula_convention=case['specula_convention'], in_place=True))

    def test_same_as_explicit_computation(self):
        for case in self._cases():
            with self.subTest(**{k: v for k, v in case.items() if k != 'idx_valid_sa'}):
                im = synim.cpuArray(synim_core.interaction_matrix(
                    8.0, self.pup_mask, self.dm_array, self.dm_mask, wfs_fov_arcsec=3.0,
                    **case))
                expected = self._explicit(case)
                self.assertEqual(im.shape, expected.shape)
                np.testing.assert_allclose(im, expected, rtol=0,
                                           atol=RTOL * np.abs(expected).max())

    def test_interpolation_matrix(self):
        """affine_interpolation_matrix gives the same result as rotshiftzoom_array."""
        rng = np.random.default_rng(2)
        cube = rng.standard_normal((70, 64, 3))
        transform = dict(dm_translation=(2.0, -1.5), dm_rotation=12.0,
                         dm_magnification=(1.05, 0.98), wfs_translation=(0.4, 0.7),
                         wfs_rotation=-33.0, wfs_magnification=(1.01, 0.99),
                         wfs_anamorphosis_45=1.02)
        for transposed in (False, True):
            with self.subTest(transposed=transposed):
                source = np.transpose(cube, (1, 0, 2)) if transposed else cube
                expected = synim.cpuArray(synim_utils.rotshiftzoom_array(
                    synim.xp.asarray(source, dtype=synim.float_dtype),
                    output_size=(60, 60), **transform))
                matrix, offset = synim_utils.rotshiftzoom_matrix(
                    source.shape[:2], (60, 60), **transform)
                operator = synim_utils.affine_interpolation_matrix(
                    source.shape[:2], (60, 60), matrix, offset, transposed=transposed)
                result = (operator @ cube.reshape(-1, 3)).reshape(60, 60, 3)
                np.testing.assert_allclose(result, expected, rtol=0,
                                           atol=RTOL * np.abs(expected).max())


if __name__ == '__main__':
    unittest.main()
