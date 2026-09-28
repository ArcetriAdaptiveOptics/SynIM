"""
Tests for data type handling and for small bug fixes:
- telescoping sum keeps the data precision and works without a mask
- apply_mask / rebin keep floating point dtypes
- dm2d_to_3d does not modify its input
- MMSE reconstructor built from per-WFS noise variances
"""
import unittest

import numpy as np

import synim
import synim.synim as synim_core
import synim.params_utils as params_utils
from synim.utils import apply_mask, rebin, dm2d_to_3d


def _circular_mask(n, radius):
    yy, xx = np.mgrid[:n, :n]
    return (np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5) < radius).astype(np.float32)


class TestTelescopingSum(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(0)
        self.n, self.nsa = 48, 8
        self.mask = _circular_mask(self.n, self.n / 2)
        self.data = rng.standard_normal((self.n, self.n, 3))

    def test_keeps_float32(self):
        tx, ty = synim_core.compute_telsum_with_extrapolation(
            self.data.astype(np.float32), mask=self.mask, wfs_nsubaps=self.nsa)
        self.assertEqual(tx.dtype, np.float32)
        self.assertEqual(ty.dtype, np.float32)

    def test_keeps_float64(self):
        tx, ty = synim_core.compute_telsum_with_extrapolation(
            self.data, mask=self.mask.astype(np.float64), wfs_nsubaps=self.nsa)
        self.assertEqual(tx.dtype, np.float64)
        self.assertEqual(ty.dtype, np.float64)

    def test_without_mask(self):
        data = self.data.astype(np.float32)
        tx, ty = synim_core.compute_telsum_with_extrapolation(
            data, mask=None, wfs_nsubaps=self.nsa)
        tx_ref, ty_ref = synim_core.compute_telsum_with_extrapolation(
            data, mask=np.ones((self.n, self.n), dtype=np.float32), wfs_nsubaps=self.nsa)
        np.testing.assert_allclose(synim.cpuArray(tx), synim.cpuArray(tx_ref), rtol=0, atol=1e-6)
        np.testing.assert_allclose(synim.cpuArray(ty), synim.cpuArray(ty_ref), rtol=0, atol=1e-6)

    def test_interaction_matrix_dtype(self):
        modes = self.data.astype(np.float32)
        for slope_method in ('derivatives', 'telsum'):
            with self.subTest(slope_method=slope_method):
                im = synim_core.interaction_matrix(
                    8.0, self.mask, modes, np.ones((self.n, self.n), dtype=np.float32),
                    0.0, 0.0, self.nsa, 4.0, (0.0, 0.0), np.inf, 0.0, (0.0, 0.0), 1.0,
                    slope_method=slope_method)
                self.assertEqual(im.dtype, synim.float_dtype)


class TestDoubleInputsSinglePrecision(unittest.TestCase):
    """
    With single precision (as in these tests), float64 inputs given by the
    user must not make the pipeline run in float64 (memory and time).
    """

    def setUp(self):
        n = 48
        self.n = n
        self.pup_mask = _circular_mask(n, n / 2).astype(np.float64)
        self.dm_mask = np.ones((n, n))
        self.dm_array = np.random.default_rng(3).standard_normal((n, n, 4))  # float64
        self.assertEqual(synim.float_dtype, np.float32)

    def test_interaction_matrix(self):
        # record the dtype of every array transformed by rotshiftzoom_array
        dtypes = []
        original = synim_core.rotshiftzoom_array

        def spy(array, *args, **kwargs):
            dtypes.append(array.dtype)
            return original(array, *args, **kwargs)

        synim_core.rotshiftzoom_array = spy
        try:
            for slope_method in ('derivatives', 'telsum'):
                im = synim_core.interaction_matrix(
                    8.0, self.pup_mask, self.dm_array, self.dm_mask, 0.0, 0.0, 8, 4.0,
                    (10.0, 0.0), 90e3, 5.0, (0.0, 0.0), 1.0, slope_method=slope_method)
                self.assertEqual(im.dtype, np.float32)
        finally:
            synim_core.rotshiftzoom_array = original
        self.assertTrue(dtypes)
        self.assertTrue(all(d == np.float32 for d in dtypes), dtypes)

    def test_interaction_matrices_multi_wfs(self):
        configs = [dict(name='a', nsubaps=8, fov_arcsec=4.0, gs_pol_coo=(10.0, 0.0),
                        gs_height=90e3, rotation=5.0)]
        im_dict, _ = synim_core.interaction_matrices_multi_wfs(
            8.0, self.pup_mask, self.dm_array, self.dm_mask, 0.0, 0.0, configs)
        self.assertEqual(im_dict['a'].dtype, np.float32)

    def test_projection_matrix(self):
        import synim.synpm as synpm
        idx = np.where(self.pup_mask > 0.5)
        base_inv = np.linalg.pinv(self.dm_array[idx[0], idx[1], :3])   # float64
        projection = synpm.projection_matrix(
            8.0, self.pup_mask, self.dm_array, self.dm_mask, base_inv,
            0.0, 0.0, 0.0, (0, 0), (1, 1), (0, 0), np.inf)
        self.assertEqual(projection.dtype, np.float32)

    def test_subaperture_illumination(self):
        illumination = synim_core.compute_subaperture_illumination(self.pup_mask, 8)
        self.assertEqual(illumination.dtype, np.float32)


class TestDtypePreservation(unittest.TestCase):

    def test_apply_mask(self):
        mask = _circular_mask(16, 6)
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=dtype.__name__):
                out = apply_mask(np.ones((16, 16, 2), dtype=dtype), mask)
                self.assertEqual(out.dtype, dtype)
        out = apply_mask(np.ones((16, 16), dtype=bool), mask)
        self.assertEqual(out.dtype, synim.float_dtype)

    def test_rebin_compression(self):
        for dtype in (np.float32, np.float64):
            for method in ('sum', 'average', 'nanmean'):
                with self.subTest(dtype=dtype.__name__, method=method):
                    out = rebin(np.ones((16, 16), dtype=dtype), (4, 4), method=method)
                    self.assertEqual(out.dtype, dtype)
        out = rebin(np.ones((16, 16), dtype=bool), (4, 4), method='sum')
        self.assertEqual(out.dtype, synim.float_dtype)
        np.testing.assert_array_equal(synim.cpuArray(out), 16)

    def test_rebin_non_integer_factor(self):
        out = rebin(np.ones((15, 15), dtype=np.float64), (4, 4), method='average')
        self.assertEqual(out.dtype, np.float64)

    def test_rebin_expansion(self):
        for dtype in (np.float32, np.float64, bool):
            with self.subTest(dtype=np.dtype(dtype).name):
                out = rebin(np.ones((4, 4), dtype=dtype), (8, 8))
                self.assertEqual(out.dtype, dtype)


class TestDm2dTo3d(unittest.TestCase):

    def test_input_not_modified(self):
        mask = _circular_mask(20, 9)
        n_valid = int(mask.sum())
        dm_2d = np.random.default_rng(1).standard_normal((4, n_valid)).astype(np.float32)
        original = dm_2d.copy()
        out = dm2d_to_3d(dm_2d, mask, normalize=True, xp_local=np, float_dtype_local=np.float32)
        np.testing.assert_array_equal(dm_2d, original)
        # normalized modes: unit RMS before piston removal, zero mean after
        values = out[mask > 0]
        np.testing.assert_allclose(values.mean(axis=0), 0, atol=1e-6)
        self.assertTrue((out[mask == 0] == 0).all())

    def test_without_normalization(self):
        mask = _circular_mask(20, 9)
        dm_2d = np.arange(3 * int(mask.sum()), dtype=np.float32).reshape(3, -1)
        out = dm2d_to_3d(dm_2d, mask, normalize=False, xp_local=np, float_dtype_local=np.float32)
        np.testing.assert_array_equal(out[mask > 0], dm_2d.T)


class TestMmseNoiseVariance(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(2)
        self.n_modes, self.n_slopes_per_wfs = 10, 40
        self.im = rng.standard_normal((2 * self.n_slopes_per_wfs, self.n_modes))
        self.c_atm = np.diag(1.0 / np.arange(1, self.n_modes + 1))

    def test_noise_variance_per_wfs(self):
        variances = [1e-2, 4e-2]
        c_noise = np.diag(np.repeat(variances, self.n_slopes_per_wfs))
        expected = params_utils.compute_mmse_reconstructor(
            self.im, self.c_atm, C_noise=c_noise, xp=np, dtype=np.float64)
        rec = params_utils.compute_mmse_reconstructor(
            self.im, self.c_atm, noise_variance=variances, xp=np, dtype=np.float64)
        np.testing.assert_allclose(rec, expected, rtol=1e-10, atol=1e-12)

    def test_default_noise(self):
        expected = params_utils.compute_mmse_reconstructor(
            self.im, self.c_atm, C_noise=np.eye(self.im.shape[0]), xp=np, dtype=np.float64)
        rec = params_utils.compute_mmse_reconstructor(self.im, self.c_atm, xp=np,
                                                      dtype=np.float64)
        np.testing.assert_allclose(rec, expected, rtol=1e-10, atol=1e-12)


class TestSubapertureIlluminationVerbose(unittest.TestCase):

    def test_verbose(self):
        illumination = synim_core.compute_subaperture_illumination(
            _circular_mask(48, 24), 8, verbose=True)
        self.assertEqual(illumination.shape, (64,))


if __name__ == '__main__':
    unittest.main()
