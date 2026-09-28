"""
Subapertures with a non-integer number of pixels (e.g. 480 pixels and 68
subapertures):
- rebin with exact area weighting (rebin_matrix), checked against an exact
  integer upsampling
- slopes (derivatives and telescoping sum) without NaN and without loss at
  the pupil edge, and close to a finely sampled (integer ratio) reference
- pupil and DM masks not modified, NaN warning on the interaction matrix
"""
import unittest
import warnings

import numpy as np

import synim
import synim.synim as synim_core
from synim.utils import rebin, rebin_matrix


def _to_xp(array):
    return synim.xp.asarray(array)


def _pupil(n, obstruction=0.25):
    c = (n - 1) / 2
    yy, xx = np.mgrid[:n, :n]
    r = np.hypot(xx - c, yy - c) / (n / 2)
    return ((r <= 1) & (r >= obstruction)).astype(np.float32)


def _modes(n):
    """Tip, tilt and smooth modes on [-1, 1] (pixel centres)."""
    c = (n - 1) / 2
    yy, xx = np.mgrid[:n, :n]
    x = (xx - c) / (n / 2)
    y = (yy - c) / (n / 2)
    modes = [x, y, x * y, x ** 2 - y ** 2, x ** 3 - 2 * x * y ** 2,
             np.sin(3 * x) * np.cos(2 * y)]
    return np.stack(modes, axis=2).astype(np.float32)


def _im(n, nsa, slope_method, pup=None, modes=None, dm_mask=None):
    pup = _pupil(n) if pup is None else pup
    modes = _modes(n) if modes is None else modes
    dm_mask = np.ones((n, n), dtype=np.float32) if dm_mask is None else dm_mask
    return synim.cpuArray(synim_core.interaction_matrix(
        8.0, pup, modes, dm_mask, 0.0, 0.0, nsa, 4.0,
        (0.0, 0.0), np.inf, 0.0, (0.0, 0.0), 1.0, specula_convention=False,
        slope_method=slope_method))


def _upsampled(array, k):
    """Each pixel replaced by k x k equal pixels (exact for area sums)."""
    return np.repeat(np.repeat(array, k, axis=0), k, axis=1)


class TestRebinMatrix(unittest.TestCase):

    def test_sums(self):
        for m, M in ((480, 68), (100, 12), (5, 2), (96, 12)):
            with self.subTest(m=m, M=M):
                matrix = synim.cpuArray(rebin_matrix(m, M, dtype=np.float64))
                self.assertEqual(matrix.shape, (M, m))
                np.testing.assert_allclose(matrix.sum(axis=0), 1.0, rtol=0, atol=1e-12)
                np.testing.assert_allclose(matrix.sum(axis=1), m / M, rtol=1e-12)
                self.assertTrue((matrix >= 0).all())

    def test_integer_ratio(self):
        matrix = synim.cpuArray(rebin_matrix(12, 3, dtype=np.float64))
        np.testing.assert_array_equal(matrix, np.kron(np.eye(3), np.ones((1, 4))))

    def test_boundaries(self):
        # 5 pixels in 2 bins: pixel 2 is split in half
        matrix = synim.cpuArray(rebin_matrix(5, 2, dtype=np.float64))
        np.testing.assert_array_equal(matrix, [[1, 1, 0.5, 0, 0], [0, 0, 0.5, 1, 1]])


class TestRebinFractional(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(0)
        # 30 pixels in 4 bins (7.5 pixels each): exact with an upsampling by 2
        self.m, self.M, self.k = 30, 4, 2
        self.data = rng.standard_normal((self.m, self.m, 3))

    def _reference(self, data, method):
        fine = _upsampled(data, self.k)
        return synim.cpuArray(rebin(fine, (self.M, self.M), method=method)) \
            / (self.k ** 2 if method == 'sum' else 1)

    def test_sum_and_average(self):
        for method in ('sum', 'average'):
            with self.subTest(method=method):
                out = synim.cpuArray(rebin(_to_xp(self.data), (self.M, self.M), method=method))
                np.testing.assert_allclose(out, self._reference(self.data, method),
                                           rtol=1e-12, atol=1e-12)

    def test_nanmean(self):
        data = self.data.copy()
        data[:7, :] = np.nan              # first bin row: only half of pixel 7 valid
        data[20:, 23:, 1] = np.nan        # NaN pattern different for each slice
        data[:, 27:, 2] = np.nan
        out = synim.cpuArray(rebin(_to_xp(data), (self.M, self.M), method='nanmean'))
        np.testing.assert_allclose(out, self._reference(data, 'nanmean'), rtol=1e-12,
                                   atol=1e-12)
        self.assertFalse(np.isnan(out).any())

    def test_nanmean_empty_bin(self):
        # First bin: pixels 0-6 and half of pixel 7
        data = np.ones((self.m, self.m))
        data[:7, :7] = np.nan
        out = synim.cpuArray(rebin(_to_xp(data), (self.M, self.M), method='nanmean'))
        self.assertFalse(np.isnan(out).any())
        np.testing.assert_allclose(out, 1.0, rtol=1e-12)
        data[:8, :8] = np.nan
        out = synim.cpuArray(rebin(_to_xp(data), (self.M, self.M), method='nanmean'))
        self.assertTrue(np.isnan(out[0, 0]))
        self.assertEqual(int(np.isnan(out).sum()), 1)

    def test_constant(self):
        # No loss at the edges and no shift: a constant stays constant
        for m, M in ((480, 68), (100, 12), (5, 2)):
            with self.subTest(m=m, M=M):
                out = synim.cpuArray(rebin(_to_xp(np.ones((m, m), dtype=np.float32)),
                                           (M, M), method='average'))
                np.testing.assert_allclose(out, 1.0, rtol=1e-5)

    def test_input_not_modified_and_dtype(self):
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=dtype.__name__):
                data = self.data.astype(dtype)
                data[0, 0, 0] = np.nan
                data_xp = _to_xp(data)
                out = rebin(data_xp, (self.M, self.M), method='nanmean')
                self.assertEqual(out.dtype, dtype)
                np.testing.assert_array_equal(synim.cpuArray(data_xp), data)

    def test_3d_equals_2d(self):
        out_3d = synim.cpuArray(rebin(_to_xp(self.data), (self.M, self.M), method='sum'))
        for i in range(self.data.shape[2]):
            out_2d = synim.cpuArray(rebin(_to_xp(self.data[:, :, i]), (self.M, self.M),
                                          method='sum'))
            np.testing.assert_allclose(out_3d[:, :, i], out_2d, rtol=1e-12, atol=1e-12)


class TestFractionalSlopes(unittest.TestCase):
    """Interaction matrices with 90 pixels and 12 subapertures (7.5 px/subap)."""

    n, nsa = 90, 12

    def test_tip_tilt_constant(self):
        # A pure tip/tilt gives the same slope in every subaperture with
        # pupil pixels, equal to the one of an integer ratio (96 pixels)
        nsa2 = self.nsa ** 2
        for slope_method in ('derivatives', 'telsum'):
            with self.subTest(slope_method=slope_method):
                im = _im(self.n, self.nsa, slope_method)
                im_ref = _im(96, self.nsa, slope_method)
                self.assertFalse(np.isnan(im).any())
                expected = np.median(im_ref[:nsa2, 0][im_ref[:nsa2, 0] != 0])
                for values in (im[:nsa2, 0], im[nsa2:, 1]):
                    nonzero = values[values != 0]
                    np.testing.assert_allclose(nonzero, expected, rtol=1e-4)
                    # cross terms are zero
                self.assertLess(np.abs(im[nsa2:, 0]).max(), 1e-4 * abs(expected))
                # every subaperture with pupil pixels has a slope (derivatives);
                # the telescoping sum needs at least one pair of adjacent pixels
                illuminated = rebin(_pupil(self.n), (self.nsa, self.nsa), method='sum')
                n_illuminated = int((synim.cpuArray(illuminated) > 0).sum())
                n_slopes = int((im[:nsa2, 0] != 0).sum())
                if slope_method == 'derivatives':
                    self.assertEqual(n_slopes, n_illuminated)
                else:
                    self.assertGreaterEqual(n_slopes, n_illuminated - 8)

    def test_vs_fine_sampling(self):
        # Reference: same modes sampled with 30 pixels per subaperture
        # (integer ratio). Maximum relative rms difference over the modes,
        # fully illuminated / partially illuminated (> 0.5) subapertures:
        #   derivatives 1.8e-3 / 2.6e-3, telsum 2.2e-3 / 2.6e-3
        # with the previous upscaling and nanmean (48 NaN with derivatives):
        #   derivatives 3.3e-3 / 1.4e-1, telsum 3.8e-3 / 2.0e-2
        nsa2 = self.nsa ** 2
        illumination = synim.cpuArray(rebin(_pupil(360), (self.nsa, self.nsa),
                                            method='sum')).ravel()
        illumination = np.tile(illumination / illumination.max(), 2)
        full = illumination > 0.999
        partial = (illumination > 0.5) & ~full
        for slope_method, tol_full, tol_partial in (('derivatives', 5e-3, 1e-2),
                                                     ('telsum', 5e-3, 1e-2)):
            with self.subTest(slope_method=slope_method):
                im = _im(self.n, self.nsa, slope_method)
                ref = _im(360, self.nsa, slope_method)
                self.assertEqual(im.shape, (2 * nsa2, 6))
                norm = np.sqrt(np.mean(ref[full] ** 2, axis=0))
                err_full = np.sqrt(np.mean((im - ref)[full] ** 2, axis=0)) / norm
                err_partial = np.sqrt(np.mean((im - ref)[partial] ** 2, axis=0)) / norm
                self.assertLess(err_full.max(), tol_full, err_full)
                self.assertLess(err_partial.max(), tol_partial, err_partial)

    def test_telsum_integer_ratio_equivalence(self):
        # The fractional weights reduce to the block computation
        rng = np.random.default_rng(1)
        data = rng.standard_normal((96, 96, 3)).astype(np.float64)
        mask = _pupil(96).astype(np.float64)
        block = synim_core.compute_telsum_with_extrapolation(
            _to_xp(data.copy()), mask=_to_xp(mask), wfs_nsubaps=12)
        fractional = synim_core._telsum_fractional(_to_xp(data), _to_xp(mask), 12)
        for a, b in zip(block, fractional):
            np.testing.assert_allclose(synim.cpuArray(b), synim.cpuArray(a),
                                       rtol=1e-10, atol=1e-12)

    def test_telsum_needs_two_pixels(self):
        with self.assertRaises(ValueError):
            synim_core.compute_telsum_with_extrapolation(
                _to_xp(np.zeros((20, 20), dtype=np.float32)), wfs_nsubaps=12)


class TestTelsumOverPupil(unittest.TestCase):
    """The telescoping sum is averaged over the pupil, as the derivatives."""

    def test_independent_of_dm_mask_outside_pupil(self):
        # Only pixel pairs inside the pupil are used: a DM mask larger than
        # the pupil (here the whole array or the pupil without obstruction,
        # dilated by 3 pixels) gives the same slopes
        for n in (96, 90):
            with self.subTest(n=n):
                c = (n - 1) / 2
                yy, xx = np.mgrid[:n, :n]
                r = np.hypot(xx - c, yy - c)
                dilated = (r <= n / 2 + 3).astype(np.float32)
                im_ones = _im(n, 12, 'telsum')
                im_dilated = _im(n, 12, 'telsum', dm_mask=dilated)
                np.testing.assert_allclose(im_dilated, im_ones, rtol=0,
                                           atol=1e-5 * np.abs(im_ones).max())

    def test_close_to_derivatives(self):
        # Smooth modes: telescoping sum and derivatives agree also in the
        # partially illuminated subapertures (illumination > 0.5). Measured
        # relative rms difference: 0.005 (0.06 - 0.15 when the telescoping
        # sum was averaged over the DM mask)
        for n in (96, 90):
            with self.subTest(n=n):
                illumination = synim.cpuArray(rebin(_pupil(n), (12, 12), method='sum')).ravel()
                illumination = np.tile(illumination / illumination.max(), 2)
                full = illumination > 0.999
                partial = (illumination > 0.5) & ~full
                der, tel = _im(n, 12, 'derivatives'), _im(n, 12, 'telsum')
                norm = np.sqrt(np.mean(der[full] ** 2, axis=0))
                for mask, tol in ((full, 0.01), (partial, 0.02)):
                    err = np.sqrt(np.mean((tel - der)[mask] ** 2, axis=0)) / norm
                    self.assertLess(err.max(), tol, err)


class TestMasksAndNaN(unittest.TestCase):

    def test_masks_not_modified(self):
        n, nsa = 90, 12
        pup = _pupil(n)
        pup[0, 0] = np.nan
        dm_mask = np.ones((n, n), dtype=np.float32)
        dm_mask[1, 1] = np.nan
        pup_xp, dm_xp = _to_xp(pup), _to_xp(dm_mask)
        rng = np.random.default_rng(2)
        dx = _to_xp(rng.standard_normal((n, n, 2)).astype(np.float32))
        dy = _to_xp(rng.standard_normal((n, n, 2)).astype(np.float32))
        im = synim_core.apply_wfs_transformations_combined(
            dx, dy, pup_xp, dm_xp, nsa, 4.0, 8.0, slope_method='derivatives')
        self.assertFalse(np.isnan(synim.cpuArray(im)).any())
        np.testing.assert_array_equal(synim.cpuArray(pup_xp), pup)
        np.testing.assert_array_equal(synim.cpuArray(dm_xp), dm_mask)
        tx, ty = synim_core.compute_telsum_with_extrapolation(
            dx, mask=_to_xp(_pupil(n)), wfs_nsubaps=nsa)
        synim_core.apply_wfs_transformations_combined(
            tx, ty, pup_xp, dm_xp, nsa, 4.0, 8.0, slope_method='telsum')
        np.testing.assert_array_equal(synim.cpuArray(pup_xp), pup)
        np.testing.assert_array_equal(synim.cpuArray(dm_xp), dm_mask)

    def test_warning(self):
        im = np.ones((6, 2), dtype=np.float32)
        im[1, 0] = im[4, 1] = im[4, 0] = np.nan
        with self.assertWarnsRegex(RuntimeWarning, '3 NaN values in 2 of 6 slopes'):
            synim_core._warn_if_nan(_to_xp(im))

    def test_no_warning(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            _im(90, 12, 'derivatives')
            configs = [dict(name='a', nsubaps=12, fov_arcsec=4.0, gs_pol_coo=(0.0, 0.0),
                            gs_height=np.inf)]
            synim_core.interaction_matrices_multi_wfs(
                8.0, _pupil(90), _modes(90), np.ones((90, 90), dtype=np.float32), 0.0, 0.0,
                configs)
        self.assertFalse([w for w in caught if 'NaN values' in str(w.message)])


if __name__ == '__main__':
    unittest.main()
