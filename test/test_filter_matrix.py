"""
Filter matrix (filtmat) generated from an interaction matrix:
- the filtered modes (tip/tilt) are removed from the slopes, the others are left unchanged
- the filtmat file written by generate_filter_matrix_from_intmat_file can be reloaded and has the same property
- the unsupported w_vec argument is rejected
"""
import os
import tempfile
import unittest

import numpy as np

from specula.data_objects.intmat import Intmat

from synim.params_utils import (
    generate_filter_matrix_from_intmat,
    generate_filter_matrix_from_intmat_file,
    load_filter_matrix,
)


class TestFilterMatrix(unittest.TestCase):

    N_SUB = 50
    N_MODES_TOTAL = 30
    N_MODES = 20
    N_FILT = 2
    # SPECULA may run in single precision (see test/__init__.py)
    RTOL = 1e-4

    def setUp(self):
        rng = np.random.default_rng(0)
        n_slopes = 2 * self.N_SUB
        im = rng.standard_normal((n_slopes, self.N_MODES_TOTAL))
        # Exact tip (x slopes only) and tilt (y slopes only), so that
        # the median smoothing of tip/tilt is a no-op
        im[:, 0] = 0.0
        im[:self.N_SUB, 0] = 0.7
        im[:, 1] = 0.0
        im[self.N_SUB:, 1] = 0.7
        self.im = im
        self.n_slopes = n_slopes

    def _check_filtering(self, im_sub, rec_sub):
        self.assertEqual(im_sub.shape, (self.n_slopes, self.N_FILT))
        self.assertEqual(rec_sub.shape, (self.N_FILT, self.n_slopes))
        im_sub = np.asarray(im_sub, dtype=np.float64)
        rec_sub = np.asarray(rec_sub, dtype=np.float64)

        for k in range(self.N_MODES):
            s = self.im[:, k]
            s_f = s - im_sub @ (rec_sub @ s)
            atol = self.RTOL * np.linalg.norm(s)
            if k < self.N_FILT:
                # filtered modes are removed
                np.testing.assert_allclose(s_f, 0.0, atol=atol)
            else:
                # the other modes are untouched
                np.testing.assert_allclose(s_f, s, atol=atol)

        # rec_sub picks the amplitudes of the filtered modes: [I | 0]
        expected = np.zeros((self.N_FILT, self.N_MODES))
        expected[:, :self.N_FILT] = np.eye(self.N_FILT)
        np.testing.assert_allclose(rec_sub @ self.im[:, :self.N_MODES], expected, atol=self.RTOL)

    def test_tip_tilt_removed_others_kept(self):
        im_sub, rec_sub = generate_filter_matrix_from_intmat(
            self.im, n_modes=self.N_MODES, n_modes_filtered=self.N_FILT)
        self._check_filtering(im_sub, rec_sub)

    def test_w_vec_not_supported(self):
        with self.assertRaises(NotImplementedError):
            generate_filter_matrix_from_intmat(
                self.im, n_modes=self.N_MODES, n_modes_filtered=self.N_FILT,
                w_vec=np.ones(self.N_MODES))

    def test_file_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            intmat_path = os.path.join(tmp, 'intmat.fits')
            filtmat_path = os.path.join(tmp, 'filtmat.fits')
            Intmat(intmat=self.im).save(intmat_path)

            generate_filter_matrix_from_intmat_file(
                intmat_path, n_modes=self.N_MODES, n_modes_filtered=self.N_FILT,
                output_filename=filtmat_path)
            im_sub, rec_sub = load_filter_matrix(filtmat_path)

        self._check_filtering(im_sub, rec_sub)


if __name__ == '__main__':
    unittest.main()
