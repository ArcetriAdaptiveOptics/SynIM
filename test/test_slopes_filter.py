"""
Offline slopes filter of the interaction matrices (_apply_slopes_filter):
- same result as the explicit IDL formula
- the filter file is read once per ParamsManager
- on GPU, the result is the same as on CPU
"""
import os
import tempfile
import unittest
from unittest import mock

import numpy as np
from astropy.io import fits

import synim.params_manager as params_manager


class _FakeCalibManager:
    def __init__(self, root_dir):
        self.root_dir = root_dir


class TestSlopesFilter(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(0)
        self.n_slopes, self.n_modes, self.n_filt = 40, 12, 3
        # big-endian float64, as the files written by IDL
        self.filtmat = rng.standard_normal((2, self.n_slopes, self.n_filt)).astype('>f8')
        self.im = rng.standard_normal((self.n_slopes, self.n_modes)).astype(np.float32)

        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        os.makedirs(os.path.join(self.tmp.name, 'data'))
        fits.writeto(os.path.join(self.tmp.name, 'data', 'filt_test.fits'), self.filtmat)

        self.pm = params_manager.ParamsManager.__new__(params_manager.ParamsManager)
        self.pm.params = {'slopec_lgs1': {'filtmat_data': 'filt_test'}}
        self.pm.cm = _FakeCalibManager(self.tmp.name)
        self.pm._filtmat_cache = {}  # set by __init__, bypassed here

    def _expected(self):
        f = self.filtmat.astype(np.float64)
        return self.im - f[0] @ (f[1].T @ self.im)

    def test_same_as_formula(self):
        out = self.pm._apply_slopes_filter(self.im, 'lgs', 1)
        np.testing.assert_allclose(out, self._expected(), rtol=0, atol=1e-12)

    def test_file_read_once(self):
        with mock.patch.object(params_manager.fits, 'open', wraps=params_manager.fits.open) as spy:
            self.pm._apply_slopes_filter(self.im, 'lgs', 1)
            self.pm._apply_slopes_filter(self.im, 'lgs', 1)
        self.assertEqual(spy.call_count, 1)

    def test_gpu_same_as_cpu(self):
        xp = params_manager.xp
        if xp is np:
            self.skipTest('GPU not enabled')
        out = self.pm._apply_slopes_filter(xp.asarray(self.im), 'lgs', 1)
        self.assertIsInstance(out, xp.ndarray)
        np.testing.assert_allclose(out.get(), self._expected(), rtol=0, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
