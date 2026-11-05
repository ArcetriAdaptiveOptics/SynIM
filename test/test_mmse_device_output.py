"""
compute_mmse_reconstructor(to_host=False) returns the reconstructor on the
compute device (xp array), equal to the host result.
"""
import unittest

import numpy as np

import synim
from synim.params_utils import compute_mmse_reconstructor


class TestMMSEDeviceOutput(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(3)
        self.A = rng.standard_normal((80, 30)).astype(np.float32)
        self.cx_inv = np.eye(30, dtype=np.float32) * 2.0
        self.cz_inv = np.eye(80, dtype=np.float32) * 0.5

    def test_device_output(self):
        xp = synim.xp
        args = (xp.asarray(self.A), xp.asarray(self.cx_inv))
        kwargs = dict(C_noise=xp.asarray(self.cz_inv), cinverse=True, xp=xp, dtype=np.float32)
        on_host = compute_mmse_reconstructor(*args, **kwargs)
        on_device = compute_mmse_reconstructor(*args, to_host=False, **kwargs)
        self.assertIsInstance(on_host, np.ndarray)
        self.assertIsInstance(on_device, xp.ndarray)
        np.testing.assert_array_equal(synim.cpuArray(on_device), on_host)


if __name__ == '__main__':
    unittest.main()
