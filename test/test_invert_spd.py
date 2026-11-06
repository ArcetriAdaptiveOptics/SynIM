"""
invert_spd: Cholesky inverse of symmetric positive definite matrices, with the
pseudo-inverse as fallback; and its use in assemble_covariance_matrix
(return_inverse=True), which inverts each covariance block.
"""
import unittest
from unittest import mock

import numpy as np

import synim
import synim.params_manager as params_manager
from synim.params_manager import ParamsManager
from synim.params_utils import invert_spd

CONVERSION = (500 / 2 / np.pi) ** 2


def _graded_spd(n, seed):
    """SPD matrix with variances spread over 4 decades, like a modal covariance."""
    rng = np.random.default_rng(seed)
    b = rng.standard_normal((n, n))
    c = b @ b.T / n + np.eye(n)
    d = np.logspace(0, -2, n)
    return c * d[:, None] * d[None, :]


class TestInvertSpd(unittest.TestCase):

    def test_spd(self):
        xp = synim.xp
        a = _graded_spd(60, 1)
        a_dev = xp.asarray(a)
        inv = invert_spd(a_dev, xp=xp)
        self.assertIsInstance(inv, xp.ndarray)
        np.testing.assert_allclose(synim.cpuArray(inv), np.linalg.inv(a), rtol=1e-9,
                                   atol=1e-9 * np.abs(np.linalg.inv(a)).max())
        # the input is not modified
        np.testing.assert_array_equal(synim.cpuArray(a_dev), a)

    def test_numpy_input_on_device(self):
        a = _graded_spd(20, 2)
        inv = invert_spd(a, xp=synim.xp)
        np.testing.assert_allclose(synim.cpuArray(inv) @ a, np.eye(20), atol=1e-9)

    def test_singular_falls_back_to_pinv(self):
        xp = synim.xp
        rng = np.random.default_rng(3)
        b = rng.standard_normal((30, 10))
        a = b @ b.T  # rank 10: not positive definite
        inv = synim.cpuArray(invert_spd(xp.asarray(a), xp=xp))
        self.assertTrue(np.all(np.isfinite(inv)))
        np.testing.assert_allclose(inv, np.linalg.pinv(a), rtol=1e-6,
                                   atol=1e-6 * np.abs(np.linalg.pinv(a)).max())


class TestAssembleCovarianceInverse(unittest.TestCase):

    def test_block_inverses(self):
        blocks = [_graded_spd(12, 4), _graded_spd(9, 5)]
        modes = [list(range(2, 12)), list(range(9))]
        weights = [0.7, 0.3]
        pm = ParamsManager.__new__(ParamsManager)
        pm.verbose = False
        pm.params = {}
        with mock.patch.object(params_manager, 'compute_layer_weights_from_turbulence',
                               return_value=weights):
            inv = pm.assemble_covariance_matrix(blocks, [1, 2], mode_indices=modes,
                                                return_inverse=True)
        self.assertIsInstance(inv, np.ndarray)
        expected = np.zeros((19, 19))
        o = 0
        for block, m, w in zip(blocks, modes, weights):
            n = len(m)
            expected[o:o + n, o:o + n] = np.linalg.inv(block[np.ix_(m, m)] * w * CONVERSION)
            o += n
        np.testing.assert_allclose(inv, expected, rtol=1e-6, atol=1e-6 * np.abs(expected).max())


if __name__ == '__main__':
    unittest.main()
