"""
MMSE reconstructor W = (A' Cz^-1 A + Cx^-1)^-1 A' Cz^-1, computed with a
Cholesky solve and, for a sparse Cz^-1, without dense products: same result
as the explicit formula (float64 reference), on CPU and, if enabled, on GPU.
"""
import unittest
from unittest import mock

import numpy as np

import synim
import synim.params_utils as params_utils
from synim.params_utils import (compute_mmse_reconstructor, _sparse_rows, _sparse_rows_matmul,
                                 _symmetric_matmul)


def _elongated_noise_inv(n_sa, rng):
    """x/y slopes of each subaperture coupled, as the elongated LGS model."""
    n = 2 * n_sa
    c = np.zeros((n, n))
    a = rng.uniform(0.5, 2.0, n_sa)
    b = rng.uniform(0.5, 2.0, n_sa)
    d = rng.uniform(-0.3, 0.3, n_sa) * np.sqrt(a * b)
    i = np.arange(n_sa)
    c[i, i] = a
    c[i + n_sa, i + n_sa] = b
    c[i, i + n_sa] = d
    c[i + n_sa, i] = d
    return c


class TestMMSESolve(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(1)
        self.n_sa, self.n_modes = 60, 40
        self.A = rng.standard_normal((2 * self.n_sa, self.n_modes))
        q = rng.standard_normal((self.n_modes, self.n_modes))
        self.cx_inv = q @ q.T / self.n_modes + np.eye(self.n_modes)
        self.cz_inv = _elongated_noise_inv(self.n_sa, rng)
        self.modules = [np]
        if synim.xp is not np:
            self.modules.append(synim.xp)

    def _reference(self, cz_inv):
        A = self.A
        B = A.T @ cz_inv
        return np.linalg.solve(B @ A + self.cx_inv, B)

    def _check(self, cz_inv, rtol=2e-4, **kwargs):
        expected = self._reference(cz_inv)
        for xp in self.modules:
            with self.subTest(xp=xp.__name__):
                W = compute_mmse_reconstructor(
                    xp.asarray(self.A), xp.asarray(self.cx_inv), C_noise=xp.asarray(cz_inv),
                    cinverse=True, xp=xp, dtype=np.float32, **kwargs)
                self.assertIsInstance(W, np.ndarray)
                self.assertEqual(W.shape, (self.n_modes, 2 * self.n_sa))
                np.testing.assert_allclose(W, expected, rtol=0,
                                           atol=rtol * np.abs(expected).max())

    def test_sparse_noise(self):
        self._check(self.cz_inv)

    def test_dense_noise(self):
        rng = np.random.default_rng(2)
        q = rng.standard_normal((2 * self.n_sa, 2 * self.n_sa))
        self._check(q @ q.T / (2 * self.n_sa) + np.eye(2 * self.n_sa))

    def test_non_symmetric_noise(self):
        cz_inv = self.cz_inv.copy()
        cz_inv[0, self.n_sa] *= 1.5
        self._check(cz_inv)

    def test_fallback_when_factorization_fails(self):
        # A failed factorization falls back to the (pseudo-)inverse of H
        with mock.patch.object(params_utils, '_solve_inplace', return_value=None):
            self._check(self.cz_inv)

    def test_sparse_rows(self):
        for xp in self.modules:
            with self.subTest(xp=xp.__name__):
                cz = xp.asarray(self.cz_inv)
                cols, vals, symmetric = _sparse_rows(cz, xp)
                self.assertEqual(cols.shape, (2 * self.n_sa, 2))
                self.assertTrue(symmetric)
                A = xp.asarray(self.A)
                out = _sparse_rows_matmul(cols, vals, A, xp, block_elements=7 * self.n_modes)
                np.testing.assert_allclose(synim.cpuArray(out), self.cz_inv @ self.A, rtol=1e-12)
                self.assertIsNone(_sparse_rows(xp.ones((10, 10)), xp))



class TestSymmetricMatmul(unittest.TestCase):

    def test_tiles(self):
        """X' A (symmetric) from the lower tiles equals the full product."""
        xp = synim.xp
        rng = np.random.default_rng(12)
        for n in (100, 1100, 4100):  # 1 (single product), 2 and 8 tiles
            a = rng.standard_normal((200, n))
            x = a * rng.uniform(0.5, 2.0, 200)[:, None]  # diagonal Cz^-1
            out = synim.cpuArray(_symmetric_matmul(xp.asarray(x).T, xp.asarray(a), xp))
            np.testing.assert_allclose(out, x.T @ a, rtol=1e-12, atol=1e-12 * 200)
            if n >= 1024:  # tiles: exactly symmetric
                np.testing.assert_array_equal(out, out.T)


if __name__ == '__main__':
    unittest.main()
