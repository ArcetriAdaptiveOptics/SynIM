"""
Physical checks of the WFS transformations in interaction_matrix.

For a linear phase (tilt), the slopes of the transformed WFS are
    g' = A^T g
where g is the gradient of the phase and A is the (inverse mapping) matrix
of the transformation. Bilinear interpolation is exact for a linear phase,
so in the central, fully illuminated subapertures the relation holds up to
rounding.

The tests estimate the 2x2 matrix C that maps the slopes of the reference
WFS (no transformation) to the slopes of the transformed WFS, and check its
invariants (eigenvalues, orthogonality, rotation angle). These do not depend
on axis ordering or sign conventions.
"""
import unittest

import numpy as np

import synim
import synim.synim as synim_core

N_PIX = 120
N_SUBAPS = 12
PUP_DIAM_M = 8.0
TOL = 1e-4


def _tilt_modes(n):
    yy, xx = np.mgrid[:n, :n].astype(np.float32)
    return np.stack([xx - n / 2, yy - n / 2], axis=2)


class TestWfsTransformationsOnTilt(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        n = N_PIX
        yy, xx = np.mgrid[:n, :n]
        r = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5)
        cls.pup_mask = (r < n / 2 - 2).astype(np.float32)
        cls.dm_mask = np.ones((n, n), dtype=np.float32)
        cls.modes = _tilt_modes(n)
        # Central subapertures, far from the pupil edge for all the tested
        # transformations
        sa = (np.arange(N_SUBAPS) + 0.5) * n / N_SUBAPS - n / 2
        sa_y, sa_x = np.meshgrid(sa, sa, indexing='ij')
        cls.central = (np.hypot(sa_x, sa_y) < 0.3 * n).ravel()
        cls.reference = cls._slope_matrix(cls, {})

    def _slope_matrix(self, wfs_params, slope_method='derivatives'):
        """2x2 matrix [slope axis, mode] averaged over central subapertures."""
        params = dict(dm_height=0.0, dm_rotation=0.0, gs_pol_coo=(0.0, 0.0),
                      gs_height=np.inf, wfs_rotation=0.0, wfs_translation=(0.0, 0.0),
                      wfs_mag_global=1.0)
        params.update(wfs_params)
        im = synim.cpuArray(synim_core.interaction_matrix(
            PUP_DIAM_M, self.pup_mask, self.modes, self.dm_mask,
            wfs_nsubaps=N_SUBAPS, wfs_fov_arcsec=4.0, slope_method=slope_method,
            **params)).astype(np.float64)
        n_sa = N_SUBAPS * N_SUBAPS
        first, second = im[:n_sa][self.central], im[n_sa:][self.central]
        # Uniform gradient: all central subapertures must give the same slopes
        for block in (first, second):
            spread = np.max(np.abs(block - block.mean(axis=0)))
            assert spread <= TOL * np.max(np.abs(block)) + 1e-12, spread
        return np.array([first.mean(axis=0), second.mean(axis=0)])

    def _transfer_matrix(self, wfs_params, slope_method='derivatives'):
        reference = self.reference if slope_method == 'derivatives' \
            else self._slope_matrix({}, slope_method)
        return self._slope_matrix(wfs_params, slope_method) @ np.linalg.inv(reference)

    def assert_eigenvalues(self, c, expected):
        np.testing.assert_allclose(np.sort(np.linalg.eigvals(c).real), np.sort(expected),
                                   rtol=0, atol=TOL)

    def test_translation(self):
        c = self._transfer_matrix(dict(wfs_translation=(1.3, -0.7)))
        np.testing.assert_allclose(c, np.eye(2), rtol=0, atol=TOL)

    def test_magnification(self):
        for mag in (0.95, 1.05):
            with self.subTest(magnification=mag):
                c = self._transfer_matrix(dict(wfs_mag_global=mag))
                np.testing.assert_allclose(c, np.eye(2) / mag, rtol=0, atol=TOL)

    def test_anamorphosis_90(self):
        anam = 1.05
        c = self._transfer_matrix(dict(wfs_anamorphosis_90=anam))
        np.testing.assert_allclose(c - np.diag(np.diag(c)), 0, rtol=0, atol=TOL)
        self.assert_eigenvalues(c, [1.0, 1.0 / anam])

    def test_anamorphosis_45(self):
        anam = 1.05
        c = self._transfer_matrix(dict(wfs_anamorphosis_45=anam))
        np.testing.assert_allclose(c, c.T, rtol=0, atol=TOL)
        self.assert_eigenvalues(c, [1.0, anam])
        # Eigenvectors along the diagonals
        _, vectors = np.linalg.eigh(c)
        np.testing.assert_allclose(np.abs(vectors), np.sqrt(0.5), rtol=0, atol=1e-3)

    def test_rotation(self):
        for angle in (4.0, 30.0):
            with self.subTest(rotation=angle):
                c = self._transfer_matrix(dict(wfs_rotation=angle))
                np.testing.assert_allclose(c @ c.T, np.eye(2), rtol=0, atol=TOL)
                self.assertAlmostEqual(np.linalg.det(c), 1.0, delta=TOL)
                measured = np.degrees(np.arccos(np.clip(c[0, 0], -1, 1)))
                self.assertAlmostEqual(measured, angle, delta=1e-2)

    def test_telsum_magnification(self):
        mag = 1.05
        c = self._transfer_matrix(dict(wfs_mag_global=mag), slope_method='telsum')
        np.testing.assert_allclose(c, np.eye(2) / mag, rtol=0, atol=TOL)


if __name__ == '__main__':
    unittest.main()
