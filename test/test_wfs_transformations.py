"""
Physical tests of the interaction matrix: results predicted from the
physics, independently of the implementation.

They run on the reference implementation (method='reference',
synim.reference), where the physics is defined, and on the optimized one
(method='operator'). Unlike the equivalence test (test_sparse_operator.py),
they also check the geometry shared by the two implementations
(synim.geometry): an error there (e.g. the sign of a rotation, or the
masks) would not show up as a difference between the two.

For a linear phase (tilt), the slopes of the transformed WFS are
    g' = A^T g
where g is the gradient of the phase and A is the (inverse mapping) matrix
of the transformation. Bilinear interpolation is exact for a linear phase,
so in the central, fully illuminated subapertures the relation holds up to
rounding.

The tests estimate the 2x2 matrix C that maps the slopes of the reference
WFS (no transformation) to the slopes of the transformed WFS, and check its
invariants (eigenvalues, orthogonality, rotation angle). These do not depend
on axis ordering or sign conventions, except for the 90 deg rotation test,
which pins the current sign convention of the rotation.
The absolute value of the slopes of a tilt is checked against the
definition of the slope units.
"""
import unittest

import numpy as np

import synim
import synim.synim as synim_core

N_PIX = 120
N_SUBAPS = 12
PUP_DIAM_M = 8.0
FOV_ARCSEC = 4.0
TOL = 1e-4


def _tilt_modes(n):
    yy, xx = np.mgrid[:n, :n].astype(np.float32)
    return np.stack([xx - n / 2, yy - n / 2], axis=2)


class TestWfsTransformationsOnTilt(unittest.TestCase):

    method = 'reference'

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
            wfs_nsubaps=N_SUBAPS, wfs_fov_arcsec=FOV_ARCSEC, slope_method=slope_method,
            method=self.method, **params)).astype(np.float64)
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

    def test_tilt_slope_value(self):
        """A tilt of 1 nm per pixel gives a phase difference of N nm across a
        subaperture of N pixels and size d, i.e. a tilt of N * 1e-9 / d rad;
        the slopes are in units of half the WFS field of view."""
        n_pix_sa = N_PIX / N_SUBAPS
        d = PUP_DIAM_M / N_SUBAPS
        expected = n_pix_sa * 1e-9 / d * 206265 / (FOV_ARCSEC / 2)
        for slope_method in ('derivatives', 'telsum'):
            for specula_convention in (True, False):
                with self.subTest(slope_method=slope_method,
                                  specula_convention=specula_convention):
                    s = self._slope_matrix(dict(specula_convention=specula_convention),
                                           slope_method)
                    # the first slope is along the second axis of the DM array
                    # (the first mode) in both conventions
                    np.testing.assert_allclose(s, expected * np.eye(2), rtol=0,
                                               atol=TOL * expected)

    def test_rotation_90(self):
        """A 90 deg rotation swaps the x and y slopes (with the current sign
        convention, x -> y and y -> -x), as three 30 deg rotations."""
        c = self._transfer_matrix(dict(wfs_rotation=90.0))
        np.testing.assert_allclose(c, [[0, -1], [1, 0]], rtol=0, atol=TOL)
        c_30 = self._transfer_matrix(dict(wfs_rotation=30.0))
        np.testing.assert_allclose(c, np.linalg.matrix_power(c_30, 3), rtol=0, atol=10 * TOL)

    def test_lgs_cone(self):
        """With a guide star at height H_lgs, the footprint on a DM at height
        h_dm is smaller by 1 - h_dm / H_lgs, and the slopes of a tilt on the DM
        scale by the same factor."""
        h_dm, H_lgs = 8000.0, 90e3
        c = self._transfer_matrix(dict(dm_height=h_dm, gs_height=H_lgs))
        np.testing.assert_allclose(c, (1 - h_dm / H_lgs) * np.eye(2), rtol=0, atol=TOL)

    def test_off_axis_shift(self):
        """An off-axis guide star shifts the footprint on an altitude DM:
        the slopes of a tilt do not change."""
        c = self._transfer_matrix(dict(dm_height=8000.0, gs_pol_coo=(30.0, 40.0)))
        np.testing.assert_allclose(c, np.eye(2), rtol=0, atol=TOL)

    def test_dm_rotation(self):
        """A DM rotation rotates the slopes of a tilt. dm_rotation and
        wfs_rotation rotate the phase map in the same direction (as the
        SPECULA rotInDeg of a DM or layer and rotAnglePhInDeg of a SH WFS),
        so the same angle gives the same transfer matrix (this pins the sign
        of the DM rotation, as test_rotation_90 does for the WFS)."""
        for angle in (20.0, -20.0, 90.0):
            with self.subTest(dm_rotation=angle):
                c = self._transfer_matrix(dict(dm_rotation=angle))
                np.testing.assert_allclose(c @ c.T, np.eye(2), rtol=0, atol=TOL)
                self.assertAlmostEqual(np.linalg.det(c), 1.0, delta=TOL)
                measured = np.degrees(np.arccos(np.clip(c[0, 0], -1, 1)))
                self.assertAlmostEqual(measured, abs(angle), delta=1e-2)
                np.testing.assert_allclose(
                    c, self._transfer_matrix(dict(wfs_rotation=angle)), rtol=0, atol=TOL)
        np.testing.assert_allclose(self._transfer_matrix(dict(dm_rotation=90.0)),
                                   [[0, -1], [1, 0]], rtol=0, atol=TOL)


class TestWfsTransformationsOnTiltOperator(TestWfsTransformationsOnTilt):
    """The same physical tests on the optimized implementation."""

    method = 'operator'


if __name__ == '__main__':
    unittest.main()
