"""
x/y axis consistency between the elongated LGS noise model and the
interaction matrices.

A WFS has two slope blocks (all the subapertures of one axis, then all of the
other). The interaction matrices (synim.sparse_operator) and the elongated
noise precision (ParamsManager._build_elongated_noise_covariance, which calls
SPECULA calc_noise_cov_elong) must agree on which block is which. If they did
not, the elongation axes (beta1/beta2) would be swapped.

Physics: the LGS spot is elongated along the line joining the subaperture to
the launcher, and the elongation grows with the distance from the launcher.
With the launcher on the x axis, the subaperture at the far edge along x
(centred in y) has a spot elongated along x: the noise variance is largest on
the slope block that measures the tilt along x. With the launcher on the y
axis it is the other block.

The test does not assume which block is x. Everything is read from the SynIM
interaction matrix of a few analytic DM modes (the phase arrays are indexed
[row, column], x is the column axis, as in SPECULA):
  - ramp along the column axis: the block that responds is the "x tilt" block
    (the block that responds to the row ramp is the "y tilt" block);
  - parabola along the column (row) axis: in the x (y) block the slope of a
    subaperture is proportional to its x (y) position, which locates every
    subaperture of the interaction matrix.
The 2x2 noise covariance of the far-edge subaperture, the inverse of the block
[[P_ii, P_i,i+n], [P_i,i+n, P_i+n,i+n]] of the precision matrix, must have its
major axis within MAX_ANGLE_DEG of the slope block of the launcher axis.

The valid subapertures are given as SPECULA SubapData does ([column, row]
pairs, listed row by row), for an asymmetric mask so that a wrong ordering of
the slopes would also show up. Everything is on CPU-sized arrays, no SPECULA
simulation is run.
"""
import unittest

import numpy as np

import synim
import synim.synim as synim_core
from synim.params_manager import ParamsManager
from test import on_backend

interaction_matrix = on_backend(synim_core.interaction_matrix)

N_PIX = 48
N_SUBAPS = 8
DIAMETER_M = 8.0
LAUNCHER_OFFSET_M = DIAMETER_M / 2
MAX_ANGLE_DEG = 20.0
MIN_AXIS_RATIO = 2.0  # major / minor variance, the elongation is large here


def _valid_subaps():
    """Asymmetric round mask, valid subapertures as [column, row] pairs row by row."""
    i, j = np.mgrid[:N_SUBAPS, :N_SUBAPS]  # i = row, j = column
    mask = np.hypot(i - 3.5, j - 3.5) < 4.2
    mask[0, 2] = False
    return np.array([(jj, ii) for ii in range(N_SUBAPS) for jj in range(N_SUBAPS)
                     if mask[ii, jj]])


def _interaction_matrix(idx_valid_sa):
    """IM (numpy) of 4 modes: row ramp, column ramp, row parabola, column parabola."""
    c = (np.arange(N_PIX) - N_PIX / 2 + 0.5) / (N_PIX / 2)
    row = np.broadcast_to(c[:, None], (N_PIX, N_PIX))
    col = np.broadcast_to(c[None, :], (N_PIX, N_PIX))
    yy, xx = np.mgrid[:N_PIX, :N_PIX]
    pup = (np.hypot(xx - N_PIX / 2 + 0.5, yy - N_PIX / 2 + 0.5) < N_PIX / 2 - 1)
    pup = pup.astype(np.float32)
    modes = (np.stack([row, col, row**2, col**2], axis=-1) * pup[..., None]).astype(np.float32)
    im = interaction_matrix(
        DIAMETER_M, pup, modes, pup, 0.0, 0.0, N_SUBAPS, 3.0, (0.0, 0.0), np.inf,
        0.0, (0.0, 0.0), 1.0, idx_valid_sa=idx_valid_sa, specula_convention=True)
    return np.asarray(im, dtype=float)


def _make_pm(idx_valid_sa, launcher_xy_m):
    """Minimal ParamsManager (no __init__) for one LGS WFS."""
    pm = ParamsManager.__new__(ParamsManager)
    pm.params = {
        'main': {'pixel_pupil': N_PIX, 'pixel_pitch': DIAMETER_M / N_PIX},
        'sh_lgs1': {'subap_on_diameter': N_SUBAPS, 'subap_npx': 8,
                    'sensor_pxscale': 0.5, 'laser_launch_tel_ref': 'launcher1'},
        'launcher1': {'tel_position': [launcher_xy_m[0], launcher_xy_m[1], 0.0],
                      'beacon_focus': 90000.0, 'spot_size': 1.0},
    }
    # tGparameter = 0: no truncated subapertures (that branch divides by the
    # elongation, zero for a subaperture aligned with the launcher)
    pm.lgs_recon_params = {'sigma2_in_nm2': 1e6, 'naThicknessInM': 10000.0,
                           'tGparameter': 0.0}
    pm.get_wfs_params = lambda wfs_type, idx, xp_local=np: {
        'idx_valid_sa': idx_valid_sa, 'wfs_nsubaps': N_SUBAPS, 'wfs_fov_arcsec': 4.0}
    pm._count_wfs = lambda wfs_type: 1
    return pm


class TestNoiseCovAxes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.idx_valid_sa = _valid_subaps()
        cls.n = len(cls.idx_valid_sa)
        im = _interaction_matrix(cls.idx_valid_sa)
        n = cls.n
        blocks = (slice(0, n), slice(n, 2 * n))

        def responding_block(mode):
            energy = [np.linalg.norm(im[b, mode]) for b in blocks]
            responding = int(np.argmax(energy))
            # the other block must be silent, not just smaller
            assert energy[1 - responding] < 1e-3 * energy[responding], energy
            return responding

        # x = column axis, y = row axis
        cls.x_block = responding_block(1)
        cls.y_block = responding_block(0)
        assert cls.x_block != cls.y_block
        # positions of the subapertures along x and y, in IM order (arbitrary units)
        cls.pos_x = im[blocks[cls.x_block], 3]
        cls.pos_y = im[blocks[cls.y_block], 2]

    def _far_subaperture(self, axis):
        """Index of the subaperture at the minimum of `axis`, centred on the other axis."""
        along, across = (self.pos_x, self.pos_y) if axis == 'x' else (self.pos_y, self.pos_x)
        scale = np.abs(along).max()
        return int(np.argmin(along / scale + 10 * np.abs(across) / scale))

    def _noise_cov_2x2(self, launcher_xy_m, k):
        pm = _make_pm(self.idx_valid_sa, launcher_xy_m)
        precision = pm._build_elongated_noise_covariance('lgs')
        self.assertEqual(precision.shape, (2 * self.n, 2 * self.n))
        self.assertTrue(np.all(np.isfinite(precision)))
        block = precision[np.ix_([k, k + self.n], [k, k + self.n])].astype(float)
        return np.linalg.inv(block)

    def _check(self, axis, launcher_xy_m):
        k = self._far_subaperture(axis)
        # the far edge is at the minimum, i.e. on the opposite side of the launcher
        along = self.pos_x if axis == 'x' else self.pos_y
        self.assertLess(along[k], 0)
        cov = self._noise_cov_2x2(launcher_xy_m, k)  # indices [block 0, block 1]
        w, v = np.linalg.eigh(cov)
        major = v[:, 1]
        launcher_block = self.x_block if axis == 'x' else self.y_block
        angle = np.degrees(np.arccos(min(1.0, abs(major[launcher_block]))))
        self.assertLess(angle, MAX_ANGLE_DEG, f'angle={angle:.1f} deg, cov={cov}')
        self.assertGreater(w[1] / w[0], MIN_AXIS_RATIO)
        # the variance on the launcher-axis block is the larger one
        self.assertGreater(cov[launcher_block, launcher_block],
                           cov[1 - launcher_block, 1 - launcher_block])

    def test_assembled_on_device(self):
        pm = _make_pm(self.idx_valid_sa, (LAUNCHER_OFFSET_M, 0.0))
        expected = pm._build_elongated_noise_covariance('lgs')
        on_device = pm._build_elongated_noise_covariance('lgs', to_host=False)
        self.assertIsInstance(on_device, synim.xp.ndarray)
        np.testing.assert_array_equal(synim.cpuArray(on_device), expected)

    def test_launcher_on_x(self):
        self._check('x', (LAUNCHER_OFFSET_M, 0.0))

    def test_launcher_on_y(self):
        self._check('y', (0.0, LAUNCHER_OFFSET_M))

    def test_axes_swap_with_launcher(self):
        """Moving the launcher from x to y moves the large variance to the other block."""
        kx = self._far_subaperture('x')
        ky = self._far_subaperture('y')
        cov_x = self._noise_cov_2x2((LAUNCHER_OFFSET_M, 0.0), kx)
        cov_y = self._noise_cov_2x2((0.0, LAUNCHER_OFFSET_M), ky)
        self.assertEqual(np.argmax(np.diag(cov_x)), self.x_block)
        self.assertEqual(np.argmax(np.diag(cov_y)), self.y_block)


if __name__ == '__main__':
    unittest.main()
