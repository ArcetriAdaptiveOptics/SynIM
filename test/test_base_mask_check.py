"""
Check that a 2D projection basis was defined on the same pupil as the pupil mask:
- check_base_mask_matches_pupil compares the valid-pixel patterns cropped to
  their bounding boxes (position in the frame is irrelevant)
- why it matters: dm2d_to_3d maps the basis rows onto the valid pixels in
  row-major order, so a moved pixel silently corrupts the map
- ParamsManager._load_base_inv_array calls the check for 2D bases

The file I/O is replaced by in-memory objects.
"""
import unittest
from unittest import mock

import numpy as np

import synim.params_manager as params_manager
from synim.params_utils import check_base_mask_matches_pupil
from synim.utils import dm2d_to_3d


def _pupil(n=24, radius=11, obstruction=4):
    """Disk with a central obstruction, as a bool array."""
    yy, xx = np.mgrid[:n, :n]
    r = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5)
    return (r < radius) & (r >= obstruction)


def _pad(mask, p):
    return np.pad(mask, p)


def _move_one_pixel(mask):
    """Switch off one valid pixel and on one invalid pixel inside the obstruction."""
    out = mask.copy()
    n = mask.shape[0]
    centre = n // 2
    off = tuple(np.argwhere(mask)[0])
    assert not out[centre, centre]
    out[off] = False
    out[centre, centre] = True
    return out


class TestCheckBaseMaskMatchesPupil(unittest.TestCase):

    def setUp(self):
        self.mask = _pupil()

    def test_identical(self):
        self.assertIsNone(check_base_mask_matches_pupil(self.mask, self.mask.copy()))

    def test_padded_both_orders(self):
        padded = _pad(self.mask, 3)
        self.assertEqual(padded.shape, (30, 30))
        self.assertIsNone(check_base_mask_matches_pupil(self.mask, padded))
        self.assertIsNone(check_base_mask_matches_pupil(padded, self.mask))

    def test_shifted_in_same_frame(self):
        big = _pad(self.mask, 4)
        shifted = np.roll(big, (2, -3), axis=(0, 1))
        self.assertIsNone(check_base_mask_matches_pupil(big, shifted))

    def test_one_pixel_moved(self):
        moved = _move_one_pixel(self.mask)
        self.assertEqual(moved.sum(), self.mask.sum())
        with self.assertRaisesRegex(ValueError, '2 pixels differ'):
            check_base_mask_matches_pupil(self.mask, moved)

    def test_different_size(self):
        other = _pupil(radius=9)
        with self.assertRaises(ValueError):
            check_base_mask_matches_pupil(self.mask, other)

    def test_float_and_bool_masks(self):
        as_float = self.mask.astype(np.float32)
        self.assertIsNone(check_base_mask_matches_pupil(as_float, self.mask))
        self.assertIsNone(check_base_mask_matches_pupil(self.mask, as_float))
        self.assertIsNone(check_base_mask_matches_pupil(as_float, _pad(as_float, 3)))
        with self.assertRaises(ValueError):
            check_base_mask_matches_pupil(as_float, _move_one_pixel(self.mask).astype(float))

    def test_empty_mask_raises(self):
        with self.assertRaises(ValueError):
            check_base_mask_matches_pupil(self.mask, np.zeros_like(self.mask))


class TestWhyTheCheckMatters(unittest.TestCase):
    """Tip mode (x coordinate on the valid pixels) mapped with dm2d_to_3d."""

    def setUp(self):
        self.mask = _pupil()
        self.xx = np.mgrid[:24, :24][1].astype(float)
        self.basis = self.xx[self.mask][None, :]  # (1, n_valid), row-major order

    def test_same_mask_padded(self):
        padded = _pad(self.mask, 3)
        cube = dm2d_to_3d(self.basis, padded, normalize=False,
                          xp_local=np, float_dtype_local=np.float64)
        self.assertEqual(cube.shape, (30, 30, 1))
        # x coordinate in the padded frame is shifted by the padding
        expected = np.where(padded, np.mgrid[:30, :30][1] - 3, 0.0)
        np.testing.assert_allclose(cube[:, :, 0], expected)
        self.assertIsNone(check_base_mask_matches_pupil(self.mask, padded))

    def test_moved_pixel_corrupts_map(self):
        moved = _move_one_pixel(self.mask)
        cube = dm2d_to_3d(self.basis, moved, normalize=False,
                          xp_local=np, float_dtype_local=np.float64)
        expected = np.where(moved, self.xx, 0.0)
        wrong = (cube[:, :, 0] != expected) & moved
        self.assertGreater(wrong.sum(), 0)
        with self.assertRaises(ValueError):
            check_base_mask_matches_pupil(self.mask, moved)


class TestLoadBaseInvArray(unittest.TestCase):

    def setUp(self):
        self.pm = params_manager.ParamsManager.__new__(params_manager.ParamsManager)
        self.pm.params = {}
        self.pm.projection_params = {'ifunc_inverse_tag': 'test_inv', 'nmodes': 3}
        self.pm.cm = object()
        self.pm.pixel_pupil = 24
        self.pm.pup_mask = _pupil().astype(np.float32)
        self.n_valid = int(_pupil().sum())

    def _load(self, basis, basis_mask):
        with mock.patch.object(params_manager, 'load_influence_functions',
                               return_value=(basis, basis_mask)):
            return self.pm._load_base_inv_array(verbose=False)

    def test_matching_mask_returns_array(self):
        basis = np.ones((3, self.n_valid), dtype=np.float32)
        out = self._load(basis, _pad(_pupil(), 3).astype(np.float32))
        self.assertEqual(out.shape, basis.shape)

    def test_mismatching_mask_raises(self):
        basis = np.ones((3, self.n_valid), dtype=np.float32)
        with self.assertRaises(ValueError):
            self._load(basis, _move_one_pixel(_pupil()).astype(np.float32))

    def test_3d_basis_not_checked(self):
        basis = np.ones((24, 24, 3), dtype=np.float32)
        out = self._load(basis, np.zeros((24, 24), dtype=np.float32))
        self.assertEqual(out.shape, basis.shape)


if __name__ == '__main__':
    unittest.main()
