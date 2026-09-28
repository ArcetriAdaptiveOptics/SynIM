"""
apply_wfs_transformations_combined must not change the result of later
calls on the same derivatives: by default it works on copies (in_place=False).
The pipeline uses in_place=True, which is valid only when all the calls on
the same derivatives use the same pupil mask.
"""
import unittest

import numpy as np

import synim
import synim.synim as synim_core


def _circular(n, radius, center=None):
    yy, xx = np.mgrid[:n, :n]
    cy, cx = center if center is not None else (n / 2 - 0.5, n / 2 - 0.5)
    return (np.hypot(xx - cx, yy - cy) < radius).astype(np.float32)


class TestDerivativesNotModified(unittest.TestCase):

    def setUp(self):
        xp = synim.xp
        rng = np.random.default_rng(0)
        # 100 px and 12 subapertures: non-integer ratio (upscaled rebin)
        self.n, self.nsa = 100, 12
        self.dm_mask = xp.asarray(np.ones((self.n, self.n), dtype=np.float32))
        self.pup_a = xp.asarray(_circular(self.n, 40))
        self.pup_b = xp.asarray(_circular(self.n, 46, center=(52.0, 47.0)))
        self.dx = xp.asarray(rng.standard_normal((self.n, self.n, 3)).astype(np.float32))
        self.dy = xp.asarray(rng.standard_normal((self.n, self.n, 3)).astype(np.float32))

    def _slopes(self, dx, dy, pup_mask, **kwargs):
        return synim.cpuArray(synim_core.apply_wfs_transformations_combined(
            dx, dy, pup_mask, self.dm_mask, self.nsa, 4.0, 8.0, **kwargs))

    def test_default_does_not_modify_inputs(self):
        dx, dy = self.dx.copy(), self.dy.copy()
        self._slopes(dx, dy, self.pup_a)
        np.testing.assert_array_equal(synim.cpuArray(dx), synim.cpuArray(self.dx))
        np.testing.assert_array_equal(synim.cpuArray(dy), synim.cpuArray(self.dy))

    def test_two_pupil_masks(self):
        # Same derivatives, two different pupil masks: the second result must
        # be the one obtained from untouched derivatives
        dx, dy = self.dx.copy(), self.dy.copy()
        self._slopes(dx, dy, self.pup_a)
        second = self._slopes(dx, dy, self.pup_b)
        expected = self._slopes(self.dx.copy(), self.dy.copy(), self.pup_b)
        np.testing.assert_array_equal(second, expected)

    def test_in_place_same_mask_is_idempotent(self):
        # What the pipeline relies on (multi-WFS groups share the pupil mask)
        dx, dy = self.dx.copy(), self.dy.copy()
        first = self._slopes(dx, dy, self.pup_a, in_place=True)
        second = self._slopes(dx, dy, self.pup_a, in_place=True)
        np.testing.assert_array_equal(first, second)
        expected = self._slopes(self.dx.copy(), self.dy.copy(), self.pup_a)
        np.testing.assert_array_equal(first, expected)


if __name__ == '__main__':
    unittest.main()
