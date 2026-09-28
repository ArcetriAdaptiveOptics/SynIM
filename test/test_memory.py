"""
Peak memory of interaction_matrix, measured on CPU with tracemalloc and
expressed in units of the input DM cube (n x n x n_modes).

Reference values (numpy, float32), before -> after the memory optimizations:
    derivatives, integer ratio pixels/subapertures      6.6 -> 3.2
    derivatives, non-integer ratio (upscaled rebin)     7.7 -> 3.7
    telsum, non-integer ratio                           5.8 -> 3.1
The thresholds leave a margin over the current values and fail if the
memory use goes back towards the previous values.
"""
import contextlib
import io
import tracemalloc
import unittest

import numpy as np

import synim
import synim.synim as synim_core

CASES = [
    # (n_pix, n_subaps, slope_method, max peak in cubes)
    (160, 20, 'derivatives', 3.8),
    (160, 23, 'derivatives', 4.3),
    (160, 23, 'telsum', 3.8),
]


@unittest.skipUnless(synim.xp is np, 'memory is measured with tracemalloc on CPU')
class TestInteractionMatrixMemory(unittest.TestCase):

    def _peak_in_cubes(self, n_pix, n_subaps, slope_method, n_modes=40):
        yy, xx = np.mgrid[:n_pix, :n_pix]
        r = np.hypot(xx - n_pix / 2 + 0.5, yy - n_pix / 2 + 0.5)
        pup_mask = ((r < n_pix / 2 - 2) & (r > 0.28 * n_pix / 2)).astype(np.float32)
        dm_mask = (r < n_pix / 2).astype(np.float32)
        dm_array = np.random.default_rng(0).standard_normal(
            (n_pix, n_pix, n_modes)).astype(np.float32)
        tracemalloc.start()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                synim_core.interaction_matrix(
                    8.0, pup_mask, dm_array, dm_mask, 0.0, 0.0, n_subaps, 4.0,
                    (20.0, 30.0), 90e3, 15.0, (0.2, -0.1), 1.0, slope_method=slope_method)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        return peak / dm_array.nbytes

    def test_peak_memory(self):
        for n_pix, n_subaps, slope_method, max_cubes in CASES:
            with self.subTest(n_pix=n_pix, n_subaps=n_subaps, slope_method=slope_method):
                peak = self._peak_in_cubes(n_pix, n_subaps, slope_method)
                self.assertLessEqual(peak, max_cubes,
                                     msg=f'peak memory {peak:.2f} cubes > {max_cubes}')


if __name__ == '__main__':
    unittest.main()
