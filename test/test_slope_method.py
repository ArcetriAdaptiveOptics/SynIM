import unittest
import numpy as np

# Initialize specula (needed to set CPU/GPU backend for synim)
import specula
specula.init(device_idx=-1, precision=1)

from synim.synim import (
    compute_gtilt_with_extrapolation,
    compute_derivatives_with_extrapolation,
    _compute_slopes_from_derivatives,
    _compute_slopes_from_gtilt
)

class TestSlopesMethods(unittest.TestCase):
    def setUp(self):
        """Set up test parameters for a single sub-aperture."""
        self.N = 10
        self.wfs_nsubaps = 1

        # Create coordinate grid
        y, x = np.mgrid[0:self.N, 0:self.N]

        # Aberration: Pure 30-degree tilt (Gx=0.866, Gy=0.5)
        self.theta = np.deg2rad(30)
        self.Gx = np.cos(self.theta)
        self.Gy = np.sin(self.theta)

        # Phase wavefront Z(x,y)
        self.phase = self.Gx * x + self.Gy * y

        # Physical mocked constants to compute slopes
        self.pup_diam_m = 1.0
        self.wfs_fov_arcsec = 2.0

        # The scale factor applied inside _compute_slopes_from_derivatives
        # coeff = 1e-9 / (pup_diam_m / wfs_nsubaps) * 206265 * (1 / (0.5 * wfs_fov_arcsec))
        self.coeff = 1e-9 / (self.pup_diam_m / self.wfs_nsubaps) * 206265
        self.coeff *= 1 / (0.5 * self.wfs_fov_arcsec)

        # Expected final slope values
        # In synim, the theoretical gradient is scaled by N (pixels per subaperture)
        # during the rebin process, and then by the physical coefficient.
        self.exp_x = self.Gx * self.N * self.coeff
        self.exp_y = self.Gy * self.N * self.coeff

    def _simulate_pipeline(self, phase, mask, method):
        """
        Runs the exact data flow used in synim.py to compute final slopes.
        """
        # 1. Compute derivatives using the selected engine
        if method == 'gtilt':
            raw_gtilt_x, raw_gtilt_y = compute_gtilt_with_extrapolation(
                data=phase, mask=mask, wfs_nsubaps=self.wfs_nsubaps
            )
            slopes = _compute_slopes_from_gtilt(
                raw_gtilt_x, raw_gtilt_y, mask, mask,
                self.wfs_nsubaps, self.wfs_fov_arcsec, self.pup_diam_m,
                None, False, False  # specula_convention=False
            )
        else:
            dx, dy = compute_derivatives_with_extrapolation(phase, mask)
            slopes = _compute_slopes_from_derivatives(
                dx, dy, mask, mask,
                self.wfs_nsubaps, self.wfs_fov_arcsec, self.pup_diam_m, None,
                False, False  # specula_convention=False
            )

        # 2. Extract final physical slopes
        # FIX: Flatten the output array and cast to standard Python floats
        slopes_flat = np.array(slopes).flatten()
        tilt_x = float(slopes_flat[0])
        tilt_y = float(slopes_flat[1])

        return tilt_x, tilt_y

    def run_comparison(self, mask_name, mask):
        """Helper to run both methods and print the comparison."""
        print(f"\n--- {mask_name} ---")
        tx_g, ty_g = self._simulate_pipeline(self.phase, mask, 'gtilt')
        tx_d, ty_d = self._simulate_pipeline(self.phase, mask, 'derivatives')

        print(f"Expected   : X={self.exp_x:.3e}, Y={self.exp_y:.3e}")
        print(f"G-Tilt     : X={tx_g:.3e}, Y={ty_g:.3e}")
        print(f"Derivatives: X={tx_d:.3e}, Y={ty_d:.3e}")

        # Assertions for G-Tilt (must be mathematically exact if illuminated)
        if tx_g != 0.0:
            np.testing.assert_allclose(tx_g, self.exp_x, err_msg=f"{mask_name} X failed")
        if ty_g != 0.0:
            np.testing.assert_allclose(ty_g, self.exp_y, err_msg=f"{mask_name} Y failed")

        return tx_g, ty_g, tx_d, ty_d

    def test_all_cases(self):
        """Test different sub-aperture partial illumination scenarios."""

        # 1. Full Mask
        self.run_comparison("Full Mask", np.ones((self.N, self.N)))

        # 2. Half Mask (e.g., edge of the pupil)
        m2 = np.zeros((self.N, self.N))
        m2[:, :self.N//2] = 1
        self.run_comparison("Half Cut Mask", m2)

        # 3. Central Spider Cross
        m3 = np.ones((self.N, self.N))
        m3[4:6, :] = 0
        m3[:, 4:6] = 0
        self.run_comparison("Spider Cross", m3)

        # 4. Diagonal Spider
        m4 = np.ones((self.N, self.N))
        np.fill_diagonal(m4, 0)
        self.run_comparison("Diagonal Spider", m4)

        # 5. Single Horizontal Line
        m5 = np.zeros((self.N, self.N))
        m5[5, :] = 1
        tx_g, ty_g, tx_d, ty_d = self.run_comparison("Single Horizontal Line", m5)
        self.assertEqual(ty_g, 0.0, "Y slope must be 0 for a single horizontal line")

        # 6. Two Isolated Lines
        m6 = np.zeros((self.N, self.N))
        m6[2, :] = 1
        m6[8, :] = 1
        tx_g, ty_g, tx_d, ty_d = self.run_comparison("Two Isolated Lines", m6)
        self.assertEqual(ty_g, 0.0, "Y slope must be 0 for isolated horizontal lines")

        # 7. Single Pixel
        m7 = np.zeros((self.N, self.N))
        m7[5, 5] = 1
        tx_g, ty_g, tx_d, ty_d = self.run_comparison("Single Pixel", m7)
        self.assertEqual(tx_g, 0.0, "X slope must be 0 for a single pixel")
        self.assertEqual(ty_g, 0.0, "Y slope must be 0 for a single pixel")
