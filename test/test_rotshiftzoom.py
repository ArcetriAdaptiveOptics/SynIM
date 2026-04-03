import unittest
import numpy as np
from synim.utils import rotshiftzoom_array
import matplotlib.pyplot as plt

import specula
specula.init(device_idx=-1, precision=1)
from specula.lib.make_mask import make_mask

class TestRotShiftZoomArray(unittest.TestCase):
    def setUp(self):
        # Create a circular mask
        self.size = 100
        self.mask = make_mask(self.size, obsratio=0.0, diaratio=1.0)
        self.mask_half = make_mask(self.size, obsratio=0.0, diaratio=0.5)
        # Create a test array: 2D Gaussian centered
        x = np.arange(self.size)
        y = np.arange(self.size)
        xx, yy = np.meshgrid(x, y)
        self.gauss = np.exp(-((xx - self.size//2)**2 + (yy - self.size//2)**2) / (2*10**2))
        self.gauss_masked = self.gauss * self.mask
        # create a stripes pattern
        self.stripes = np.zeros((self.size, self.size), dtype=np.float32)
        self.stripes[:, ::10] = 1.0
        self.stripes_masked = self.stripes * self.mask
        # create a rectangular pattern
        self.rectangular = np.zeros((self.size, self.size), dtype=np.float32)
        self.rectangular[30:70, 40:60] = 1.0

    def test_shift(self):
        # Apply a shift
        shifted = rotshiftzoom_array(
            self.gauss_masked, dm_translation=(10, 5), output_size=(self.size, self.size)
        )
        diff = shifted - np.roll(self.gauss_masked, shift=(10, 5), axis=(0, 1))

        plot_debug = False
        if plot_debug:
            fig, axs = plt.subplots(1, 3, figsize=(12, 4))
            axs[0].imshow(self.gauss_masked, cmap='viridis')
            axs[0].set_title('Original Masked Gaussian')
            axs[1].imshow(shifted, cmap='viridis')
            axs[1].set_title('Shifted Masked Gaussian')
            axs[2].imshow(diff, cmap='bwr')
            axs[2].set_title('Difference')
            plt.show()

        rel_error = np.sum(np.abs(diff)) / np.sum(np.abs(self.gauss_masked))
        self.assertTrue(rel_error < 1e-2)

    def test_magnification(self):
        # Apply magnification
        mag = rotshiftzoom_array(
            self.mask_half, dm_magnification=(2.0, 2.0), output_size=(self.size, self.size)
        )

        plot_debug = False
        if plot_debug:
            fig, axs = plt.subplots(1, 3, figsize=(12, 4))
            axs[0].imshow(self.mask_half, cmap='viridis')
            axs[0].set_title('Original Half Mask')
            axs[1].imshow(mag, cmap='viridis')
            axs[1].set_title('Magnified Mask')
            axs[2].imshow(mag - self.mask, cmap='bwr')
            axs[2].set_title('Difference')
            plt.show()

        # half mask should now cover full size
        expected_area = np.sum(self.mask)
        actual_area = np.sum(mag > 0.5)
        self.assertTrue(abs(expected_area - actual_area)/expected_area <  0.1)

    def test_magnification_and_shift_combined2(self):
        # Apply magnification and shift
        mag_trans = rotshiftzoom_array(
            self.mask,
            wfs_magnification=(0.5, 0.5),
            wfs_translation=(25, 25),
            output_size=(self.size, self.size)
        )

        # Check that the area is roughly 25% of the original
        expected_area = np.sum(self.mask) * 0.25
        actual_area = np.sum(mag_trans > 0.5)
        self.assertTrue(abs(expected_area - actual_area) / expected_area < 0.1, 
                        "Area non corretta nel test combinato!")

        # Check that the new center is shifted by EXACTLY +25, +25
        y_indices, x_indices = np.where(mag_trans > 0.5)
        self.assertAlmostEqual(np.median(y_indices), self.size // 2 + 25, delta=1.0, 
                               msg="Combined WFS shift errato in Y!")
        self.assertAlmostEqual(np.median(x_indices), self.size // 2 + 25, delta=1.0, 
                               msg="Combined WFS shift errato in X!")

    def test_anamorphosis_45(self):
        # Apply 45° anamorphosis
        anam45 = rotshiftzoom_array(
            self.mask, wfs_anamorphosis_45=1.5, output_size=(self.size, self.size)
        )
        # Check 180° rotational symmetry
        rot180 = np.rot90(anam45, 2)

        plot_debug = False
        if plot_debug:
            fig, axs = plt.subplots(1, 2, figsize=(12, 4))
            axs[0].imshow(self.mask, cmap='viridis')
            axs[0].set_title('Original Masked Gaussian')
            axs[1].imshow(anam45, cmap='viridis')
            axs[1].set_title('Anamorphosis 45°')
            fig, axs = plt.subplots(1, 3, figsize=(8, 4))
            im = axs[0].imshow(anam45, cmap='viridis')
            axs[0].set_title('Anamorphosis 45°')
            fig.colorbar(im, ax=axs[0])
            im = axs[1].imshow(rot180, cmap='viridis')
            axs[1].set_title('Rotated 180°')
            fig.colorbar(im, ax=axs[1])
            diff = anam45 - rot180
            im = axs[2].imshow(diff, cmap='bwr')
            axs[2].set_title('Difference')
            fig.colorbar(im, ax=axs[2])
            plt.show()

        self.assertTrue(np.allclose(anam45, rot180, atol=1e-2))
        # Check that horizontal and vertical flips are NOT symmetric
        flip_h = np.flip(anam45, axis=1)
        flip_v = np.flip(anam45, axis=0)
        self.assertFalse(np.allclose(anam45, flip_h, atol=1e-2))
        self.assertFalse(np.allclose(anam45, flip_v, atol=1e-2))

    def test_anamorphosis_and_rotation_combined_vs_separate(self):
        # Test parameters
        anam45_factor = 1.5
        rotation_deg = 30

        mask = self.mask.astype(np.float32)

        # Combined application
        combined = rotshiftzoom_array(
            mask,
            wfs_anamorphosis_45=anam45_factor,
            wfs_rotation=rotation_deg,
            output_size=(self.size, self.size)
        )

        # Separate application in two steps: first anamorphosis, then rotation
        step1 = rotshiftzoom_array(
            self.mask,
            wfs_anamorphosis_45=anam45_factor,
            output_size=(self.size, self.size)
        )
        separate = rotshiftzoom_array(
            step1,
            wfs_rotation=rotation_deg,
            output_size=(self.size, self.size)
        )

        # Comparison
        diff = combined - separate
        rel_error = np.sum(np.abs(diff)) / np.sum(np.abs(combined))

        plot_debug = False
        if plot_debug:
            fig, axs = plt.subplots(1, 3, figsize=(12, 4))
            im = axs[0].imshow(combined, cmap='viridis')
            axs[0].set_title('Combined')
            fig.colorbar(im, ax=axs[0])
            im = axs[1].imshow(separate, cmap='viridis')
            axs[1].set_title('Separate')
            fig.colorbar(im, ax=axs[1])
            im = axs[2].imshow(diff, cmap='bwr')
            axs[2].set_title('Difference')
            fig.colorbar(im, ax=axs[2])
            plt.show()

        # The results should be very similar
        self.assertTrue(rel_error < 2e-2)

    def test_anamorphosis_90_preserves_vertical_stripes(self):

        # Apply 90° anamorphosis to vertical stripes
        anam90 = rotshiftzoom_array(
            self.stripes, dm_magnification=(0.8, 1.0), output_size=(self.size, self.size)
        )

        # Calculate average profiles over columns to check stripe orientation
        original_profile = np.mean(self.stripes, axis=0)
        anam_profile = np.mean(anam90, axis=0) / 0.8

        plot_debug = False
        if plot_debug:
            fig, ax = plt.subplots(1, 2, figsize=(8, 4))
            im = ax[0].imshow(self.stripes, cmap='viridis')
            ax[0].set_title('Original Stripes Pattern')
            fig.colorbar(im, ax=ax[0])
            im = ax[1].imshow(anam90, cmap='viridis')
            ax[1].set_title('Anamorphosed Stripes Pattern')
            fig.colorbar(im, ax=ax[1])
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(original_profile, label='Original Stripes Profile')
            ax.plot(anam_profile, label='Anamorphosed Stripes Profile')
            ax.set_title('Profile Comparison')
            ax.set_xlabel('Column Index')
            ax.set_ylabel('Average Intensity')
            ax.legend()
            plt.show()

        # The stripe pattern should be preserved
        self.assertTrue(np.allclose(original_profile, anam_profile, atol=0.1))

    def test_rotation(self):
        # Apply a 90° rotation
        rotated = rotshiftzoom_array(
            self.rectangular, wfs_rotation=90, output_size=(self.size, self.size)
        )
        # The mask rotated by 90° should be equal to the transpose and vertical flip
        expected = np.rot90(self.rectangular, k=1)
        diff = rotated - np.roll(expected, shift=(1, 0), axis=(0, 1))

        plot_debug = False
        if plot_debug:
            fig, axs = plt.subplots(1, 3, figsize=(12, 4))
            axs[0].imshow(self.rectangular, cmap='gray')
            axs[0].set_title('Original Rectangular')
            axs[1].imshow(rotated, cmap='gray')
            axs[1].set_title('Rotated 90°')
            axs[2].imshow(diff, cmap='bwr')
            axs[2].set_title('Difference')
            plt.show()

        self.assertTrue(np.allclose(rotated,
                        np.roll(expected, shift=(1, 0), axis=(0, 1)), atol=1e-2))

    def test_wfs_translation_independence(self):
        """
        Verifies that wfs_translation represents a pure camera shift on the 
        focal plane and is not erroneously scaled by wfs_magnification.
        """
        trans_y, trans_x = 20, 0
        mag = 2.0  # Strong WFS zoom

        transformed = rotshiftzoom_array(
            self.gauss_masked,
            wfs_translation=(trans_y, trans_x),
            wfs_magnification=(mag, mag),
            output_size=(self.size, self.size)
        )

        peak_y, peak_x = np.unravel_index(np.argmax(transformed), transformed.shape)
        center_y, center_x = self.size // 2, self.size // 2

        expected_y = center_y + trans_y
        expected_x = center_x + trans_x

        self.assertEqual(peak_y, expected_y,
                         "WFS Y translation was corrupted by WFS magnification!")
        self.assertEqual(peak_x, expected_x,
                         "WFS X translation was corrupted by WFS magnification!")

    def test_dm_translation_independence_from_magnification(self):
        """
        Verifies that dm_translation (e.g., LGS off-axis shift) 
        is strictly independent of the Cone Effect (dm_magnification).
        """
        trans_y, trans_x = 20, -10
        mag = 0.5  # Strong cone effect (shrink)

        transformed = rotshiftzoom_array(
            self.gauss_masked,
            dm_translation=(trans_y, trans_x),
            dm_magnification=(mag, mag),
            output_size=(self.size, self.size)
        )

        peak_y, peak_x = np.unravel_index(np.argmax(transformed), transformed.shape)
        center_y, center_x = self.size // 2, self.size // 2

        expected_y = int(center_y + trans_y * mag)
        expected_x = int(center_x + trans_x * mag)

        self.assertEqual(peak_y, expected_y,
                         "DM Y translation was corrupted by Cone Effect!")
        self.assertEqual(peak_x, expected_x,
                         "DM X translation was corrupted by Cone Effect!")

    def test_dm_translation_invariance_under_rotation(self):
        """
        Verifies that the LGS footprint shift in the sky remains invariant 
        even if the DM rotates around the optical axis.
        """
        trans_y, trans_x = 30, 0
        rot_deg = 90.0

        # We use the asymmetric rectangular mask to verify BOTH the fixed shift and the rotation.
        transformed = rotshiftzoom_array(
            self.rectangular,
            dm_translation=(trans_y, trans_x),
            dm_rotation=rot_deg,
            output_size=(self.size, self.size)
        )

        expected_center_y = self.size // 2 + trans_y
        expected_center_x = self.size // 2 + trans_x

        y_indices, x_indices = np.where(transformed > 0.5)

        # Check that the center hasn't moved
        self.assertAlmostEqual(np.median(y_indices), expected_center_y, delta=1.0, 
                               msg="DM center shifted incorrectly due to rotation!")
        self.assertAlmostEqual(np.median(x_indices), expected_center_x, delta=1.0, 
                               msg="DM center shifted incorrectly due to rotation!")

        # Check that the shape actually rotated
        width_y = np.max(y_indices) - np.min(y_indices)
        width_x = np.max(x_indices) - np.min(x_indices)
        self.assertTrue(width_x > 35 and width_y < 25,
                        "Rectangle did not physically rotate!")

    def test_wfs_translation_independence_from_rotation(self):
        """
        Verifies that WFS translation (camera physical shift) is an absolute shift 
        on the final detector, strictly independent of the WFS optical rotation.
        """
        trans_y, trans_x = 20, 0
        rot_deg = 45.0

        transformed = rotshiftzoom_array(
            self.gauss_masked,
            wfs_translation=(trans_y, trans_x),
            wfs_rotation=rot_deg,
            output_size=(self.size, self.size)
        )

        peak_y, peak_x = np.unravel_index(np.argmax(transformed), transformed.shape)
        center_y, center_x = self.size // 2, self.size // 2

        self.assertEqual(peak_y, center_y + trans_y, "WFS Y translation corrupted by rotation!")
        self.assertEqual(peak_x, center_x + trans_x, "WFS X translation corrupted by rotation!")

    def test_3d_cube_integrity(self):
        """
        Verifies that 3D arrays (representing DM modes) are transformed 
        correctly slice by slice, without any cross-talk along the Z-axis.
        """
        cube = np.zeros((self.size, self.size, 3), dtype=np.float32)
        cube[49:52, 49:52, 0] = 1.0  # Center
        cube[19:22, 19:22, 1] = 1.0  # Top-Left
        cube[79:82, 79:82, 2] = 1.0  # Bottom-Right

        transformed_cube = rotshiftzoom_array(
            cube,
            dm_translation=(10, -10),
            dm_magnification=(0.8, 0.8),
            wfs_rotation=45.0,
            output_size=(self.size, self.size)
        )

        # 1. Z-Axis Integrity: Check that the block core survived interpolation
        for mode in range(3):
            self.assertTrue(np.max(transformed_cube[:, :, mode]) > 0.5, 
                            f"Energy heavily diluted or lost in mode {mode}!")

        # 2. Cross-talk check: The peak of mode 1 MUST NOT bleed into mode 0
        peak_y_m1, peak_x_m1 = np.unravel_index(
            np.argmax(transformed_cube[:, :, 1]), (self.size, self.size)
        )
        self.assertAlmostEqual(transformed_cube[peak_y_m1, peak_x_m1, 0], 0.0, places=5,
                               msg="Cross-talk detected between DM modes!")
