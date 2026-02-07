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
