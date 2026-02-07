import unittest
import numpy as np
import os
from astropy.io import fits

import specula
specula.init(device_idx=-1, precision=1)
from specula.data_objects.ifunc import IFunc
from specula.lib.modal_base_generator import compute_ifs_covmat


class TestCovariance(unittest.TestCase):
    """Test covariance matrix computation against reference data."""

    def setUp(self):
        """Set up test parameters."""
        self.nmodes = 100
        self.npixels = 128
        self.diameter = 8.0  # meters
        self.r0 = 0.2  # meters
        self.L0 = 50.0  # meters

        # Reference file path
        test_dir = os.path.dirname(os.path.abspath(__file__))
        self.ref_file = os.path.join(
            test_dir,
            'data',
            'cov_zern_100m_128p_8m_0.2r0_50L0.fits'
        )

    def test_zernike_covariance_vs_reference(self):
        """
        Test that computed Zernike covariance matches reference file.
        
        Computes covariance matrix for 100 Zernike modes on 128x128 pixels
        with D=8m, r0=0.2m, L0=50m and compares with reference FITS file.
        """
        print(f"\n{'='*70}")
        print(f"Testing Zernike Covariance vs Reference File")
        print(f"{'='*70}")
        print(f"Parameters:")
        print(f"  Modes: {self.nmodes}")
        print(f"  Pixels: {self.npixels}")
        print(f"  Diameter: {self.diameter} m")
        print(f"  r0: {self.r0} m")
        print(f"  L0: {self.L0} m")
        print(f"{'='*70}\n")

        # Generate Zernike influence functions
        print("Generating Zernike modes...")
        ifunc = IFunc(
            type_str='zernike',
            nmodes=self.nmodes,
            npixels=self.npixels,
            obsratio=0.0,
            diaratio=1.0,
            precision=1,
            target_device_idx=-1
        )

        save_ifunc = False
        if save_ifunc: # pragma: no cover
            ifunc.save(f"/Users/guido/GitHub/SynIM/test/calib/ifunc/"
                       f"zernike_{self.nmodes}m_{self.npixels}px.ifunc")

        # Check if reference file exists
        if not os.path.exists(self.ref_file):
            self.skipTest(f"Reference file not found: {self.ref_file}")

        pupil_mask = ifunc.mask_inf_func.astype(np.float32)
        z_if_3d = ifunc.ifunc_2d_to_3d(normalize=True)

        # Flatten inside pupil
        idx = np.where(pupil_mask.ravel() > 0.5)[0]
        npupil = idx.size
        influence_functions = np.zeros((self.nmodes, npupil), dtype=np.float32)

        for k in range(self.nmodes):
            influence_functions[k, :] = z_if_3d[:, :, k].ravel()[idx]

        print(f"  ✓ Generated {self.nmodes} Zernike modes")
        print(f"  ✓ Valid pupil pixels: {npupil}")

        # Compute covariance matrix
        print("\nComputing covariance matrix...")
        cov_computed = compute_ifs_covmat(
            pupil_mask,
            self.diameter,
            influence_functions,
            self.r0,
            self.L0,
            oversampling=2,
            xp=np,
            dtype=np.float32
        )

        print(f"  ✓ Computed covariance shape: {cov_computed.shape}")
        print(f"  ✓ Diagonal RMS (rad): {np.sqrt(np.diag(cov_computed)).mean():.6f}")

        # Load reference covariance
        print(f"\nLoading reference from: {os.path.basename(self.ref_file)}")
        with fits.open(self.ref_file) as hdul:
            cov_reference = hdul[0].data.astype(np.float32)

            # Print header info if available
            header = hdul[0].header
            if 'R0' in header:
                print(f"  Reference r0: {header['R0']} m")
            if 'L0' in header:
                print(f"  Reference L0: {header['L0']} m")
            if 'DIAMM' in header:
                print(f"  Reference diameter: {header['DIAMM']} m")

        print(f"  ✓ Loaded reference shape: {cov_reference.shape}")
        print(f"  ✓ Diagonal RMS (rad): {np.sqrt(np.diag(cov_reference)).mean():.6f}")

        # Compare shapes
        self.assertEqual(cov_computed.shape, cov_reference.shape,
                        "Covariance matrix shapes don't match")

        # Compare values
        print(f"\nComparing matrices...")

        # Absolute difference
        diff = cov_computed - cov_reference
        abs_diff = np.abs(diff)
        max_abs_diff = np.max(abs_diff)
        mean_abs_diff = np.mean(abs_diff)

        # Relative difference (avoid division by zero)
        rel_diff = np.abs((cov_computed - cov_reference) /
                         (np.abs(cov_reference) + 1e-10))
        max_rel_diff = np.max(rel_diff)
        mean_rel_diff = np.mean(rel_diff)

        print(f"  Maximum absolute difference: {max_abs_diff:.2e} rad²")
        print(f"  Mean absolute difference: {mean_abs_diff:.2e} rad²")
        print(f"  Maximum relative difference: {max_rel_diff*100:.2f}%")
        print(f"  Mean relative difference: {mean_rel_diff*100:.2f}%")

        # RMS
        rms_diff = np.sqrt(np.mean(diff**2))
        rms_rel_diff = rms_diff/np.sqrt(np.mean(cov_reference**2))

        print(f"  RMS absolute difference: {rms_diff:.2e} rad²")
        print(f"  RMS relative difference: {rms_rel_diff*100:.2f}%")

        # Compare diagonal elements specifically
        diag_computed = np.diag(cov_computed)
        diag_reference = np.diag(cov_reference)
        diag_diff = diag_computed - diag_reference
        diag_rel_diff = np.abs(diag_diff / diag_reference)

        print(f"\n  Diagonal comparison:")
        print(f"    Max relative difference: {np.max(diag_rel_diff)*100:.2f}%")
        print(f"    Mean relative difference: {np.mean(diag_rel_diff)*100:.2f}%")

        # RMS
        rms_diag_diff = np.sqrt(np.mean(diag_diff**2))
        rms_diag_rel_diff = rms_diag_diff / np.sqrt(np.mean(diag_reference**2))

        print(f"    RMS absolute difference: {rms_diag_diff:.2e} rad²")
        print(f"    RMS relative difference: {rms_diag_rel_diff*100:.2f}%")

        # Visual comparison (optional)
        plot_debug = False
        if plot_debug: # pragma: no cover
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(2, 3, figsize=(15, 10))

            # Computed covariance
            im0 = axes[0, 0].imshow(cov_computed, cmap='viridis')
            axes[0, 0].set_title('Computed Covariance')
            plt.colorbar(im0, ax=axes[0, 0])

            # Reference covariance
            im1 = axes[0, 1].imshow(cov_reference, cmap='viridis')
            axes[0, 1].set_title('Reference Covariance')
            plt.colorbar(im1, ax=axes[0, 1])

            # Absolute difference
            im2 = axes[0, 2].imshow(abs_diff, cmap='hot')
            axes[0, 2].set_title(f'Abs Difference (max={max_abs_diff:.2e})')
            plt.colorbar(im2, ax=axes[0, 2])

            # Diagonal comparison
            axes[1, 0].semilogy(diag_computed, 'b-', label='Computed', alpha=0.7)
            axes[1, 0].semilogy(diag_reference, 'r--', label='Reference', alpha=0.7)
            axes[1, 0].set_xlabel('Mode index')
            axes[1, 0].set_ylabel('Variance (rad²)')
            axes[1, 0].set_title('Diagonal Elements')
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)

            # Relative difference
            im3 = axes[1, 1].imshow(rel_diff * 100, cmap='hot', vmin=0, vmax=5)
            axes[1, 1].set_title('Relative Difference (%)')
            plt.colorbar(im3, ax=axes[1, 1], label='%')

            # Histogram of relative differences
            axes[1, 2].hist(rel_diff.ravel() * 100, bins=50, edgecolor='black')
            axes[1, 2].set_xlabel('Relative Difference (%)')
            axes[1, 2].set_ylabel('Count')
            axes[1, 2].set_title('Distribution of Relative Errors')
            axes[1, 2].axvline(mean_rel_diff * 100, color='r',
                             linestyle='--', label=f'Mean={mean_rel_diff*100:.2f}%')
            axes[1, 2].legend()
            axes[1, 2].grid(True, alpha=0.3)

            plt.tight_layout()
            plt.show()

        # Assertions with reasonable tolerances
        # Allow up to 1% on the diagonal elements
        self.assertLess(np.max(diag_rel_diff), 0.01,
             f"Max diagonal relative difference {np.max(diag_rel_diff)*100:.2f}% exceeds 1%")

        # Allow up to 1/1000 RMS relative error due to numerical differences
        self.assertLess(rms_rel_diff, 0.001,
             f"RMS relative difference {rms_rel_diff*100:.2f}% exceeds 0.1%")

        # Check that most elements are very close
        close_elements = np.sum(abs_diff < 0.01*cov_computed.mean())
        total_elements = rel_diff.size
        close_fraction = close_elements / total_elements

        print(f"\n  Elements within 1% tolerance: {close_fraction*100:.1f}%")

        self.assertGreater(close_fraction, 0.95,
                          f"Only {close_fraction*100:.1f}% of elements within 1% tolerance")

        print(f"\n{'='*70}")
        print(f"✓ Covariance matrix matches reference within tolerances")
        print(f"{'='*70}\n")


if __name__ == '__main__':
    unittest.main()
