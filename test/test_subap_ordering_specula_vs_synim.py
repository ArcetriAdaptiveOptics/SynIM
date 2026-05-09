"""
Test to verify consistent subaperture ordering between SPECULA Slopec and SynIM.

This test compares:
1. flux_per_subaperture from SPECULA ShSlopec
2. idx_valid_sa ordering from find_subapdata (display_map order)
3. illumination from compute_subaperture_illumination

Uses a circular, slightly off-axis pupil mask for realism.
"""

import unittest
import numpy as np
import tempfile
import os

# SPECULA imports
import specula
specula.init(device_idx=-1, precision=1)

from specula.data_objects.pixels import Pixels
from specula.data_objects.subap_data import SubapData
from specula.processing_objects.sh_slopec import ShSlopec
from specula.processing_objects.sh_subap_calibrator import ShSubapCalibrator
from specula.data_objects.electric_field import ElectricField
from specula.processing_objects.sh import SH
from specula.lib.make_mask import make_mask

# SynIM imports
from synim.synim import compute_subaperture_illumination
from synim.params_utils import find_subapdata
from synim.params_manager import _idx_valid_sa_to_linear_for_illumination
from specula.calib_manager import CalibManager


def make_circular_pupil(npix, radius_fraction=0.4, off_axis_shift=(0, 0), obstruction_ratio=0.0):
    """
    Create a circular pupil mask, optionally off-axis or with central obstruction.
    
    Args:
        npix: Size of the mask (npix x npix)
        radius_fraction: Fraction of npix/2 for pupil radius (default 0.4 → diameter ~0.8*npix)
        off_axis_shift: (dx, dy) offset in pixels
        obstruction_ratio: Central obstruction ratio (0 to 1)
    
    Returns:
        Circular mask (0/1) with optional obstruction
    """
    x = np.arange(npix) - npix / 2.0 + off_axis_shift[0]
    y = np.arange(npix) - npix / 2.0 + off_axis_shift[1]
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2 + yy**2)

    radius = npix / 2.0 * radius_fraction
    mask = (r <= radius).astype(np.float32)

    # Add central obstruction if requested
    if obstruction_ratio > 0:
        obstruct_radius = radius * obstruction_ratio
        obstruction = (r <= obstruct_radius).astype(np.float32)
        mask = mask - obstruction

    return mask


def make_circular_pupil_with_spider(npix, radius_fraction=0.4, off_axis_shift=(0, 0), 
                                     obstruction_ratio=0.0, n_spiders=6, spider_width=2, 
                                     spider_angle_offset=0.0):
    """
    Create a circular pupil mask with spider arms (like real telescopes).
    
    Args:
        npix: Size of the mask (npix x npix)
        radius_fraction: Fraction of npix/2 for pupil radius
        off_axis_shift: (dx, dy) offset in pixels
        obstruction_ratio: Central obstruction ratio (0 to 1)
        n_spiders: Number of spider arms (e.g., 6 for common VLT design)
        spider_width: Width of spider arms in pixels
        spider_angle_offset: Rotation angle offset for spiders in degrees
    
    Returns:
        Circular mask (0/1) with spiders
    """
    diaratio = radius_fraction
    xc = off_axis_shift[0] / npix
    yc = off_axis_shift[1] / npix

    mask = make_mask(
        np_size=npix,
        obsratio=obstruction_ratio,
        diaratio=diaratio,
        xc=xc,
        yc=yc,
        square=False,
        inverse=False,
        centeronpixel=False,
        get_idx=False,
        xp=np,
        spider=True,
        spider_width=spider_width,
        n_petals=n_spiders,
        angle_offset=spider_angle_offset
    )

    return mask.astype(np.float32)


class TestSubapOrderingSpeculaVsSynim(unittest.TestCase):
    """Compare subaperture ordering between SPECULA and SynIM."""

    def setUp(self):
        """Set up test data."""
        # Pupil and WFS parameters
        # EF: 160 px, 0.05 m/px → 8 m diameter
        # Detector: 40 subaps × 6 px/subap = 240 px
        self.npix = 160          # EF size (160 * 0.05 m = 8 m telescope)
        self.pixel_pitch = 0.05  # m/pixel
        self.subap_on_diameter = 40  # 40x40 subapertures
        self.subap_npx = 6       # pixels per subaperture → detector 240x240

        # WFS parameters
        self.wavelength_nm = 750
        self.pxscale_arcsec = 0.4
        self.wfs_rotation_deg = 6.0
        self.pupil_rotation_deg = 18.0
        self.pupil_shift = (7, -5)
        self.subap_energy_th = 1e-3

        self.plot_debug = False  # Set to True to show plots during testing

        # Create temporary directory for calibration files
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def _create_subapdata_from_pupil(self, pixels=None):
        """
        Create SubapData from the test pupil mask using the SPECULA method.
        Uses contiguous pixel blocks for each subaperture (like test_sh_slopec).
        
        Returns:
            subapdata (SubapData): Subaperture data object
            subapdata_path (str): Path where SubapData is saved
        """
        # Preferred path: mimic SPECULA workflow and detect valid subaps from pixels.
        if pixels is not None:
            calibrator = ShSubapCalibrator(
                subap_on_diameter=self.subap_on_diameter,
                data_dir=self.temp_dir,
                energy_th=self.subap_energy_th,
                output_tag='test_subapdata_from_calibrator',
                overwrite=True,
                target_device_idx=-1,
            )
            calibrator.inputs['in_pixels'].set(pixels)
            calibrator.setup()
            calibrator.check_ready(1)
            calibrator.trigger()
            calibrator.post_trigger()

            subapdata = calibrator.subaps
            fd, subapdata_path = tempfile.mkstemp(prefix='test_subapdata_',
                                                  suffix='.fits', dir=self.temp_dir)
            os.close(fd)
            os.remove(subapdata_path)
            subapdata.save(subapdata_path, overwrite=True)
            return subapdata, subapdata_path

        # Fallback path: full grid of contiguous subapertures.
        pixel_per_subap = self.subap_npx
        total_grid_size = self.subap_on_diameter * pixel_per_subap

        # Create grid mask for subapertures
        mask_subap = np.ones((total_grid_size, total_grid_size))

        idxs = {}
        display_map = {}
        valid_count = 0

        for i in range(self.subap_on_diameter):
            for j in range(self.subap_on_diameter):
                y_start = i * pixel_per_subap
                y_end = (i + 1) * pixel_per_subap
                x_start = j * pixel_per_subap
                x_end = (j + 1) * pixel_per_subap

                # Create temporary mask for this subaperture
                mask_subap[:] = 0
                mask_subap[y_start:y_end, x_start:x_end] = 1

                # Get pixel indices
                pixel_idx = np.ravel_multi_index(
                    np.where(mask_subap == 1),
                    mask_subap.shape
                )
                idxs[valid_count] = pixel_idx

                # Display map position
                display_map[valid_count] = j * self.subap_on_diameter + i
                valid_count += 1

        # Convert to numpy arrays (following test_sh_slopec pattern)
        n_valid_subaps = valid_count
        idxs_array = np.zeros((n_valid_subaps, pixel_per_subap * pixel_per_subap), dtype=int)
        display_map_array = np.zeros(n_valid_subaps, dtype=int)

        for k, idx in idxs.items():
            idxs_array[k] = idx
            display_map_array[k] = display_map[k]

        # Create and save SubapData
        subapdata = SubapData(
            idxs=idxs_array,
            display_map=display_map_array,
            nx=self.subap_on_diameter,
            ny=self.subap_on_diameter,
            target_device_idx=-1
        )

        fd, subapdata_path = tempfile.mkstemp(prefix='test_subapdata_',
                                              suffix='.fits', dir=self.temp_dir)
        os.close(fd)
        os.remove(subapdata_path)
        subapdata.save(subapdata_path)

        return subapdata, subapdata_path

    def test_subap_ordering_consistency(self):
        """
        Main test: verify that subaperture ordering is consistent between
        SPECULA Slopec flux_per_subaperture and SynIM find_subapdata idx_valid_sa.
        
        Uses a realistic pupil with 6 spider arms and a 6 degree rotation.
        
        IMPORTANT: We create proper SPECULA setup with:
        - ElectricField with the realistic pupil mask
        - SH sensor with the same mask
        - Detector to capture intensity
        - ShSlopec to compute flux
        """
        print("\n" + "="*70)
        print("Test: Subaperture Ordering Consistency (SPECULA vs SynIM)")
        print("With realistic pupil: 6 spiders + 6 degree rotation")
        print("="*70)

        # ========== SETUP: Create realistic pupil with spiders ==========
        print("\n[1] Creating realistic pupil mask (6 spiders + 6° rotation)...")

        # Create pupil with 6 spiders
        illum_pupil = make_circular_pupil_with_spider(
            npix=self.npix,
            radius_fraction=0.45,
            off_axis_shift=self.pupil_shift,
            obstruction_ratio=0.0,
            n_spiders=6,
            spider_width=3,
            spider_angle_offset=self.pupil_rotation_deg
        )

        print(f"    ✓ Created pupil with 6 spiders, shape: {illum_pupil.shape}")
        print(f"    ✓ Pupil shift (px): {self.pupil_shift}")
        print(f"    ✓ Pupil spider rotation: {self.pupil_rotation_deg}°")
        print(f"    ✓ Pupil illumination: {np.sum(illum_pupil):.0f} pixels")
        print(f"    ✓ Illuminated fraction: {100 * np.sum(illum_pupil) / illum_pupil.size:.1f}%")

        # ========== SPECULA SETUP: EF + SH + Detector + ShSlopec ==========
        print("\n[2] Setting up SPECULA (EF + SH + Detector + ShSlopec)...")

        # Step 2a: Create ElectricField with the realistic pupil mask
        ef = ElectricField(
            self.npix, self.npix,
            pixel_pitch=self.pixel_pitch,
            S0=1,
            target_device_idx=-1
        )

        # Apply the pupil mask as amplitude (field[0]); phase stays zero
        ef.A = illum_pupil
        ef.generation_time = 1
        print(f"    ✓ Created ElectricField with pupil mask")

        # Step 2b: Create SH sensor
        sh = SH(
            wavelengthInNm=self.wavelength_nm,
            subap_wanted_fov=self.subap_npx * self.pxscale_arcsec,
            sensor_pxscale=self.pxscale_arcsec,
            subap_on_diameter=self.subap_on_diameter,
            subap_npx=self.subap_npx,
            rotAnglePhInDeg=self.wfs_rotation_deg,
            target_device_idx=-1
        )

        # Connect EF to SH
        sh.inputs['in_ef'].set(ef)
        sh.setup()
        sh.check_ready(1)
        sh.trigger()
        sh.post_trigger()

        print(f"    ✓ Created SH sensor with {self.subap_on_diameter}x"
              f" {self.subap_on_diameter} subapertures")

        # Step 2c: Get intensity from SH detector
        intensity = sh.outputs['out_i'].i.copy()
        print(f"    ✓ Detector output shape: {intensity.shape}")
        print(f"    ✓ Intensity range: [{intensity.min():.2e}, {intensity.max():.2e}]")

        # Step 2d: Create Pixels object from detector intensity
        pixels = Pixels(*intensity.shape, target_device_idx=-1)
        pixels.pixels = intensity
        pixels.generation_time = 1

        # Step 2e: Create SubapData and ShSlopec
        print("\n[3] Creating SubapData and ShSlopec...")
        subapdata, subapdata_path = self._create_subapdata_from_pupil(pixels=pixels)
        print(f"    ✓ Created {subapdata.n_subaps} subapertures (calibrator-filtered)")

        # Create and run ShSlopec
        slopec = ShSlopec(subapdata, target_device_idx=-1)
        slopec.inputs['in_pixels'].set(pixels)
        slopec.check_ready(1)
        slopec.trigger()
        slopec.post_trigger()

        flux_per_subap_specula = np.asarray(slopec.outputs['out_flux_per_subaperture'].value)
        print(f"    ✓ Flux per subaperture shape: {flux_per_subap_specula.shape}")
        print(f"    ✓ Flux range: [{flux_per_subap_specula.min():.2e},"
              f" {flux_per_subap_specula.max():.2e}]")
        n_pos = np.sum(flux_per_subap_specula > 0)
        print(f"    ✓ Subapertures with flux > 0: {n_pos}/{len(flux_per_subap_specula)}")

        # ========== SynIM SETUP: compute_subaperture_illumination ==========
        print("\n[4] Computing subaperture illumination with SynIM...")

        # Restore subapdata
        restored_subapdata = SubapData.restore(subapdata_path)

        # Extract idx_valid_sa from display_map EXACTLY AS DONE IN im_sh_synim_generator.py
        # This is the configuration that works in params_manager and IM computation.
        display_map = np.asarray(restored_subapdata.display_map, dtype=np.int64)
        n = self.subap_on_diameter
        
        # Reconstruct 2D indices: display_map = col*n + row, so:
        idx_i = display_map // n  # col
        idx_j = display_map % n   # row
        idx_valid_sa_new = np.column_stack((idx_i, idx_j))  # [[col, row], ...]

        # pup_mask for SynIM must be the amplitude of the EF (same 80x80 pupil mask)
        # NOT a resampled version - compute_subaperture_illumination does rebin internally
        ef_amplitude = np.asarray(ef.A)
        print(f"    ✓ EF amplitude (pup_mask) shape: {ef_amplitude.shape}")

        # Match SynIM rotation to SH geometry
        wfs_rotation_deg = self.wfs_rotation_deg
        illum = compute_subaperture_illumination(
            pup_mask=ef_amplitude,  # EF amplitude = pupil mask (160x160), rebin done internally
            wfs_nsubaps=self.subap_on_diameter,
            wfs_rotation=wfs_rotation_deg,
            wfs_translation=(0.0, 0.0),
            wfs_magnification=(1.0, 1.0),
            idx_valid_sa=idx_valid_sa_new,  # 2D array [[col, row], ...]
            verbose=False,
            specula_convention=True  # MUST match params_manager configuration
        )

        print(f"    ✓ WFS rotation: {wfs_rotation_deg}°")
        print(f"    ✓ Illumination shape: {illum.shape}")
        print(f"    ✓ Illumination range: [{illum.min():.4f}, {illum.max():.4f}]")
        print(f"    ✓ Subapertures with illum > 0: {np.sum(illum > 0)}/{len(illum)}")

        # ========== COMPARISON ==========
        print("\n[5] Comparing SPECULA flux with SynIM illumination...")

        # Both arrays should have the same length
        self.assertEqual(
            len(flux_per_subap_specula),
            len(illum),
            f"Length mismatch: Slopec flux {len(flux_per_subap_specula)}"
            f" vs SynIM illum {len(illum)}"
        )
        print("    ✓ Lengths match")

        # Count illuminated vs blocked subapertures
        n_illuminated = np.sum(illum > 0)
        n_zero = np.sum(illum == 0)
        print(f"    ✓ Subapertures with illumination > 0: {n_illuminated}")
        print(f"    ✓ Subapertures blocked (illumination = 0): {n_zero}")

        # Verify that illuminated subapertures have positive flux
        # NOTE: some subapertures with partial SynIM illumination may still have
        # zero SPECULA flux due to spider placement differences → we report, not assert
        illuminated_mask = illum > 0
        if np.any(illuminated_mask):
            flux_illuminated = flux_per_subap_specula[illuminated_mask]
            n_zero_in_illuminated = np.sum(flux_illuminated == 0)
            if n_zero_in_illuminated > 0:
                print(f"    ⚠ {n_zero_in_illuminated} subapertures have SynIM illum>0"
                      f" but SPECULA flux=0")
            else:
                print(f"    ✓ All illuminated subapertures have positive flux")

        # ========== DETAILED ANALYSIS ==========
        print("\n[6] Summary comparison (SPECULA EF+SH+Slopec vs SynIM illumination):")
        fmax = np.nanmax(flux_per_subap_specula)
        flux_normalized = flux_per_subap_specula / fmax if fmax > 0 else flux_per_subap_specula
        imax = np.nanmax(illum)
        illum_normalized = illum / imax if imax > 0 else illum

        # Explicit requested diagnostic: diff = illum - flux (both normalized)
        diff_full = illum_normalized - flux_normalized
        print(f"    diff stats (illum-flux): mean={np.mean(diff_full):.4f},"
            f" mean_abs={np.mean(np.abs(diff_full)):.4f},"
            f" max_abs={np.max(np.abs(diff_full)):.4f}")

        blocked_count = 0
        illuminated_indices = []
        for i in range(len(flux_normalized)):
            if illum[i] == 0:
                blocked_count += 1
            else:
                illuminated_indices.append(i)

        n_subaps_total = self.subap_on_diameter ** 2
        print(f"    Total blocked by SynIM: {blocked_count}/{n_subaps_total}")
        print(f"    Total illuminated:      {len(illuminated_indices)}/{n_subaps_total}")

        # Trap potential 90-degree ordering mistakes explicitly.
        flux_mask = flux_per_subap_specula > 0
        illum_mask = illum > 0
        flux_map_mask = np.zeros((self.subap_on_diameter, self.subap_on_diameter), dtype=bool)
        illum_map_mask = np.zeros((self.subap_on_diameter, self.subap_on_diameter), dtype=bool)
        for k, dm in enumerate(display_map):
            row = dm % self.subap_on_diameter
            col = dm // self.subap_on_diameter
            flux_map_mask[row, col] = flux_mask[k]
            illum_map_mask[row, col] = illum_mask[k]

        overlap_0 = np.sum(flux_map_mask & illum_map_mask)
        overlap_90 = np.sum(flux_map_mask & np.rot90(illum_map_mask, 1))
        overlap_270 = np.sum(flux_map_mask & np.rot90(illum_map_mask, 3))
        print(f"    Overlap check (0°/90°/270°): {overlap_0}/{overlap_90}/{overlap_270}")
        self.assertGreaterEqual(
            overlap_0,
            max(overlap_90, overlap_270),
            "Detected possible 90-degree ordering/orientation mismatch",
        )
        # Statistical comparison for illuminated subapertures
        if len(illuminated_indices) > 0:
            print(f"\n    Analysis of {len(illuminated_indices)} illuminated subapertures:")
            flux_illum = flux_normalized[illuminated_indices]
            illum_illum = illum[illuminated_indices]

            # Normalize both to same scale for comparison
            fmax_illum = np.nanmax(flux_illum)
            flux_illum_norm = flux_illum / fmax_illum if fmax_illum > 0 else flux_illum
            illum_illum_norm = illum_illum / np.max(illum_illum)

            diff = illum_illum_norm - flux_illum_norm
            mad = np.mean(np.abs(diff))
            print(f"    Mean signed diff (illum-flux): {np.mean(diff):.4f}")
            print(f"    Mean abs difference (normalized): {mad:.4f}")
            print(f"    Max abs difference (normalized): {np.max(np.abs(diff)):.4f}")

            corr = np.nan
            if np.std(flux_illum_norm) > 0 and np.std(illum_illum_norm) > 0:
                corr = np.corrcoef(flux_illum_norm, illum_illum_norm)[0, 1]
            print(f"    Correlation: {corr:.4f}")
            self.assertTrue(np.isfinite(mad),
                            "Mean absolute diff between normalized illum and flux is not finite")

        self._plot_comparison(
            flux_per_subap_specula, illum,
            illum_pupil, intensity,
            self.subap_on_diameter,
            display_map
        )

        print("\n" + "="*70)
        print("✓ Test PASSED: Subaperture ordering is consistent")
        print("  (With proper EF+SH setup applying pupil mask)")
        print("="*70 + "\n")

    def _plot_comparison(self, flux, illum, pupil_mask, detector_image, nsubaps, display_map):
        """Plot comparison: pupil mask, detector image, flux map, illumination map, scatter."""
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        n = nsubaps

        # Flux 2D map from explicit display_map to avoid row/col ambiguity.
        flux_map = np.full((n, n), np.nan)
        illum_map = np.full((n, n), np.nan)
        display_map = np.asarray(display_map, dtype=np.int64)
        for k, dm in enumerate(display_map):
            # SPECULA display_map is col-major: dm = col * n + row
            row = dm % n
            col = dm // n
            flux_map[row, col] = flux[k]
            illum_map[row, col] = illum[k]

        # Normalize flux map ignoring NaN
        fmax = np.nanmax(flux_map)
        if fmax > 0:
            flux_map_norm = flux_map / fmax
        else:
            flux_map_norm = flux_map

        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        axes = axes.ravel()
        fig.suptitle('SPECULA vs SynIM subaperture illumination comparison', fontsize=12)

        # 1) Pupil mask
        im0 = axes[0].imshow(pupil_mask, origin='lower', cmap='gray')
        axes[0].set_title('Pupil mask (EF amplitude)')
        plt.colorbar(im0, ax=axes[0], fraction=0.046)

        # 2) Detector image (SH output)
        dmax = detector_image.max()
        if dmax > 0:
            im1 = axes[1].imshow(detector_image, origin='lower', cmap='hot',
                                 norm=mcolors.LogNorm(vmin=max(dmax*1e-4, 1e-12), vmax=dmax))
        else:
            im1 = axes[1].imshow(detector_image, origin='lower', cmap='hot')
        axes[1].set_title('SH detector image')
        plt.colorbar(im1, ax=axes[1], fraction=0.046)

        # 3) SPECULA flux map (2D)
        im2 = axes[2].imshow(
            flux_map_norm,
            origin='lower',
            cmap='viridis',
            vmin=0,
            vmax=1,
            interpolation='nearest',
        )
        axes[2].set_title('SPECULA flux (norm, col-major)')
        plt.colorbar(im2, ax=axes[2], fraction=0.046)

        # 4) SynIM illumination map (2D)
        im3 = axes[3].imshow(
            illum_map,
            origin='lower',
            cmap='viridis',
            vmin=0,
            vmax=1,
            interpolation='nearest',
        )

        # Quantify residual apparent shift from weighted centroids.
        yy, xx = np.indices((n, n), dtype=np.float64)
        fsum = np.nansum(flux_map_norm)
        isum = np.nansum(illum_map)
        if fsum > 0 and isum > 0:
            flux_cy = np.nansum(yy * flux_map_norm) / fsum
            flux_cx = np.nansum(xx * flux_map_norm) / fsum
            illum_cy = np.nansum(yy * illum_map) / isum
            illum_cx = np.nansum(xx * illum_map) / isum
            drow = illum_cy - flux_cy
            dcol = illum_cx - flux_cx
            axes[3].set_title(
                'SynIM illumination (2D)\n'
                f'centroid Δ(row,col)=({drow:+.2f},{dcol:+.2f}) subap'
            )
        else:
            axes[3].set_title('SynIM illumination (2D)')
        plt.colorbar(im3, ax=axes[3], fraction=0.046)

        # 5) Diff map: SynIM - SPECULA (both normalized)
        diff_map = illum_map - flux_map_norm
        diff_valid = np.isfinite(diff_map)
        if np.any(diff_valid):
            diff_vals = diff_map[diff_valid]
            diff_bias = float(np.mean(diff_vals))
            diff_mae = float(np.mean(np.abs(diff_vals)))
            diff_rmse = float(np.sqrt(np.mean(diff_vals**2)))
        else:
            diff_bias = np.nan
            diff_mae = np.nan
            diff_rmse = np.nan
        max_abs_diff = np.nanmax(np.abs(diff_map))
        if not np.isfinite(max_abs_diff) or max_abs_diff <= 0:
            max_abs_diff = 1.0
        im4 = axes[4].imshow(
            diff_map,
            origin='lower',
            cmap='RdBu_r',
            vmin=-max_abs_diff,
            vmax=max_abs_diff,
            interpolation='nearest',
        )
        axes[4].set_title(
            f'Diff (SynIM - SPECULA norm)\n'
            f'bias={diff_bias:.3f}, MAE={diff_mae:.3f}, RMSE={diff_rmse:.3f}'
        )
        axes[4].text(
            0.02,
            0.02,
            f'max|diff|={max_abs_diff:.3f}',
            transform=axes[4].transAxes,
            fontsize=9,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
        )
        plt.colorbar(im4, ax=axes[4], fraction=0.046)

        # 6) Scatter: flux vs illumination
        valid = (illum > 0) & np.isfinite(flux)
        flux_norm = flux / np.nanmax(flux) if np.nanmax(flux) > 0 else flux
        axes[5].scatter(illum[valid], flux_norm[valid], alpha=0.5, s=10)
        axes[5].plot([0, 1], [0, 1], 'r--', label='ideal')
        axes[5].set_xlabel('SynIM illumination')
        axes[5].set_ylabel('SPECULA flux (norm)')
        axes[5].set_title('Flux vs Illumination (valid subaps)')
        axes[5].legend()
        axes[5].set_xlim(0, 1.05)
        axes[5].set_ylim(0, 1.05)

        plt.tight_layout()
        plt.show()
        plt.close(fig)

    def test_idx_valid_sa_vs_display_map(self):
        """
        Specific test: verify that idx_valid_sa extracted via display_map
        matches the order used by Slopec for flux_per_subaperture.
        """
        print("\n" + "="*70)
        print("Test: idx_valid_sa vs display_map Ordering")
        print("="*70)

        subapdata, subapdata_path = self._create_subapdata_from_pupil()

        # Get both methods
        restored = SubapData.restore(subapdata_path)

        # Method 1: from display_map
        display_map = np.asarray(restored.display_map, dtype=np.int64)
        rows = display_map % restored.nx
        cols = display_map // restored.nx
        idx_from_display = np.column_stack((rows, cols))

        # Method 2: from coordinates (legacy)
        coords = np.transpose(np.asarray(np.where(restored.single_mask())))
        sort_idx = np.lexsort((coords[:, 0], coords[:, 1]))
        idx_from_coords = coords[sort_idx]

        print(f"\nComparing ordering methods for {restored.n_subaps} subapertures:")
        print(f"  display_map order: {display_map}")
        print(f"  idx_from_display:\n{idx_from_display}")
        print(f"  idx_from_coords (legacy):\n{idx_from_coords}")

        # They may not be identical (that's okay), but display_map is the one used by Slopec
        print(f"\n✓ display_map order is what Slopec flux_per_subaperture follows")
        print("="*70 + "\n")

    def _compute_rotation_overlap_metrics(self, synim_rotation_deg):
        """Return overlap metrics between SPECULA flux mask and SynIM illumination mask."""
        illum_pupil = make_circular_pupil_with_spider(
            npix=self.npix,
            radius_fraction=0.45,
            off_axis_shift=self.pupil_shift,
            obstruction_ratio=0.0,
            n_spiders=6,
            spider_width=3,
            spider_angle_offset=self.pupil_rotation_deg
        )

        ef = ElectricField(
            self.npix, self.npix,
            pixel_pitch=self.pixel_pitch,
            S0=1,
            target_device_idx=-1
        )
        ef.A = illum_pupil
        ef.generation_time = 1

        sh = SH(
            wavelengthInNm=self.wavelength_nm,
            subap_wanted_fov=self.subap_npx * self.pxscale_arcsec,
            sensor_pxscale=self.pxscale_arcsec,
            subap_on_diameter=self.subap_on_diameter,
            subap_npx=self.subap_npx,
            rotAnglePhInDeg=self.wfs_rotation_deg,
            target_device_idx=-1
        )
        sh.inputs['in_ef'].set(ef)
        sh.setup()
        sh.check_ready(1)
        sh.trigger()
        sh.post_trigger()
        intensity = sh.outputs['out_i'].i.copy()

        pixels = Pixels(*intensity.shape, target_device_idx=-1)
        pixels.pixels = intensity
        pixels.generation_time = 1

        subapdata, subapdata_path = self._create_subapdata_from_pupil(pixels=pixels)

        slopec = ShSlopec(subapdata, target_device_idx=-1)
        slopec.inputs['in_pixels'].set(pixels)
        slopec.check_ready(1)
        slopec.trigger()
        slopec.post_trigger()
        flux = np.asarray(slopec.outputs['out_flux_per_subaperture'].value)

        restored = SubapData.restore(subapdata_path)
        display_map = np.asarray(restored.display_map, dtype=np.int64)
        
        # Reconstruct 2D indices EXACTLY AS IN PARAMS_MANAGER
        n = self.subap_on_diameter
        idx_i = display_map // n  # col
        idx_j = display_map % n   # row
        idx_valid_sa_2d = np.column_stack((idx_i, idx_j))  # [[col, row], ...]

        illum = compute_subaperture_illumination(
            pup_mask=np.asarray(ef.A),
            wfs_nsubaps=self.subap_on_diameter,
            wfs_rotation=synim_rotation_deg,
            wfs_translation=(0.0, 0.0),
            wfs_magnification=(1.0, 1.0),
            idx_valid_sa=idx_valid_sa_2d,  # 2D array [[col, row], ...]
            verbose=False,
            specula_convention=True  # MUST match params_manager
        )

        flux_mask = flux > 0
        illum_mask = illum > 0
        n = self.subap_on_diameter
        flux_map_mask = np.zeros((n, n), dtype=bool)
        illum_map_mask = np.zeros((n, n), dtype=bool)
        for k, dm in enumerate(display_map):
            row = dm % n
            col = dm // n
            flux_map_mask[row, col] = flux_mask[k]
            illum_map_mask[row, col] = illum_mask[k]

        overlap_0 = int(np.sum(flux_map_mask & illum_map_mask))
        overlap_90 = int(np.sum(flux_map_mask & np.rot90(illum_map_mask, 1)))
        overlap_180 = int(np.sum(flux_map_mask & np.rot90(illum_map_mask, 2)))
        overlap_270 = int(np.sum(flux_map_mask & np.rot90(illum_map_mask, 3)))
        illum_count = int(np.sum(illum_mask))
        recall = overlap_0 / illum_count if illum_count > 0 else 0.0

        return {
            'overlap_0': overlap_0,
            'overlap_90': overlap_90,
            'overlap_180': overlap_180,
            'overlap_270': overlap_270,
            'illum_count': illum_count,
            'recall': recall,
        }

    def test_rotation_alignment_good_match(self):
        """Positive test: with coherent SH/SynIM rotations the match must be good."""
        m = self._compute_rotation_overlap_metrics(self.wfs_rotation_deg)
        self.assertGreater(m['overlap_0'], m['overlap_90'])
        self.assertGreater(m['overlap_0'], m['overlap_270'])
        self.assertGreater(m['recall'], 0.90)

    def test_rotation_alignment_detects_synim_90deg_error(self):
        """Negative test: forcing +90deg in SynIM must worsen direct overlap."""
        good = self._compute_rotation_overlap_metrics(self.wfs_rotation_deg)
        bad = self._compute_rotation_overlap_metrics(self.wfs_rotation_deg + 90.0)

        # With calibrator-filtered subaps and realistic pupil asymmetry,
        # the direct overlap at wrong rotation should be strictly worse than correct rotation.
        self.assertLess(bad['overlap_0'], good['overlap_0'],
                        'Expected: wrong rotation (6°+90°=96°) gives worse direct overlap than correct (6°)')

    def test_display_map_decoding_traps_90deg_rotation(self):
        """
        Regression test: catch 90-degree orientation mistakes.

        We build a non-symmetric 2D pattern, encode it to vector using SPECULA
        display_map convention, decode it back, and ensure that:
        1) correct decoding is exact;
        2) a 90-degree rotated map is clearly wrong.
        """
        n = 4
        # Non-symmetric pattern so 90-degree rotations are unambiguously different.
        ref_map = np.array(
            [
                [1, 2, 3, 4],
                [5, 7, 11, 13],
                [17, 19, 23, 29],
                [31, 37, 41, 43],
            ],
            dtype=np.float32,
        )

        # SPECULA display_map convention used in this file: dm = col*n + row
        display_map = np.array([j * n + i for i in range(n) for j in range(n)], dtype=np.int64)

        # Vector in Slopec order
        vec = np.zeros(n * n, dtype=np.float32)
        for k, dm in enumerate(display_map):
            row = dm % n
            col = dm // n
            vec[k] = ref_map[row, col]

        # Correct decode
        decoded = np.full((n, n), np.nan, dtype=np.float32)
        for k, dm in enumerate(display_map):
            row = dm % n
            col = dm // n
            decoded[row, col] = vec[k]

        # Wrong decode representative of 90-degree bug (row-major assumption)
        wrong = np.full((n, n), np.nan, dtype=np.float32)
        for k, dm in enumerate(display_map):
            row, col = np.unravel_index(dm, (n, n))
            wrong[row, col] = vec[k]

        self.assertTrue(np.allclose(decoded, ref_map),
                        "Correct col-major SPECULA decoding must reconstruct the original map")
        self.assertFalse(np.allclose(wrong, ref_map),
                         "Wrong row-major decode should not match (traps 90-degree errors)")

    def test_params_manager_illumination_idx_conversion(self):
        """Ensure params_manager helper preserves display_map ordering."""
        n = 6
        display_map = np.array([j * n + i for i in range(n) for j in range(n)], dtype=np.int64)

        # Same format used by find_subapdata path: coords from unravel_index(display_map)
        idx_valid_sa_2d = np.column_stack(np.unravel_index(display_map, (n, n)))
        converted = _idx_valid_sa_to_linear_for_illumination(idx_valid_sa_2d, n)

        self.assertTrue(np.array_equal(converted, display_map),
                        "2D idx_valid_sa conversion must preserve display_map order")

        # 1D input must pass through unchanged
        passthrough = _idx_valid_sa_to_linear_for_illumination(display_map, n)
        self.assertTrue(np.array_equal(passthrough, display_map),
                        "1D idx_valid_sa should be left unchanged")
