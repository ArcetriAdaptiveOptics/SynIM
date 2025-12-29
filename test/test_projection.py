import unittest
import numpy as np

from synim.synpm import projection_matrix
from synim.utils import (
    rotshiftzoom_array,
    shiftzoom_from_source_dm_params,
    apply_mask
)

import specula
specula.init(device_idx=-1, precision=1)
from specula.data_objects.ifunc import IFunc
from specula.lib.make_mask import make_mask

# ============================================================================
# LEGACY FUNCTION - USED ONLY FOR TESTING
# ============================================================================

def update_dm_pup(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                  wfs_rotation, wfs_translation, wfs_magnification,
                  gs_pol_coo, gs_height, verbose=False, specula_convention=True):
    """
    Legacy function - used only for testing the new implementation.
    Update the DM and pupil array to be used in the computation of projection matrix.
    
    NOTE: specula_convention is now handled in the caller (projection_matrix_former)
    """

    # *** TRANSPOSE ONLY IF REQUESTED ***
    if specula_convention:
        dm_array = np.transpose(dm_array, (1, 0, 2))
        dm_mask = np.transpose(dm_mask)
        pup_mask = np.transpose(pup_mask)

    pup_diam_pix = pup_mask.shape[0]
    pixel_pitch = pup_diam_m / pup_diam_pix

    if dm_mask.shape[0] != dm_array.shape[0]:
        raise ValueError('Error in input data, '
                         'the dm and mask array must have the same dimensions.')

    dm_translation, dm_magnification = shiftzoom_from_source_dm_params(
        gs_pol_coo, gs_height, dm_height, pixel_pitch
    )
    output_size = (pup_diam_pix, pup_diam_pix)

    trans_dm_array = rotshiftzoom_array(
        dm_array,
        dm_translation=dm_translation,
        dm_rotation=dm_rotation,
        dm_magnification=dm_magnification,
        wfs_translation=wfs_translation,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        output_size=output_size
    )

    trans_dm_mask = rotshiftzoom_array(
        dm_mask,
        dm_translation=dm_translation,
        dm_rotation=dm_rotation,
        dm_magnification=dm_magnification,
        wfs_translation=(0, 0),
        wfs_rotation=0,
        wfs_magnification=(1, 1),
        output_size=output_size
    )
    trans_dm_mask[trans_dm_mask < 0.5] = 0

    if np.max(trans_dm_mask) <= 0:
        raise ValueError('Error in input data, the rotated dm mask is empty.')

    trans_pup_mask = rotshiftzoom_array(
        pup_mask,
        dm_translation=(0, 0),
        dm_rotation=0,
        dm_magnification=(1, 1),
        wfs_translation=wfs_translation,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        output_size=output_size
    )
    trans_pup_mask[trans_pup_mask < 0.5] = 0

    if np.max(trans_pup_mask) <= 0:
        raise ValueError('Error in input data, the rotated pup mask is empty.')

    trans_dm_array = apply_mask(trans_dm_array, trans_dm_mask)

    if np.max(trans_dm_array) <= 0:
        raise ValueError('Error in input data, the rotated dm array is empty.')

    return trans_dm_array, trans_dm_mask, trans_pup_mask


def projection_matrix_former(pup_diam_m, pup_mask,
                             dm_array, dm_mask,
                             base_inv_array, dm_height,
                             dm_rotation, base_rotation,
                             base_translation, base_magnification,
                             gs_pol_coo, gs_height,
                             verbose=False,
                             specula_convention=True):
    """
    Legacy function - used only for testing the new implementation.
    Computes a projection matrix using the old method.
    """

    # *** IMPORT NEEDED FUNCTION ***
    from synim.synpm import transpose_base_array_for_specula

    # *** SPECULA CONVENTION: Save original mask FIRST ***
    if specula_convention:
        # Save ORIGINAL mask BEFORE transposing
        pup_mask_original = pup_mask.copy()

        # Transpose arrays
        dm_array = np.transpose(dm_array, (1, 0, 2))
        dm_mask = np.transpose(dm_mask)
        pup_mask = np.transpose(pup_mask)

        # *** USE HELPER FUNCTION WITH ORIGINAL MASK ***
        base_inv_array = transpose_base_array_for_specula(
            base_inv_array,
            pup_mask_original,  # Pass ORIGINAL (non-transposed) mask
            verbose=False
        )

    trans_dm_array, trans_dm_mask, trans_pup_mask = update_dm_pup(
        pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
        base_rotation, base_translation, base_magnification,
        gs_pol_coo, gs_height, verbose=verbose,
        specula_convention=False  # Already transposed above
    )

    # Create mask for valid pixels (both in DM and pupil)
    valid_mask = trans_pup_mask.copy() #trans_dm_mask * trans_pup_mask

    # *** USE INDICES INSTEAD OF BOOLEAN MASK ***
    idx_valid = np.where(valid_mask > 0.5)
    n_valid_pixels = len(idx_valid[0])

    # *** EXTRACT DM VALUES USING INDICES ***
    height, width, n_modes = trans_dm_array.shape
    dm_valid_values = trans_dm_array[idx_valid[0], idx_valid[1], :]
    # Shape: (n_valid_pixels, n_modes)

    # *** EXTRACT BASE VALUES USING INDICES ***
    if base_inv_array.ndim == 3:
        # 3D: Extract using same indices
        base_valid_values = base_inv_array[idx_valid[0], idx_valid[1], :]
        # Shape: (n_valid_pixels, n_modes_base)
    elif base_inv_array.ndim == 2:
        # 2D: Handle IFunc or IFuncInv format
        n_rows, n_cols = base_inv_array.shape

        if n_cols == n_valid_pixels:
            # IFunc: (nmodes, npixels_valid)
            base_valid_values = base_inv_array.T
        elif n_rows == n_valid_pixels:
            # IFuncInv: (npixels_valid, nmodes)
            base_valid_values = base_inv_array
        else:
            raise ValueError(
                f"Base shape {base_inv_array.shape} doesn't match {n_valid_pixels} valid pixels"
            )
    else:
        raise ValueError(f"base_inv_array must be 2D or 3D, got {base_inv_array.ndim}D")

    # Perform matrix multiplication to get projection coefficients
    projection = np.dot(dm_valid_values.T, base_valid_values)

    return projection

# ============================================================================
# TESTS
# ============================================================================

class TestProjection(unittest.TestCase):

    def setUp(self):
        """Set up common test parameters"""
        # Pupil parameters
        self.pixel_pupil = 100
        self.dm_meta_pupil = 120
        self.pixel_pitch = 0.01  # 1cm per pixel -> 1m pupil
        self.pup_diam_m = self.pixel_pupil * self.pixel_pitch

        # Create circular pupil mask
        self.pup_mask = make_mask(self.pixel_pupil, obsratio=0.0, diaratio=1.0)

        # DM parameters
        self.dm_height = 1000.0  # meters
        self.dm_rotation = 10.0  # degrees
        self.nmodes_dm = 50  # DM Zernike modes
        self.nmodes_base = 30  # Basis modes (e.g., KL modes)

        # Create DM with Zernike modes using specula
        self.dm_ifunc = IFunc(
            type_str='zern',
            npixels=self.dm_meta_pupil,
            nmodes=self.nmodes_dm,
            obsratio=0.0,
            diaratio=1.0,
            target_device_idx=-1
        )

        # Get DM array and mask from ifunc
        self.dm_array = self.dm_ifunc.ifunc_2d_to_3d(normalize=True)
        self.dm_mask = self.dm_ifunc.mask_inf_func

        self.dm0_ifunc = IFunc(
            type_str='zern',
            npixels=self.pixel_pupil,
            nmodes=self.nmodes_dm,
            obsratio=0.0,
            diaratio=1.0,
            target_device_idx=-1
        )
        self.dm0_array = self.dm0_ifunc.ifunc_2d_to_3d(normalize=True)
        self.dm0_mask = self.dm0_ifunc.mask_inf_func

        # Create basis (e.g., KL modes or another Zernike set)
        # For testing, we use another set of Zernike modes
        self.base_ifunc = IFunc(
            type_str='zern',
            npixels=self.pixel_pupil,
            nmodes=self.nmodes_base,
            obsratio=0.0,
            diaratio=1.0,
            target_device_idx=-1
        )

        # Invert self.base_ifunc and use as base_inv_array
        base_ifunc_inv = self.base_ifunc.inverse()
        self.base_inv_array = np.zeros(
            (self.pixel_pupil, self.pixel_pupil, self.nmodes_base), dtype=self.dm_array.dtype
        )
        idx = np.where(self.pup_mask > 0)
        self.base_inv_array[idx[0], idx[1], :] = base_ifunc_inv.ifunc_inv

    def test_compare_former_vs_new_on_axis(self):
        """Compare projection_matrix_former vs projection_matrix for on-axis case"""
        # On-axis guide star
        gs_pol_coo = (0.0, 0.0)
        gs_height = np.inf

        base_rotation = 0.0
        base_translation = (0.0, 0.0)
        base_magnification = (1.0, 1.0)

        # Compute with former method
        pm_former = projection_matrix_former(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True
        )

        # Compute with new method
        pm_new = projection_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True,
            specula_convention_inv=True
        )

        plot_debug = False
        if plot_debug:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(12,5))
            plt.subplot(1,2,1)
            plt.title('Projection Matrix - Former Method')
            plt.imshow(pm_former, cmap='viridis', aspect='auto')
            plt.colorbar()
            plt.subplot(1,2,2)
            plt.title('Projection Matrix - New Method')
            plt.imshow(pm_new, cmap='viridis', aspect='auto')
            plt.colorbar()
            plt.show()

        # Compare results
        np.testing.assert_allclose(pm_former, pm_new, rtol=1e-6, atol=1e-8,
                                   err_msg="Former and new methods differ for on-axis case")

        # Verify shapes
        self.assertEqual(pm_new.shape[0], self.nmodes_dm,
                        "Wrong number of DM modes")
        self.assertEqual(pm_new.shape[1], self.nmodes_base,
                        "Wrong number of basis modes")

    def test_base_rotation_only(self):
        """Test with only basis rotation applied"""
        gs_pol_coo = (0.0, 0.0)
        gs_height = np.inf

        base_rotation = 30.0
        base_translation = (0.0, 0.0)
        base_magnification = (1.0, 1.0)

        pm_former = projection_matrix_former(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True
        )

        pm_new = projection_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True,
            specula_convention_inv=True
        )

        np.testing.assert_allclose(pm_former, pm_new, rtol=1e-6, atol=1e-8,
                                   err_msg="Former and new methods differ for basis rotation only")

    def test_off_axis_guide_star(self):
        """Test with off-axis guide star (no basis transformations)"""
        gs_pol_coo = (15.0, 0.0)  # 15 arcsec off-axis
        gs_height = np.inf

        base_rotation = 0.0
        base_translation = (0.0, 0.0)
        base_magnification = (1.0, 1.0)

        pm_former = projection_matrix_former(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True
        )

        pm_new = projection_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True,
            specula_convention_inv=True
        )

        np.testing.assert_allclose(pm_former, pm_new, rtol=1e-6, atol=1e-8,
                                   err_msg="Former and new methods differ for off-axis GS")

    def test_combined_transformations(self):
        """Test with both DM and basis transformations"""
        gs_pol_coo = (10.0, 45.0)
        gs_height = np.inf

        base_rotation = 15.0
        base_translation = (0.5, 0.3)
        base_magnification = (1.0, 1.0)

        pm_former = projection_matrix_former(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True
        )

        pm_new = projection_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            self.dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True,
            specula_convention_inv=True
        )

        np.testing.assert_allclose(pm_former, pm_new, rtol=1e-6, atol=1e-8,
                                   err_msg="Former and new methods differ for combined transforms")

    def test_dm_at_ground_level(self):
        """Test with DM at ground level (height=0)"""
        gs_pol_coo = (0.0, 0.0)
        gs_height = np.inf

        base_rotation = 0.0
        base_translation = (0.0, 0.0)
        base_magnification = (1.0, 1.0)

        dm_height = 0.0  # Ground level DM

        pm_former = projection_matrix_former(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True
        )

        pm_new = projection_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.base_inv_array,
            dm_height, self.dm_rotation,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True,
            specula_convention_inv=True
        )

        np.testing.assert_allclose(pm_former, pm_new, rtol=1e-6, atol=1e-8,
                                   err_msg="Former and new methods differ for ground-level DM")

    def test_projection_identity(self):
        """Test that the projection between two identical bases has a dominant diagonal"""
        gs_pol_coo = (0.0, 0.0)
        gs_height = np.inf
        base_rotation = 0.0
        base_translation = (0.0, 0.0)
        base_magnification = (1.0, 1.0)

        pm = projection_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm0_array, self.dm0_mask,
            self.base_inv_array,
            0.0, 0.0,
            base_rotation, base_translation, base_magnification,
            gs_pol_coo, gs_height,
            verbose=False, specula_convention=True,
            specula_convention_inv=True
        )
        # Consider only the square part for diagonal dominance
        n_diag = min(pm.shape[0], pm.shape[1])
        diag = np.abs(np.array([pm[i, i] for i in range(n_diag)]))
        off_diag = np.array([np.delete(np.abs(pm[i, :]), i) for i in range(n_diag)])
        max_off_diag = off_diag.max(axis=1)
        ratio = diag / (max_off_diag + 1e-12)
        self.assertTrue(np.all(ratio > 5),
                        "The diagonal is not sufficiently dominant in the"
                        " projection matrix between identical bases")

    def test_projection_zero_mask(self):
        """Test that the projection with a zero mask raises ValueError"""
        gs_pol_coo = (0.0, 0.0)
        gs_height = np.inf
        base_rotation = 0.0
        base_translation = (0.0, 0.0)
        base_magnification = (1.0, 1.0)
        zero_mask = np.zeros_like(self.pup_mask)

        with self.assertRaises(ValueError):
            projection_matrix(
                self.pup_diam_m, zero_mask,
                self.dm_array, self.dm_mask,
                self.base_inv_array,
                self.dm_height, self.dm_rotation,
                base_rotation, base_translation, base_magnification,
                gs_pol_coo, gs_height,
                verbose=False, specula_convention=True,
                specula_convention_inv=True
            )
