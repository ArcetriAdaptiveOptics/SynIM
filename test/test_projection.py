import unittest
import numpy as np
from unittest.mock import patch, MagicMock

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


# ============================================================================
# TOMOGRAPHIC PROJECTION ASSEMBLY TESTS
# ============================================================================

class TestTomographicProjectionAssembly(unittest.TestCase):
    """
    Test suite to verify the mode truncation and zero-padding behavior 
    in the tomographic projection matrix assembly.
    """

    @patch('synim.params_manager.ParamsManager.compute_projection_matrix')
    @patch('synim.params_manager.extract_dm_list')
    @patch('synim.params_manager.extract_layer_list')
    def test_mode_truncation_and_ttf_padding(self, mock_layer_list, mock_dm_list, mock_compute_pm):
        from synim.params_manager import ParamsManager

        # 1. Setup a dummy ParamsManager without calling __init__
        pm = ParamsManager.__new__(ParamsManager)
        pm.verbose = False
        pm.proj_reg_factor = 1e-4

        # Mock configuration dictionaries indicating total available modes per component
        pm.params = {
            'dm1': {'nmodes': 100, 'start_mode': 0},
            'dm2': {'nmodes': 100, 'start_mode': 0},
            'layer1': {'nmodes': 200, 'start_mode': 0},
            'layer2': {'nmodes': 200, 'start_mode': 0}
        }

        mock_dm_list.return_value = [{'index': '1'}, {'index': '2'}]
        mock_layer_list.return_value = [{'index': '1'}, {'index': '2'}]

        pm._resolve_component_key = lambda ctype, idx: f"{ctype}{idx}"

        # We want to skip 3 low modes (Tip, Tilt, Focus)
        pm._get_recon_params = MagicMock(return_value={'n_low_modes_zero_filt': 3})

        # Set target truncation configurations (e.g. from YAML file)
        pm.ref_n_modes_dm = [50, 40]      # Target DM modes: 90 total
        pm.ref_n_modes_layer = [60, 70]   # Target Layer modes: 130 total

        # Single optical source with weight 1.0
        pm.projection_params = {
            'opt_sources': [{'weight': 1.0}]
        }

        # 2. Create dummy full projection matrices (as if loaded from disk)
        # Shapes: (n_opt_sources, n_total_modes, n_pupil_modes)
        # We use n_pupil_modes = 500 for testing
        dummy_full_dm = np.random.randn(1, 200, 500).astype(np.float32)
        dummy_full_layer = np.random.randn(1, 400, 500).astype(np.float32)

        mock_compute_pm.return_value = (None, dummy_full_dm, dummy_full_layer)

        # 3. Call the target function
        p_opt, _, _, info = pm.compute_tomographic_projection_matrix(
            output_dir=None, save=False, wfs_type='ref', verbose=False
        )

        # 4. Assertions

        # A. Check if the final shape strictly matches the target truncation limits
        expected_dm_modes = 50 + 40
        expected_layer_modes = 60 + 70
        self.assertEqual(
            p_opt.shape,
            (expected_dm_modes, expected_layer_modes),
            "The projection matrix shape does not match the truncated targets."
        )

        # B. Check if the low modes (Tip, Tilt, Focus) were correctly zero-padded
        # Layer 1 target block: cols 0 to 59. Low modes are 0, 1, 2.
        self.assertTrue(np.all(p_opt[:, 0:3] == 0.0), "Layer 1 low modes are not zeroed.")
        # Verify the rest of Layer 1 is not zeroed out
        self.assertFalse(np.all(p_opt[:, 3:60] == 0.0), "Layer 1 high modes should not be zero.")

        # Layer 2 target block: cols 60 to 129. Low modes are 60, 61, 62.
        self.assertTrue(np.all(p_opt[:, 60:63] == 0.0), "Layer 2 low modes are not zeroed.")
        # Verify the rest of Layer 2 is not zeroed out
        self.assertFalse(np.all(p_opt[:, 63:130] == 0.0), "Layer 2 high modes should not be zero.")

        # C. Check info metadata
        self.assertEqual(info['n_dm_modes_target'], 90)
        self.assertEqual(info['n_layer_modes_target'], 130)
        self.assertEqual(info['n_low_modes_zero_filt'], 3)
