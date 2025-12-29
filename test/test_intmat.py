import unittest
import numpy as np

from synim.synim import (
    interaction_matrix,
    interaction_matrices_multi_wfs,
    rotshiftzoom_array,
    shiftzoom_from_source_dm_params,
    compute_derivatives_with_extrapolation
)
from synim.utils import apply_mask, rebin

import specula
specula.init(device_idx=-1, precision=1)
from specula.lib.make_mask import make_mask

# ============================================================================
# LEGACY FUNCTION - USED ONLY FOR TESTING
# ============================================================================

def update_dm_pup(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                  wfs_rotation, wfs_translation, wfs_magnification,
                  gs_pol_coo, gs_height, verbose=False, specula_convention=True):
    """
    Legacy function - used only for testing the new implementation.
    Update the DM and pupil array to be used in the computation of interaction matrix.
    """

    if specula_convention:
        dm_array = np.transpose(dm_array, (1, 0, 2))
        dm_mask = np.transpose(dm_mask)
        pup_mask = np.transpose(pup_mask)

    pup_diam_pix = pup_mask.shape[0]
    pixel_pitch = pup_diam_m/pup_diam_pix

    if dm_mask.shape[0] != dm_array.shape[0]:
        raise ValueError('Error in input data,'
                         ' the dm and mask array must have the same dimensions.')

    pixel_pitch = pup_diam_m / pup_diam_pix

    dm_translation, dm_magnification = shiftzoom_from_source_dm_params(
        gs_pol_coo, gs_height, dm_height, pixel_pitch
    )
    output_size = (pup_diam_pix,pup_diam_pix)

    trans_dm_array = rotshiftzoom_array(
        dm_array, dm_translation=dm_translation,
        dm_rotation=dm_rotation, dm_magnification=dm_magnification,
        wfs_translation=wfs_translation, wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification, output_size=output_size
    )

    trans_dm_mask  = rotshiftzoom_array(
        dm_mask, dm_translation=dm_translation,
        dm_rotation=dm_rotation, dm_magnification=dm_magnification,
        wfs_translation=(0,0), wfs_rotation=0,
        wfs_magnification=(1,1), output_size=output_size
    )
    trans_dm_mask[trans_dm_mask<0.5] = 0
    if np.max(trans_dm_mask) <= 0:
        raise ValueError('Error in input data, the rotated dm mask is empty.')

    trans_pup_mask  = rotshiftzoom_array(
        pup_mask, dm_translation=(0,0),
        dm_rotation=0, dm_magnification=(1,1),
        wfs_translation=wfs_translation, wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification, output_size=output_size
    )
    trans_pup_mask[trans_pup_mask<0.5] = 0

    if np.max(trans_pup_mask) <= 0:
        raise ValueError('Error in input data, the rotated pup mask is empty.')

    trans_dm_array = apply_mask(trans_dm_array,trans_dm_mask)
    if np.max(trans_dm_array) <= 0:
        raise ValueError('Error in input data, the rotated dm array is empty.')

    return trans_dm_array, trans_dm_mask, trans_pup_mask


def interaction_matrix_former(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                       wfs_nsubaps, wfs_rotation, wfs_translation, wfs_magnification,
                       wfs_fov_arcsec, gs_pol_coo, gs_height, idx_valid_sa=None,
                       verbose=False, display=False, specula_convention=True):
    """
    Legacy function - used only for testing the new implementation.
    Computes a single interaction matrix using the old method.
    """

    trans_dm_array, trans_dm_mask, trans_pup_mask = update_dm_pup(
                  pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                  wfs_rotation, wfs_translation, wfs_magnification,
                  gs_pol_coo, gs_height, verbose=verbose, specula_convention=specula_convention)

    der_dx, der_dy = compute_derivatives_with_extrapolation(trans_dm_array,mask=trans_dm_mask)

    if np.isnan(trans_pup_mask).any():
        np.nan_to_num(trans_pup_mask, copy=False, nan=0.0, posinf=None, neginf=None)
    if np.isnan(trans_dm_mask).any():
        np.nan_to_num(trans_dm_mask, copy=False, nan=0.0, posinf=None, neginf=None)

    pup_mask_sa = rebin(trans_pup_mask, (wfs_nsubaps,wfs_nsubaps), method='sum')
    pup_mask_sa = pup_mask_sa * 1/np.max(pup_mask_sa)

    dm_mask_sa = rebin(trans_dm_mask, (wfs_nsubaps,wfs_nsubaps), method='sum')
    if np.max(dm_mask_sa) <= 0:
        raise ValueError('Error in input data, the dm mask is empty.')
    dm_mask_sa = dm_mask_sa * 1/np.max(dm_mask_sa)

    if np.isnan(der_dx).any():
        np.nan_to_num(der_dx, copy=False, nan=0.0, posinf=None, neginf=None)
    if np.isnan(der_dy).any():
        np.nan_to_num(der_dy, copy=False, nan=0.0, posinf=None, neginf=None)

    der_dx = apply_mask(der_dx, trans_pup_mask, fill_value=np.nan)
    der_dy = apply_mask(der_dy, trans_pup_mask, fill_value=np.nan)

    scale_factor = (der_dx.shape[0]/wfs_nsubaps) \
        / np.median(rebin(trans_pup_mask, (wfs_nsubaps,wfs_nsubaps), method='average'))
    wfs_signal_x = rebin(der_dx, (wfs_nsubaps,wfs_nsubaps), method='nanmean') * scale_factor
    wfs_signal_y = rebin(der_dy, (wfs_nsubaps,wfs_nsubaps), method='nanmean') * scale_factor

    combined_mask_sa = (dm_mask_sa > 0.0) & (pup_mask_sa > 0.0)

    wfs_signal_x = apply_mask(wfs_signal_x, combined_mask_sa, fill_value=0)
    wfs_signal_y = apply_mask(wfs_signal_y, combined_mask_sa, fill_value=0)

    wfs_signal_x_2D = wfs_signal_x.reshape((-1,wfs_signal_x.shape[2]))
    wfs_signal_y_2D = wfs_signal_y.reshape((-1,wfs_signal_y.shape[2]))

    if idx_valid_sa is not None:
        if specula_convention:
            sa2D = np.zeros((wfs_nsubaps,wfs_nsubaps))
            sa2D[idx_valid_sa[:,0], idx_valid_sa[:,1]] = 1
            sa2D = np.transpose(sa2D)
            idx_temp = np.where(sa2D>0)
            idx_valid_sa_new = idx_valid_sa*0.
            idx_valid_sa_new[:,0] = idx_temp[0]
            idx_valid_sa_new[:,1] = idx_temp[1]
        else:
            idx_valid_sa_new = idx_valid_sa

        if len(idx_valid_sa_new.shape) > 1 and idx_valid_sa_new.shape[1] == 2:
            width = wfs_nsubaps
            linear_indices = idx_valid_sa_new[:,0] * width + idx_valid_sa_new[:,1]
            wfs_signal_x_2D = wfs_signal_x_2D[linear_indices.astype(int),:]
            wfs_signal_y_2D = wfs_signal_y_2D[linear_indices.astype(int),:]
        else:
            wfs_signal_x_2D = wfs_signal_x_2D[idx_valid_sa_new.astype(int),:]
            wfs_signal_y_2D = wfs_signal_y_2D[idx_valid_sa_new.astype(int),:]

    if specula_convention:
        im = np.concatenate((wfs_signal_y_2D, wfs_signal_x_2D))
    else:
        im = np.concatenate((wfs_signal_x_2D, wfs_signal_y_2D))

    coeff = 1e-9/(pup_diam_m/wfs_nsubaps) * 206265
    coeff *= 1/(0.5 * wfs_fov_arcsec)
    im = im * coeff

    return im


# ============================================================================
# TESTS
# ============================================================================

class TestIntmat(unittest.TestCase):

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
        self.nmodes = 50  # Zernike modes

        # WFS parameters
        self.wfs_nsubaps = 10
        self.wfs_fov_arcsec = 4.0

        # Create DM with Zernike modes using specula
        from specula.data_objects.ifunc import IFunc

        self.dm_ifunc = IFunc(
            type_str='zern',
            npixels=self.dm_meta_pupil,
            nmodes=self.nmodes,
            obsratio=0.0,
            diaratio=1.0,
            target_device_idx=-1
        )

        # Get DM array and mask from ifunc
        self.dm_array = self.dm_ifunc.ifunc_2d_to_3d(normalize=True)
        self.dm_mask = self.dm_ifunc.mask_inf_func

        # All subapertures valid for simplicity
        self.idx_valid_sa = None

    def test_compare_former_vs_new_single_wfs(self):
        """Compare interaction_matrix_former vs interaction_matrix for on-axis case"""
        # On-axis guide star
        gs_pol_coo = (0.0, 0.0)  # (radius_arcsec, angle_deg)
        gs_height = np.inf

        wfs_rotation = 0.0
        wfs_translation = (0.5, 0.0)
        wfs_magnification = (1.0, 1.0)

        # Compute with former method
        im_former = interaction_matrix_former(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.dm_height, self.dm_rotation,
            self.wfs_nsubaps, wfs_rotation, wfs_translation, wfs_magnification,
            self.wfs_fov_arcsec, gs_pol_coo, gs_height,
            idx_valid_sa=self.idx_valid_sa,
            verbose=False, display=False, specula_convention=True
        )

        # Compute with new method
        im_new = interaction_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.dm_height, self.dm_rotation,
            self.wfs_nsubaps, wfs_rotation, wfs_translation, wfs_magnification,
            self.wfs_fov_arcsec, gs_pol_coo, gs_height,
            idx_valid_sa=self.idx_valid_sa,
            verbose=False, display=False, specula_convention=True
        )

        plot_debug = False
        if plot_debug:
            import matplotlib.pyplot as plt
            modes_to_plot = [0, 2, 10, 25, 49]
            # WFS 0 -- on axis
            for mode in modes_to_plot:
                plt.figure(figsize=(8, 4))
                plt.suptitle(f"WFS 0 (on-axis) - Mode {mode}")
                plt.subplot(1, 2, 1)
                plt.title("Former")
                plt.plot(im_former[:, mode])
                plt.subplot(1, 2, 2)
                plt.title("New")
                plt.plot(im_new[:, mode])
                plt.tight_layout()
                # difference plot new and multi with former
                plt.figure(figsize=(4, 4))
                plt.title("New - Former")
                plt.plot(im_new[:, mode] - im_former[:, mode])
                plt.tight_layout()
                plt.show()
            plt.figure(figsize=(12, 10))
            diff = np.abs(im_former - im_new)
            print("Max difference:", np.max(diff))
            print("Mean difference:", np.mean(diff))
            print("Standard deviation of difference:", np.std(diff))
            plt.imshow(diff, cmap='bwr', aspect='auto')
            plt.colorbar()
            plt.title("Difference between Former and New Interaction Matrices")
            plt.tight_layout()
            plt.show()

        # Compare results
        np.testing.assert_allclose(im_former, im_new, rtol=1e-6, atol=1e-8,
                                   err_msg="Former and new methods differ for on-axis case")

    def test_compare_three_directions(self):
        """Test interaction matrices for 3 different WFS directions"""

        # Define 3 WFS configurations with their own gs_pol_coo
        wfs_configs = [
            {
                'nsubaps': self.wfs_nsubaps,
                'rotation': 0.0,
                'translation': (0.0, 0.0),
                'magnification': (1.0, 1.0),
                'fov_arcsec': self.wfs_fov_arcsec,
                'idx_valid_sa': self.idx_valid_sa,
                'gs_pol_coo': (0.0, 0.0),  # On-axis
                'gs_height': np.inf,
                'name': 'wfs_on_axis'
            },
            {
                'nsubaps': self.wfs_nsubaps,
                'rotation': 0.0,
                'translation': (0.5, 0.0),
                'magnification': (1.0, 1.0),
                'fov_arcsec': self.wfs_fov_arcsec,
                'idx_valid_sa': self.idx_valid_sa,
                'gs_pol_coo': (15.0, 0.0),  # 15" at 0°
                'gs_height': np.inf,
                'name': 'wfs_15arcsec_0deg'
            },
            {
                'nsubaps': self.wfs_nsubaps,
                'rotation': 30.0,
                'translation': (0.0, 0.0),
                'magnification': (1.0, 1.0),
                'fov_arcsec': self.wfs_fov_arcsec,
                'idx_valid_sa': self.idx_valid_sa,
                'gs_pol_coo': (120.0, 135.0),  # 120" at 135°
                'gs_height': np.inf,
                'name': 'wfs_120arcsec_135deg'
            }
        ]

        # Compute with former method (one at a time)
        im_former_list = []
        for wfs_config in wfs_configs:
            im = interaction_matrix_former(
                self.pup_diam_m, self.pup_mask,
                self.dm_array, self.dm_mask,
                self.dm_height, self.dm_rotation,
                wfs_config['nsubaps'],
                wfs_config['rotation'],
                wfs_config['translation'],
                wfs_config['magnification'],
                wfs_config['fov_arcsec'],
                wfs_config['gs_pol_coo'],
                wfs_config['gs_height'],
                idx_valid_sa=self.idx_valid_sa,
                verbose=False, display=False, specula_convention=True
            )
            im_former_list.append(im)

        # Compute with new method (one at a time)
        im_new_list = []
        for wfs_config in wfs_configs:
            im = interaction_matrix(
                self.pup_diam_m, self.pup_mask,
                self.dm_array, self.dm_mask,
                self.dm_height, self.dm_rotation,
                wfs_config['nsubaps'],
                wfs_config['rotation'],
                wfs_config['translation'],
                wfs_config['magnification'],
                wfs_config['fov_arcsec'],
                wfs_config['gs_pol_coo'],
                wfs_config['gs_height'],
                idx_valid_sa=self.idx_valid_sa,
                verbose=False, display=False, specula_convention=True
            )
            im_new_list.append(im)

        # Compute with multi-WFS method (no global gs_pol_coo/gs_height)
        im_multi_dict, info = interaction_matrices_multi_wfs(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.dm_height, self.dm_rotation,
            wfs_configs,  # gs_pol_coo and gs_height are in each config
            verbose=True, specula_convention=True
        )

        print(f"\nMulti-WFS workflow: {info['workflow']}")

        plot_debug = False
        if plot_debug:
            import matplotlib.pyplot as plt
            modes_to_plot = [0, 2, 10, 25, 49]
            idx_gs = 1
            for mode in modes_to_plot:
                plt.figure(figsize=(12, 4))
                plt.suptitle(f"WFS {idx_gs} - Mode {mode}")
                plt.subplot(1, 3, 1)
                plt.title("Former")
                plt.plot(im_former_list[idx_gs][:, mode])
                plt.subplot(1, 3, 2)
                plt.title("New")
                plt.plot(im_new_list[idx_gs][:, mode])
                plt.subplot(1, 3, 3)
                plt.title("Multi-WFS")
                plt.plot(list(im_multi_dict.values())[idx_gs][:, mode])
                plt.tight_layout()
                #difference plots
                plt.figure(figsize=(8, 4))
                plt.suptitle(f"WFS {idx_gs} - Mode {mode} Differences")
                plt.subplot(1, 2, 1)
                plt.title("New - Former")
                plt.plot(im_new_list[idx_gs][:, mode] - im_former_list[idx_gs][:, mode])
                plt.subplot(1, 2, 2)
                plt.title("Multi-WFS - Former")
                plt.plot(list(im_multi_dict.values())[idx_gs][:, mode] \
                    - im_former_list[idx_gs][:, mode])
                plt.tight_layout()
                plt.show()

        # Compare former vs new for each WFS
        for i, (im_former, im_new) in enumerate(zip(im_former_list, im_new_list)):
            np.testing.assert_allclose(
                im_former, im_new, rtol=1e-6, atol=1e-8,
                err_msg=f"Former and new methods differ for WFS {i}"
            )

        # Compare new vs multi-WFS for each configuration
        for i, (im_new, wfs_config) in enumerate(zip(im_new_list, wfs_configs)):
            im_multi = im_multi_dict[wfs_config['name']]
            np.testing.assert_allclose(
                im_new, im_multi, rtol=1e-6, atol=1e-8,
                err_msg=f"New single and multi-WFS methods differ for WFS {i}"
            )

        # Verify shapes
        expected_nslopes = 2 * self.wfs_nsubaps * self.wfs_nsubaps  # x and y slopes
        for i, im in enumerate(im_multi_dict.values()):
            self.assertEqual(im.shape[0], expected_nslopes,
                           f"WFS {i}: Wrong number of slopes")
            self.assertEqual(im.shape[1], self.nmodes,
                           f"WFS {i}: Wrong number of modes")

    def test_off_axis_lgs(self):
        """Test with off-axis LGS configuration"""
        # LGS at 15 arcsec, 90km altitude
        gs_pol_coo = (15.0, 45.0)
        gs_height = 90000.0  # meters

        wfs_rotation = 15.0
        wfs_translation = (0.2, 0.3)
        wfs_magnification = (1.0, 1.0)

        # Compute with both methods
        im_former = interaction_matrix_former(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.dm_height, self.dm_rotation,
            self.wfs_nsubaps, wfs_rotation, wfs_translation, wfs_magnification,
            self.wfs_fov_arcsec, gs_pol_coo, gs_height,
            idx_valid_sa=self.idx_valid_sa,
            verbose=False, display=False, specula_convention=True
        )

        im_new = interaction_matrix(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.dm_height, self.dm_rotation,
            self.wfs_nsubaps, wfs_rotation, wfs_translation, wfs_magnification,
            self.wfs_fov_arcsec, gs_pol_coo, gs_height,
            idx_valid_sa=self.idx_valid_sa,
            verbose=False, display=False, specula_convention=True
        )

        # Uses more relaxed tolerances for the LGS case (complex geometric transformations)
        np.testing.assert_allclose(im_former, im_new, rtol=1e-6, atol=1e-8,
                                   err_msg="Former and new methods differ for LGS case")

    def test_consistency_separated_vs_combined_workflow(self):
        """Verify that separated and combined workflows
        produce identical results without transformations"""
        # Identical configuration for 2 WFS
        wfs_configs = []
        for i in range(2):
            wfs_configs.append({
                'nsubaps': self.wfs_nsubaps,
                'rotation': 0.0,
                'translation': (0.0, 0.0),
                'magnification': (1.0, 1.0),
                'fov_arcsec': self.wfs_fov_arcsec,
                'idx_valid_sa': self.idx_valid_sa,
                'gs_pol_coo': (0.0, 0.0),
                'gs_height': np.inf,
                'name': f'wfs_{i}'
            })

        # Separate calculation
        im_separated = []
        for wfs_config in wfs_configs:
            im = interaction_matrix(
                self.pup_diam_m, self.pup_mask,
                self.dm_array, self.dm_mask,
                self.dm_height, self.dm_rotation,
                wfs_config['nsubaps'],
                wfs_config['rotation'],
                wfs_config['translation'],
                wfs_config['magnification'],
                wfs_config['fov_arcsec'],
                wfs_config['gs_pol_coo'],
                wfs_config['gs_height'],
                idx_valid_sa=self.idx_valid_sa,
                verbose=False, display=False, specula_convention=True
            )
            im_separated.append(im)

        # Combined calculation
        im_multi_dict, info = interaction_matrices_multi_wfs(
            self.pup_diam_m, self.pup_mask,
            self.dm_array, self.dm_mask,
            self.dm_height, self.dm_rotation,
            wfs_configs,
            verbose=False, specula_convention=True
        )

        # Comparison
        for i, wfs_config in enumerate(wfs_configs):
            im_sep = im_separated[i]
            im_comb = im_multi_dict[wfs_config['name']]
            np.testing.assert_allclose(
                im_sep, im_comb, rtol=1e-6, atol=1e-8,
                err_msg=f"Separated and combined workflow differ for WFS {i}"
            )
