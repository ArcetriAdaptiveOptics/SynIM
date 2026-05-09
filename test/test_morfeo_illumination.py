
import unittest
import os
import numpy as np
import matplotlib.pyplot as plt

import specula
specula.init(device_idx=-1, precision=1)

from specula.data_objects.electric_field import ElectricField
from specula.processing_objects.sh import SH
from specula.data_objects.pixels import Pixels
from specula.data_objects.subap_data import SubapData
from specula.processing_objects.sh_slopec import ShSlopec

from synim.params_manager import ParamsManager, _idx_valid_sa_to_linear_for_illumination
import synim.synim as synim


def map_1d_to_2d(data_1d, display_map, n_subaps):
    """
    Maps a 1D array of valid subapertures back to a 2D grid using the display_map.
    SPECULA display_map is column-major: linear_idx = col * n_subaps + row.
    """
    map_2d = np.full((n_subaps, n_subaps), np.nan)
    for k, linear_idx in enumerate(display_map):
        row = linear_idx % n_subaps
        col = linear_idx // n_subaps
        map_2d[row, col] = data_1d[k]
    return map_2d


class TestMorfeoIlluminationSpecula(unittest.TestCase):
    """
    Test to verify consistency between the analytical subaperture illumination
    computed by SynIM and the actual optical flux computed by SPECULA primitives.
    
    Evaluates the 6 MORFEO LGS WFS, which feature different rotations.
    """

    def setUp(self):
        """Initialize the ParamsManager and load the MORFEO configuration."""
        self.yaml_file = '/home/guido/pythonLib/SPECULA_scripts/morfeo/params_morfeo_calib.yml'
        self.root_dir = '/raid1/guido/PASSATA/MAORYC'
        
        # Check if files exist to avoid test failure on machines without the data
        if not os.path.exists(self.yaml_file) or not os.path.exists(self.root_dir):
            self.skipTest("MORFEO configuration files or root directory not found.")
            
        self.pm = ParamsManager(self.yaml_file, root_dir=self.root_dir, verbose=False)
        self.main_params = self.pm.params['main']
        
        # The pupil mask loaded by ParamsManager
        self.pupil_mask = np.asarray(self.pm.pup_mask)
        self.npix = self.pupil_mask.shape[0]
        self.pixel_pitch = self.main_params['pixel_pitch']

    def test_morfeo_lgs_illumination(self):
        """Iterate over all 6 LGS and compare SynIM vs SPECULA flux."""
        print("\n" + "="*70)
        print("Testing SynIM Illumination vs SPECULA Primitives (MORFEO LGS)")
        print("="*70)

        # Iterate over the 6 LGS
        for wfs_idx in range(1, 7):
            print(f"\n--- Analyzing LGS {wfs_idx} ---")
            
            wfs_key = f'sh_lgs{wfs_idx}'
            slopec_key = f'slopec_lgs{wfs_idx}'
            
            # Fetch WFS parameters from YAML
            sh_cfg = self.pm.params[wfs_key]
            slopec_cfg = self.pm.params[slopec_key]
            
            # 1. Setup SPECULA ElectricField
            ef = ElectricField(
                self.npix, self.npix,
                pixel_pitch=self.pixel_pitch,
                S0=1,
                target_device_idx=-1
            )
            # Inject the pupil mask purely as amplitude. 
            # Flat phase guarantees we only measure geometric footprint.
            ef.A = self.pupil_mask
            ef.generation_time = 1

            # 2. Setup SPECULA SH Sensor
            sh = SH(
                wavelengthInNm=sh_cfg['wavelengthInNm'],
                subap_wanted_fov=sh_cfg['subap_wanted_fov'],
                sensor_pxscale=sh_cfg['sensor_pxscale'],
                subap_on_diameter=sh_cfg['subap_on_diameter'],
                subap_npx=sh_cfg['subap_npx'],
                rotAnglePhInDeg=sh_cfg.get('rotAnglePhInDeg', 0.0),
                target_device_idx=-1
            )
            sh.inputs['in_ef'].set(ef)
            sh.setup()
            sh.check_ready(1)
            sh.trigger()
            sh.post_trigger()
            
            # Extract CCD intensity
            intensity = sh.outputs['out_i'].i.copy()
            pixels = Pixels(*intensity.shape, target_device_idx=-1)
            pixels.pixels = intensity
            pixels.generation_time = 1

            # 3. Setup SPECULA ShSlopec
            subap_tag = slopec_cfg['subapdata_object']
            subap_path = self.pm.cm.filename('subapdata', subap_tag)
            subapdata = SubapData.restore(subap_path)

            slopec = ShSlopec(subapdata, target_device_idx=-1)
            slopec.inputs['in_pixels'].set(pixels)
            slopec.check_ready(1)
            slopec.trigger()
            slopec.post_trigger()

            # Extract SPECULA Flux
            flux_specula = np.asarray(slopec.outputs['out_flux_per_subaperture'].value)
            max_flux = np.nanmax(flux_specula)
            flux_specula_norm = flux_specula / max_flux if max_flux > 0 else flux_specula

            # 4. Compute SynIM Illumination
            wfs_params = self.pm.get_wfs_params('lgs', wfs_idx, xp_local=np)
            n_subaps = wfs_params['wfs_nsubaps']

            idx_valid_sa_illum = _idx_valid_sa_to_linear_for_illumination(
                wfs_params['idx_valid_sa'], n_subaps
            )

            illum_synim = synim.compute_subaperture_illumination(
                pup_mask=self.pupil_mask,
                wfs_nsubaps=n_subaps,
                wfs_rotation=wfs_params['wfs_rotation'],
                wfs_translation=wfs_params['wfs_translation'],
                wfs_magnification=wfs_params['wfs_magnification'],
                idx_valid_sa=wfs_params['idx_valid_sa'],
                verbose=False,
                specula_convention=True
            )

            # 5. Compare Results
            diff = illum_synim - flux_specula_norm
            mean_abs_err = float(np.mean(np.abs(diff)))
            max_abs_err = float(np.max(np.abs(diff)))

            print(f"  Rotation angle: {wfs_params['wfs_rotation']}°")
            print(f"  Valid subaps:   {len(illum_synim)}")
            print(f"  Mean Abs Error: {mean_abs_err:.4f}")
            print(f"  Max Abs Error:  {max_abs_err:.4f}")

            # Optional: Generate plots for visual inspection
            self._plot_comparison(
                wfs_idx, wfs_params['wfs_rotation'],
                flux_specula_norm, illum_synim, diff,
                subapdata.display_map, n_subaps, max_abs_err
            )

            # Assert that the maximum absolute error is below an acceptable threshold
            # (Tolerance is usually ~0.05 due to minor differences in resampling techniques)
            self.assertLess(
                max_abs_err, 0.05,
                f"LGS {wfs_idx}: SynIM illumination deviates too much from SPECULA flux."
            )

        print("\n✓ All MORFEO LGS WFS matched successfully.")

    def _plot_comparison(self, wfs_idx, rotation, specula_flux, synim_illum, diff, display_map, n_subaps, max_abs_err):
        """Generates a 2D map comparison and saves it to disk."""
        display_map = np.asarray(display_map, dtype=np.int64)

        specula_2d = map_1d_to_2d(specula_flux, display_map, n_subaps)
        synim_2d = map_1d_to_2d(synim_illum, display_map, n_subaps)
        diff_2d = map_1d_to_2d(diff, display_map, n_subaps)

        fig, axs = plt.subplots(1, 4, figsize=(22, 5))
        fig.suptitle(f"MORFEO LGS {wfs_idx} Illumination (Rotation: {rotation}°)", fontsize=14)

        im0 = axs[0].imshow(specula_2d, cmap='viridis', vmin=0, vmax=1, origin='lower')
        axs[0].set_title('SPECULA Pipeline Flux (Normalized)')
        fig.colorbar(im0, ax=axs[0], fraction=0.046)

        im1 = axs[1].imshow(synim_2d, cmap='viridis', vmin=0, vmax=1, origin='lower')
        axs[1].set_title('SynIM Analytic Illumination')
        fig.colorbar(im1, ax=axs[1], fraction=0.046)

        limit = max(0.05, np.nanmax(np.abs(diff_2d)))
        im2 = axs[2].imshow(diff_2d, cmap='RdBu_r', vmin=-limit, vmax=limit, origin='lower')
        axs[2].set_title(f'Difference (SynIM - SPECULA)\nMax Err: {max_abs_err:.3f}')
        fig.colorbar(im2, ax=axs[2], fraction=0.046)

        axs[3].scatter(specula_flux, synim_illum, alpha=0.5, s=10)
        axs[3].plot([0, 1], [0, 1], 'r--', label='Ideal Match')
        axs[3].set_xlabel('SPECULA Flux')
        axs[3].set_ylabel('SynIM Illumination')
        axs[3].set_title('Correlation')
        axs[3].grid(True, alpha=0.3)
        axs[3].legend()

        plt.tight_layout()
        plt.savefig(f'illumination_comparison_lgs{wfs_idx}.png', dpi=150)
        plt.close(fig)

if __name__ == '__main__':
    unittest.main()
