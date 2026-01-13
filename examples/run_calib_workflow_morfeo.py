"""
Automated calibration workflow for MORFEO MCAO system.

This script orchestrates the complete calibration pipeline:
0. Influence functions and modal bases computation (SPECULA)
1. Subaperture geometry calibration (SPECULA)
2. Interaction matrices computation (SynIM)
3. Filter matrices generation for slope filtering
4. LO (NGS) and reference (REF) interaction matrices
5. Focus reconstruction matrix
6. Tomographic reconstructor computation
7. NGS reconstructor
8. Projection matrices
9. Final YAML configuration update
10. Workflow summary generation

The workflow automatically updates YAML files with generated filenames
and manages dependencies between calibration steps.
"""
import os
import sys
import yaml
import shutil
from datetime import datetime
from pathlib import Path
import subprocess
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import re
        
import numpy as np

import specula
specula.init(device_idx=-1, precision=1)

import synim
synim.init(device_idx=0, precision=1)
from synim.params_manager import ParamsManager
from synim.params_utils import compute_influence_functions_and_modalbases
from synim.params_utils import generate_filter_matrix_from_intmat_file

class MORFEOCalibrationWorkflow:
    """Manages the complete MORFEO calibration pipeline."""

    def __init__(self,
                 base_yml_path: str,
                 root_dir: str,
                 output_base_dir: str = None,
                 device_idx: int = 0,
                 overwrite: bool = False,
                 verbose: bool = True):
        """
        Initialize the calibration workflow.
        
        Parameters
        ----------
        base_yml_path : str
            Path to the initial MORFEO YAML configuration file (template)
        root_dir : str
            Root directory for calibration files (SPECULA CalibManager root)
        output_base_dir : str, optional
            Base directory for workflow outputs. If None, uses root_dir/workflow/
        device_idx : int
            GPU device index (-1 for CPU)
        overwrite : bool
            Whether to overwrite existing calibration files
        verbose : bool
            Verbose output
        """
        self.base_yml_path = Path(base_yml_path)
        self.root_dir = Path(root_dir)
        self.device_idx = device_idx
        self.overwrite = overwrite
        self.verbose = verbose

        # Create timestamped workflow directory
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if output_base_dir is None:
            self.workflow_dir = self.root_dir / "workflow" / self.timestamp
        else:
            self.workflow_dir = Path(output_base_dir) / self.timestamp

        self.workflow_dir.mkdir(parents=True, exist_ok=True)

        # Initialize directory structure
        self.yml_dir = self.workflow_dir / "yml"
        self.yml_dir.mkdir(exist_ok=True)

        # Copy and load initial YAML
        self.initial_yml_path = self.yml_dir / f"params_morfeo_initial_{self.timestamp}.yml"
        shutil.copy(self.base_yml_path, self.initial_yml_path)

        with open(self.initial_yml_path) as f:
            self.config = yaml.safe_load(f)

        # Store generated filenames for tracking
        self.generated_files = {}

        # Extract calibration parameters from YAML
        self._extract_calibration_params()

        self._log(f"Initialized MORFEO calibration workflow")
        self._log(f"  Base YAML: {self.base_yml_path}")
        self._log(f"  Root dir: {self.root_dir}")
        self._log(f"  Workflow dir: {self.workflow_dir}")
        self._log(f"  Timestamp: {self.timestamp}")

    def _extract_calibration_params(self):
        """Extract calibration parameters from YAML configuration."""
        # Get main parameters
        main_config = self.config.get('main', {})
        self.pixel_pupil = main_config.get('pixel_pupil', 160)
        self.pixel_pitch = main_config.get('pixel_pitch', 0.2406)
        self.telescope_diameter = self.pixel_pupil * self.pixel_pitch

        # LGS filter parameters
        self.lgs_filter_n_modes = 1000
        self.lgs_filter_n_modes_filtered = 3  # tip, tilt, focus

        # REF focus parameters
        self.ref_focus_n_modes = 30

        # Get reconstruction parameters
        if 'reconstruction' in self.config:
            rec_config = self.config['reconstruction']
            self.ngs_n_modes_dm = rec_config.get('ngs_n_modes_dm', [2, 0, 3])
            self.ref_n_modes_dm = rec_config.get('ref_n_modes_dm', [50, 45, 48])
        else:
            self._log("⚠ No reconstruction section in YAML, using defaults")
            self.ngs_n_modes_dm = [2, 0, 3]
            self.ref_n_modes_dm = [50, 45, 48]

        # Atmospheric parameters from YAML
        atmo_config = self.config.get('atmo', {})
        self.r0 = 0.15  # Will be set at runtime
        self.L0 = atmo_config.get('L0', 25.0)

        # DM and layer configuration
        self.dm_altitudes = [600, 6500, 17500]  # From dm1, dm2, dm3 heights in YAML
        self.layer_altitudes = [0, 2000, 4500, 11000, 15000, 18000, 22000]  # From layer1-7

        # Number of actuators per component (from morfeo_compute_influence_functions.py)
        self.dm_nacts = [41, 37, 37]  # dm1, dm2, dm3
        self.layer_nacts = [41, 41, 41, 37, 37, 37, 37]  # layer1-7

        # Get projection parameters
        if 'projection' in self.config:
            proj_config = self.config['projection']
            self.proj_reg_factor = proj_config.get('reg_factor', 1e-4)
        else:
            self._log("⚠ No projection section in YAML, using defaults")
            self.proj_reg_factor = 1e-4

        # Pupil mask parameters
        self.n_petals = 6
        self.obsratio = 0.283

    def _log(self, message: str):
        """Print log message if verbose."""
        if self.verbose:
            print(f"[MORFEO Workflow] {message}")

    def _save_yaml(self, config: dict, filename: str) -> Path:
        """Save YAML configuration to workflow directory."""
        output_path = self.yml_dir / filename
        with open(output_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        self._log(f"  Saved YAML: {output_path}")
        return output_path

    def _update_config_field(self, key_path: str, value):
        """
        Update a field in the configuration using dot notation.
        
        Example: key_path = 'slopec_lgs1.filtmat_data'
        """
        keys = key_path.split('.')
        config = self.config

        # Navigate to the parent key
        for key in keys[:-1]:
            if key not in config:
                config[key] = {}
            config = config[key]

        # Set the final value
        config[keys[-1]] = value

    def _run_specula(self, *yml_files):
        """Run SPECULA simulation with given YAML files."""
        cmd = ['specula'] + [str(f) for f in yml_files]
        self._log(f"  Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"SPECULA failed:\n{result.stderr}")
        return result

    # =======================================================================
    # STEP 0: Influence Functions and Modal Bases
    # =======================================================================

    def step0_influence_functions(self):
        """
        Step 0: Compute influence functions and modal bases for DMs and layers.
        
        Generates:
        - Pupil mask with petals
        - Influence functions for 3 DMs at different altitudes
        - Influence functions for 7 atmospheric layers
        - Modal bases (KL modes) and M2C matrices
        - Inverse modal base for ground layer (for modal_analysis)
        """
        self._log("\n" + "="*70)
        self._log("STEP 0: INFLUENCE FUNCTIONS AND MODAL BASES")
        self._log("="*70)

        # Call the library function
        result = compute_influence_functions_and_modalbases(
            root_dir=self.root_dir,
            pixel_pupil=self.pixel_pupil,
            pixel_pitch=self.pixel_pitch,
            dm_altitudes=self.dm_altitudes,
            dm_nacts=self.dm_nacts,
            layer_altitudes=self.layer_altitudes,
            layer_nacts=self.layer_nacts,
            fov_arcsec=self.config.get('atmo', {}).get('fov', 160),
            n_petals=self.n_petals,
            obsratio=self.obsratio,
            r0=self.r0,
            L0=self.L0,
            overwrite=self.overwrite,
            verbose=self.verbose,
            return_pupil_mask=False  # Don't need to return it
        )

        # Store results
        self.generated_files['pupil_mask_tag'] = result['pupil_mask_tag']
        self.generated_files['ifunc_tags'] = result['ifunc_tags']
        self.generated_files['m2c_tags'] = result['m2c_tags']
        self.generated_files['ifunc_inv_tag'] = result['ifunc_inv_tag']

        # Update config with generated tags
        self._update_config_field('pupilstop.tag', result['pupil_mask_tag'])

        # Update DMs
        for i in result['dm_indices']:
            dm_key = f"dm{i}"
            if dm_key in result['ifunc_tags']:
                self._update_config_field(f'{dm_key}.ifunc_object',
                                          result['ifunc_tags'][dm_key])
                self._update_config_field(f'{dm_key}.m2c_object',
                                          result['m2c_tags'][dm_key])
                # Extract n_modes from m2c tag
                n_modes = int(result['m2c_tags'][dm_key].split('_')[-1].replace('modes', ''))
                self._update_config_field(f'{dm_key}.nmodes', n_modes)

        # Update layers
        for i in result['layer_indices']:
            layer_key = f"layer{i}"
            if layer_key in result['ifunc_tags']:
                self._update_config_field(f'{layer_key}.ifunc_object',
                                          result['ifunc_tags'][layer_key])
                self._update_config_field(f'{layer_key}.m2c_object',
                                          result['m2c_tags'][layer_key])
                n_modes = int(result['m2c_tags'][layer_key].split('_')[-1].replace('modes', ''))
                self._update_config_field(f'{layer_key}.nmodes', n_modes)

        # Update modal_analysis with inverse ifunc
        if result['ifunc_inv_tag']:
            self._update_config_field('modal_analysis.ifunc_inv_object',
                                      result['ifunc_inv_tag'])

        self._log(f"\n✓ Generated influence functions for {len(result['ifunc_tags'])} components")

        return result

    # =======================================================================
    # STEP 1: Subaperture Calibration
    # =======================================================================

    def step1_subap_calibration(self):
        """
        Step 1: Calibrate subaperture geometry for all WFS types.
        
        Generates:
        - LGS subaperture data (3 unique rotations)
        - NGS subaperture data (1x1)
        - REF subaperture data (8x8)
        """
        self._log("\n" + "="*70)
        self._log("STEP 1: SUBAPERTURE CALIBRATION")
        self._log("="*70)

        # Find the calib YAML file
        calib_yml = self.base_yml_path.parent / "calib_morfeo_tiny_subaps.yml"
        if not calib_yml.exists():
            raise FileNotFoundError(f"Subaperture calibration YAML not found: {calib_yml}")

        # Run SPECULA calibration
        self._log("Running subaperture calibration with SPECULA...")
        self._run_specula(self.initial_yml_path, calib_yml)

        # Extract generated subapdata tags from the calibration YAML
        with open(calib_yml) as f:
            calib_config = yaml.safe_load(f)

        subap_tags = {}
        for key, value in calib_config.items():
            if key.startswith('sh_subaps_'):
                wfs_name = key.replace('sh_subaps_', '')
                subap_tags[wfs_name] = value['output_tag']

        self.generated_files['subap_tags'] = subap_tags

        # Update config with subap tags
        for wfs_name, tag in subap_tags.items():
            slopec_key = f"slopec_{wfs_name}"
            if slopec_key in self.config:
                self._update_config_field(f"{slopec_key}.subapdata_object", tag)

        self._log(f"✓ Generated subaperture data for {len(subap_tags)} WFS configurations")

        for wfs_name, tag in subap_tags.items():
            self._log(f"  {wfs_name}: {tag}")

        return subap_tags

    # =======================================================================
    # STEP 2: LGS Interaction Matrices
    # =======================================================================

    def step2_lgs_interaction_matrices(self):
        """
        Step 2: Compute LGS interaction matrices for all layers.
        
        Uses SynIM to compute:
        - Layer 0 (ground) interaction matrices for 3 unique LGS rotations
        - Layer interaction matrices for tomography (all 7 layers)
        """
        self._log("\n" + "="*70)
        self._log("STEP 2: LGS INTERACTION MATRICES")
        self._log("="*70)

        # Save current config state
        current_yml = self._save_yaml(
            self.config,
            f"params_morfeo_step2_{self.timestamp}.yml"
        )

        # Initialize ParamsManager
        params_mgr = ParamsManager(
            str(current_yml),
            root_dir=str(self.root_dir),
            verbose=self.verbose
        )

        # Compute interaction matrices for LGS
        self._log("Computing LGS interaction matrices...")
        im_paths = params_mgr.compute_interaction_matrices(
            wfs_type='lgs',
            output_im_dir=str(self.root_dir / "synim"),
            output_rec_dir=str(self.root_dir / "synrec"),
            overwrite=self.overwrite,
            verbose=self.verbose,
            display=False
        )

        self.generated_files['lgs_im_paths'] = im_paths
        self._log(f"✓ Generated {len(im_paths)} LGS interaction matrices")

        # Log details
        layer0_count = sum(1 for p in im_paths if 'layH0.0' in str(p))
        self._log(f"  Layer 0 (ground): {layer0_count} matrices")
        self._log(f"  Other layers: {len(im_paths) - layer0_count} matrices")

        return im_paths

    # =======================================================================
    # STEP 3: Filter Matrices
    # =======================================================================

    def step3_filter_matrices(self):
        """
        Step 3: Generate filter matrices for slope filtering.
        
        Creates filtmat files from layer 0 interaction matrices.
        Automatically detects unique WFS configurations (rotations, positions, etc.)
        from the interaction matrix filenames.
        
        These matrices filter tip, tilt, and focus modes from the slopes.
        """
        self._log("\n" + "="*70)
        self._log("STEP 3: FILTER MATRICES")
        self._log("="*70)

        lgs_im_paths = self.generated_files.get('lgs_im_paths', [])
        if not lgs_im_paths:
            raise RuntimeError("No LGS interaction matrices found. Run step2 first.")

        # Filter only layer 0 IMs (ground level, where filtering is needed)
        layer0_ims = [p for p in lgs_im_paths if 'layH0.0' in str(p)]

        if not layer0_ims:
            raise RuntimeError("No layer 0 interaction matrices found")

        filtmat_paths = []
        filtmat_mapping = {}  # unique_id -> filtmat info
        wfs_to_filtmat = {}   # wfs_name -> filtmat_name

        # Group IMs by unique configuration (rotation/position)
        config_groups = defaultdict(list)

        for im_path in layer0_ims:
            filename = im_path.stem

            # Extract WFS identifier and configuration from filename
            # Expected patterns like: "IM_syn_..._lgs1_rot6.2_..." or "IM_syn_..._lgs1_..."

            # Find WFS name
            wfs_match = re.search(r'_(lgs\d+|ngs\d+|ref\d+)_', filename)
            if not wfs_match:
                self._log(f"  ⚠ Could not extract WFS name from: {im_path.name}")
                continue

            wfs_name = wfs_match.group(1)

            # Extract rotation angle (if present)
            rot_match = re.search(r'rot(-?\d+\.?\d*)', filename)
            rot_str = f"rot{rot_match.group(1)}" if rot_match else "rot0.0"

            # Extract position coordinates (if present)
            pos_match = re.search(r'pos\[(-?\d+\.?\d*),(-?\d+\.?\d*)\]', filename)
            pos_str = f"pos[{pos_match.group(1)},{pos_match.group(2)}]" if pos_match else ""

            # Create unique configuration identifier
            if pos_str:
                config_id = f"{rot_str}_{pos_str}"
            else:
                config_id = rot_str

            config_groups[config_id].append((im_path, wfs_name))

        self._log(f"Found {len(config_groups)} unique WFS configurations:")
        for config_id, wfs_list in config_groups.items():
            wfs_names = [wfs for _, wfs in wfs_list]
            self._log(f"  {config_id}: {wfs_names}")

        # Generate filtmat for each unique configuration
        for config_id, wfs_list in config_groups.items():
            # Use the first IM of each group (they should all be equivalent)
            im_path, _ = wfs_list[0]

            # Create filtmat name based on configuration
            # Clean up config_id for filename (remove special characters)
            config_str = config_id.replace('[', '').replace(']', '').replace(',', '_')

            output_name = (f"filtmat_IM_syn_{config_str}_"
                        f"mn{self.lgs_filter_n_modes}_{self.lgs_filter_n_modes_filtered}")
            output_path = self.root_dir / "data" / f"{output_name}.fits"
            output_path.parent.mkdir(parents=True, exist_ok=True)

            self._log(f"\nGenerating filter matrix for {config_id}...")
            self._log(f"  Using IM: {im_path.name}")
            self._log(f"  Output: {output_name}")

            # Use library function
            generate_filter_matrix_from_intmat_file(
                intmat_filename=str(im_path),
                n_modes=self.lgs_filter_n_modes,
                n_modes_filtered=self.lgs_filter_n_modes_filtered,
                output_filename=str(output_path),
                cut_modes=0,
                smooth_tt=True,
                overwrite=self.overwrite,
                verbose=False  # Don't clutter output
            )

            filtmat_paths.append(output_path)
            filtmat_mapping[config_id] = {
                'filename': output_name,
                'path': output_path,
                'wfs_names': [wfs for _, wfs in wfs_list],
                'im_path': im_path
            }

            # Map each WFS to its filtmat
            for _, wfs_name in wfs_list:
                wfs_to_filtmat[wfs_name] = output_name

        self.generated_files['filtmat_paths'] = filtmat_paths
        self.generated_files['filtmat_mapping'] = filtmat_mapping
        self.generated_files['wfs_to_filtmat'] = wfs_to_filtmat

        # Update config with filtmat filenames
        for wfs_name, filtmat_name in wfs_to_filtmat.items():
            slopec_key = f"slopec_{wfs_name}"
            if slopec_key in self.config:
                self._update_config_field(f"{slopec_key}.filtmat_data", filtmat_name)
                self._log(f"  Updated {slopec_key}.filtmat_data = {filtmat_name}")

        self._log(f"\n✓ Generated {len(filtmat_paths)} filter matrices")
        for config_id, info in filtmat_mapping.items():
            self._log(f"  {config_id}: {info['filename']}")
            self._log(f"    → Applied to: {', '.join(info['wfs_names'])}")

        return filtmat_paths

    # =======================================================================
    # STEP 4: NGS and REF Interaction Matrices
    # =======================================================================

    def step4_ngs_ref_interaction_matrices(self):
        """
        Step 4: Compute NGS (LO) and REF interaction matrices.
        
        Uses SynIM to compute:
        - NGS interaction matrices for 3 DMs (tip-tilt control)
        - REF interaction matrices for layer 0 (focus sensing)
        """
        self._log("\n" + "="*70)
        self._log("STEP 4: NGS AND REF INTERACTION MATRICES")
        self._log("="*70)

        # Save current config state
        current_yml = self._save_yaml(
            self.config,
            f"params_morfeo_step4_{self.timestamp}.yml"
        )

        params_mgr = ParamsManager(
            str(current_yml),
            root_dir=str(self.root_dir),
            verbose=self.verbose
        )

        # Compute NGS interaction matrices
        self._log("Computing NGS interaction matrices...")
        ngs_paths = params_mgr.compute_interaction_matrices(
            wfs_type='ngs',
            output_im_dir=str(self.root_dir / "synim"),
            output_rec_dir=str(self.root_dir / "synrec"),
            overwrite=self.overwrite,
            verbose=self.verbose,
            display=False
        )

        # Compute REF interaction matrices
        self._log("Computing REF interaction matrices...")
        ref_paths = params_mgr.compute_interaction_matrices(
            wfs_type='ref',
            output_im_dir=str(self.root_dir / "synim"),
            output_rec_dir=str(self.root_dir / "synrec"),
            overwrite=self.overwrite,
            verbose=self.verbose,
            display=False
        )

        self.generated_files['ngs_im_paths'] = ngs_paths
        self.generated_files['ref_im_paths'] = ref_paths

        self._log(f"✓ Generated {len(ngs_paths)} NGS interaction matrices")
        self._log(f"✓ Generated {len(ref_paths)} REF interaction matrices")

        return ngs_paths, ref_paths

    # =======================================================================
    # STEP 5: Focus Reconstruction Matrix
    # =======================================================================

    def step5_focus_reconstruction_matrix(self):
        """
        Step 5: Generate focus reconstruction matrix.
        
        Creates a reconstruction matrix for focus sensing from the 3 REF WFS,
        by stacking their interaction matrices and computing the pseudoinverse.
        """
        self._log("\n" + "="*70)
        self._log("STEP 5: FOCUS RECONSTRUCTION MATRIX")
        self._log("="*70)

        ref_im_paths = self.generated_files.get('ref_im_paths', [])
        if not ref_im_paths:
            raise RuntimeError("No REF interaction matrices found. Run step4 first.")

        # Find the layer 0 REF IM (should be only one)
        layer0_ref_im = [p for p in ref_im_paths if 'layH0.0' in str(p)]
        if not layer0_ref_im:
            raise RuntimeError("No layer 0 REF interaction matrix found")

        im_path = layer0_ref_im[0]

        output_name = f"focus_recmat_{self.timestamp}_nmodes{self.ref_focus_n_modes}"
        output_path = self.root_dir / "synrec" / f"{output_name}.fits"

        self._log(f"Generating focus reconstruction matrix (n_modes={self.ref_focus_n_modes})...")
        self._log(f"  From IM: {im_path.name}")

        # Import required classes
        from specula.data_objects.intmat import Intmat
        from specula.data_objects.recmat import Recmat

        # Load IM and stack 3 copies (for 3 REF WFS)
        intmat_obj = Intmat.restore(str(im_path))
        intmat_obj.intmat = np.vstack([
            intmat_obj.intmat,
            intmat_obj.intmat,
            intmat_obj.intmat
        ])

        self._log(f"  Stacked IM shape: {intmat_obj.intmat.shape}")

        # Generate reconstruction matrix
        recmat_obj = intmat_obj.generate_rec(
            nmodes=self.ref_focus_n_modes,
            cut_modes=0,
            interactive=False
        )

        self._log(f"  Recmat shape: {recmat_obj.recmat.shape}")

        # Save
        focus_recmat = Recmat(
            recmat=recmat_obj.recmat,
            norm_factor=recmat_obj.norm_factor
        )
        focus_recmat.save(str(output_path), overwrite=self.overwrite)

        self.generated_files['focus_recmat_path'] = output_path
        self.generated_files['focus_recmat_name'] = output_name

        # Update config
        self._update_config_field('rec_focus.recmat_object', f"{self.timestamp}/{output_name}")

        self._log(f"✓ Generated focus reconstruction matrix: {output_name}")

        return output_path

    # =======================================================================
    # STEP 6: LGS Tomographic Reconstructor
    # =======================================================================

    def step6_lgs_tomographic_reconstructor(self):
        """
        Step 6: Compute LGS tomographic reconstructor.
        
        Uses SynIM to compute MMSE reconstructor with:
        - Atmospheric turbulence modeling
        - Noise modeling (photon + readout)
        - Regularization
        """
        self._log("\n" + "="*70)
        self._log("STEP 6: LGS TOMOGRAPHIC RECONSTRUCTOR")
        self._log("="*70)

        # Save current config state
        current_yml = self._save_yaml(
            self.config,
            f"params_morfeo_step6_{self.timestamp}.yml"
        )

        params_mgr = ParamsManager(
            str(current_yml),
            root_dir=str(self.root_dir),
            verbose=self.verbose
        )

        output_rec_dir = self.root_dir / "synrec"

        self._log(f"Computing LGS tomographic reconstructor...")
        self._log(f"  r0 = {self.r0} m")
        self._log(f"  L0 = {self.L0} m")

        result = params_mgr.compute_tomographic_reconstructor(
            r0=self.r0,
            L0=self.L0,
            wfs_type='lgs',
            component_type='layer',  # Reconstruct on atmospheric layers
            noise_variance=None,  # Auto-compute from magnitude
            C_noise=None,
            output_dir=str(output_rec_dir),
            save=True,
            verbose=self.verbose
        )

        rec_filename = result.get('rec_filename')
        if rec_filename:
            rec_path = Path(rec_filename)
            rec_name = rec_path.stem

            self.generated_files['lgs_rec_path'] = rec_path
            self.generated_files['lgs_rec_name'] = rec_name

            # Update config
            self._update_config_field('tomo_polc_lgs.recmat_object', f"{self.timestamp}/{rec_name}")

        # Also save the assembled interaction matrix
        self._log("Saving assembled LGS interaction matrix...")
        im_path = params_mgr.save_assembled_interaction_matrix(
            wfs_type='lgs',
            component_type='layer',
            output_dir=str(self.root_dir / "synim"),
            overwrite=self.overwrite,
            apply_filter=True,  # Apply tip/tilt/focus filtering
            verbose=self.verbose
        )

        if im_path:
            im_path = Path(im_path)
            im_name = im_path.stem

            self.generated_files['lgs_im_assembled_path'] = im_path
            self.generated_files['lgs_im_assembled_name'] = im_name

            # Update config
            self._update_config_field('tomo_polc_lgs.intmat_object', f"{self.timestamp}/{im_name}")

        self._log(f"✓ Generated LGS tomographic reconstructor")
        self._log(f"  Reconstructor: {rec_name if rec_filename else 'N/A'}")
        self._log(f"  Interaction matrix: {im_name if im_path else 'N/A'}")
        self._log(f"  Reconstructor shape: {result['reconstructor'].shape}")

        return result

    # =======================================================================
    # STEP 7: NGS Reconstructor
    # =======================================================================

    def step7_ngs_reconstructor(self):
        """
        Step 7: Compute NGS (LO) reconstructor for tip-tilt control.
        """
        self._log("\n" + "="*70)
        self._log("STEP 7: NGS RECONSTRUCTOR")
        self._log("="*70)

        # Save current config state
        current_yml = self._save_yaml(
            self.config,
            f"params_morfeo_step7_{self.timestamp}.yml"
        )

        params_mgr = ParamsManager(
            str(current_yml),
            root_dir=str(self.root_dir),
            verbose=self.verbose
        )

        output_rec_dir = self.root_dir / "synrec"

        self._log(f"Computing NGS reconstructor...")

        result = params_mgr.compute_tomographic_reconstructor(
            r0=self.r0,
            L0=self.L0,
            wfs_type='ngs',
            component_type='dm',  # Reconstruct on DMs
            noise_variance=None,
            C_noise=None,
            output_dir=str(output_rec_dir),
            save=True,
            verbose=self.verbose
        )

        rec_filename = result.get('rec_filename')
        if rec_filename:
            rec_path = Path(rec_filename)
            rec_name = rec_path.stem

            self.generated_files['ngs_rec_path'] = rec_path
            self.generated_files['ngs_rec_name'] = rec_name

            # Update config
            self._update_config_field('tomo_ngs.recmat_object', f"{self.timestamp}/{rec_name}")

        self._log(f"✓ Generated NGS reconstructor")
        self._log(f"  Reconstructor: {rec_name if rec_filename else 'N/A'}")
        self._log(f"  Shape: {result['reconstructor'].shape}")

        return result

    # =======================================================================
    # STEP 8: Projection Matrices
    # =======================================================================

    def step8_projection_matrices(self):
        """
        Step 8: Compute projection matrices for tomography.
        
        Uses SynIM to compute:
        - Individual projection matrices for each DM-layer combination
        - Tomographic projection matrix (optimal MMSE combination)
        """
        self._log("\n" + "="*70)
        self._log("STEP 8: PROJECTION MATRICES")
        self._log("="*70)

        # Save current config state
        current_yml = self._save_yaml(
            self.config,
            f"params_morfeo_step8_{self.timestamp}.yml"
        )

        params_mgr = ParamsManager(
            str(current_yml),
            root_dir=str(self.root_dir),
            verbose=self.verbose
        )

        output_pm_dir = self.root_dir / "synpm"

        # Compute individual projection matrices
        self._log("Computing projection matrices...")
        pm_paths = params_mgr.compute_projection_matrices(
            output_dir=str(output_pm_dir),
            overwrite=self.overwrite
        )

        # Compute tomographic projection matrix
        self._log("Computing tomographic projection matrix...")
        p_opt, pm_full_dm, pm_full_layer, info = params_mgr.compute_tomographic_projection_matrix(
            output_dir=str(output_pm_dir),
            save=True,
            verbose=self.verbose
        )

        self.generated_files['pm_paths'] = pm_paths
        self.generated_files['p_opt_info'] = info

        # Update config with projection matrix
        proj_name = f"p{self.timestamp}/projmatAll"
        self._update_config_field('tomo_polc_lgs.projmat_object', proj_name)

        self._log(f"✓ Generated {len(pm_paths)} projection matrices")
        self._log(f"✓ Generated tomographic projection matrix")
        self._log(f"  Shape: {p_opt.shape}")
        self._log(f"  Optical sources: {info['n_opt_sources']}")
        self._log(f"  Regularization: {info['reg_factor']}")

        return pm_paths, info

    # =======================================================================
    # STEP 9: Update Final YAML Configuration
    # =======================================================================

    def step9_update_yaml_configuration(self):
        """
        Step 9: Update YAML configuration with all generated filenames.
        
        Creates final YAML files:
        - Complete configuration with all calibration references
        - Override to remove calibration objects for simulation
        """
        self._log("\n" + "="*70)
        self._log("STEP 9: UPDATE YAML CONFIGURATION")
        self._log("="*70)

        # Save final calibrated YAML
        final_yml = self._save_yaml(
            self.config,
            f"params_morfeo_calibrated_{self.timestamp}.yml"
        )

        # Create "remove_calib" override to clean up calibration objects
        remove_config = {
            'remove': ['reconstruction', 'projection', 
                      'layer1', 'layer2', 'layer3', 'layer4',
                      'layer5', 'layer6', 'layer7']
        }
        remove_yml = self._save_yaml(
            remove_config,
            f"override_remove_calib_{self.timestamp}.yml"
        )

        self.generated_files['final_yml'] = final_yml
        self.generated_files['remove_yml'] = remove_yml

        self._log(f"✓ Generated calibrated YAML: {final_yml.name}")
        self._log(f"✓ Generated removal override: {remove_yml.name}")

        return final_yml, remove_yml

    # =======================================================================
    # STEP 10: Generate Workflow Summary
    # =======================================================================

    def step10_generate_summary(self):
        """
        Step 10: Generate workflow summary and instructions.
        
        Creates a comprehensive summary file with:
        - All generated files and their locations
        - Commands to run final SPECULA simulation
        - Validation checks
        - Timing information
        """
        self._log("\n" + "="*70)
        self._log("STEP 10: WORKFLOW SUMMARY")
        self._log("="*70)

        summary_path = self.workflow_dir / "workflow_summary.txt"

        with open(summary_path, 'w') as f:
            f.write("MORFEO CALIBRATION WORKFLOW SUMMARY\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Workflow timestamp: {self.timestamp}\n")
            f.write(f"Base YAML: {self.base_yml_path}\n")
            f.write(f"Root directory: {self.root_dir}\n")
            f.write(f"Workflow directory: {self.workflow_dir}\n\n")

            f.write("CALIBRATION PARAMETERS:\n")
            f.write("-" * 70 + "\n")
            f.write(f"  r0: {self.r0} m\n")
            f.write(f"  L0: {self.L0} m\n")
            f.write(f"  LGS filter modes: {self.lgs_filter_n_modes} "
                   f"(filtered: {self.lgs_filter_n_modes_filtered})\n")
            f.write(f"  REF focus modes: {self.ref_focus_n_modes}\n")
            f.write(f"  NGS modes per DM: {self.ngs_n_modes_dm}\n")
            f.write(f"  REF modes per DM: {self.ref_n_modes_dm}\n")
            f.write(f"  Projection regularization: {self.proj_reg_factor}\n\n")

            f.write("GENERATED FILES:\n")
            f.write("-" * 70 + "\n\n")

            # Subaperture data
            if 'subap_tags' in self.generated_files:
                f.write("1. Subaperture Data:\n")
                for wfs, tag in self.generated_files['subap_tags'].items():
                    f.write(f"     {wfs}: {tag}\n")
                f.write("\n")
 
            # Interaction matrices
            if 'lgs_im_paths' in self.generated_files:
                f.write("2. LGS Interaction Matrices:\n")
                f.write(f"     Count: {len(self.generated_files['lgs_im_paths'])}\n")
                layer0_count = sum(1 for p in self.generated_files['lgs_im_paths'] 
                                  if 'layH0.0' in str(p))
                f.write(f"     Layer 0: {layer0_count}\n")
                f.write(f"     Other layers: {len(self.generated_files['lgs_im_paths']) - layer0_count}\n\n")
 
            # Filter matrices
            if 'filtmat_mapping' in self.generated_files:
                f.write("3. Filter Matrices:\n")
                for rot, fname in self.generated_files['filtmat_mapping'].items():
                    f.write(f"     {rot}: {fname}\n")
                f.write("\n")
 
            if 'ngs_im_paths' in self.generated_files:
                f.write("4. NGS Interaction Matrices:\n")
                f.write(f"     Count: {len(self.generated_files['ngs_im_paths'])}\n\n")

            if 'ref_im_paths' in self.generated_files:
                f.write("5. REF Interaction Matrices:\n")
                f.write(f"     Count: {len(self.generated_files['ref_im_paths'])}\n\n")

            # Focus reconstruction
            if 'focus_recmat_name' in self.generated_files:
                f.write("6. Focus Reconstruction Matrix:\n")
                f.write(f"     {self.generated_files['focus_recmat_name']}\n\n")

            # LGS reconstructor
            if 'lgs_rec_name' in self.generated_files:
                f.write("7. LGS Tomographic Reconstructor:\n")
                f.write(f"     Recmat: {self.generated_files['lgs_rec_name']}\n")
                if 'lgs_im_assembled_name' in self.generated_files:
                    f.write(f"     Intmat: {self.generated_files['lgs_im_assembled_name']}\n")
                f.write("\n")

            # NGS reconstructor
            if 'ngs_rec_name' in self.generated_files:
                f.write("8. NGS Reconstructor:\n")
                f.write(f"     {self.generated_files['ngs_rec_name']}\n\n")

            # Projection matrices
            if 'pm_paths' in self.generated_files:
                f.write("9. Projection Matrices:\n")
                f.write(f"     Individual: {len(self.generated_files['pm_paths'])}\n")
                if 'p_opt_info' in self.generated_files:
                    info = self.generated_files['p_opt_info']
                    f.write(f"     Optical sources: {info['n_opt_sources']}\n")
                    f.write(f"     Regularization: {info['reg_factor']}\n")
                f.write("\n")

            # Final YAML
            if 'final_yml' in self.generated_files:
                f.write("10. Final Configuration:\n")
                f.write(f"      Main: {self.generated_files['final_yml'].name}\n")
                if 'remove_yml' in self.generated_files:
                    f.write(f"      Override: {self.generated_files['remove_yml'].name}\n")
                f.write("\n")

            f.write("TO RUN FINAL SPECULA SIMULATION:\n")
            f.write("-" * 70 + "\n\n")
            if 'final_yml' in self.generated_files and 'remove_yml' in self.generated_files:
                f.write(f"cd {self.yml_dir}\n")
                f.write(f"specula {self.generated_files['final_yml'].name} "
                       f"{self.generated_files['remove_yml'].name}\n\n")
            else:
                f.write("(Workflow incomplete - missing final YAML)\n\n")
 
            f.write("VALIDATION CHECKS:\n")
            f.write("-" * 70 + "\n")
            checks_passed = 0
            checks_total = 0

            check_list = [
                ('subap_tags', 'Subaperture calibration'),
                ('lgs_im_paths', 'LGS interaction matrices'),
                ('filtmat_paths', 'Filter matrices'),
                ('ngs_im_paths', 'NGS interaction matrices'),
                ('ref_im_paths', 'REF interaction matrices'),
                ('focus_recmat_path', 'Focus reconstruction matrix'),
                ('lgs_rec_path', 'LGS tomographic reconstructor'),
                ('ngs_rec_path', 'NGS reconstructor'),
                ('pm_paths', 'Projection matrices'),
                ('final_yml', 'Final YAML configuration')
            ]

            for key, desc in check_list:
                checks_total += 1
                if key in self.generated_files:
                    f.write(f"  ✓ {desc}\n")
                    checks_passed += 1
                else:
                    f.write(f"  ✗ {desc}\n")
    
            f.write(f"\nPassed: {checks_passed}/{checks_total}\n\n")

        self._log(f"✓ Generated workflow summary: {summary_path}")

        # Print summary to console
        print("\n" + "="*70)
        with open(summary_path) as f:
            print(f.read())
        print("="*70 + "\n")
 
        return summary_path

    # =======================================================================
    # Main Workflow Execution
    # =======================================================================

    def run_full_workflow(self, steps: list = None):
        """
        Run the complete calibration workflow.
        
        Parameters
        ----------
        steps : list of int, optional
            Which steps to run. If None, runs all steps 0-10.
            Example: [0, 1, 2, 3] to run only first four steps
        """
        if steps is None:
            steps = list(range(1, 11))

        self._log("\n" + "="*70)
        self._log("STARTING MORFEO CALIBRATION WORKFLOW")
        self._log("="*70)
        self._log(f"Steps to run: {steps}")

        start_time = datetime.now()

        try:
            if 0 in steps:
                self.step0_influence_functions()

            if 1 in steps:
                self.step1_subap_calibration()

            if 2 in steps:
                self.step2_lgs_interaction_matrices()

            if 3 in steps:
                self.step3_filter_matrices()

            if 4 in steps:
                self.step4_ngs_ref_interaction_matrices()

            if 5 in steps:
                self.step5_focus_reconstruction_matrix()

            if 6 in steps:
                self.step6_lgs_tomographic_reconstructor()

            if 7 in steps:
                self.step7_ngs_reconstructor()

            if 8 in steps:
                self.step8_projection_matrices()

            if 9 in steps:
                self.step9_update_yaml_configuration()

            if 10 in steps:
                self.step10_generate_summary()

            end_time = datetime.now()
            duration = end_time - start_time

            self._log("\n" + "="*70)
            self._log("WORKFLOW COMPLETED SUCCESSFULLY")
            self._log("="*70)
            self._log(f"Total duration: {duration}")
            self._log(f"Workflow directory: {self.workflow_dir}")
            self._log(f"Timestamp: {self.timestamp}")
            
        except Exception as e:
            self._log("\n" + "="*70)
            self._log("WORKFLOW FAILED")
            self._log("="*70)
            self._log(f"Error: {str(e)}")
            import traceback
            traceback.print_exc()
            raise


# =======================================================================
# Example Usage
# =======================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run MORFEO calibration workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete workflow
  python run_calib_workflow_morfeo.py params_morfeo_ref_focus_tiny_initial.yml
  
  # Run only first 3 steps
  python run_calib_workflow_morfeo.py params_morfeo_ref_focus_tiny_initial.yml --steps 1 2 3
  
  # Run with overwrite
  python run_calib_workflow_morfeo.py params_morfeo_ref_focus_tiny_initial.yml --overwrite
  
  # Run on CPU
  python run_calib_workflow_morfeo.py params_morfeo_ref_focus_tiny_initial.yml --device -1
        """
    )
    parser.add_argument(
        "yaml_file",
        help="Path to initial MORFEO YAML configuration file (template)"
    )
    parser.add_argument(
        "--root-dir",
        default="/raid1/guido/PASSATA/MORFEO/",
        help="Root directory for calibration files (default: /raid1/guido/PASSATA/MORFEO/)"
    )
    parser.add_argument(
        "--output-dir",
        help="Base directory for workflow outputs (default: root_dir/workflow/)"
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        type=int,
        help="Steps to run (default: all). Example: --steps 1 2 3"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing calibration files"
    )
    parser.add_argument(
        "--device",
        type=int,
        default=0,
        help="GPU device index, -1 for CPU (default: 0)"
    )

    args = parser.parse_args()

    # Create and run workflow
    workflow = MORFEOCalibrationWorkflow(
        base_yml_path=args.yaml_file,
        root_dir=args.root_dir,
        output_base_dir=args.output_dir,
        device_idx=args.device,
        overwrite=args.overwrite,
        verbose=True
    )

    workflow.run_full_workflow(steps=args.steps)
