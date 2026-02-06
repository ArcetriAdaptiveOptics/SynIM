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
import glob
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
from specula.calib_manager import CalibManager
from specula.data_objects.intmat import Intmat
from specula.data_objects.recmat import Recmat

import synim
synim.init(device_idx=0, precision=1)
from synim.params_manager import ParamsManager
from synim.params_utils import compute_influence_functions_and_modalbases
from synim.params_utils import generate_filter_matrix_from_intmat_file

class MORFEOCalibrationWorkflow:

    def __init__(self,
                 workflow_yml_path: str,
                 base_yml_path: str = None,
                 params_dir: str = None,
                 root_dir: str = None,  # Only for override
                 output_base_dir: str = None,
                 device_idx: int = None,
                 overwrite: bool = None,
                 verbose: bool = None):
        """
        Initialize the calibration workflow.
        
        Parameters
        ----------
        workflow_yml_path : str
            Path to workflow configuration YAML (contains all parameters)
        base_yml_path : str, optional
            Path to the initial MORFEO YAML configuration file (template).
            If None, uses base_simulation_yml from workflow_yml_path (relative to params_dir).
        params_dir : str, optional
            Directory containing parameter YAML files.
            If None, uses the directory of workflow_yml_path.
        root_dir : str, optional
            Root directory for SPECULA calibration files (CalibManager root).
            If None, reads from base_simulation_yml's main.root_dir.
            Command-line override takes precedence.
        output_base_dir : str, optional
            Base directory for workflow outputs. If None, uses root_dir/workflow/
        device_idx : int, optional
            GPU device index (-1 for CPU). Overrides workflow YAML if provided.
        overwrite : bool, optional
            Whether to overwrite existing calibration files. Overrides workflow YAML.
        verbose : bool, optional
            Verbose output. Overrides workflow YAML if provided.
        """
        # Load workflow configuration first
        self.workflow_yml_path = Path(workflow_yml_path)

        if not self.workflow_yml_path.exists():
            raise FileNotFoundError(f"Workflow configuration not found: {self.workflow_yml_path}")

        with open(self.workflow_yml_path) as f:
            self.workflow_config = yaml.safe_load(f)

        # Get input files from workflow config
        input_config = self.workflow_config.get('input_files', {})

        # Resolve params directory (YAML parameter files)
        if params_dir is not None:
            # Command-line override
            self.params_dir = Path(params_dir)
        else:
            # Use workflow config
            params_dir_from_config = input_config.get('params_dir')
            if params_dir_from_config is None:
                # Default: same as workflow YAML directory
                self.params_dir = self.workflow_yml_path.parent
            else:
                self.params_dir = Path(params_dir_from_config)

        # Resolve base YAML path
        if base_yml_path is not None:
            # Command-line override
            self.base_yml_path = Path(base_yml_path)
        else:
            # Use workflow config
            base_yml_from_config = input_config.get('base_simulation_yml')
            if base_yml_from_config is None:
                raise ValueError(
                    "base_simulation_yml not specified in workflow config"
                    " and not provided via command line"
                )

            base_yml_path_obj = Path(base_yml_from_config)
            if base_yml_path_obj.is_absolute():
                self.base_yml_path = base_yml_path_obj
            else:
                # Relative to params_dir
                self.base_yml_path = self.params_dir / base_yml_path_obj

        if not self.base_yml_path.exists():
            raise FileNotFoundError(f"Base simulation YAML not found: {self.base_yml_path}")

        # Load base simulation config to extract root_dir
        with open(self.base_yml_path) as f:
            self.config = yaml.safe_load(f)

        # Resolve root directory
        if root_dir is not None:
            # Command-line override
            self.root_dir = Path(root_dir)
        else:
            # Extract from base simulation YAML
            main_config = self.config.get('main', {})
            root_dir_from_sim = main_config.get('root_dir')
            if root_dir_from_sim is None:
                raise ValueError(
                    "root_dir not found in base simulation YAML (main.root_dir)"
                    " and not provided via command line"
                )
            self.root_dir = Path(root_dir_from_sim)

        # Extract parameters from workflow config (with command-line overrides)
        exec_config = self.workflow_config['execution']
        self.device_idx = device_idx if device_idx is not None else exec_config['device_idx']
        self.overwrite = overwrite if overwrite is not None else exec_config['overwrite']
        self.verbose = verbose if verbose is not None else exec_config['verbose']

        # Create timestamped workflow directory
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if output_base_dir is None:
            output_base = self.root_dir / exec_config['output_dirs']['workflow']
            self.workflow_dir = output_base / self.timestamp
        else:
            self.workflow_dir = Path(output_base_dir) / self.timestamp

        self.workflow_dir.mkdir(parents=True, exist_ok=True)

        # Initialize directory structure
        self.yml_dir = self.workflow_dir / "yml"
        self.yml_dir.mkdir(exist_ok=True)

        # Copy and load initial YAML (already loaded above, but keep for tracking)
        self.initial_yml_path = self.yml_dir / f"params_morfeo_initial_{self.timestamp}.yml"
        shutil.copy(self.base_yml_path, self.initial_yml_path)

        # Store generated filenames for tracking
        self.generated_files = {}

        # Extract calibration parameters from workflow YAML
        self._extract_workflow_parameters()

        self._log(f"Initialized MORFEO calibration workflow")
        self._log(f"  Workflow YAML: {self.workflow_yml_path}")
        self._log(f"  Base simulation YAML: {self.base_yml_path}")
        self._log(f"  Root dir (from simulation YAML): {self.root_dir}")
        self._log(f"  Params dir: {self.params_dir}")
        self._log(f"  Workflow dir: {self.workflow_dir}")
        self._log(f"  Timestamp: {self.timestamp}")

    def _extract_workflow_parameters(self):
        """Extract calibration parameters from workflow YAML configuration."""
        wf = self.workflow_config

        # System parameters
        sys_config = wf['system']
        self.pixel_pupil = sys_config['pixel_pupil']
        self.pixel_pitch = sys_config['pixel_pitch']
        self.telescope_diameter = self.pixel_pupil * self.pixel_pitch
        self.n_petals = sys_config['n_petals']
        self.obsratio = sys_config['obsratio']
        self.r0 = sys_config['r0']
        self.L0 = sys_config['L0']

        # Extract DM parameters
        dms_config = wf['dms']
        self.dm_altitudes = []
        self.dm_nacts = []
        for dm_key in sorted(dms_config.keys()):  # dm1, dm2, dm3
            dm = dms_config[dm_key]
            self.dm_altitudes.append(dm['altitude'])
            self.dm_nacts.append(dm['n_actuators'])

        # Extract layer parameters
        layers_config = wf['layers']
        self.layer_altitudes = []
        self.layer_nacts = []
        for layer_key in sorted(layers_config.keys()):  # layer1, layer2, ...
            layer = layers_config[layer_key]
            self.layer_altitudes.append(layer['altitude'])
            self.layer_nacts.append(layer['n_actuators'])

        # WFS parameters
        wfs_config = wf['wfs']
        self.lgs_filter_n_modes = wfs_config['lgs']['filter_n_modes']
        self.lgs_filter_n_modes_filtered = wfs_config['lgs']['filter_n_modes_filtered']
        self.ref_focus_n_modes = wfs_config['ref']['focus_n_modes']
        self.ngs_n_modes_dm = wfs_config['ngs']['n_modes_dm']
        self.ref_n_modes_dm = wfs_config['ref']['n_modes_dm']

        # Reconstruction parameters
        rec_config = wf['reconstruction']
        self.sigma2_in_nm2 = rec_config['sigma2_in_nm2']
        self.noise_elong_model = rec_config['noise_elong_model']
        self.naThicknessInM = rec_config['na_thickness_in_m']
        self.tGparameter = rec_config['tG_parameter']

        # Projection parameters
        proj_config = wf['projection']
        self.proj_reg_factor = proj_config['reg_factor']

        # Calibration files
        calib_files_config = wf.get('calibration_files', {})
        self.subap_calib_yml = calib_files_config.get('subap_calib_yml',
                                                       'calib_morfeo_simplified_subaps.yml')

        if self.verbose:
            self._log("\nWorkflow parameters loaded:")
            self._log(f"  Pupil: {self.pixel_pupil} px × {self.pixel_pitch} m")
            self._log(f"  DMs: {len(self.dm_altitudes)} at altitudes {self.dm_altitudes}")
            self._log(f"  Layers: {len(self.layer_altitudes)} at altitudes {self.layer_altitudes}")
            self._log(f"  LGS filter: {self.lgs_filter_n_modes} modes "
                     f"({self.lgs_filter_n_modes_filtered} filtered)")
            self._log(f"  Projection reg_factor: {self.proj_reg_factor}")

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

        # Helper function to normalize key (remove underscores)
        def normalize_key(key):
            """Convert 'dm_1' or 'layer_1' to 'dm1' or 'layer1'"""
            return key.replace('_', '')

        # Update DMs - iterate over result keys and match to config
        for result_key in result['ifunc_tags'].keys():
            if result_key.startswith('dm'):
                config_key = normalize_key(result_key)  # dm_1 -> dm1

                if config_key in self.config:
                    self._update_config_field(f'{config_key}.ifunc_object',
                                            result['ifunc_tags'][result_key])
                    self._update_config_field(f'{config_key}.m2c_object',
                                            result['m2c_tags'][result_key])
                    # Extract n_modes from m2c tag
                    m2c_tag = result['m2c_tags'][result_key]
                    n_modes = int(m2c_tag.split('_')[-1].replace('modes', ''))
                    self._update_config_field(f'{config_key}.nmodes', n_modes)

                    if self.verbose:
                        self._log(f"  Updated {config_key}:")
                        self._log(f"    ifunc: {result['ifunc_tags'][result_key]}")
                        self._log(f"    m2c: {m2c_tag}")
                        self._log(f"    nmodes: {n_modes}")

        # Update layers - iterate over result keys and match to config
        for result_key in result['ifunc_tags'].keys():
            if result_key.startswith('layer'):
                config_key = normalize_key(result_key)  # layer_1 -> layer1

                if config_key in self.config:
                    self._update_config_field(f'{config_key}.ifunc_object',
                                            result['ifunc_tags'][result_key])
                    self._update_config_field(f'{config_key}.m2c_object',
                                            result['m2c_tags'][result_key])
                    m2c_tag = result['m2c_tags'][result_key]
                    n_modes = int(m2c_tag.split('_')[-1].replace('modes', ''))
                    self._update_config_field(f'{config_key}.nmodes', n_modes)

                    if self.verbose:
                        self._log(f"  Updated {config_key}:")
                        self._log(f"    ifunc: {result['ifunc_tags'][result_key]}")
                        self._log(f"    m2c: {m2c_tag}")
                        self._log(f"    nmodes: {n_modes}")

        # Update modal_analysis with inverse ifunc
        if result['ifunc_inv_tag']:
            self._update_config_field('modal_analysis.ifunc_inv_object',
                                    result['ifunc_inv_tag'])
            self._update_config_field('projection.ifunc_inverse_tag',
                                    result['ifunc_inv_tag'])

            if self.verbose:
                self._log(f"  Updated modal_analysis.ifunc_inv_object: {result['ifunc_inv_tag']}")

        # Count components
        n_dms = len([k for k in result['ifunc_tags'] if k.startswith('dm')])
        n_layers = len([k for k in result['ifunc_tags'] if k.startswith('layer')])

        self._log(f"\n✓ Generated influence functions for {len(result['ifunc_tags'])} components")
        self._log(f"  DMs: {n_dms}")
        self._log(f"  Layers: {n_layers}")

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

        # Resolve the calib YAML path (now relative to params_dir)
        if self.subap_calib_yml is None:
            # Auto-detect: look for calib_*_subaps.yml in params_dir
            calib_pattern = self.params_dir / "calib_morfeo_*_subaps.yml"
            matches = glob.glob(str(calib_pattern))
            if not matches:
                raise FileNotFoundError(
                    f"No subaperture calibration YAML found matching {calib_pattern}. "
                    f"Please specify 'subap_calib_yml' in workflow config."
                )
            calib_yml = Path(matches[0])
            self._log(f"Auto-detected subaperture calibration YAML: {calib_yml.name}")
        else:
            # Check if it's an absolute path
            calib_yml = Path(self.subap_calib_yml)
            if not calib_yml.is_absolute():
                # Relative to params_dir
                calib_yml = self.params_dir / self.subap_calib_yml

            if not calib_yml.exists():
                raise FileNotFoundError(
                    f"Subaperture calibration YAML not found: {calib_yml}"
                )

        self._log(f"Using subaperture calibration YAML: {calib_yml}")
        self._log(f"  (from params_dir: {self.params_dir})")

        # Load calibration YAML to extract expected tags
        with open(calib_yml) as f:
            calib_config = yaml.safe_load(f)

        # Extract expected subapdata tags
        expected_tags = {}
        for key, value in calib_config.items():
            if key.startswith('sh_subaps_'):
                wfs_name = key.replace('sh_subaps_', '')
                if 'output_tag' in value:
                    expected_tags[wfs_name] = value['output_tag']

        if not expected_tags:
            raise ValueError(
                f"No 'sh_subaps_*' entries with 'output_tag' found in {calib_yml}"
            )

        self._log(f"Expected {len(expected_tags)} subaperture configurations")

        # Check which subapdata files already exist
        cm = CalibManager(str(self.root_dir))

        existing_tags = {}
        missing_tags = {}

        for wfs_name, tag in expected_tags.items():
            # Loo for the subapdata
            subap_data_path = Path(cm.filename('subapdata', tag))
            if subap_data_path.exists():
                existing_tags[wfs_name] = tag
                self._log(f"  ✓ {wfs_name}: {tag} (already exists)")
            else:
                missing_tags[wfs_name] = tag
                self._log(f"  ✗ {wfs_name}: {tag} (needs computation)")

        # Determine if we need to run SPECULA
        need_computation = bool(missing_tags) or self.overwrite

        if need_computation:
            if self.overwrite:
                self._log("\nOverwrite mode: recomputing all subapertures")
            else:
                self._log(f"\nMissing {len(missing_tags)} subaperture files, running SPECULA...")

            # Run SPECULA calibration
            self._log("Running subaperture calibration with SPECULA...")
            self._run_specula(self.initial_yml_path, calib_yml)

            # Verify all files were created
            for wfs_name, tag in expected_tags.items():
                subap_data_path = Path(cm.filename('subapdata', tag))
                if subap_data_path.exists():
                    self._log(f"  ✓ Generated: {wfs_name} → {tag}")
                else:
                    self._log(f"  ✗ Failed to generate: {wfs_name} → {tag}")
                    raise RuntimeError(
                        f"SPECULA failed to generate subapdata for {wfs_name} (tag: {tag})"
                    )
        else:
            self._log("\n✓ All subaperture files already exist, skipping computation")

        # Store all tags (existing + newly generated)
        subap_tags = expected_tags
        self.generated_files['subap_tags'] = subap_tags
        self.generated_files['subap_calib_yml_used'] = calib_yml

        # Update config with subap tags
        for wfs_name, tag in subap_tags.items():
            slopec_key = f"slopec_{wfs_name}"
            if slopec_key in self.config:
                self._update_config_field(f"{slopec_key}.subapdata_object", tag)
                if self.verbose:
                    self._log(f"  Updated {slopec_key}.subapdata_object = {tag}")

        self._log(f"\n✓ Subaperture calibration complete")
        self._log(f"  Total configurations: {len(subap_tags)}")
        if not need_computation:
            self._log(f"  Existing: {len(existing_tags)}")
        else:
            self._log(f"  Existing: {len(existing_tags)}")
            self._log(f"  Generated: {len(missing_tags)}")

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

        lgs_im_dict = self.generated_files.get('lgs_im_paths', [])
        if not lgs_im_dict:
            raise RuntimeError("No LGS interaction matrices found. Run step2 first.")

        # Extract Path objects from dictionary values
        lgs_im_paths = [Path(p) for p in lgs_im_dict.values()]

        # Filter only layer 0 IMs (ground level, where filtering is needed)
        layer0_ims = [p for p in lgs_im_paths if 'layH0.0' in p.name]

        if not layer0_ims:
            print("Available LGS IMs:")
            for p in lgs_im_paths:
                print(f"  {p.name}")
            raise RuntimeError("No layer 0 interaction matrices found")

        self._log(f"Found {len(layer0_ims)} layer 0 interaction matrices")

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

        ref_im_dict = self.generated_files.get('ref_im_paths', {})
        if not ref_im_dict:
            raise RuntimeError("No REF interaction matrices found. Run step2 first.")

        # Extract Path objects from dictionary values
        ref_im_paths = [Path(p) for p in ref_im_dict.values()]

        # Filter only layer 0 IMs (ground level, where filtering is needed)
        layer0_ref_im = [p for p in ref_im_paths if 'layH0.0' in p.name]

        if not layer0_ref_im:
            raise RuntimeError("No layer 0 REF interaction matrix found")

        im_path = layer0_ref_im[0]

        output_name = f"focus_recmat_{self.timestamp}_nmodes{self.ref_focus_n_modes}"
        output_path = self.root_dir / "synrec" / f"{output_name}.fits"

        self._log(f"Generating focus reconstruction matrix (n_modes={self.ref_focus_n_modes})...")
        self._log(f"  From IM: {im_path.name}")

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
            f.write(f"Workflow YAML: {self.workflow_yml_path}\n")
            f.write(f"Base simulation YAML: {self.base_yml_path}\n\n")
            
            f.write("CONFIGURATION HIERARCHY:\n")
            f.write("-" * 70 + "\n")
            f.write(f"1. Workflow YAML: {self.workflow_yml_path.name}\n")
            f.write(f"   - Provides: calibration parameters, execution settings\n")
            f.write(f"2. Base simulation YAML: {self.base_yml_path.name}\n")
            f.write(f"   - Provides: root_dir, system configuration, sources, WFS setup\n")
            f.write(f"3. Command-line: (if used, overrides both YAMLs)\n\n")
            
            f.write("DIRECTORY STRUCTURE:\n")
            f.write("-" * 70 + "\n")
            f.write(f"Root directory (from simulation YAML main.root_dir):\n")
            f.write(f"  {self.root_dir}\n")
            f.write(f"  ├── pupilstop/  (pupil masks)\n")
            f.write(f"  ├── ifunc/      (influence functions)\n")
            f.write(f"  ├── m2c/        (modal-to-commands matrices)\n")
            f.write(f"  ├── synim/      (interaction matrices)\n")
            f.write(f"  ├── synrec/     (reconstructors)\n")
            f.write(f"  ├── synpm/      (projection matrices)\n")
            f.write(f"  └── workflow/   (workflow outputs)\n\n")
            
            f.write(f"Params directory:\n")
            f.write(f"  {self.params_dir}\n")
            f.write(f"  Contains: simulation YAMLs, calibration YAMLs\n\n")
            
            f.write(f"Workflow directory (this run):\n")
            f.write(f"  {self.workflow_dir}\n")
            f.write(f"  ├── yml/  (generated YAML configurations)\n")
            f.write(f"  └── workflow_summary.txt (this file)\n\n")

            f.write("INPUT FILES:\n")
            f.write("-" * 70 + "\n")
            if 'subap_calib_yml_used' in self.generated_files:
                f.write(f"  Subaperture calib: {self.generated_files['subap_calib_yml_used']}\n")
            f.write("\n")

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
                f.write(f"     Other layers:"
                        f" {len(self.generated_files['lgs_im_paths']) - layer0_count}\n\n")

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
            steps = list(range(0, 11))

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
Configuration hierarchy (from lowest to highest priority):
  1. Workflow YAML (params_workflow_morfeo.yml)
  2. Base simulation YAML (params_morfeo_*.yml - includes root_dir)
  3. Command-line arguments (override both YAMLs)

Directory structure:
  root_dir/           # From base simulation YAML (main.root_dir)
    ├── pupilstop/    # Pupil masks
    ├── ifunc/        # Influence functions
    ├── m2c/          # Modal-to-commands matrices
    ├── synim/        # Interaction matrices
    ├── synrec/       # Reconstructors
    ├── synpm/        # Projection matrices
    └── workflow/     # Workflow outputs (timestamped)
        └── YYYYMMDD_HHMMSS/
            ├── yml/  # YAML configurations
            └── workflow_summary.txt

  params_dir/         # Parameter YAML files (default: workflow YAML directory)
    ├── params_morfeo_*.yml  (contains root_dir in main.root_dir)
    ├── calib_morfeo_*_subaps.yml
    └── params_workflow_morfeo.yml

Examples:
  # Run with all defaults (root_dir from base simulation YAML)
  python run_calib_workflow_morfeo.py params_workflow_morfeo.yml
  
  # Override root_dir (ignore value in simulation YAML)
  python run_calib_workflow_morfeo.py params_workflow_morfeo.yml \\
      --root-dir /raid1/guido/PASSATA/MORFEOtest/
  
  # Override params directory
  python run_calib_workflow_morfeo.py params_workflow_morfeo.yml \\
      --params-dir /home/guido/pythonLib/SPECULA_scripts/params_morfeo_test/
  
  # Override base simulation YAML
  python run_calib_workflow_morfeo.py params_workflow_morfeo.yml \\
      --base-yml /path/to/custom_params_morfeo.yml
  
  # Run only specific steps
  python run_calib_workflow_morfeo.py params_workflow_morfeo.yml --steps 0 1 2
  
  # Run with overwrite on CPU
  python run_calib_workflow_morfeo.py params_workflow_morfeo.yml \\
      --overwrite --device -1
        """
    )
    parser.add_argument(
        "workflow_config",
        help="Path to workflow configuration YAML"
    )
    parser.add_argument(
        "--base-yml",
        help="Path to base MORFEO simulation YAML (overrides workflow config). "
             "The root_dir is extracted from main.root_dir in this file."
    )
    parser.add_argument(
        "--root-dir",
        help="Root directory for SPECULA calibration files "
             "(overrides main.root_dir from base simulation YAML)"
    )
    parser.add_argument(
        "--params-dir",
        help="Directory containing parameter YAML files (overrides workflow config)"
    )
    parser.add_argument(
        "--output-dir",
        help="Base directory for workflow outputs (default: root_dir/workflow/)"
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        type=int,
        help="Steps to run (default: all). Example: --steps 0 1 2"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing calibration files (overrides workflow config)"
    )
    parser.add_argument(
        "--device",
        type=int,
        help="GPU device index, -1 for CPU (overrides workflow config)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output (overrides workflow config)"
    )

    args = parser.parse_args()

    # Create and run workflow
    workflow = MORFEOCalibrationWorkflow(
        workflow_yml_path=args.workflow_config,
        base_yml_path=args.base_yml,
        params_dir=args.params_dir,
        root_dir=args.root_dir,
        output_base_dir=args.output_dir,
        device_idx=args.device,
        overwrite=args.overwrite,
        verbose=args.verbose
    )

    workflow.run_full_workflow(steps=args.steps)
