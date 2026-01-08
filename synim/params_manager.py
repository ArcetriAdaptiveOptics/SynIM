import os
import re
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# *** Import xp, cpuArray, to_xp, float_dtype ***
from synim import (
    xp, cpuArray, to_xp, float_dtype, cpu_float_dtype,
    default_target_device_idx, global_precision
)

import synim.synim as synim
import synim.synpm as synpm

# Import all utility functions from params_utils
from synim.utils import *
from synim.params_utils import *

import specula
specula.init(device_idx=-1, precision=1)

from specula.calib_manager import CalibManager
from specula.data_objects.intmat import Intmat
from specula.data_objects.recmat import Recmat
from specula.lib.modal_base_generator import compute_ifs_covmat

verbose_gpu_flag = False  # Global flag for verbose GPU memory printing

class ParamsManager:
    """
    Class for managing parameters needed to compute interaction matrices 
    for all combinations of DMs and WFSs without redundant loading.
    """

    def __init__(self, params_file, root_dir=None, verbose=False):
        """
        Initialize the manager and load all common parameters.
        
        Args:
            params_file (str or dict): Path to YAML/PRO configuration file or dictionary
            root_dir (str, optional): Root directory to override in params
            verbose (bool): Whether to print detailed information
        """
        # Load configuration
        self.params_file = params_file
        if isinstance(params_file, str):
            self.params = parse_params_file(params_file)
        else:
            self.params = params_file

        self.verbose = verbose

        self.target_device_idx = default_target_device_idx
        self.precision = global_precision

        # Set root_dir if provided
        if root_dir:
            if 'main' in self.params:
                self.params['main']['root_dir'] = root_dir
                if verbose:
                    print(f"Root directory set to: {self.params['main']['root_dir']}")

        # initialize im_dir, pm_dir, rec_dir, cov_dir
        self.im_dir = self.params['main']['root_dir'] + '/synim/'
        self.pm_dir = self.params['main']['root_dir'] + '/synpm/'
        self.rec_dir = self.params['main']['root_dir'] + '/synrec/'
        self.cov_dir = self.params['main']['root_dir'] + '/covariance/'

        # Initialize the CalibManager
        self.main_params = self.params['main']
        self.cm = CalibManager(self.main_params['root_dir'])

        # Extract common parameters
        self.pixel_pupil = self.main_params['pixel_pupil']
        self.pixel_pitch = self.main_params['pixel_pitch']
        self.pup_diam_m = self.pixel_pupil * self.pixel_pitch

        # Determine configuration type
        self.is_simple_config = is_simple_config(self.params)

        # Load all WFS and DM configurations
        self.wfs_list = extract_wfs_list(self.params)
        self.dm_list = extract_dm_list(self.params)

        # *** Extract and cache projection/reconstructor params ***
        self._extract_advanced_params()

        if self.verbose:
            print(f"Found {len(self.wfs_list)} WFS(s) and {len(self.dm_list)} DM(s)")
            for wfs in self.wfs_list:
                print(f"  WFS: {wfs['name']} (type: {wfs['type']}, index: {wfs['index']})")
            for dm in self.dm_list:
                print(f"  DM: {dm['name']} (index: {dm['index']})")
           # Validate optical sources if present
            validate_opt_sources(self.params, verbose=True)

        # Pre-load pupil mask
        self.pup_mask = self._load_pupil_mask()

        # Cache for DM and WFS parameters (loaded on demand)
        self.dm_cache = {}  # Key: dm_index, Value: dict with dm_array, dm_mask, etc.
        self.wfs_cache = {}  # Key: (wfs_type, wfs_index), Value: dict with wfs params

    def _load_pupil_mask(self):
        """
        Load or create the pupil mask.
        
        Returns:
            numpy.ndarray: Pupil mask array
        """
        pup_mask = None
        pupilstop_params = None
        if 'pupilstop' in self.params:
            pupilstop_params = self.params['pupilstop']
            if self.verbose:
                print("Found 'pupilstop' in params")
        elif 'pupil_stop' in self.params:
            pupilstop_params = self.params['pupil_stop']
            if self.verbose:
                print("Found 'pupil_stop' in params")
        if pupilstop_params is None:
            raise ValueError("Pupilstop parameters not found in configuration")
        pup_mask = load_pupilstop(self.cm, pupilstop_params, self.pixel_pupil, self.pixel_pitch,
                                verbose=self.verbose)

        # If no pupilstop defined, create a default circular pupil
        if pup_mask is None:
            simul_params = SimulParams(pixel_pupil=self.pixel_pupil, pixel_pitch=self.pixel_pitch)
            pupilstop = Pupilstop(
                simul_params,
                mask_diam=1.0,
                obs_diam=0.0,
                target_device_idx=-1,
                precision=0
            )
            pup_mask = pupilstop.A

        # *** Convert to xp with float_dtype ***
        pup_mask = to_xp(xp, pup_mask, dtype=float_dtype)

        print('---> valid pixels: ', int(np.sum(pup_mask > 0.5)))

        return pup_mask


    def _extract_advanced_params(self):
        """
        Extract and cache advanced parameters for projection and reconstruction.
        
        This method consolidates all parameter extraction logic in one place,
        determining the source of parameters based on priority:
        1. YAML modern format (projection/reconstructor sections)
        2. IDL legacy format (modalrec1 section)
        3. Defaults
        
        Sets the following attributes:
        - self.projection_params: dict or None
        - self.proj_reg_factor: float
        - self.reconstructor_params: dict or None (for future use)
        """

        # ==================== PROJECTION PARAMETERS ====================
        # Priority 1: Try new YAML 'projection' section
        self.projection_params = extract_projection_params(self.params)

        if self.projection_params is not None:
            # Use reg_factor from projection section
            self.proj_reg_factor = float(self.projection_params['reg_factor'])

            if self.verbose:
                print(f"\n✓ Found 'projection' section:")
                print(f"  Optical sources: {len(self.projection_params['opt_sources'])}")
                print(f"  reg_factor: {self.proj_reg_factor}")
                if 'ifunc_inverse_tag' in self.projection_params:
                    print(f"  ifunc_inverse_tag: {self.projection_params['ifunc_inverse_tag']}")
        else:
            # Priority 2: Try IDL-style 'modalrec1.proj_regFactor'
            modalrec1 = self.params.get('modalrec1', {})
            proj_reg_factor_val = modalrec1.get('proj_regFactor', None)
            if proj_reg_factor_val is not None:
                self.proj_reg_factor = float(proj_reg_factor_val)
            else:
                self.proj_reg_factor = None

            if self.proj_reg_factor is None:
                # Priority 3: Default
                self.proj_reg_factor = 1e-8

                if self.verbose:
                    print(f"\n⚠ No projection parameters found, using defaults:")
                    print(f"  reg_factor: {self.proj_reg_factor}")
            else:
                if self.verbose:
                    print(f"\n✓ Using IDL-style projection parameters:")
                    print(f"  reg_factor from modalrec1.proj_regFactor: {self.proj_reg_factor}")

        # ==================== RECONSTRUCTOR PARAMETERS ====================
        self.reconstructor_params = {}

        # Try YAML format first (reconstruction section)
        if 'reconstruction' in self.params:
            recon_section = self.params['reconstruction']
            self.sigma2_in_nm2 = float(recon_section.get('sigma2inNm2', 2e4))

            # Extract NGS and REF mode configurations ***
            self.ngs_n_modes_dm = recon_section.get('ngs_n_modes_dm', [])
            self.ngs_n_modes_layer = recon_section.get('ngs_n_modes_layer', [])
            self.ref_n_modes_dm = recon_section.get('ref_n_modes_dm', [])
            self.ref_n_modes_layer = recon_section.get('ref_n_modes_layer', [])

            if self.verbose:
                print(f"  Using YAML 'reconstruction' section")
                print(f"  sigma2_in_nm2: {self.sigma2_in_nm2}")

            self.reconstructor_params = {
                'r0': recon_section.get('r0', None),
                'L0': recon_section.get('L0', None),
                'reg_factor': recon_section.get('reg_factor', 1e-8),
                'wfs_type': recon_section.get('wfs_type', 'lgs'),
                'component_type': recon_section.get('component_type', 'layer'),
                'sigma2_in_nm2': self.sigma2_in_nm2,
                'ngs_n_modes_dm': self.ngs_n_modes_dm,
                'ngs_n_modes_layer': self.ngs_n_modes_layer,
                'ref_n_modes_dm': self.ref_n_modes_dm,
                'ref_n_modes_layer': self.ref_n_modes_layer
            }

        # Fallback to IDL-style modalrec1
        elif 'modalrec1' in self.params:
            modalrec1 = self.params['modalrec1']
            self.sigma2_in_nm2 = float(modalrec1.get('sigma2inNm2', 2e4))

            if self.verbose:
                print(f"  Using IDL 'modalrec1' section")
                print(f"  sigma2_in_nm2: {self.sigma2_in_nm2}")

            self.reconstructor_params = {
                'sigma2_in_nm2': self.sigma2_in_nm2
            }

            self.ngs_n_modes_dm = []
            self.ngs_n_modes_layer = []
            self.ref_n_modes_dm = []
            self.ref_n_modes_layer = []

        # Default fallback
        else:
            self.sigma2_in_nm2 = 2e4
            self.ngs_n_modes_dm = []
            self.ngs_n_modes_layer = []
            self.ref_n_modes_dm = []
            self.ref_n_modes_layer = []

            if self.verbose:
                print(f"  Using default sigma2_in_nm2: {self.sigma2_in_nm2}")

            self.reconstructor_params = {
                'sigma2_in_nm2': self.sigma2_in_nm2
            }


    def _print_gpu_memory(self, label="", verbose=True):
        """Print current GPU memory usage."""
        if not verbose or not hasattr(xp, 'cuda'):
            return

        try:
            import cupy as cp
            mempool = cp.get_default_memory_pool()
            print(f"\n  --- GPU Memory {label}:")
            print(f"        Used:  {mempool.used_bytes() / 1024**3:.2f} GB")
            print(f"        Total: {mempool.total_bytes() / 1024**3:.2f} GB")
        except:
            pass


    def count_mcao_stars(self):
        """
        Count the number of LGS, NGS, reference stars, DMs, optimisation optics,
        science stars and layers in the parameter configuration, similar to
        count_mcao_stars of IDL.

            Returns:
        dict: Dictionary with counts.
        """
        def count_keys_with_prefix(params, prefix):
            # Count keys starting with prefix and followed by a number
            return len([k for k in params.keys() if re.match(rf"{prefix}\d+$", k)])

        params = self.params

        out = {}
        out['n_lgs'] = count_keys_with_prefix(params, 'source_lgs')
        out['n_ngs'] = count_keys_with_prefix(params, 'source_ngs')
        out['n_ref'] = count_keys_with_prefix(params, 'source_ref')

        if 'source_opt1' in params:
            out['n_opt'] = count_keys_with_prefix(params, 'source_opt')
            out['n_gs'] = out['n_lgs'] + out['n_ngs'] + out['n_opt']
        else:
            out['n_opt'] = 0
            out['n_gs'] = out['n_lgs'] + out['n_ngs']

        if 'dm1' in params:
            out['n_dm'] = count_keys_with_prefix(params, 'dm')
        else:
            out['n_dm'] = 1

        if 'science_source1' in params:
            out['n_star'] = count_keys_with_prefix(params, 'science_source')
        else:
            out['n_star'] = 0

        if 'layer1' in params:
            out['n_rec_layer'] = count_keys_with_prefix(params, 'layer')
        else:
            out['n_rec_layer'] = 0

        return out

    def get_component_params(self,
                             component_idx,
                             is_layer=False,
                             cut_start_mode=False,
                             n_modes_to_use=None,
                             xp_local=xp,
                             float_dtype_local=float_dtype):
        """
        Get DM or layer parameters, loading from cache if available.

        Args:
            component_idx (int): Index of the DM or layer to load
            is_layer (bool): Whether to load a layer instead of a DM
            cut_start_mode (bool): Whether to select modes from start_mode to start_mode+n_modes
            n_modes_to_use (int, optional): Number of modes to use starting from start_mode.
                                        If None, uses all available modes after start_mode.
            xp_local: array module to use (numpy or cupy)
            float_dtype_local: data type for arrays (e.g., float_dtype or cpu_float_dtype)

        Returns:
            dict: DM or layer parameters
        """
        component_type = "layer" if is_layer else "dm"

        # Create cache key that includes n_modes_to_use
        cache_key_base = f"{component_type}_{component_idx}"
        if cut_start_mode:
            cache_key = f"{cache_key_base}_cut"
        else:
            cache_key = cache_key_base
        if n_modes_to_use is not None:
            cache_key += f"_{n_modes_to_use}"

        if cache_key in self.dm_cache:
            return self.dm_cache[cache_key]

        # Try with index first (dm1, dm2, layer1, layer2...)
        component_key = f"{component_type}{component_idx}"

        if component_key not in self.params:
            if component_type == "dm" and component_idx == 1 and "dm" in self.params:
                component_key = "dm"
            else:
                raise ValueError(f"{component_type.capitalize()} {component_idx} not found")

        component_params = self.params[component_key]

        dm_array, dm_mask = load_influence_functions(
            self.cm, component_params, self.pixel_pupil, verbose=self.verbose
        )

        # *** Select modes from start_mode to start_mode+n_modes ***
        if cut_start_mode or n_modes_to_use is not None:
            if 'start_mode' in component_params:
                start_mode = component_params['start_mode']
            else:
                start_mode = 0
            total_modes = dm_array.shape[2]

            if n_modes_to_use is not None:
                # Use specified number of modes starting from start_mode
                end_mode = start_mode + n_modes_to_use
                if end_mode > total_modes:
                    end_mode = total_modes
                    if self.verbose:
                        print(f"  Warning: Requested {n_modes_to_use} modes from {start_mode},"
                            f" but only {total_modes - start_mode} available")

                dm_array = dm_array[:, :, start_mode:end_mode]

                if self.verbose:
                    print(f"  *** AFTER CUT: dm_array.shape = {dm_array.shape} ***")
                    print(f"  Selected modes [{start_mode}:{end_mode}] "
                        f"({end_mode - start_mode} modes out of {total_modes})")
            else:
                # Use all modes after start_mode (old behavior)
                dm_array = dm_array[:, :, start_mode:]

                if self.verbose:
                    print(f"  Selected modes from {start_mode} onwards "
                        f"({total_modes - start_mode} modes)")

        # Convert to xp with float_dtype
        dm_array = to_xp(xp_local, dm_array, dtype=float_dtype_local)
        dm_mask = to_xp(xp_local, dm_mask, dtype=float_dtype_local)

        self.dm_cache[cache_key] = {
            'dm_array': dm_array,
            'dm_mask': dm_mask,
            'dm_height': component_params.get('height', 0.0),
            'dm_rotation': component_params.get('rotation', 0.0),
            'component_key': component_key
        }

        return self.dm_cache[cache_key]

    def get_wfs_params(self, wfs_type=None, wfs_index=None,
                       xp_local=xp):
        """
        Get WFS parameters for a specific WFS.
        
        Args:
            wfs_type (str, optional): Type of WFS ('sh', 'pyr') or source type ('lgs', 'ngs', 'ref')
            wfs_index (int, optional): Index of the WFS (1-based)
            xp_local: array module to use (numpy or cupy)
            
        Returns:
            dict: WFS parameters
        """
        # Create a cache key
        cache_key = (wfs_type, wfs_index)

        # Check if already loaded in cache
        if cache_key in self.wfs_cache:
            return self.wfs_cache[cache_key]

        # WFS selection logic
        selected_wfs = None
        source_type = None

        if self.verbose:
            print("WFS -- Looking for WFS parameters...")

        # Simple SCAO configuration case
        if self.is_simple_config:
            if self.verbose:
                print("     Simple SCAO configuration detected")
            if len(self.wfs_list) > 0:
                selected_wfs = self.wfs_list[0]
                source_type = 'ngs'  # Default for simple configs
                if self.verbose:
                    print(f"     Using WFS: {selected_wfs['name']} of type {selected_wfs['type']}")
            else:
                raise ValueError("No WFS configuration found in the configuration file.")
        else:
            # Complex MCAO configuration
            if self.verbose:
                print("     Complex MCAO configuration detected")

            # Case 1: wfs_type specifies the sensor type ('sh', 'pyr')
            if wfs_type in ['sh', 'pyr']:
                if self.verbose:
                    print(f"     Looking for WFS of type: {wfs_type}")
                matching_wfs = [wfs for wfs in self.wfs_list if wfs['type'] == wfs_type]

                if wfs_index is not None:
                    # Try to find specific index
                    for wfs in matching_wfs:
                        if wfs['index'] == str(wfs_index):
                            selected_wfs = wfs
                            if self.verbose:
                                print(f"     Found WFS with specified index: {wfs['name']}")
                            break

                # If no specific index found, use the first one
                if selected_wfs is None and matching_wfs:
                    selected_wfs = matching_wfs[0]
                    if self.verbose:
                        print(f"     Using first WFS of type {wfs_type}: {selected_wfs['name']}")

            # Case 2: wfs_type specifies the source type ('lgs', 'ngs', 'ref')
            elif wfs_type in ['lgs', 'ngs', 'ref']:
                source_type = wfs_type
                if self.verbose:
                    print(f"     Looking for WFS associated with {wfs_type} source")

                # Pattern for WFS names corresponding to the source type
                pattern = f"sh_{source_type}"
                matching_wfs = [wfs for wfs in self.wfs_list if pattern in wfs['name']]

                if wfs_index is not None:
                    # Try to find specific index within the source type
                    target_name = f"{pattern}{wfs_index}"
                    for wfs in matching_wfs:
                        if wfs['name'] == target_name:
                            selected_wfs = wfs
                            if self.verbose:
                                print(f"     Found WFS with specified index: {wfs['name']}")
                            break

                # If no specific index found, use the first one
                if selected_wfs is None and matching_wfs:
                    selected_wfs = matching_wfs[0]
                    if self.verbose:
                        print(f"     Using first WFS for {wfs_type}: {selected_wfs['name']}")

            # Case 3: Only wfs_index is specified (no wfs_type)
            elif wfs_index is not None:
                if self.verbose:
                    print(f"     Looking for WFS with index: {wfs_index}")
                for wfs in self.wfs_list:
                    if wfs['index'] == str(wfs_index):
                        selected_wfs = wfs
                        if self.verbose:
                            print(f"     Found WFS with specified index: {wfs['name']}")
                        break

            # Case 4: No specific criteria, use the first available WFS
            if selected_wfs is None and self.wfs_list:
                selected_wfs = self.wfs_list[0]
                if self.verbose:
                    print(f"     Using first available WFS: {selected_wfs['name']}")

        # If no WFS found, raise error
        if selected_wfs is None:
            raise ValueError("No matching WFS configuration found in the configuration file.")

        wfs_key = selected_wfs['name']
        wfs_params = selected_wfs['config']
        wfs_type_detected = selected_wfs['type']

        # Determine source type from WFS name if not already set
        if source_type is None:
            source_type = determine_source_type(wfs_key)

        if self.verbose:
            print(f"     Source type determined: {source_type}")

        # Extract WFS parameters
        wfs_rotation = wfs_params.get('rotation', wfs_params.get('rotAnglePhInDeg', 0.0))
        x_shift = wfs_params.get('xShiftPhInPixel', 0.0)
        y_shift = wfs_params.get('yShiftPhInPixel', 0.0)
        wfs_translation = wfs_params.get('translation', [x_shift, y_shift])
        wfs_magnification = wfs_params.get('magnification', 1.0)
        wfs_anamorphosis_90 = wfs_params.get('anamorph90', 1.0)
        wfs_anamorphosis_45 = wfs_params.get('anamorph45', 1.0)
        wfs_fov_arcsec = wfs_fov_from_config(wfs_params)

        if wfs_type_detected == 'sh':
            # Shack-Hartmann specific parameters
            wfs_nsubaps = wfs_params.get('subap_on_diameter', 0)
        else:
            # Pyramid specific parameters
            wfs_nsubaps = wfs_params.get('pup_diam', 0)

        # Load SubapData for valid subapertures
        idx_valid_sa = find_subapdata(
            self.cm, wfs_params, wfs_key, self.params, verbose=self.verbose
        )

        # *** Convert to xp if not None ***
        if idx_valid_sa is not None:
            idx_valid_sa = to_xp(xp_local, idx_valid_sa)

        # Guide star parameters
        if source_type == 'lgs':
            # LGS is at finite height
            # Try to get height from source or use typical LGS height
            gs_height = None

            # Check if there's a specific source for this WFS and try to get height
            source_match = re.search(r'(lgs\d+)', wfs_key)
            if source_match:
                source_key = f'source_{source_match.group(1)}'
                if source_key in self.params:
                    gs_height = self.params[source_key].get('height', None)

            # If still no height, use default
            if gs_height is None:
                gs_height = 90000.0  # Default LGS height in meters
        else:
            # NGS and REF are at infinite distance
            gs_height = float('inf')

        # Get source polar coordinates
        gs_pol_coo = extract_source_coordinates(self.params, wfs_key)

        # Create dictionary with WFS parameters
        wfs_data = {
            'wfs_key': wfs_key,
            'wfs_type': wfs_type_detected,
            'wfs_nsubaps': wfs_nsubaps,
            'wfs_rotation': wfs_rotation,
            'wfs_translation': wfs_translation,
            'wfs_magnification': wfs_magnification,
            'wfs_anamorphosis_90': wfs_anamorphosis_90,
            'wfs_anamorphosis_45': wfs_anamorphosis_45,
            'wfs_fov_arcsec': wfs_fov_arcsec,
            'gs_pol_coo': gs_pol_coo,
            'gs_height': gs_height,
            'idx_valid_sa': idx_valid_sa,
            'source_type': source_type
        }

        # Cache for future use
        self.wfs_cache[cache_key] = wfs_data

        return wfs_data

    def prepare_interaction_matrix_params(self, wfs_type=None, wfs_index=None,
                                        dm_index=None, layer_index=None):
        """
        Prepare parameters for computing an interaction matrix.
        
        Args:
            wfs_type (str, optional): Type of WFS ('sh', 'pyr') or source type ('lgs', 'ngs', 'ref')
            wfs_index (int, optional): Index of the WFS (1-based)
            dm_index (int, optional): Index of the DM (1-based)
            layer_index (int, optional): Index of the Layer (1-based)
            
        Returns:
            dict: Parameters ready to be passed to synim.interaction_matrix
        """
        # Check that only one of dm_index or layer_index is specified
        if dm_index is not None and layer_index is not None:
            raise ValueError("Cannot specify both dm_index and layer_index")

        if dm_index is None and layer_index is None:
            raise ValueError("Must specify either dm_index or layer_index")

        # Get DM or Layer parameters
        if dm_index is not None:
            component_params = self.get_component_params(dm_index,
                                                         xp_local=np,
                                                         float_dtype_local=cpu_float_dtype)
            component_key = component_params['component_key']
        else:
            component_params = self.get_component_params(layer_index, is_layer=True,
                                                         xp_local=np,
                                                         float_dtype_local=cpu_float_dtype)
            component_key = component_params['component_key']

        # Get WFS parameters
        wfs_params = self.get_wfs_params(wfs_type, wfs_index,
                                         xp_local=np)
        # Combine them into a single dictionary with all needed parameters
        params = {
            'pup_diam_m': self.pup_diam_m,
            'pup_mask': self.pup_mask,
            'dm_array': component_params['dm_array'],
            'dm_mask': component_params['dm_mask'],
            'dm_height': component_params['dm_height'],
            'dm_rotation': component_params['dm_rotation'],
            'wfs_key': wfs_params['wfs_key'],
            'wfs_type': wfs_params['wfs_type'],
            'wfs_nsubaps': wfs_params['wfs_nsubaps'],
            'wfs_rotation': wfs_params['wfs_rotation'],
            'wfs_translation': wfs_params['wfs_translation'],
            'wfs_magnification': wfs_params['wfs_magnification'],
            'wfs_fov_arcsec': wfs_params['wfs_fov_arcsec'],
            'gs_pol_coo': wfs_params['gs_pol_coo'],
            'gs_height': wfs_params['gs_height'],
            'idx_valid_sa': wfs_params['idx_valid_sa'],
            'dm_key': component_key if dm_index is not None else None,
            'layer_key': component_key if layer_index is not None else None,
            'component_key': component_key,
            'source_type': wfs_params['source_type']
        }

        return params

    def compute_interaction_matrix(self, wfs_type=None, wfs_index=None, dm_index=None,
                                layer_index=None, verbose=None, display=False):
        """
        Compute an interaction matrix for a specific WFS-DM/Layer combination.

        Args:
            wfs_type (str, optional): Type of WFS ('sh', 'pyr') or source type ('lgs', 'ngs', 'ref')
            wfs_index (int, optional): Index of the WFS (1-based)
            dm_index (int, optional): Index of the DM (1-based)
            layer_index (int, optional): Index of the Layer (1-based)
            verbose (bool, optional): Override the class's verbose setting
            display (bool): Whether to display plots

        Returns:
            numpy.ndarray: Computed interaction matrix
        """
        # Check that only one of dm_index or layer_index is specified
        if dm_index is not None and layer_index is not None:
            raise ValueError("Cannot specify both dm_index and layer_index")

        if dm_index is None and layer_index is None:
            raise ValueError("Must specify either dm_index or layer_index")

        # Use class verbose setting if not overridden
        verbose_flag = self.verbose if verbose is None else verbose

        # Prepare parameters
        params = self.prepare_interaction_matrix_params(wfs_type, wfs_index, dm_index, layer_index)

        component_name = params['dm_key'] if dm_index is not None else params['layer_key']

        if verbose_flag:
            print("Computing interaction matrix with parameters:")
            print(f"      WFS: {params['wfs_key']}, Component: {component_name}")
            print(f"      WFS type: {params['wfs_type']}, nsubaps: {params['wfs_nsubaps']}")
            print(f"      Guide star: {params['gs_pol_coo']} at height {params['gs_height']} m")

        # Calculate the interaction matrix
        im = synim.interaction_matrix(
            pup_diam_m=params['pup_diam_m'],
            pup_mask=to_xp(xp, params['pup_mask'], dtype=float_dtype),
            dm_array=to_xp(xp, params['dm_array'], dtype=float_dtype),
            dm_mask=to_xp(xp, params['dm_mask'], dtype=float_dtype),
            dm_height=params['dm_height'],
            gs_pol_coo=params['gs_pol_coo'],
            gs_height=params['gs_height'],
            dm_rotation=params['dm_rotation'],
            wfs_nsubaps=params['wfs_nsubaps'],
            wfs_rotation=params['wfs_rotation'],
            wfs_translation=params['wfs_translation'],
            wfs_mag_global=params['wfs_magnification'],
            wfs_fov_arcsec=params['wfs_fov_arcsec'],
            idx_valid_sa=to_xp(xp, params['idx_valid_sa'], dtype=float_dtype),
            verbose=verbose_flag,
            display=display
        )

        # *** Convert to CPU for saving ***
        im = cpuArray(im)

        return im

    def compute_interaction_matrices(self, output_im_dir, output_rec_dir,
                                wfs_type=None, overwrite=False, verbose=None,
                                display=False):
        """
        Compute and save interaction matrices for all combinations of WFSs and DMs/Layers.
        Uses multi-WFS optimization when possible.
        
        Args:
            output_im_dir (str): Output directory for saved matrices
            output_rec_dir (str): Output directory for reconstruction matrices
            wfs_type (str, optional): Type of WFS ('ngs', 'lgs', 'ref') to use
            overwrite (bool, optional): Whether to overwrite existing files
            verbose (bool, optional): Override the class's verbose setting
            display (bool, optional): Whether to display plots
            
        Returns:
            dict: Dictionary mapping WFS-Component pairs to saved interaction matrix paths
        """
        saved_matrices = {}
        os.makedirs(output_im_dir, exist_ok=True)
        os.makedirs(output_rec_dir, exist_ok=True)

        verbose_flag = self.verbose if verbose is None else verbose

        self._print_gpu_memory("(initial)", verbose_gpu_flag)

        # Filter WFS list by type if specified
        if wfs_type is not None:
            filtered_wfs_list = [wfs for wfs in self.wfs_list if wfs_type in wfs['name']]
        else:
            filtered_wfs_list = self.wfs_list

        # Combine DM and layer in a single components list
        components = []
        for dm in self.dm_list:
            components.append({
                'type': 'dm',
                'index': int(dm['index']),
                'name': dm['name']
            })
        for layer in extract_layer_list(self.params):
            components.append({
                'type': 'layer',
                'index': int(layer['index']),
                'name': layer['name']
            })

        if verbose_flag:
            print(f"Computing interaction matrices for {len(filtered_wfs_list)} WFS(s) "
                f"and {len(components)} component(s)")
            for wfs in filtered_wfs_list:
                print(f"  WFS: {wfs['name']} (type: {wfs['type']}, index: {wfs['index']})")
            for comp in components:
                print(f"  {comp['type'].upper()}: {comp['name']} (index: {comp['index']})")

        # ==================== PROCESS EACH COMPONENT ====================
        for component in components:
            comp_idx = component['index']
            comp_name = component['name']
            comp_type = component['type']

            if verbose_flag:
                print(f"\n{'='*60}")
                print(f"Loading {comp_type.upper()} {comp_name} (index {comp_idx})")
                print(f"{'='*60}")

            # Load component parameters once
            component_params = self.get_component_params(
                comp_idx,
                is_layer=(comp_type == 'layer'),
                cut_start_mode=False,
                xp_local=np,
                float_dtype_local=cpu_float_dtype
            )

            if verbose_flag:
                print(f"  Component array shape: {component_params['dm_array'].shape}")
                print(f"  Component height: {component_params['dm_height']}")
                print(f"  Component rotation: {component_params['dm_rotation']}")

            # ========== BUILD WFS CONFIGURATIONS FOR MULTI-WFS ==========
            wfs_configs = []
            wfs_to_compute = []  # Track which WFS need computation

            for wfs in filtered_wfs_list:
                wfs_idx = int(wfs['index'])
                wfs_name = wfs['name']
                source_type = determine_source_type(wfs_name)

                # Generate filename
                im_filename = generate_im_filename(
                    self.params_file,
                    wfs_type=source_type,
                    wfs_index=wfs_idx,
                    dm_index=comp_idx if comp_type == 'dm' else None,
                    layer_index=comp_idx if comp_type == 'layer' else None
                )

                im_path = os.path.join(output_im_dir, im_filename)

                # Check if file exists
                if os.path.exists(im_path) and not overwrite:
                    if verbose_flag:
                        print(f"  File already exists: {im_filename}")
                    saved_matrices[f"{wfs_name}_{comp_name}"] = im_path
                    continue
                else:
                    if verbose_flag:
                        if os.path.exists(im_path):
                            print(f"  Overwriting existing file: {im_filename}")
                        else:
                            print(f"  Will compute and save: {im_filename}")

                # Get WFS parameters
                wfs_params = self.get_wfs_params(source_type, wfs_idx,
                                                 xp_local=np)

                # Add to configurations for multi-WFS computation
                wfs_configs.append({
                    'nsubaps': wfs_params['wfs_nsubaps'],
                    'rotation': wfs_params['wfs_rotation'],
                    'translation': wfs_params['wfs_translation'],
                    'magnification': wfs_params['wfs_magnification'],
                    'fov_arcsec': wfs_params['wfs_fov_arcsec'],
                    'idx_valid_sa': wfs_params['idx_valid_sa'],
                    'gs_pol_coo': wfs_params['gs_pol_coo'],
                    'gs_height': wfs_params['gs_height'],
                    'name': wfs_name
                })

                wfs_to_compute.append({
                    'name': wfs_name,
                    'index': wfs_idx,
                    'source_type': source_type,
                    'filename': im_filename,
                    'path': im_path,
                    'params': wfs_params
                })

            # ========== COMPUTE INTERACTION MATRICES ==========
            if len(wfs_configs) > 0:
                if verbose_flag:
                    print(f"\n  Computing {len(wfs_configs)} interaction"
                          f" matrices using multi-WFS...")

                self._print_gpu_memory("(before IM computation)", verbose_gpu_flag)

                # Use multi-WFS computation
                im_dict, derivatives_info = synim.interaction_matrices_multi_wfs(
                    pup_diam_m=self.pup_diam_m,
                    pup_mask=to_xp(xp, self.pup_mask, dtype=float_dtype),
                    dm_array=to_xp(xp, component_params['dm_array'], dtype=float_dtype),
                    dm_mask=to_xp(xp, component_params['dm_mask'], dtype=float_dtype),
                    dm_height=component_params['dm_height'],
                    dm_rotation=component_params['dm_rotation'],
                    wfs_configs=wfs_configs,
                    verbose=verbose_flag,
                    specula_convention=True,
                    im_on_cpu=True,
                    minimize_memory=True
                )

                self._print_gpu_memory("(after IM computation)", verbose_gpu_flag)

                if verbose_flag:
                    print(f"  Used {derivatives_info['workflow'].upper()} workflow")

                # Save each computed interaction matrix
                for wfs_info in wfs_to_compute:
                    wfs_name = wfs_info['name']
                    im = im_dict[wfs_name]

                    # *** Convert to CPU for transpose and saving ***
                    im = cpuArray(im)

                    if verbose_flag:
                        print(f"  Saving {wfs_name}: {wfs_info['filename']}")
                        print(f"    IM shape: {im.shape}")

                    # Display if requested
                    if display:
                        plt.figure(figsize=(10, 8))
                        plt.imshow(im, cmap='viridis')
                        plt.colorbar()
                        plt.title(f"IM: {wfs_name} - {comp_name}")
                        plt.tight_layout()
                        plt.show()

                    # Save interaction matrix
                    config_name = (os.path.basename(self.params_file).split('.')[0]
                                if isinstance(self.params_file, str) else "config")
                    pupdata_tag = f"{config_name}_{wfs_info['params']['wfs_type']}"
                    pupdata_tag += f"_{wfs_info['params']['wfs_nsubaps']}"

                    intmat_obj = Intmat(
                        im,
                        pupdata_tag=pupdata_tag,
                        norm_factor=1.0
                    )
                    intmat_obj.save(wfs_info['path'])

                    saved_matrices[f"{wfs_name}_{comp_name}"] = wfs_info['path']

        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Completed: {len(saved_matrices)} interaction matrices computed/loaded")
            print(f"{'='*60}\n")

        return saved_matrices


    def assemble_interaction_matrices(self, wfs_type='ngs', output_im_dir=None,
                                    component_type='dm', save=False,
                                    apply_filter=True):
        """
        Assemble interaction matrices for a specific type of WFS into
        a single full interaction matrix.

        NOTE: This method applies filtering AFTER loading from disk,
        ensuring consistent treatment of both computed and pre-existing matrices.

        Args:
            wfs_type (str): The type of WFS to assemble matrices for ('ngs', 'lgs', 'ref')
            output_im_dir (str, optional): Directory where IM files are stored
            component_type (str): Type of component to assemble ('dm' or 'layer')
            save (bool): Whether to save the assembled matrix to disk
            apply_filter (bool): Whether to apply filtmat_tag filtering if present
            
        Returns:
            tuple: (im_full, n_slopes_per_wfs, mode_indices, component_indices) -
                Assembled matrix and associated parameters
        """
        # Set up output directory
        if output_im_dir is None:
            raise ValueError("output_im_dir must be specified.")

        # Validate component_type
        if component_type not in ['dm', 'layer']:
            raise ValueError("component_type must be either 'dm' or 'layer'")

        # Count WFSs of the specified type in the configuration
        wfs_list = [wfs for wfs in self.wfs_list if wfs_type in wfs['name']]
        n_wfs = len(wfs_list)

        if self.verbose:
            print(f"Found {n_wfs} {wfs_type.upper()} WFSs")

        # Get the number of slopes per WFS (from idx_valid_sa)
        n_tot_slopes = 0
        n_slopes_list = []
        for wfs in wfs_list:
            wfs_params = self.get_wfs_params(wfs_type, int(wfs['index']),
                                             xp_local=np)
            if wfs_params['idx_valid_sa'] is not None:
                # Each valid subaperture produces X and Y slopes
                n_slopes_this_wfs = len(wfs_params['idx_valid_sa']) * 2
            else:
                print(f"Warning: idx_valid_sa is None for WFS {wfs['name']},"
                      f" setting n_slopes to 0")
                n_slopes_this_wfs = 0
            n_slopes_list.append(n_slopes_this_wfs)
            n_tot_slopes += n_slopes_this_wfs

        n_slopes_per_wfs = n_tot_slopes // n_wfs

        if self.verbose:
            print(f"Each WFS has {n_slopes_per_wfs} slopes")

        # Get component list based on type
        if component_type == 'dm':
            component_list = self.dm_list
        else:
            component_list = extract_layer_list(self.params)

        # Component indices and start modes based on config
        component_indices = []
        component_start_modes = []
        mode_indices = []
        total_modes = 0

        # Check for reconstruction section based on WFS type ***
        use_reconstruction_modes = False
        reconstruction_n_modes = []

        if wfs_type in ['ngs', 'ref']:
            # Try to get modes from reconstruction section
            if wfs_type == 'ngs':
                if component_type == 'dm':
                    reconstruction_n_modes = self.ngs_n_modes_dm
                else:  # layer
                    reconstruction_n_modes = self.ngs_n_modes_layer
            else:  # ref
                if component_type == 'dm':
                    reconstruction_n_modes = self.ref_n_modes_dm
                else:  # layer
                    reconstruction_n_modes = self.ref_n_modes_layer

            if reconstruction_n_modes:
                use_reconstruction_modes = True
                if self.verbose:
                    print(f"Using {wfs_type.upper()} mode configuration"
                          f" from reconstruction section")
                    print(f"  {component_type.upper()} modes: {reconstruction_n_modes}")

        # Check if modal_combination exists for this WFS type and component type
        modal_key = f'modes_{wfs_type}_{component_type}'
        has_modal_combination = ('modal_combination' in self.params and
                                modal_key in self.params['modal_combination'])

        if use_reconstruction_modes:
            # Use reconstruction section for NGS/REF
            for i, n_modes in enumerate(reconstruction_n_modes):
                # *** SKIP components with 0 modes ***
                if n_modes == 0:
                    if self.verbose:
                        print(f"  Skipping component {i+1}: 0 modes requested")
                    continue

                if i < len(component_list):
                    comp = component_list[i]
                    comp_idx = int(comp['index'])

                    comp_config = self.params[comp['name']]
                    start_mode = comp_config.get('start_mode', 0)
                    component_start_modes.append(start_mode)
                    component_indices.append(comp_idx)

                    # Modes from start_mode to n_modes
                    modes = list(range(start_mode, start_mode + n_modes))
                    mode_indices.append(modes)
                    total_modes += len(modes)

                    if self.verbose:
                        print(f"  Component {comp_idx}: {n_modes} modes from {start_mode}")

        elif has_modal_combination:
            if self.verbose:
                print(f"Using modal_combination for {component_type}s")

            modes_config = self.params['modal_combination'][modal_key]
            for i, n_modes in enumerate(modes_config):
                # *** SKIP components with 0 modes ***
                if n_modes == 0:
                    if self.verbose:
                        print(f"  Skipping component {i+1}: 0 modes in modal_combination")
                    continue

                if n_modes > 0:
                    # Check for start_mode in component config
                    comp_key = f'{component_type}{i+1}'
                    if comp_key in self.params and 'start_mode' in self.params[comp_key]:
                        comp_start_mode = self.params[comp_key]['start_mode']
                    else:
                        comp_start_mode = 0

                    n_modes_effective = n_modes - comp_start_mode
                    component_start_modes.append(comp_start_mode)
                    component_indices.append(i + 1)
                    mode_indices.append(
                        list(range(comp_start_mode, comp_start_mode + n_modes_effective))
                    )
                    total_modes += n_modes_effective
        else:
            # Default: use all components with all modes
            if self.verbose:
                print(f"No modal_combination found, using all {component_type}s with all modes")

            for comp in component_list:
                comp_idx = int(comp['index'])

                comp_key = f'{component_type}{comp_idx}'
                if comp_key not in self.params:
                    # Fallback: if dm1 does not exist but dm does and comp_idx==1
                    if component_type == 'dm' and comp_idx == 1 and 'dm' in self.params:
                        comp_key = 'dm'
                    else:
                        raise ValueError(f"Component {comp_key} not found in configuration")

                comp_config = self.params[comp_key]

                # Get total modes from config
                n_modes = comp_config.get('nmodes', None)
                if n_modes is None:
                    raise ValueError(f"nmodes not specified for {comp_key}")

                # Check for start_mode
                if 'start_mode' in comp_config:
                    comp_start_mode = comp_config['start_mode']
                else:
                    comp_start_mode = 0

                n_modes_effective = n_modes - comp_start_mode

                component_start_modes.append(comp_start_mode)
                component_indices.append(comp_idx)
                mode_indices.append(
                    list(range(comp_start_mode, comp_start_mode + n_modes_effective))
                )
                total_modes += n_modes_effective

        n_tot_modes = total_modes  # Total number of modes

        if self.verbose:
            print(f"Total modes: {n_tot_modes}, Total slopes: {n_tot_slopes}")
            print(f"{component_type.upper()} indices for {wfs_type}: {component_indices}")
            print(f"{component_type.upper()} start modes: {component_start_modes}")

        # *** DETERMINE DTYPE FROM FIRST IM ***
        # Load first IM to get dtype
        first_wfs_idx = int(wfs_list[0]['index'])
        first_comp_idx = component_indices[0]
        first_im_filename = generate_im_filename(
            self.params_file,
            wfs_type=wfs_type,
            wfs_index=first_wfs_idx,
            dm_index=first_comp_idx if component_type == 'dm' else None,
            layer_index=first_comp_idx if component_type == 'layer' else None
        )
        first_im_path = os.path.join(output_im_dir, first_im_filename)
        first_intmat_obj = Intmat.restore(first_im_path)
        im_dtype = first_intmat_obj.intmat.dtype

        # Create the full interaction matrix
        im_full = np.zeros((n_tot_slopes, n_tot_modes), dtype=im_dtype)

        # ==================== ASSEMBLE MATRICES ====================
        # Load and assemble the interaction matrices
        for ii in range(n_wfs):
            for jj, comp_ind in enumerate(component_indices):
                # Get the appropriate mode indices for this component
                mode_idx = mode_indices[jj] # absolute mode indices
                comp_start_mode = component_start_modes[jj]  # start_mode offset

                # Generate and load the interaction matrix file
                im_filename = generate_im_filename(
                    self.params_file,
                    wfs_type=wfs_type,
                    wfs_index=ii+1,
                    dm_index=comp_ind if component_type == 'dm' else None,
                    layer_index=comp_ind if component_type == 'layer' else None
                )
                im_path = os.path.join(output_im_dir, im_filename)

                if self.verbose:
                    print(f"--> Indices: WFS {ii+1}, Component {comp_ind}")
                    print(f"--> Loading IM: {im_filename}")

                # Load the interaction matrix
                intmat_obj = Intmat.restore(im_path)

                if self.verbose:
                    print(f"    IM shape: {intmat_obj.intmat.shape}")

                # Use intmat_obj.intmat as base
                im = intmat_obj.intmat

                # *** APPLY FILTER HERE (after loading) ***
                if apply_filter:
                    im = self._apply_slopes_filter(
                        im,
                        wfs_type,
                        ii+1,  # WFS index (1-based)
                        verbose=self.verbose
                    )

                # Convert absolute mode indices to relative
                n_modes_available = im.shape[1]

                # mode_idx contains absolute indices [start_mode, start_mode+1, ...]
                # IM contains all modes [0, 1, 2, ...]
                # We need to extract [start_mode:start_mode+n_eff] from IM

                # Calculate relative mode indices (which start from 0)
                relative_mode_idx = [mi for mi in mode_idx if mi < n_modes_available]

                if len(relative_mode_idx) == 0:
                    if self.verbose:
                        print(f"    Warning: No valid modes for WFS {ii+1}, Component {comp_ind}")
                    continue

                # Calculate starting column in full IM (absolute indices)
                # but extract from im using relative indices
                start_col_in_full = sum(len(mode_indices[k]) for k in range(jj))
                n_modes_this_comp = len(relative_mode_idx)

                # Fill the appropriate section of the full interaction matrix
                if ii == 0:
                    idx_start = 0
                else:
                    idx_start = sum(n_slopes_list[:ii])

                try:
                    im_full[idx_start:idx_start+n_slopes_list[ii],
                        start_col_in_full:start_col_in_full+n_modes_this_comp] = \
                        im[:, relative_mode_idx]
                except Exception as e:
                    print(f"Error assembling IM for WFS {ii+1}, Component {comp_ind}: {e}")
                    print(f"  IM shape: {im.shape}")
                    print(f"  relative_mode_idx: {relative_mode_idx}")
                    print(f"  Filling im_full[{idx_start}:{idx_start+n_slopes_list[ii]}, "
                        f"{start_col_in_full}:{start_col_in_full+n_modes_this_comp}]")
                    raise

        # Display summary
        if self.verbose:
            print(f"\nAssembled interaction matrix shape: {im_full.shape}")

        # Save the full interaction matrix if requested
        if save:
            filter_suffix = "_filtered" if apply_filter else ""
            output_filename = f"im_full_{wfs_type}_{component_type}{filter_suffix}.npy"
            np.save(os.path.join(output_im_dir, output_filename), cpuArray(im_full))
            if self.verbose:
                print(f"Saved full interaction matrix to {output_filename}")

        return im_full, n_slopes_per_wfs, mode_indices, component_indices


    def _apply_slopes_filter(self, im, wfs_type, wfs_index, verbose=False):
        """
        Apply slopes filtering to interaction matrix if filtmat_tag is present.
        
        Implements the offline TT filtering as in IDL:
        1. Load filtmat with shape (n_modes, n_slopes, 2)
        2. filt_int_mat = filtmat[:,:,0].T
        3. filt_rec_mat = filtmat[:,:,1]
        4. m = im @ filt_rec_mat
        5. im0 = m @ filt_int_mat
        6. im_filtered = im - im0
        
        Args:
            im (np.ndarray): Interaction matrix (n_slopes, n_modes)
            wfs_type (str): WFS type ('lgs', 'ngs', 'ref')
            wfs_index (int): WFS index (1-based)
            verbose (bool): Whether to print information
            
        Returns:
            np.ndarray: Filtered interaction matrix
        """
        # Determine slopec key
        if wfs_type == 'lgs':
            slopec_key = f'slopec{wfs_index}'
        elif wfs_type == 'ngs':
            slopec_key = f'slopec_ngs{wfs_index}'
        elif wfs_type == 'ref':
            slopec_key = f'slopec_ref{wfs_index}'
        else:
            if verbose:
                print(f"  Unknown WFS type: {wfs_type}, cannot apply filter.")
            return im  # No filtering

        # Check if filtmat_tag exists
        if slopec_key not in self.params:
            slopec_key = f'slopec_lgs{wfs_index}'  # Fallback key
            if slopec_key not in self.params:
                if verbose:
                    print(f"  No slopec parameters for key: {slopec_key}")
                return im

        slopec_params = self.params[slopec_key]

        if 'filtmat_tag' not in slopec_params:
            if 'filtmat_data' not in slopec_params:
                if verbose:
                    print(f"  No filtmat_tag or filtmat_data in slopec parameters for key:"
                          f" {slopec_key}")
                return im
            else:
                filtmat_key = 'filtmat_data'
        else:
            filtmat_key = 'filtmat_tag'

        # Check if filtName is present (inline filtering, already applied)
        if 'filtName' in slopec_params:
            if verbose:
                print(f"  Filter already applied inline (filtName present)")
            return im

        filtmat_tag = slopec_params[filtmat_key]

        if verbose:
            print(f"  Applying offline slopes filter: {filtmat_tag}")

        # Load filter matrix from CalibManager
        # if filtmat directory exists, use it, otherwise use default data directory
        filtmat_dir = os.path.join(self.cm.root_dir, 'filtmat')
        if os.path.exists(filtmat_dir):
            filtmat_path = os.path.join(filtmat_dir, f'{filtmat_tag}.fits')
        else:
            filtmat_path = os.path.join(self.cm.root_dir, 'data', f'{filtmat_tag}.fits')

        if not os.path.exists(filtmat_path):
            if verbose:
                print(f"    Warning: Filter file not found: {filtmat_path}")
            return im

        with fits.open(filtmat_path) as hdul:
            filtmat = hdul[0].data  # Shape: (2, n_slopes, n_modes)

        debug_verbose = False
        if debug_verbose:
            print(f"    Loaded filtmat shape: {filtmat.shape}")

        # Extract filter components
        # IDL convention: filtmat[*,*,0] is intmat, filtmat[*,*,1] is recmat
        filt_int_mat = filtmat[0, :, :]  # (n_slopes, n_modes)
        filt_rec_mat = filtmat[1, :, :].T  # (n_modes, n_slopes)

        if debug_verbose:
            print(f"    filt_int_mat shape: {filt_int_mat.shape}")
            print(f"    filt_rec_mat shape: {filt_rec_mat.shape}")
            print(f"    IM shape before filtering: {im.shape}")


        m = filt_rec_mat @ im  # (n_filter_modes, n_dm_modes)
        im0 = filt_int_mat @ m   # (n_slopes, n_dm_modes)
        im_filtered = im - im0

        if debug_verbose:
            print(f"    IM shape after filtering: {im_filtered.shape}")
            rms_before = np.sqrt(np.mean(im**2))
            rms_after = np.sqrt(np.mean(im_filtered**2))
            print(f"    RMS before: {rms_before:.4e}, after: {rms_after:.4e}")
            print(f"    Filtered power: {100*(1-rms_after/rms_before):.2f}%")

        return im_filtered


    def _load_base_inv_array(self, verbose):
        """Extract base_inv_array loading logic."""
        base_inv_array = None

        # Priority 1: projection.ifunc_inverse_tag (YAML style)**
        if self.projection_params is not None and 'ifunc_inverse_tag' in self.projection_params:
            if verbose:
                print(f"Loading from projection.ifunc_inverse_tag")
            ifunc_inv_tag = self.projection_params['ifunc_inverse_tag']
            dm_inv_params = {'ifunc_tag': ifunc_inv_tag}
            base_inv_array, _ = load_influence_functions(
                self.cm, dm_inv_params, self.pixel_pupil,
                verbose=verbose, is_inverse_basis=True
            )
            return to_xp(xp, base_inv_array, dtype=float_dtype)

        # Priority 2: modal_analysis
        if 'modal_analysis' in self.params:
            if verbose:
                print("Loading from modal_analysis")
            base_inv_array, _ = load_influence_functions(
                self.cm, self.params['modal_analysis'], self.pixel_pupil,
                verbose=verbose, is_inverse_basis=True
            )

        # Priority 3: modalrec1.tag_ifunc4proj
        elif 'modalrec1' in self.params and 'tag_ifunc4proj' in self.params['modalrec1']:
            if verbose:
                print("Loading from modalrec1.tag_ifunc4proj")
            inv_tag = f"{self.params['modalrec1']['tag_ifunc4proj']}_inv"
            dm_inv_params = {'ifunc_tag': inv_tag}
            if 'tag_m2c4proj' in self.params['modalrec1']:
                dm_inv_params['m2c_tag'] = self.params['modalrec1']['tag_m2c4proj']
            base_inv_array, _ = load_influence_functions(
                self.cm, dm_inv_params, self.pixel_pupil,
                verbose=verbose, is_inverse_basis=True
            )

        # Priority 4: modalrec.tag_ifunc4proj
        elif 'modalrec' in self.params and 'tag_ifunc4proj' in self.params['modalrec']:
            if verbose:
                print("Loading from modalrec.tag_ifunc4proj")
            inv_tag = f"{self.params['modalrec']['tag_ifunc4proj']}_inv"
            dm_inv_params = {'ifunc_tag': inv_tag}
            if 'tag_m2c4proj' in self.params['modalrec']:
                dm_inv_params['m2c_tag'] = self.params['modalrec']['tag_m2c4proj']
            base_inv_array, _ = load_influence_functions(
                self.cm, dm_inv_params, self.pixel_pupil,
                verbose=verbose, is_inverse_basis=True
            )

        if base_inv_array is None:
            raise ValueError("No valid base_inv_array found in configuration")

        # *** Convert to xp with float_dtype ***
        base_inv_array = to_xp(xp, base_inv_array, dtype=float_dtype)

        return base_inv_array


    def _save_projection_matrix(self, pm, source_info, comp_name, saved_matrices, verbose):
        """Save projection matrix helper."""
        if verbose:
            print(f"    Saving {source_info['name']}: {source_info['filename']}")
            print(f"    PM shape: {pm.shape}")

        plot_debug = False
        if plot_debug:
            plt.figure(figsize=(10, 8))
            plt.imshow(cpuArray(pm), cmap='viridis')
            plt.colorbar()
            plt.title(f"PM: {source_info['name']} - {comp_name}")
            plt.tight_layout()
            plt.show()

        config_name = (os.path.basename(self.params_file).split('.')[0]
                    if isinstance(self.params_file, str) else "config")
        pupdata_tag = f"{config_name}_opt{source_info['index']}"

        pm_obj = Intmat(pm, pupdata_tag=pupdata_tag, norm_factor=1.0,
                        target_device_idx=None, precision=None)
        pm_obj.save(source_info['path'])
        saved_matrices[f"{source_info['name']}_{comp_name}"] = source_info['path']


    def compute_projection_matrices(self, output_dir=None, overwrite=False,
                                    verbose=None):
        """
        Compute and save projection matrices for all combinations of optical sources and DMs/layers.
        Uses multi-base optimization when possible.
        
        Args:
            output_dir (str, optional): Output directory for saved matrices
            overwrite (bool, optional): Whether to overwrite existing files
            verbose (bool, optional): Override the class's verbose setting

        Returns:
            dict: Dictionary mapping Source-Component pairs to saved projection matrix paths
        """
        saved_matrices = {}

        # Set up directories
        if output_dir is None:
            raise ValueError("output_dir must be specified.")

        os.makedirs(output_dir, exist_ok=True)

        # Use verbose flag from instance if not overridden
        verbose_flag = self.verbose if verbose is None else verbose

        # ==================== GET COMPONENTS ====================
        components = []
        for dm in self.dm_list:
            components.append({'type': 'dm',
                               'index': int(dm['index']),
                               'name': dm['name']})
        for layer in extract_layer_list(self.params):
            components.append({'type': 'layer',
                               'index': int(layer['index']),
                               'name': layer['name']})

        # ==================== GET OPTICAL SOURCES ====================
        opt_sources = extract_opt_list(self.params)

        if verbose_flag:
            print(f"Computing PMs for {len(opt_sources)} sources × {len(components)} components")

        # ==================== CHECK WHICH FILES NEED COMPUTATION ====================
        # First pass: check which files exist to avoid loading base unnecessarily
        all_files_exist = True
        sources_needing_computation = []

        for component in components:
            comp_idx = component['index']
            comp_type = component['type']

            for source in opt_sources:
                source_idx = source['index']

                pm_filename = generate_pm_filename(
                    self.params_file, opt_index=source_idx,
                    dm_index=comp_idx if comp_type == 'dm' else None,
                    layer_index=comp_idx if comp_type == 'layer' else None
                )
                pm_path = os.path.join(output_dir, pm_filename)

                if not os.path.exists(pm_path) or overwrite:
                    all_files_exist = False
                    sources_needing_computation.append({
                        'component': component,
                        'source': source
                    })

                    if verbose_flag:
                        status = "will overwrite" if os.path.exists(pm_path) else "missing"
                        print(f"  {pm_filename}: {status}")
                else:
                    # File exists and not overwriting - add to saved_matrices
                    saved_matrices[f"{source['name']}_{component['name']}"] = pm_path

        # ==================== EARLY EXIT IF ALL FILES EXIST ====================
        if all_files_exist and not overwrite:
            if verbose_flag:
                print(f"\n✓ All {len(saved_matrices)} projection matrices already exist")
                print(f"  Set overwrite=True to recompute")
            return saved_matrices

        # ==================== LOAD BASE (ONLY IF NEEDED) ====================
        if verbose_flag:
            print(f"\n  Loading base inverse array (needed for"
                  f" {len(sources_needing_computation)} computations)...")

        base_inv_array = self._load_base_inv_array(verbose_flag)

        # ==================== PROCESS EACH COMPONENT ====================
        for component in components:
            comp_idx = component['index']
            comp_name = component['name']
            comp_type = component['type']

            if verbose_flag:
                print(f"\n{'='*60}")
                print(f"Processing {comp_type.upper()} {comp_name} (index {comp_idx})")
                print(f"{'='*60}")

            # Load component ONCE
            component_params = self.get_component_params(
                comp_idx, is_layer=(comp_type == 'layer'), cut_start_mode=True,
                xp_local=np, float_dtype_local=cpu_float_dtype
            )

            # ==================== CHECK ALL SOURCES ====================
            sources_to_compute = []
            for source in opt_sources:
                source_name = source['name']
                source_idx = source['index']
                source_config = source['config']

                pm_filename = generate_pm_filename(
                    self.params_file, opt_index=source_idx,
                    dm_index=comp_idx if comp_type == 'dm' else None,
                    layer_index=comp_idx if comp_type == 'layer' else None
                )
                pm_path = os.path.join(output_dir, pm_filename)

                if os.path.exists(pm_path) and not overwrite:
                    if verbose_flag:
                        print(f"  Exists: {pm_filename}")
                    continue

                if 'polar_coordinate' in source_config:
                    gs_pol_coo = source_config.get('polar_coordinate')
                elif 'polar_coordinates' in source_config:
                    gs_pol_coo = source_config.get('polar_coordinates')
                else:
                    raise ValueError(f"No polar coordinates found for source {source_name}")

                sources_to_compute.append({
                    'name': source_name,
                    'index': source_idx,
                    'filename': pm_filename,
                    'path': pm_path,
                    'gs_pol_coo': gs_pol_coo,
                    'gs_height': source_config.get('height', float('inf'))
                })

            if not sources_to_compute:
                continue

            # ==================== COMPUTE PROJECTION MATRICES ====================
            # *** REMOVE the all_gs_same optimization - not applicable here ***

            if verbose_flag:
                print(f"\n  Computing {len(sources_to_compute)} PMs individually")

            for source_info in sources_to_compute:
                if verbose_flag:
                    print(f"\n  Processing {source_info['name']}:")

                # Transpose ONCE using ORIGINAL pupil mask
                base_inv_array_transposed = synpm.transpose_base_array_for_specula(
                    to_xp(xp, base_inv_array, dtype=float_dtype),
                    to_xp(xp, self.pup_mask, dtype=float_dtype),
                    verbose=verbose_flag
                )

                pm = synpm.projection_matrix(
                    pup_diam_m=self.pup_diam_m,
                    pup_mask=to_xp(xp, self.pup_mask, dtype=float_dtype),
                    dm_array=to_xp(xp, component_params['dm_array'], dtype=float_dtype),
                    dm_mask=to_xp(xp, component_params['dm_mask'], dtype=float_dtype),
                    base_inv_array=base_inv_array_transposed,
                    dm_height=component_params['dm_height'],
                    dm_rotation=component_params['dm_rotation'],
                    base_rotation=0.0,
                    base_translation=(0.0, 0.0),
                    base_magnification=(1.0, 1.0),
                    gs_pol_coo=source_info['gs_pol_coo'],
                    gs_height=source_info['gs_height'],
                    verbose=verbose_flag,
                    specula_convention=True
                )

                pm = cpuArray(pm)

                self._save_projection_matrix(
                    pm, source_info, comp_name, saved_matrices, verbose_flag
                )

        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Completed: {len(saved_matrices)} projection matrices")
            print(f"{'='*60}\n")

        return saved_matrices


    def compute_projection_matrix(self, reg_factor=1e-8, output_dir=None, save=False):
        """
        Assemble 4D projection matrices from individual PM files and
        calculate the final projection matrix using the full DM and layer matrices.
        
        Optimized version: loads each DM/layer only once and iterates over optical sources.
        
        Args:
            reg_factor (float, optional): Regularization factor for the pseudoinverse calculation.
                Default is 1e-8.
            output_dir (str, optional): Directory where PM files are stored and where
                assembled matrices will be saved.
            save (bool): Whether to save the assembled matrices to disk
        
        Returns:
            popt (numpy.ndarray): Final projection matrix (n_dm_modes, n_layer_modes)
            pm_full_dm (numpy.ndarray): Full DM projection matrix 
                (n_opt_sources, n_dm_modes, n_pupil_modes)
            pm_full_layer (numpy.ndarray): Full Layer projection matrix 
                (n_opt_sources, n_layer_modes, n_pupil_modes)
        """

        # Set up output directory
        if output_dir is None:
            raise ValueError("output_dir must be specified.")

        os.makedirs(output_dir, exist_ok=True)

        # Extract all necessary lists
        opt_sources = extract_opt_list(self.params)
        dm_list = extract_dm_list(self.params)
        layer_list = extract_layer_list(self.params)

        # check if list are None or empty
        if opt_sources is None or len(opt_sources) == 0:
            raise ValueError("No optical sources found in configuration.")
        if (dm_list is None or len(dm_list) == 0):
            raise ValueError("No DMs found in configuration.")
        if (layer_list is None or len(layer_list) == 0):
            raise ValueError("No layers found in configuration.")

        n_opt = len(opt_sources)
        n_dm = len(dm_list)
        n_layer = len(layer_list)

        if self.verbose:
            print(f"Found {n_opt} optical sources, {n_dm} DMs, and {n_layer} layers")

        # Extract weights
        weights_array = np.array([src['config'].get('weight', 1.0) for src in opt_sources])

        # Initialize output arrays (we'll determine size from first PM loaded)
        pm_full_dm = None
        pm_full_layer = None

        # ==================== PROCESS DMs ====================
        if n_dm > 0:
            if self.verbose:
                print("\n=== Processing DM Projection Matrices ===")

            for jj, dm in enumerate(dm_list):
                dm_index = dm['index']
                dm_name = dm['name']

                if self.verbose:
                    print(f"\nLoading PMs for {dm_name} (index {dm_index})...")

                # Load PMs for all optical sources for this DM
                dm_pms_list = []

                for ii, opt_source in enumerate(opt_sources):
                    opt_index = opt_source['index']

                    pm_filename = generate_pm_filename(
                        self.params_file, opt_index=opt_index, dm_index=dm_index
                    )

                    if pm_filename is None:
                        raise ValueError(f"Could not generate filename"
                                         f" for opt{opt_index}, dm{dm_index}")

                    pm_path = os.path.join(output_dir, pm_filename)

                    if not os.path.exists(pm_path):
                        raise FileNotFoundError(f"File {pm_path} does not exist")

                    if self.verbose:
                        print(f"  Loading opt{opt_index}: {pm_filename}")

                    intmat_obj = Intmat.restore(pm_path)
                    dm_pms_list.append(intmat_obj.intmat)

                # Stack all optical sources for this DM
                # dm_pms_list[i] has shape (n_dm_modes_i, n_pupil_modes)
                # We stack along axis 0 to get (n_opt, n_dm_modes_i, n_pupil_modes)
                dm_stack = np.stack(dm_pms_list, axis=0)

                if self.verbose:
                    print(f"  Stacked PM shape for {dm_name}: {dm_stack.shape}")

                # Concatenate along DM modes (axis 1)
                if pm_full_dm is None:
                    pm_full_dm = dm_stack
                else:
                    pm_full_dm = np.concatenate((pm_full_dm, dm_stack), axis=1)

                if self.verbose:
                    print(f"  Current pm_full_dm shape: {pm_full_dm.shape}")

        # ==================== PROCESS LAYERS ====================
        if n_layer > 0:
            if self.verbose:
                print("\n=== Processing Layer Projection Matrices ===")

            for jj, layer in enumerate(layer_list):
                layer_index = layer['index']
                layer_name = layer['name']

                if self.verbose:
                    print(f"\nLoading PMs for {layer_name} (index {layer_index})...")

                # Load PMs for all optical sources for this layer
                layer_pms_list = []

                for ii, opt_source in enumerate(opt_sources):
                    opt_index = opt_source['index']

                    pm_filename = generate_pm_filename(
                        self.params_file, opt_index=opt_index, layer_index=layer_index
                    )

                    if pm_filename is None:
                        if self.verbose:
                            print(f"  Warning: Could not generate filename"
                                  f" for opt{opt_index}, layer{layer_index}")
                        # Create zero matrix as placeholder
                        if layer_pms_list:
                            layer_pms_list.append(np.zeros_like(layer_pms_list[0]))
                        continue

                    pm_path = os.path.join(output_dir, pm_filename)

                    if not os.path.exists(pm_path):
                        if self.verbose:
                            print(f"  Warning: File {pm_path} does not exist, using zeros")
                        if layer_pms_list:
                            layer_pms_list.append(np.zeros_like(layer_pms_list[0]))
                        continue

                    if self.verbose:
                        print(f"  Loading opt{opt_index}: {pm_filename}")

                    intmat_obj = Intmat.restore(pm_path)
                    layer_pms_list.append(intmat_obj.intmat)

                # Stack all optical sources for this layer
                layer_stack = np.stack(layer_pms_list, axis=0)

                if self.verbose:
                    print(f"  Stacked PM shape for {layer_name}: {layer_stack.shape}")

                # Concatenate along layer modes (axis 1)
                if pm_full_layer is None:
                    pm_full_layer = layer_stack
                else:
                    pm_full_layer = np.concatenate((pm_full_layer, layer_stack), axis=1)

                if self.verbose:
                    print(f"  Current pm_full_layer shape: {pm_full_layer.shape}")

        # Display summary information
        if self.verbose:
            print("\n=== Final 3D Projection Matrices ===")
            if pm_full_dm is not None:
                print(f"DM projection matrix shape: {pm_full_dm.shape} "
                    f"(n_opt_sources, n_dm_modes, n_pupil_modes)")
            if pm_full_layer is not None:
                print(f"Layer projection matrix shape: {pm_full_layer.shape} "
                    f"(n_opt_sources, n_layer_modes, n_pupil_modes)")

        # Save the matrices if requested
        if save:
            if pm_full_dm is not None:
                dm_output_filename = "pm_full_dm.npy"
                np.save(os.path.join(output_dir, dm_output_filename), cpuArray(pm_full_dm))
                if self.verbose:
                    print(f"Saved DM projection matrix to {dm_output_filename}")

            if pm_full_layer is not None:
                layer_output_filename = "pm_full_layer.npy"
                np.save(os.path.join(output_dir, layer_output_filename), cpuArray(pm_full_layer))
                if self.verbose:
                    print(f"Saved Layer projection matrix to {layer_output_filename}")

        # ==================== COMPUTE OPTIMAL PROJECTION ====================
        if self.verbose:
            print("\n=== Computing Optimal Projection Matrix ===")

        # Weighted combination
        tpdm_pdm = np.zeros((pm_full_dm.shape[1], pm_full_dm.shape[1]), dtype=cpu_float_dtype)
        tpdm_pl = np.zeros((pm_full_dm.shape[1], pm_full_layer.shape[1]), dtype=cpu_float_dtype)

        total_weight = np.sum(weights_array)

        for i in range(n_opt):
            pdm_i = pm_full_dm[i, :, :]      # shape: (n_dm_modes, n_pupil_modes)
            pl_i = pm_full_layer[i, :, :]    # shape: (n_layer_modes, n_pupil_modes)
            w = weights_array[i] / total_weight

            tpdm_pdm += pdm_i @ pdm_i.T * w
            tpdm_pl += pdm_i @ pl_i.T * w

            if self.verbose:
                print(f"  Processed opt{opt_sources[i]['index']} with weight {w:.3f}")

        # Pseudoinverse with regularization
        if self.verbose:
            print(f"\nApplying regularization (reg_factor={reg_factor})")
            print(f"Shape of tpdm_pdm: {tpdm_pdm.shape}")
            print(f"Shape of tpdm_pl: {tpdm_pl.shape}")

        eps = 1e-14
        tpdm_pdm_inv = np.linalg.pinv(
            tpdm_pdm + reg_factor * np.eye(tpdm_pdm.shape[0]),
            rcond=eps
        )
        p_opt = tpdm_pdm_inv @ tpdm_pl

        if self.verbose:
            print(f"Final projection matrix shape: {p_opt.shape} "
                f"(n_dm_modes, n_layer_modes)")
            print("=== Computation Complete ===\n")

        return p_opt, pm_full_dm, pm_full_layer


    def save_assembled_interaction_matrix(self, wfs_type='lgs', component_type='dm',
                                        output_dir=None, overwrite=False,
                                        apply_filter=True, verbose=None):
        """
        Assemble and save the full interaction matrix for a specific WFS type and component type.
        
        This is useful for:
        - Debugging reconstructor computation
        - Reusing assembled IMs without recomputing
        - Comparing filtered vs unfiltered IMs
        
        Args:
            wfs_type (str): Type of WFS ('lgs', 'ngs', 'ref')
            component_type (str): Type of component ('dm' or 'layer')
            output_dir (str, optional): Directory to save the assembled IM
            overwrite (bool): Whether to overwrite existing file
            apply_filter (bool): Whether to apply filtmat_tag filtering
            verbose (bool, optional): Override the class's verbose setting
            
        Returns:
            str: Path to saved file
        """
        verbose_flag = self.verbose if verbose is None else verbose

        if output_dir is None:
            output_dir = self.im_dir

        os.makedirs(output_dir, exist_ok=True)

        # Generate filename
        config_name = (os.path.basename(self.params_file).split('.')[0]
                    if isinstance(self.params_file, str) else "config")
        filter_suffix = "_filtered" if apply_filter else ""
        output_filename = f"im_full_{config_name}_{wfs_type}_to_{component_type}" \
                          f"{filter_suffix}.fits"
        output_path = os.path.join(output_dir, output_filename)

        # Check if file exists
        if os.path.exists(output_path) and not overwrite:
            if verbose_flag:
                print(f"✓ Assembled IM already exists: {output_filename}")
            return output_path

        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Assembling and Saving Interaction Matrix")
            print(f"{'='*60}")
            print(f"  WFS type: {wfs_type}")
            print(f"  Component type: {component_type}")
            print(f"  Apply filter: {apply_filter}")
            print(f"  Output: {output_filename}")
            print(f"{'='*60}\n")

        # Assemble the IM
        im_full, n_slopes_per_wfs, mode_indices, component_indices = \
            self.assemble_interaction_matrices(
                wfs_type=wfs_type,
                output_im_dir=self.im_dir,
                component_type=component_type,
                save=False,  # Don't save as .npy, we'll save as FITS
                apply_filter=apply_filter
            )

        if verbose_flag:
            print(f"\n  Assembled IM shape: {im_full.shape}")
            print(f"  Components: {component_indices}")
            print(f"  Modes per component: {[len(mi) for mi in mode_indices]}")

        # Save as Intmat (SPECULA format)
        pupdata_tag = f"{config_name}_{wfs_type}_{component_type}"

        intmat_obj = Intmat(
            intmat=im_full,
            pupdata_tag=pupdata_tag,
            norm_factor=1.0,
            target_device_idx=self.target_device_idx \
                if hasattr(self, 'target_device_idx') else None,
            precision=self.precision if hasattr(self, 'precision') else None
        )

        intmat_obj.save(output_path, overwrite=True)

        if verbose_flag:
            print(f"\n  ✓ Saved to: {output_filename}")
            print(f"{'='*60}\n")

        return output_path


    def compute_tomographic_reconstructor(self, r0, L0,
                                        wfs_type='lgs', component_type='layer',
                                        noise_variance=None,
                                        C_noise=None, output_dir=None,
                                        save=False, overwrite=False,
                                        skip_gpu_covariance=False,
                                        verbose=None):
        """
        Compute full tomographic reconstructor from interaction matrices and covariances.
        
        This method integrates:
        1. Interaction matrix assembly (computed on-the-fly, not saved)
        2. Covariance matrix computation/loading (cached to disk)
        3. MMSE reconstructor calculation
        
        Args:
            r0 (float): Fried parameter in meters
            L0 (float): Outer scale in meters
            wfs_type (str): Type of WFS ('lgs', 'ngs', 'ref')
            component_type (str): Type of component ('dm' or 'layer')
            noise_variance (float or array, optional): Noise variance per WFS
            C_noise (np.ndarray, optional): Full noise covariance matrix
            output_dir (str, optional): Directory for saving results
            save (bool): Whether to save the reconstructor
            overwrite (bool): Whether to overwrite existing files
            skip_gpu_covariance (bool): Whether to skip GPU acceleration for covariance computation
            verbose (bool, optional): Override the class's verbose setting
            
        Returns:
            dict: Dictionary with:
                - 'reconstructor': MMSE reconstructor matrix
                - 'im_full': Full interaction matrix
                - 'C_atm_full': Full atmospheric covariance matrix
                - 'C_noise': Noise covariance matrix
                - 'mode_indices': Mode indices per component
                - 'component_indices': Component indices
                - 'n_slopes_per_wfs': Number of slopes per WFS
                - 'rec_filename': Filename if saved
        """
        verbose_flag = self.verbose if verbose is None else verbose

        if verbose_flag:
            print(f"\n{'='*70}")
            print(f"COMPUTING TOMOGRAPHIC RECONSTRUCTOR")
            print(f"{'='*70}")
            print(f"  WFS type: {wfs_type}")
            print(f"  Component type: {component_type}")
            print(f"  r0: {r0} m, L0: {L0} m")
            print(f"{'='*70}\n")

        # ==================== STEP 1: Compute IMs on-the-fly ====================
        if verbose_flag:
            print(f"STEP 1: Computing Interaction Matrices (on-the-fly)")
            print(f"-" * 70)

        # Compute IMs without saving permanently
        _ = self.compute_interaction_matrices(
            output_im_dir=self.im_dir,
            output_rec_dir=self.rec_dir,
            wfs_type=wfs_type,
            overwrite=overwrite,
            verbose=verbose_flag,
            display=False
        )

        # Assemble full IM
        im_full, n_slopes_per_wfs, mode_indices, component_indices = \
            self.assemble_interaction_matrices(
                wfs_type=wfs_type,
                output_im_dir=self.im_dir,
                component_type=component_type,
                save=False  # Don't save assembled IM
            )

        # Temp dir is automatically cleaned up here

        if verbose_flag:
            print(f"  ✓ IM assembled: {im_full.shape}")
            print(f"  Components: {component_indices}")
            print(f"  Modes: {[len(mi) for mi in mode_indices]}")
            print()

        # ==================== STEP 2: Compute/Load Covariances ====================
        if verbose_flag:
            print(f"STEP 2: Computing/Loading Covariance Matrices")
            print(f"-" * 70)

        cov_result = self.compute_covariance_matrices(
            r0=r0,
            L0=L0,
            component_type=component_type,
            output_dir=self.cov_dir,
            overwrite=overwrite,
            full_modes=True,
            skip_gpu_covariance=skip_gpu_covariance,
            verbose=verbose_flag
        )

        # Assemble covariance for selected modes
        C_atm_full = self.assemble_covariance_matrix(
            C_atm_blocks=cov_result['C_atm_blocks'],
            component_indices=cov_result['component_indices'],
            mode_indices=mode_indices,
            verbose=verbose_flag,
            return_inverse=True
        )

        if verbose_flag:
            print(f"  ✓ Covariance assembled: {C_atm_full.shape}")
            print()

        # ==================== STEP 3: Build Noise Covariance ====================
        if verbose_flag:
            print(f"STEP 3: Building Noise Covariance Matrix")
            print(f"-" * 70)

        if C_noise is None:
            # Count WFSs
            out = self.count_mcao_stars()
            if wfs_type == 'lgs':
                n_wfs = out['n_lgs']
            elif wfs_type == 'ngs':
                n_wfs = out['n_ngs']
            else:
                n_wfs = out['n_ref']

            # Compute noise variance if not provided
            if noise_variance is None:
                # ***  Use sigma2_in_nm2 from params ***
                params = self.params
                wfs_params = params[f'sh_{wfs_type}1']
                sa_side_in_m = (params['main']['pixel_pupil'] *
                            params['main']['pixel_pitch'] /
                            wfs_params['subap_on_diameter'])
                subap_npx = wfs_params.get('subap_npx', wfs_params.get('sensor_npx'))
                if subap_npx is None:
                    raise KeyError(f"Neither 'subap_npx' nor 'sensor_npx' found"
                                f" in params['sh_{wfs_type}1']")
                sensor_fov = (wfs_params['sensor_pxscale'] * subap_npx)

                rad2arcsec = 3600. * 180. / np.pi

                # *** Use extracted sigma2_in_nm2 ***
                # convert it to cpu_float_dtype for consistency
                if np.isscalar(self.sigma2_in_nm2):
                    sigma2_in_nm2 = cpu_float_dtype(self.sigma2_in_nm2)
                else:
                    sigma2_in_nm2 = np.array(self.sigma2_in_nm2, dtype=cpu_float_dtype)

                sigma2_in_arcsec2 = sigma2_in_nm2 / (1./rad2arcsec * sa_side_in_m / 4. * 1e9)**2.
                sigma2_in_slope = sigma2_in_arcsec2 * 1./(sensor_fov/2.)**2.

                # *** This is the BASE noise variance (for fully illuminated subaperture) ***
                base_noise_variance = sigma2_in_slope

                if verbose_flag:
                    print(f"  Base noise variance (fully illuminated): {base_noise_variance:.2e}")
                    print(f"  (from sigma2_in_nm2 = {sigma2_in_nm2:.2e})")

                # *** Compute illumination for each WFS ***
                n_slopes_total = im_full.shape[0]
                C_noise = np.zeros((n_slopes_total, n_slopes_total), dtype=cpu_float_dtype)

                n_slopes_list = []
                illumination_list = []

                for i in range(n_wfs):
                    wfs_params = self.get_wfs_params(wfs_type, i+1,
                                                     xp_local=np)

                    if wfs_params['idx_valid_sa'] is not None:
                        n_slopes_this_wfs = len(wfs_params['idx_valid_sa']) * 2
                    else:
                        n_slopes_this_wfs = 0

                    n_slopes_list.append(n_slopes_this_wfs)

                    # *** Compute subaperture illumination for this WFS ***
                    if n_slopes_this_wfs > 0:
                        illumination = synim.compute_subaperture_illumination(
                            pup_mask=self.pup_mask,
                            wfs_nsubaps=wfs_params['wfs_nsubaps'],
                            wfs_rotation=wfs_params['wfs_rotation'],
                            wfs_translation=wfs_params['wfs_translation'],
                            wfs_magnification=wfs_params['wfs_magnification'],
                            idx_valid_sa=wfs_params['idx_valid_sa'],
                            verbose=verbose_flag
                        )

                        # *** Noise variance is INVERSELY proportional to illumination ***
                        # (less light = more noise)
                        # Each subaperture has noise variance = base_variance / illumination
                        # We repeat for X and Y slopes
                        sa_noise_variance = base_noise_variance / \
                            (illumination + 1e-10)  # Avoid division by zero
                        wfs_noise_variance = np.repeat(sa_noise_variance, 2).astype(cpu_float_dtype)  # X and Y slopes

                        illumination_list.append(illumination)
                    else:
                        wfs_noise_variance = np.array([], dtype=cpu_float_dtype)
                        illumination_list.append(np.array([], dtype=cpu_float_dtype))

                    if verbose_flag:
                        if len(wfs_noise_variance) > 0:
                            print(f"  WFS {i+1} noise variance:")
                            print(f"    Min: {np.min(wfs_noise_variance):.2e}")
                            print(f"    Max: {np.max(wfs_noise_variance):.2e}")
                            print(f"    Mean: {np.mean(wfs_noise_variance):.2e}")

                # Build diagonal noise covariance with illumination weighting
                for i in range(n_wfs):
                    start_idx = sum(n_slopes_list[:i])
                    end_idx = sum(n_slopes_list[:i+1])

                    if n_slopes_list[i] > 0:
                        # Use illumination-weighted noise variance
                        sa_noise_variance = base_noise_variance / (illumination_list[i] + 1e-10)
                        wfs_noise_variance = np.repeat(sa_noise_variance, 2)

                        # C_noise diagonal elements are 1/variance (precision)
                        C_noise[start_idx:end_idx, start_idx:end_idx] = np.diag(
                            1.0 / wfs_noise_variance
                        )

                if verbose_flag:
                    print(f"  ✓ Noise covariance built with illumination weighting:"
                          f"{C_noise.shape}")
            else:
                # User provided noise variance
                n_slopes_total = im_full.shape[0]
                if np.isscalar(noise_variance):
                    C_noise = np.eye(n_slopes_total, dtype=cpu_float_dtype) * (1.0 / noise_variance)
                else:
                    if len(noise_variance) != n_slopes_total:
                        raise ValueError("Length of noise_variance does not match"
                                         " number of slopes in IM")
                    C_noise = np.diag(1.0 / noise_variance)

                if verbose_flag:
                    print(f"  ✓ Noise covariance built from provided variance:"
                          f"{C_noise.shape}")
            C_noise_from_input = False
        else:
            if verbose_flag:
                print(f"  Using provided C_noise: {C_noise.shape}")
            C_noise = C_noise.astype(cpu_float_dtype)
            C_noise_from_input = True

        print()

        # ==================== STEP 4: Compute MMSE Reconstructor ====================
        if verbose_flag:
            print(f"STEP 4: Computing MMSE Reconstructor")
            print(f"-" * 70)

        # Check: all must be cpu_float_dtype (regardless of endianness)
        expected_dtype = cpu_float_dtype
        for name, arr in [('im_full', im_full), ('C_atm_full', C_atm_full), ('C_noise', C_noise)]:
            arr_dtype = np.dtype(arr.dtype)
            # Check kind and itemsize (float32: kind='f', itemsize=4)
            if not (arr_dtype.kind == expected_dtype.kind and
                    arr_dtype.itemsize == expected_dtype.itemsize):
                msg = (f"Data type mismatch for {name}: got {arr_dtype}, "
                    f"expected {expected_dtype} (any endianness)")
                if verbose_flag:
                    print("WARNING:", msg)
                else:
                    raise TypeError(msg)

        reconstructor = compute_mmse_reconstructor(
            im_full,  # Transpose for SPECULA convention
            C_atm_full,
            noise_variance=None,  # Already in C_noise
            C_noise=C_noise,
            cinverse=True,
            verbose=verbose_flag
        )

        if verbose_flag:
            print(f"  ✓ Reconstructor computed: {reconstructor.shape}")
            print()

        # ==================== STEP 5: Save if requested ====================
        rec_filename = None
        rec_path = None

        if save:
            if output_dir is None:
                raise ValueError("output_dir must be specified when save=True")

            os.makedirs(output_dir, exist_ok=True)

            if verbose_flag:
                print(f"STEP 5: Saving Reconstructor")
                print(f"-" * 70)

            # Generate filename
            config_name = (os.path.basename(self.params_file).split('.')[0]
                        if isinstance(self.params_file, str) else "config")

            rec_filename = (f"rec_{config_name}_{wfs_type}_{component_type}_"
                        f"r0{r0:.3f}_L0{L0:.1f}")
            if C_noise_from_input:
                rec_filename += f"_Cnoise"
            else:
                if noise_variance is None:
                    rec_filename += f"_var{self.sigma2_in_nm2:.3f}nm2"
                else:
                    rec_filename += f"_var{noise_variance:.3f}uSl2"
            rec_filename += ".fits"
            rec_path = os.path.join(output_dir, rec_filename)

            # Save as Recmat (SPECULA format)
            recmat_obj = Recmat(
                recmat=reconstructor,
                norm_factor=1.0,
                target_device_idx=self.target_device_idx \
                    if hasattr(self, 'target_device_idx') else None,
                precision=self.precision if hasattr(self, 'precision') else None
            )
            recmat_obj.save(rec_path, overwrite=True)

            if verbose_flag:
                print(f"  ✓ Saved: {rec_filename}")
                print()

        # ==================== Summary ====================
        if verbose_flag:
            print(f"{'='*70}")
            print(f"TOMOGRAPHIC RECONSTRUCTOR COMPLETE")
            print(f"{'='*70}")
            print(f"  Reconstructor shape: {reconstructor.shape}")
            print(f"  (n_modes, n_slopes) = ({reconstructor.shape[0]}, {reconstructor.shape[1]})")
            if save:
                print(f"  Saved to: {rec_filename}")
            print(f"{'='*70}\n")

        return {
            'reconstructor': reconstructor,
            'im_full': im_full,
            'C_atm_full': C_atm_full,
            'C_noise': C_noise,
            'mode_indices': mode_indices,
            'component_indices': component_indices,
            'n_slopes_per_wfs': n_slopes_per_wfs,
            'n_wfs': n_wfs if 'n_wfs' in locals() else None,
            'rec_filename': rec_filename,
            'rec_path': rec_path,
            'r0': r0,
            'L0': L0
        }


    def compute_tomographic_projection_matrix(self, output_dir=None, save=False,
                                            verbose=None):
        """
        Compute tomographic projection matrix following IDL compute_mcao_popt logic.
        
        Uses parameters extracted during initialization (from YAML 'projection' 
        or IDL 'modalrec1' sections).
        
        Implements the standard MCAO projection:
            p_opt = (P_DM^T @ P_DM + reg_factor*I)^(-1) @ P_DM^T @ P_Layer
        
        where P_DM and P_Layer are weighted combinations of projection matrices
        from multiple optical sources.
        
        Args:
            output_dir (str, optional): Directory where PM files are stored
            save (bool): Whether to save the projection matrix to disk
            verbose (bool, optional): Override the class's verbose setting
        
        Returns:
            tuple: (p_opt, pm_full_dm, pm_full_layer, info)
                - p_opt: Tomographic projection matrix (n_dm_modes, n_layer_modes)
                - pm_full_dm: Full DM projection matrix (n_opt, n_dm_modes, n_pupil_modes)
                - pm_full_layer: Full Layer projection matrix (n_opt, n_layer_modes, n_pupil_modes)
                - info: dict with computation metadata
        """

        verbose_flag = self.verbose if verbose is None else verbose

        # Use cached reg_factor from initialization
        reg_factor = self.proj_reg_factor

        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Computing Tomographic Projection Matrix")
            print(f"{'='*60}")
            if self.projection_params is not None:
                print(f"  Using YAML 'projection' section")
            else:
                print(f"  Using IDL-style or default parameters")
            print(f"  Regularization factor: {reg_factor}")
            print(f"{'='*60}\n")

        # ==================== LOAD FULL PM MATRICES ====================
        # Reuse existing function to load/compute all projection matrices
        _, pm_full_dm, pm_full_layer = self.compute_projection_matrix(
            reg_factor=reg_factor,
            output_dir=output_dir,
            save=False  # Don't save intermediate matrices
        )

        if verbose_flag:
            print(f"\n  Loaded PM matrices:")
            print(f"    pm_full_dm shape: {pm_full_dm.shape}")
            print(f"      (n_opt_sources, n_dm_modes, n_pupil_modes)")
            print(f"    pm_full_layer shape: {pm_full_layer.shape}")
            print(f"      (n_opt_sources, n_layer_modes, n_pupil_modes)")

        # ==================== GET OPTICAL SOURCE WEIGHTS ====================            
        if self.projection_params is not None:
            # Use from projection section
            opt_sources = self.projection_params['opt_sources']
            weights_array = np.array([src['weight'] for src in opt_sources])
        else:
            # Fallback to extract_opt_list (which now handles both formats)
            opt_sources_list = extract_opt_list(self.params)
            opt_sources = [{'index': src['index'], **src['config']} for src in opt_sources_list]
            weights_array = np.array([src.get('weight', 1.0) for src in opt_sources])

        n_opt = len(opt_sources)
        total_weight = np.sum(weights_array)

        if n_opt == 0:
            raise ValueError(
                "No optical sources (source_optX) found in configuration. "
                "Cannot compute tomographic projection matrix."
            )

        if verbose_flag:
            print(f"\n  Optical sources: {n_opt}")
            print(f"  Weights: {weights_array}")
            print(f"  Total weight: {total_weight}")

        # ==================== COMPUTE WEIGHTED MATRICES ====================
        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Computing Weighted Combination")
            print(f"{'='*60}")

        # Get dimensions
        n_dm_modes = pm_full_dm.shape[1]
        n_layer_modes = pm_full_layer.shape[1]
        n_pupil_modes = pm_full_dm.shape[2]

        # Initialize accumulation matrices
        tpdm_pdm = np.zeros((n_dm_modes, n_dm_modes))
        tpdm_pl = np.zeros((n_dm_modes, n_layer_modes))

        # Accumulate weighted contributions from each optical source
        for i in range(n_opt):
            pdm_i = pm_full_dm[i, :, :]      # (n_dm_modes, n_pupil_modes)
            pl_i = pm_full_layer[i, :, :]    # (n_layer_modes, n_pupil_modes)
            w = weights_array[i] / total_weight

            # Accumulate: P_DM^T @ P_DM and P_DM^T @ P_Layer
            tpdm_pdm += pdm_i @ pdm_i.T * w
            tpdm_pl += pdm_i @ pl_i.T * w

            if verbose_flag:
                print(f"  Processed opt{i+1} (weight={w:.3f})")

        if verbose_flag:
            print(f"\n  tpdm_pdm shape: {tpdm_pdm.shape} (P_DM^T @ P_DM)")
            print(f"  tpdm_pl shape: {tpdm_pl.shape} (P_DM^T @ P_Layer)")

        # ==================== TIKHONOV REGULARIZATION ====================
        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Applying Tikhonov Regularization")
            print(f"{'='*60}")
            print(f"  Adding reg_factor * I to P_DM^T @ P_DM")

        tpdm_pdm_reg = tpdm_pdm + reg_factor * np.eye(n_dm_modes)

        # ==================== PSEUDOINVERSE ====================
        if verbose_flag:
            print(f"  Computing pseudoinverse...")

        # Condition number check
        cond_number_check = False
        if cond_number_check:
            cond_number = np.linalg.cond(tpdm_pdm_reg)
            print(f"  Condition number: {cond_number:.2e}")

        # Compute pseudoinverse
        rcond = 1e-14  # Same as IDL default
        tpdm_pdm_inv = np.linalg.pinv(tpdm_pdm_reg, rcond=rcond)

        if verbose_flag:
            print(f"  ✓ Pseudoinverse computed (rcond={rcond})")

        # ==================== FINAL PROJECTION MATRIX ====================
        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Computing Final Projection Matrix")
            print(f"{'='*60}")
            print(f"  p_opt = (P_DM^T @ P_DM + λI)^(-1) @ P_DM^T @ P_Layer")

        # p_opt = (tpdm_pdm + reg_factor*I)^(-1) @ tpdm_pl
        p_opt = tpdm_pdm_inv @ tpdm_pl

        if verbose_flag:
            print(f"\n  ✓ Tomographic projection matrix computed: {p_opt.shape}")
            print(f"    (n_dm_modes, n_layer_modes) = ({n_dm_modes}, {n_layer_modes})")

        # ==================== SAVE MATRICES ====================
        rec_filename = None
        rec_path = None

        if save and output_dir is not None:
            os.makedirs(output_dir, exist_ok=True)

            if verbose_flag:
                print(f"\n{'='*60}")
                print(f"Saving Results")
                print(f"{'='*60}")

            # Create Recmat object for p_opt
            recmat_obj = Recmat(
                recmat=p_opt,
                norm_factor=1.0,
                target_device_idx=self.target_device_idx \
                    if hasattr(self, 'target_device_idx') else None,
                precision=self.precision if hasattr(self, 'precision') else None
            )

            # Generate filename
            config_name = (os.path.basename(self.params_file).split('.')[0]
                        if isinstance(self.params_file, str) else "config")

            rec_filename = f"rec_{config_name}_tomographic_r{reg_factor:.0e}.fits"
            rec_path = os.path.join(output_dir, rec_filename)

            # Save as FITS (SPECULA format)
            recmat_obj.save(rec_path, overwrite=True)

            if verbose_flag:
                print(f"  ✓ Saved tomographic reconstruction matrix (SPECULA format):")
                print(f"    {rec_filename}")
                print(f"    Shape: {p_opt.shape} (n_dm_modes, n_layer_modes)")

            # Also save intermediate matrices as numpy arrays for debugging
            np.save(os.path.join(output_dir, "tpdm_pdm.npy"), cpuArray(tpdm_pdm))
            np.save(os.path.join(output_dir, "tpdm_pl.npy"), cpuArray(tpdm_pl))
            np.save(os.path.join(output_dir, "tpdm_pdm_reg.npy"), cpuArray(tpdm_pdm_reg))

            if verbose_flag:
                print(f"\n  ✓ Also saved debug matrices (NumPy format):")
                print(f"    - tpdm_pdm.npy (P_DM^T @ P_DM)")
                print(f"    - tpdm_pl.npy (P_DM^T @ P_Layer)")
                print(f"    - tpdm_pdm_reg.npy (with regularization)")

        # ==================== METADATA ====================
        info = {
            'n_opt_sources': n_opt,
            'n_dm_modes': n_dm_modes,
            'n_layer_modes': n_layer_modes,
            'n_pupil_modes': n_pupil_modes,
            'weights': weights_array,
            'reg_factor': reg_factor,
            'rcond': rcond,
            'rec_filename': rec_filename,
            'rec_path': rec_path
        }

        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Tomographic Projection Computation Complete")
            print(f"{'='*60}\n")

        return p_opt, pm_full_dm, pm_full_layer, info


    def list_wfs(self):
        """Return a list of all WFS names and types."""
        return [(wfs['name'], wfs['type'], wfs['index']) for wfs in self.wfs_list]


    def list_dm(self):
        """Return a list of all DM names and indices."""
        return [(dm['name'], dm['index']) for dm in self.dm_list]


    def get_source_info(self, wfs_name):
        """Return source info for a given WFS."""
        return extract_source_info(self.params, wfs_name)


    def generate_im_filename(self, wfs_type=None, wfs_index=None,
                            dm_index=None, layer_index=None, timestamp=False, verbose=False):
        """Generate the interaction matrix filename for a given WFS-DM/Layer combination."""
        return generate_im_filename(
            self.params_file,
            wfs_type=wfs_type,
            wfs_index=wfs_index,
            dm_index=dm_index,
            layer_index=layer_index,  # Aggiungi esplicitamente questo parametro
            timestamp=timestamp,
            verbose=verbose
        )

    def compute_covariance_matrices(self, r0, L0, component_type='layer',
                                    output_dir=None, overwrite=False,
                                    full_modes=True, skip_gpu_covariance=False,
                                    verbose=None):
        """
        Compute atmospheric covariance matrices for all components (DMs or layers).
        
        Implements smart caching similar to IDL's ifs_covmat:
        - Generates unique filename based on component parameters, r0, L0
        - Loads from disk if file exists (unless overwrite=True)
        - Computes and saves if file doesn't exist
        - Saves in FITS format (like IDL) for compatibility
        
        Args:
            r0 (float): Fried parameter in meters
            L0 (float): Outer scale in meters
            component_type (str): Type of component ('dm' or 'layer')
            output_dir (str, optional): Directory to save covariance matrices
            overwrite (bool): Whether to overwrite existing files
            full_modes (bool): If True, compute covariance for ALL modes available.
                            If False, use only modes from modal_combination.
            skip_gpu_covariance (bool): If True, skip GPU-accelerated covariance computation
            verbose (bool, optional): Override the class's verbose setting
            
        Returns:
            dict: Dictionary with:
                - 'C_atm_blocks': List of covariance matrices for each component
                - 'component_indices': List of component indices
                - 'n_modes_per_component': List of number of modes for each component
                - 'start_modes': List of start mode indices for each component
                - 'r0': Fried parameter used
                - 'L0': Outer scale used
                - 'files': List of file paths (saved or loaded)
        """

        verbose_flag = self.verbose if verbose is None else verbose

        # Validate component_type
        if component_type not in ['dm', 'layer']:
            raise ValueError("component_type must be either 'dm' or 'layer'")

        # Get wavelength for conversion nm -> rad^2
        wavelengthInNm = 500.0  # Default
        if 'pyramid' in self.params:
            wavelengthInNm = self.params['pyramid'].get('wavelengthInNm', 500.0)
        elif 'sh' in self.params:
            wavelengthInNm = self.params['sh'].get('wavelengthInNm', 500.0)
        elif 'sh1' in self.params:
            wavelengthInNm = self.params['sh1'].get('wavelengthInNm', 500.0)

        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Computing Covariance Matrices")
            print(f"{'='*60}")
            print(f"  r0 = {r0} m")
            print(f"  L0 = {L0} m")
            print(f"  Component type: {component_type}")
            print(f"  Full modes: {full_modes}")
            print(f"  Wavelength: {wavelengthInNm} nm")
            print(f"{'='*60}\n")

        # Get component list
        if component_type == 'dm':
            component_list = self.dm_list
        else:
            component_list = extract_layer_list(self.params)

        # Set up output directory (use 'ifunc' subdirectory like IDL)
        if output_dir is None:
            output_dir = os.path.join(self.cm.root_dir, 'ifunc')
        os.makedirs(output_dir, exist_ok=True)

        component_indices = []
        C_atm_blocks = []
        n_modes_per_component = []
        start_modes = []
        cov_files = []

        # ==================== COMPUTE/LOAD EACH COMPONENT ====================
        for idx, comp in enumerate(component_list):
            comp_idx = int(comp['index'])
            component_indices.append(comp_idx)

            if verbose_flag:
                print(f"\n[{idx+1}/{len(component_list)}] Processing {component_type}{comp_idx}...")

            # Get component parameters
            comp_params = self.get_component_params(
                comp_idx,
                is_layer=(component_type == 'layer'),
                cut_start_mode=False,  # Load ALL modes
                xp_local=np
            )

            # Get start_mode
            comp_key = f'{component_type}{comp_idx}'
            if comp_key not in self.params and component_type == 'dm' and comp_idx == 1 and 'dm' in self.params:
                comp_key = 'dm'
            if comp_key in self.params and 'start_mode' in self.params[comp_key]:
                start_mode = self.params[comp_key]['start_mode']
            else:
                start_mode = 0
            start_modes.append(start_mode)

            # Total modes available
            total_modes = comp_params['dm_array'].shape[2]
            modes_to_use = list(range(total_modes))
            n_modes = len(modes_to_use)

            if verbose_flag:
                print(f"  Total modes available: {total_modes}")
                print(f"  Start mode: {start_mode}")
                print(f"  Computing for: {n_modes} modes")

            # ========== GENERATE FILENAME (EXACTLY LIKE IDL) ==========
            cov_filename, base_tag = generate_cov_filename(
                self.params[comp_key], self.pup_diam_m, r0, L0
            )
            cov_path = os.path.join(output_dir, cov_filename)
            cov_files.append(cov_path)

            if verbose_flag:
                print(f"  Filename: {cov_filename}")

            # ========== CHECK IF FILE EXISTS ==========
            if os.path.exists(cov_path) and not overwrite:
                if verbose_flag:
                    print(f"  ✓ Loading existing: {cov_filename}")

                # Load from FITS
                with fits.open(cov_path) as hdul:
                    C_atm = hdul[0].data.copy()

                if verbose_flag:
                    print(f"    Shape: {C_atm.shape}")

                C_atm_blocks.append(C_atm)
                continue

            # ========== COMPUTE COVARIANCE MATRIX ==========
            if verbose_flag:
                if os.path.exists(cov_path):
                    print(f"  ⚠ File exists but overwrite=True, recomputing...")
                else:
                    print(f"  File not found, computing...")

            # Convert 3D DM array to 2D
            dm2d = dm3d_to_2d(cpuArray(comp_params['dm_array']), cpuArray(comp_params['dm_mask']),
                              xp_local=np, float_dtype_local=cpu_float_dtype)

            # Select modes
            dm2d_selected = dm2d[modes_to_use, :]

            if verbose_flag:
                print(f"  dm2d_selected shape: {dm2d_selected.shape}")
                print(f"  Computing covariance matrix...")

            # Compute covariance matrix
            if skip_gpu_covariance:
                C_atm_rad2 = compute_ifs_covmat(
                    comp_params['dm_mask'],
                    self.pup_diam_m,
                    dm2d_selected,
                    r0,
                    L0,
                    xp=np,
                    dtype=cpu_float_dtype,
                    oversampling=2,
                    verbose=False
                )
            else:
                self._print_gpu_memory(label="Before covariance computation: ",
                                       verbose=verbose_gpu_flag)

                C_atm_rad2 = compute_ifs_covmat(
                    to_xp(xp, comp_params['dm_mask'], dtype=float_dtype),
                    self.pup_diam_m,
                    to_xp(xp, dm2d_selected , dtype=float_dtype),
                    r0,
                    L0,
                    xp=xp,
                    dtype=float_dtype,
                    oversampling=2,
                    verbose=False
                )

                C_atm_rad2 = cpuArray(C_atm_rad2)

            if verbose_flag:
                print(f"  ✓ Covariance computed: {C_atm_rad2.shape}")
                print(f"    RMS (nm):"
                      f" {np.sqrt(np.diag(C_atm_rad2*(500**2/2/np.pi**2))).mean():.2f}")
                print(f"    RMS (rad): {np.sqrt(np.diag(C_atm_rad2)).mean():.4f}")

            # ========== SAVE TO FITS (LIKE IDL) ==========
            # IDL: writefits, fileNameCov, turb_covmat
            hdu = fits.PrimaryHDU(cpuArray(C_atm_rad2))
            hdu.header['R0'] = (r0, 'Fried parameter [m]')
            hdu.header['L0'] = (L0, 'Outer scale [m]')
            hdu.header['UNITS'] = ('rad^2', 'Covariance units')
            hdu.header['DIAMM'] = (self.pup_diam_m, 'Pupil diameter [m]')
            hdu.header['WAVELNM'] = (wavelengthInNm, 'Wavelength [nm]')
            hdu.header['NMODES'] = (n_modes, 'Number of modes')
            hdu.header['STARTMOD'] = (0, 'Covariance includes ALL modes from 0')
            hdu.header['TOTMODES'] = (total_modes, 'Total modes in covariance')
            hdu.header['REFSTART'] = (start_mode, 'Reference start_mode (not applied)')
            hdu.header['COMPTAG'] = (base_tag, 'Component tag')
            hdu.header['COMPTYPE'] = (component_type, 'Component type (dm/layer)')
            hdu.header['COMPIDX'] = (comp_idx, 'Component index')

            hdu.writeto(cov_path, overwrite=True)

            if verbose_flag:
                print(f"  ✓ Saved to FITS: {cov_filename}")

            C_atm_blocks.append(C_atm_rad2)

        if verbose_flag:
            print(f"\n{'='*60}")
            print(f"Covariance Computation Complete")
            print(f"{'='*60}")
            print(f"  Components processed: {len(C_atm_blocks)}")
            print(f"  Files in: {output_dir}")
            print(f"{'='*60}\n")

        return {
            'C_atm_blocks': C_atm_blocks,
            'component_indices': component_indices,
            'n_modes_per_component': n_modes_per_component,
            'start_modes': start_modes,
            'r0': r0,
            'L0': L0,
            'wavelength_nm': wavelengthInNm,
            'files': cov_files
        }


    def assemble_covariance_matrix(self, C_atm_blocks, component_indices,
                                    mode_indices=None,
                                    wfs_type=None, component_type='layer',
                                    verbose=None, return_inverse=False):
        """
        Assemble the full covariance matrix from individual blocks,
        extracting only the modes specified in mode_indices.
        
        Args:
            C_atm_blocks (list): List of covariance matrices for each component
            component_indices (list): List of component indices
            mode_indices (list, optional): List of mode index arrays for each component.
                                        If None, uses modal_combination.
            wfs_type (str, optional): WFS type for modal_combination lookup
            component_type (str): Type of component ('dm' or 'layer')
            verbose (bool, optional): Override the class's verbose setting
            return_inverse (bool, optional): Whether to return the inverse of the covariance matrix
            
        Returns:
            np.ndarray: Full covariance matrix with selected modes
        """
        verbose_flag = self.verbose if verbose is None else verbose

        # If mode_indices not provided, get from modal_combination
        if mode_indices is None:
            if wfs_type is None:
                raise ValueError("Either mode_indices or wfs_type must be provided")

            modal_key = f'modes_{wfs_type}_{component_type}'
            if ('modal_combination' not in self.params or
                modal_key not in self.params['modal_combination']):
                raise ValueError(f"modal_combination.{modal_key} not found in params")

            modes_config = self.params['modal_combination'][modal_key]
            mode_indices = []

            for i, (comp_idx, n_modes_cfg) in enumerate(zip(component_indices, modes_config)):
                if n_modes_cfg > 0:
                    # Get start_mode
                    comp_key = f'{component_type}{comp_idx}'
                    if comp_key not in self.params and component_type == 'dm' and comp_idx == 1 and 'dm' in self.params:
                        comp_key = 'dm'
                    if comp_key in self.params and 'start_mode' in self.params[comp_key]:
                        start_mode = self.params[comp_key]['start_mode']
                    else:
                        start_mode = 0

                    n_total = C_atm_blocks[i].shape[0]
                    mode_indices.append(list(range(start_mode, n_total)))
                else:
                    mode_indices.append([])

            if verbose_flag:
                print(f"Using modal_combination: {modal_key}")
                print(f"  Mode counts: {modes_config}")

        # ==================== COMPUTE WEIGHTS FROM ATMOSPHERIC PROFILE ====================
        weights = compute_layer_weights_from_turbulence(
            self.params,
            component_indices,
            component_type=component_type,
            verbose=verbose_flag
        )

        # Calculate total modes
        total_modes = sum(len(mi) for mi in mode_indices)

        if verbose_flag:
            print(f"\nAssembling covariance matrix:")
            print(f"  Components: {component_indices}")
            print(f"  Modes per component: {[len(mi) for mi in mode_indices]}")
            print(f"  Total modes: {total_modes}")
            print(f"  Weights: {weights}")

        # Initialize full covariance matrix
        C_atm_full = np.zeros((total_modes, total_modes), dtype=cpu_float_dtype)

        # Conversion factor (nm to rad^2 at 500nm)
        conversion_factor = (500 / 2 / np.pi) ** 2

        # Fill the blocks
        current_idx = 0
        for i, (C_atm_block, modes, weight) in enumerate(zip(C_atm_blocks, mode_indices, weights)):
            if len(modes) == 0:
                continue

            # Extract the sub-block for selected modes
            # C_atm_block has shape (n_total_modes, n_total_modes)
            # We want to extract modes[i] × modes[j]
            valid_modes = [m for m in modes if m < C_atm_block.shape[0]]
            if len(valid_modes) == 0:
                if verbose_flag:
                    print(f"  Warning: No valid modes for component {i+1} (requested:"
                          f" {modes}, available: {C_atm_block.shape[0]}) -- skipping.")
                continue

            idx_modes = np.ix_(valid_modes, valid_modes)
            try:
                C_atm_sub = C_atm_block[idx_modes]
            except IndexError:
                print(f"Error extracting modes for component {i+1}:")
                print(f"  Available modes: {C_atm_block.shape[0]}")
                print(f"  Requested modes: {valid_modes}")
                raise

            # Place in full matrix
            idx_full = slice(current_idx, current_idx + len(valid_modes))
            if return_inverse:
                # Check if C_atm_sub is diagonal
                if np.allclose(C_atm_sub, np.diag(np.diagonal(C_atm_sub))):
                    # Invert only the diagonal elements
                    diag = np.diagonal(C_atm_sub * weight * conversion_factor)
                    C_atm_full[idx_full, idx_full] = np.diag(1.0 / diag)
                else:
                    C_atm_full[idx_full, idx_full] = \
                        np.linalg.pinv(C_atm_sub * weight * conversion_factor)
            else:
                C_atm_full[idx_full, idx_full] = C_atm_sub * weight * conversion_factor

            if verbose_flag:
                print(f"  Component {i+1}: modes {valid_modes[0]}-{valid_modes[-1]} → "
                    f"full matrix [{current_idx}:{current_idx + len(valid_modes)}]")

            current_idx += len(valid_modes)

        if verbose_flag:
            print(f"\n  ✓ Full covariance matrix assembled: {C_atm_full.shape}")

        return C_atm_full
