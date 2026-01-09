"""Utility functions for handling simulation parameters in SynIM."""

import os
import re
import yaml
import datetime
import numpy as np

from synim import cpuArray, to_xp, cpu_float_dtype

# Import all utility functions from utils
from synim.utils import *

import specula
specula.init(device_idx=-1, precision=1)

from specula.calib_manager import CalibManager
from specula.data_objects.ifunc import IFunc
from specula.data_objects.m2c import M2C
from specula.data_objects.simul_params import SimulParams
from specula.data_objects.pupilstop import Pupilstop
from specula.data_objects.subap_data import SubapData
from specula.lib.compute_zern_ifunc import compute_zern_ifunc


def is_simple_config(config):
    """
    Detect if this is a simple SCAO config or a complex MCAO config.
    
    Args:
        config (dict): Configuration dictionary
        
    Returns:
        bool: True for simple SCAO config, False for complex MCAO config
    """
    # Check for multiple DMs
    dm_count = sum(1 for key in config if key.startswith('dm') and key != 'dm')

    # Check for multiple WFSs
    wfs_count = sum(1 for key in config if
                   (key.startswith('sh_') or key.startswith('pyramid')) and key != 'pyramid')

    return dm_count == 0 and wfs_count == 0

def wfs_fov_from_config(wfs_params):
    """
    Extract field of view value from WFS parameters.
    
    Args:
        wfs_params (dict): Dictionary with WFS parameters
        
    Returns:
        float: Field of view in arcseconds
    """
    if wfs_params.get('sensor_fov') is not None:
        wfs_fov_arcsec = wfs_params['sensor_fov']
    elif wfs_params.get('fov') is not None:
        wfs_fov_arcsec = wfs_params['fov']
    elif wfs_params.get('subap_wanted_fov') is not None:
        wfs_fov_arcsec = wfs_params['subap_wanted_fov']
    else:
        wfs_fov_arcsec = 0
    return wfs_fov_arcsec

def determine_source_type(wfs_name):
    """
    Determine the source type from a WFS name.
    
    Args:
        wfs_name (str): Name of the WFS
        
    Returns:
        str: Source type ('lgs', 'ngs', or 'ref')
    """
    if 'lgs' in wfs_name:
        return 'lgs'
    elif 'ngs' in wfs_name:
        return 'ngs'
    elif 'ref' in wfs_name:
        return 'ref'
    return 'ngs'  # default

def extract_source_coordinates(config, wfs_key):
    """
    Extract polar coordinates for a given source.
    
    Args:
        config (dict): Configuration dictionary
        wfs_key (str): Key of the WFS in the config
        
    Returns:
        list: [distance, angle] polar coordinates
    """
    # First check if coordinates are in WFS parameters
    if wfs_key in config and 'gs_pol_coo' in config[wfs_key]:
        return config[wfs_key]['gs_pol_coo']

    # Try to find source corresponding to this WFS
    source_match = re.search(r'((?:lgs|ngs|ref)\d+)', wfs_key)
    if source_match:
        source_key = f'source_{source_match.group(1)}'
        if source_key in config:
            if 'polar_coordinates' in config[source_key]:
                return config[source_key]['polar_coordinates']
            elif 'polar_coordinate' in config[source_key]:
                return config[source_key]['polar_coordinate']

    # Try on_axis_source for simple configs
    if 'on_axis_source' in config:
        if 'polar_coordinates' in config['on_axis_source']:
            return config['on_axis_source']['polar_coordinates']
        elif 'polar_coordinate' in config['on_axis_source']:
            return config['on_axis_source']['polar_coordinate']

    # Default to on-axis
    return [0.0, 0.0]

def extract_source_height(config, wfs_key):
    """
    Extracts the actual height of the source from the configuration.
    If zenithAngleInDeg is present, returns height * airmass.
    If height is not present, returns np.inf.
    """
    # Compute airmass
    if 'main' in config:
        zenith_angle = config['main'].get('zenithAngleInDeg', None)
        zenith_rad = np.deg2rad(zenith_angle)
        airmass = 1.0 / np.cos(zenith_rad)
    else:
        airmass = 1.0

    # First check if coordinates are in WFS parameters
    if wfs_key in config and 'height' in config[wfs_key]:
        return config[wfs_key]['height'] * airmass

    # Try to find source corresponding to this WFS
    source_match = re.search(r'((?:lgs|ngs|ref)\d+)', wfs_key)
    if source_match:
        source_key = f'source_{source_match.group(1)}'
        if source_key in config:
            if 'height' in config[source_key]:
                return config[source_key]['height'] * airmass

    # Try on_axis_source for simple configs
    if 'on_axis_source' in config:
        if 'height' in config['on_axis_source']:
            return config['on_axis_source']['height'] * airmass

    return np.inf  # Default to infinity if no height is specified

def get_tag_or_object(params, base_key, look_for_tag=False):
    """
    Looks for and returns the value associated with '{base_key}_tag' or
    '{base_key}_object' or '{base_key}_ref' in params.
    Returns None if neither is present.
    """
    if f"{base_key}_tag" in params:
        return params[f"{base_key}_tag"]
    elif f"{base_key}_object" in params:
        return params[f"{base_key}_object"]
    elif f"{base_key}_ref" in params:
        ref_name = params[f"{base_key}_ref"]
        # Look for referenced section in params
        if ref_name in params:
            ref_section = params[ref_name]
            if 'tag' in ref_section:
                return ref_section['tag']
            else:
                raise ValueError(f"Referenced section '{ref_name}' does not contain a 'tag'")
        else:
            raise ValueError(f"Referenced section '{ref_name}' not found in config")
    elif f"tag" in params and look_for_tag:
        return params["tag"]
    return None

def load_pupilstop(cm, pupilstop_params, pixel_pupil, pixel_pitch, verbose=False):
    """
    Load or create a pupilstop.
    
    NOTE: Returns numpy array from disk - caller should convert with to_xp if needed
    """

    pupilstop_tag = get_tag_or_object(pupilstop_params, 'pupil_mask', look_for_tag=True)
    if pupilstop_tag is None:
        pupilstop_tag = get_tag_or_object(pupilstop_params, 'pupilstop')
    if pupilstop_tag is not None:
        if verbose:
            print(f"     Loading pupilstop from file, tag: {pupilstop_tag}")
        pupilstop_path = cm.filename('pupilstop', pupilstop_tag)
        pupilstop = Pupilstop.restore(pupilstop_path)
        # *** Returns numpy from FITS file ***
        return pupilstop.A
    else:
        # Create pupilstop from parameters
        mask_diam = pupilstop_params.get('mask_diam', 1.0)
        obs_diam = pupilstop_params.get('obs_diam', 0.0)

        # Create a new Pupilstop instance with the given parameters
        simul_params = SimulParams(pixel_pupil=pixel_pupil, pixel_pitch=pixel_pitch)
        pupilstop = Pupilstop(
            simul_params,
            mask_diam=mask_diam,
            obs_diam=obs_diam,
            target_device_idx=-1,
            precision=0
        )
        # *** Returns numpy ***
        return pupilstop.A

def load_influence_functions(cm, dm_params, pixel_pupil, verbose=False, is_inverse_basis=False):
    """
    Load or generate DM influence functions.
    
    Args:
        cm (CalibManager): SPECULA calibration manager
        dm_params (dict): DM parameters
        pixel_pupil (int): Number of pixels across pupil
        verbose (bool): Whether to print details
        is_inverse_basis (bool): If True, don't convert to 3D (for projection matrices)
        
    Returns:
        tuple: (dm_array, dm_mask) - For DM: 3D array, For inverse basis: 2D array
    """

    ifunc_tag = get_tag_or_object(dm_params, 'ifunc')

    if ifunc_tag is not None:
        if verbose:
            print(f"     Loading influence function from file, tag: {ifunc_tag}")
        ifunc_path = cm.filename('ifunc', ifunc_tag)
        ifunc = IFunc.restore(ifunc_path)

        m2c_tag = get_tag_or_object(dm_params, 'm2c')
        if m2c_tag is not None:
            if verbose:
                print(f"     Loading M2C from file, tag: {m2c_tag}")
            m2c_path = cm.filename('m2c', m2c_tag)
            m2c = M2C.restore(m2c_path)
            # multiply the influence function by the M2C
            ifunc.influence_function = m2c.m2c.T @ ifunc.influence_function

        if 'nmodes' in dm_params:
            if ifunc.influence_function.shape[0] > dm_params['nmodes']:
                ifunc.influence_function = ifunc.influence_function[:dm_params['nmodes'],:]

        # Check dimensions to determine if this is an inverse basis
        n_rows = ifunc.influence_function.shape[0]
        n_cols = ifunc.influence_function.shape[1]

        if ifunc.mask_inf_func is not None:
            n_valid_pixels = np.sum(ifunc.mask_inf_func > 0.5)
        else:
            raise ValueError("IFunc without mask_inf_func is not supported."
                           " Mask is required.")

        # If n_rows >> n_cols, this is likely an inverse basis (pixels x modes)
        # If n_cols >> n_rows, this is a normal basis (modes x pixels)
        is_inverse = (n_rows > n_cols * 2)  # Heuristic: if rows >> cols

        if verbose:
            print(f"     Influence function shape: {ifunc.influence_function.shape}")
            print(f"     Valid pixels in mask: {n_valid_pixels}")
            print(f"     Detected as: {'INVERSE basis' if is_inverse else 'NORMAL basis'}")

        # For inverse basis, ALWAYS return 2D
        if is_inverse or is_inverse_basis:
            if verbose:
                print(f"     Returning 2D inverse basis (optimized format)")

            # Make sure it's (n_modes, n_pixels) format
            if ifunc.influence_function.shape[0] < ifunc.influence_function.shape[1]:
                # Already correct: n_modes x n_pixels
                return ifunc.influence_function, ifunc.mask_inf_func
            else:
                # Need to transpose: n_pixels x n_modes -> n_modes x n_pixels
                return ifunc.influence_function.T, ifunc.mask_inf_func

        # For normal basis, convert to 3D
        else:
            # Convert influence function from 2D to 3D
            dm_array = dm2d_to_3d(ifunc.influence_function, ifunc.mask_inf_func,
                                  xp_local=np, float_dtype_local=cpu_float_dtype)
            if verbose:
                print(f"     DM array shape: {dm_array.shape}")
            dm_mask = ifunc.mask_inf_func.copy()
            if verbose:
                print(f"     DM mask shape: {dm_mask.shape}")
                print(f"     DM mask sum: {np.sum(dm_mask)}")
            return dm_array, dm_mask

    elif 'type_str' in dm_params:
        # ... existing code for Zernike generation ...
        if verbose:
            print(f"     Loading influence function from type_str: {dm_params['type_str']}")

        nmodes = dm_params.get('nmodes', 100)
        obsratio = dm_params.get('obsratio', 0.0)
        npixels = dm_params.get('npixels', pixel_pupil)

        mask_tag = get_tag_or_object(dm_params, 'mask')
        if mask_tag is not None:
            mask_path = cm.filename('pupilstop', mask_tag)
            if verbose:
                print(f"     Loading mask from file, tag: {mask_tag}")
            pupilstop = Pupilstop.restore(mask_path)
            mask = pupilstop.A
        else:
            mask = None
            if verbose:
                print("     No mask provided. Using default mask.")

        z_ifunc, z_mask = compute_zern_ifunc(npixels, nmodes, xp=np, dtype=float,
                                             obsratio=obsratio, diaratio=1.0,
                                             start_mode=0, mask=mask)

        # For Zernike, always return 3D
        dm_array = dm2d_to_3d(z_ifunc, z_mask,
                              xp_local=np, float_dtype_local=cpu_float_dtype)
        if verbose:
            print(f"     DM array shape: {dm_array.shape}")
        dm_mask = z_mask
        if verbose:
            print(f"     DM mask shape: {dm_mask.shape}")
            print(f"     DM mask sum: {np.sum(dm_mask)}")
        return dm_array, dm_mask
    else:
        raise ValueError("No valid influence function configuration found."
                         " Need either 'ifunc_tag', 'ifunc_object', or 'type_str'.")


def find_subapdata(cm, wfs_params, wfs_key, params, verbose=False):
    """
    Find and load SubapData for valid subapertures.
    
    NOTE: Returns numpy array from disk - caller should convert with to_xp if needed
    
    Returns:
        numpy.ndarray: Array of valid subaperture indices (numpy) or None
    """

    subap_path = None
    subap_tag = None

    # First check - Try to get subapdata from WFS params
    subap_tag = get_tag_or_object(wfs_params, 'subapdata')
    if subap_tag is not None:
        if verbose:
            print("     Loading subapdata from file, tag:", subap_tag)
        subap_path = cm.filename('subapdata', subap_tag)

    # Second check - Try to find corresponding slopec section based on WFS name
    elif wfs_key is not None:
        # Determine potential slopec key based on WFS key (e.g., sh_lgs1 -> slopec_lgs1)
        slopec_key = None
        if wfs_key.startswith('sh_'):
            potential_slopec = 'slopec_' + wfs_key[3:]
            if potential_slopec in params:
                slopec_key = potential_slopec
        elif wfs_key.startswith('pyramid_'):
            potential_slopec = 'slopec_' + wfs_key[8:]
            if potential_slopec in params:
                slopec_key = potential_slopec
        # Handle numeric indices (e.g., sh1 -> slopec1)
        elif any(char.isdigit() for char in wfs_key):
            # Extract numeric portion
            numeric_part = ''.join(char for char in wfs_key if char.isdigit())
            if numeric_part:
                potential_slopec = f'slopec{numeric_part}'
                if potential_slopec in params:
                    slopec_key = potential_slopec

        # Check standard slopec key
        if slopec_key is None and 'slopec' in params:
            slopec_key = 'slopec'

        if slopec_key:
            slopec_params = params[slopec_key]
            subap_tag = get_tag_or_object(slopec_params, 'subapdata')
            if subap_tag is not None:
                if verbose:
                    print(f"     Loading subapdata from {slopec_key}, tag:", subap_tag)
                subap_path = cm.filename('subapdata', subap_tag)

    # Third check - Try generic slopec section
    elif 'slopec' in params:
        slopec_params = params['slopec']
        subap_tag = get_tag_or_object(slopec_params, 'subapdata')
        if subap_tag is not None:
            if verbose:
                print("     Loading subapdata from slopec, tag:", subap_tag)
            subap_path = cm.filename('subapdata', subap_tag)

    if subap_path is None:
        if verbose:
            print("     No subapdata file found. Using default.")
        return None
    else:
        if verbose:
            print("     Subapdata file found:", subap_path)

    # Try to load the subapdata if a path was found
    if subap_path and os.path.exists(subap_path):
        if verbose:
            print("     Loading subapdata from file:", subap_path)
        subap_data = SubapData.restore(subap_path)
        # *** Returns numpy from FITS file ***
        return np.transpose(np.asarray(np.where(subap_data.single_mask())))

    return None

def insert_interaction_matrix_part(im_full, intmat_obj, mode_idx,
                                   slope_idx_start, slope_idx_end, verbose=False):
    """
    Insert part of an interaction matrix into a combined matrix.
    
    Args:
        im_full (numpy.ndarray): Target combined interaction matrix
        intmat_obj (Intmat): Source interaction matrix object
        mode_idx (list): Indices of modes to extract
        slope_idx_start (int): Start index for slopes in target matrix
        slope_idx_end (int): End index for slopes in target matrix
        verbose (bool): Whether to print details
        
    Returns:
        bool: True if insertion was successful
    """
    # Make sure we don't exceed matrix dimensions
    if not mode_idx or slope_idx_end > im_full.shape[1]:
        if verbose:
            print(f"  Warning: Invalid indices for matrix insertion")
        return False

    # Calculate how many modes we can actually use from this IM
    available_dm_modes = intmat_obj.intmat.shape[0]
    actual_mode_indices = [idx for idx in mode_idx if idx < available_dm_modes]

    if not actual_mode_indices:
        if verbose:
            print(f"  Warning: No valid mode indices. "
                  f"Available modes: {available_dm_modes}, requested: {mode_idx}")
        return False

    # Insert the valid modes into our combined matrix
    n_slopes = slope_idx_end - slope_idx_start
    im_full[mode_idx, slope_idx_start:slope_idx_end] = \
        intmat_obj.intmat[actual_mode_indices, :n_slopes]

    if verbose:
        print(f"  Inserted {len(actual_mode_indices)} modes at indices {actual_mode_indices}, "
              f"slopes {slope_idx_start}:{slope_idx_end}")

    return True

def build_source_filename_part(source_config, zenith_angle=None):
    """
    Build the source-specific part of the filename.
    
    Args:
        source_config (dict): Source configuration parameters
        zenith_angle (float, optional): Zenith angle override
        
    Returns:
        str: Source filename component
    """
    parts = []

    # Polar coordinates
    if 'polar_coordinate' in source_config or 'polar_coordinates' in source_config:
        pol_coords = source_config.get('polar_coordinate', source_config.get('polar_coordinates'))
        if isinstance(pol_coords, (list, tuple)) and len(pol_coords) >= 2:
            # Format: pdXXaYY where XX is angle in arcsec, YY is azimuth in degrees
            angle_arcsec = pol_coords[0]
            azimuth_deg = pol_coords[1]
            parts.append(f"pd{angle_arcsec:.1f}a{azimuth_deg:.0f}")

    # Height
    height = source_config.get('height', float('inf'))
    if zenith_angle is not None:
        zenith_rad = np.deg2rad(zenith_angle)
        airmass = 1.0 / np.cos(zenith_rad)
        height *= airmass
    if height != float('inf'):
        parts.append(f"h{height:.0f}")

    # Join all parts with underscore
    return "_".join(parts) if parts else "on_axis"


def build_pupil_filename_part(pupil_params):
    """
    Build filename part for pupil parameters.
    
    Args:
        pupil_params (dict): Pupil configuration
        
    Returns:
        list: Filename parts for pupil
    """
    parts = []

    if pupil_params:
        ps = pupil_params.get('pixel_pupil', 0)
        pp = pupil_params.get('pixel_pitch', 0)
        parts.append(f"ps{ps}p{pp:.4f}")

        if 'obsratio' in pupil_params and pupil_params['obsratio'] > 0:
            parts.append(f"o{pupil_params['obsratio']:.3f}")

    return parts

def build_wfs_filename_part(wfs_config, wfs_type=None):
    """
    Build the WFS-specific part of the filename.
    
    Args:
        wfs_config (dict): WFS configuration parameters
        wfs_type (str, optional): WFS type for context
        
    Returns:
        str: WFS filename component
    """
    parts = []

    # WFS type
    wfs_class = wfs_config.get('class', 'sh')

    # Number of subapertures or pupil diameter
    if 'subap_on_diameter' in wfs_config:
        n_subaps = wfs_config['subap_on_diameter']
        parts.append(f"{wfs_class}{n_subaps}x{n_subaps}")
    elif 'pup_diam' in wfs_config:
        pup_diam = wfs_config['pup_diam']
        parts.append(f"{wfs_class}{pup_diam}")

    # Wavelength
    if 'wavelength' in wfs_config:
        wl = wfs_config['wavelength']
        parts.append(f"wl{wl:.0f}")

    # Field of view
    if 'sensor_pxscale' in wfs_config and 'subap_npx' in wfs_config:
        fov = wfs_config['sensor_pxscale'] * wfs_config['subap_npx']
        parts.append(f"fv{fov:.1f}")

    # Number of pixels
    if 'subap_npx' in wfs_config:
        npx = wfs_config['subap_npx']
        parts.append(f"np{npx}")

    # Shifts
    if 'xShiftPhInPixel' in wfs_config or 'yShiftPhInPixel' in wfs_config:
        x_shift = wfs_config.get('xShiftPhInPixel', 0)
        y_shift = wfs_config.get('yShiftPhInPixel', 0)
        parts.append(f"shiftX{x_shift}Y{y_shift}")

    # Rotation
    if 'rotAnglePhInDeg' in wfs_config:
        parts.append(f"rot{wfs_config['rotAnglePhInDeg']}")

    # Magnification
    if 'magnification' in wfs_config:
        parts.append(f"mag{wfs_config['magnification']:.3f}")

    # Anamorphosis
    if 'anamorph90' in wfs_config:
        parts.append(f"an90d{wfs_config['anamorph90']:.3f}")
    if 'anamorph45' in wfs_config:
        parts.append(f"an45d{wfs_config['anamorph45']:.3f}")

    return "_".join(parts) if parts else "wfs"


def build_component_filename_part(component_config, component_type='dm', 
                                  include_height=True, n_modes_override=None):
    """
    Build the DM/Layer-specific part of the filename.
    
    Args:
        component_config (dict): DM or Layer configuration parameters
        component_type (str): Type of component ('dm' or 'layer')
        include_height (bool): Whether to include height in the output
        n_modes_override (int, optional): Number of modes ACTUALLY USED.
                                         Only adds suffix if < total nmodes.
        
    Returns:
        str: DM or Layer filename component
    """
    parts = []

    # Only include height if requested
    if include_height:
        prefix = 'layH' if component_type == 'layer' else 'dmH'
        parts.append(f"{prefix}{component_config.get('height', 0.0):.1f}")

    ifunc_tag = get_tag_or_object(component_config, 'ifunc')
    if ifunc_tag is not None:
        parts.append(f"ifunc_{ifunc_tag}")

    m2c_tag = get_tag_or_object(component_config, 'm2c')
    if m2c_tag is not None:
        # *** Check if m2c_tag already contains "modes" ***
        if 'modes' in m2c_tag.lower():
            # Tag already has modes info (e.g., "base_162px_41acts_1260modes")

            if n_modes_override is not None:
                # Need to add override info
                # Extract total modes from tag if possible
                import re
                match = re.search(r'(\d+)modes', m2c_tag, re.IGNORECASE)
                if match:
                    nmodes_in_tag = int(match.group(1))

                    # Only add suffix if using fewer modes
                    if n_modes_override < nmodes_in_tag:
                        parts.append(f"m2c_{m2c_tag}_{n_modes_override}modes")
                    else:
                        # Using all modes from tag
                        parts.append(f"m2c_{m2c_tag}")
                else:
                    # Can't parse, just add the tag with override
                    parts.append(f"m2c_{m2c_tag}_{n_modes_override}modes")
            else:
                # No override, just use the tag as-is
                parts.append(f"m2c_{m2c_tag}")
        else:
            # Tag doesn't contain "modes", add it from config
            nmodes_total = component_config.get('nmodes', None)

            if n_modes_override is not None and nmodes_total is not None:
                # Only add suffix if using fewer modes
                if n_modes_override < nmodes_total:
                    parts.append(f"m2c_{m2c_tag}_{n_modes_override}modes")
                else:
                    # Using all modes
                    parts.append(f"m2c_{m2c_tag}_{nmodes_total}modes")
            elif nmodes_total is not None:
                # No override, use total
                parts.append(f"m2c_{m2c_tag}_{nmodes_total}modes")
            else:
                # No modes info
                parts.append(f"m2c_{m2c_tag}")

    return "_".join(parts) if parts else component_type


def parse_pro_file(pro_file_path):
    """
    Parse a .pro file and extract its structure into a Python dictionary.
    Improved version to handle IDL-style syntax better.

    Args:
        pro_file_path (str): Path to the .pro file.

    Returns:
        dict: Parsed data as a dictionary.
    """
    data = {}
    current_section = None

    with open(pro_file_path, 'r') as file:
        for line_num, line in enumerate(file, 1):
            # Remove comments (everything after ;)
            comment_pos = line.find(';')
            if comment_pos != -1:
                line = line[:comment_pos]

            line = line.strip()
            if not line:
                continue

            try:
                # Recognize the start of a new section (e.g., {main, {dm1, etc.)
                section_match = re.match(r'^\{(\w+),?', line)
                if section_match:
                    current_section = section_match.group(1).lower()
                    data[current_section] = {}
                    continue

                # Recognize the end of a section
                if line == '}':
                    current_section = None
                    continue

                # If we're in a section, process key-value pairs
                if current_section:
                    # Match both : and = assignments
                    key_value_match = re.match(r'(\w+)\s*[:=]\s*(.+)', line)
                    if key_value_match:
                        key = key_value_match.group(1).strip()
                        value = key_value_match.group(2).strip()

                        # Remove trailing comma
                        if value.endswith(','):
                            value = value[:-1].strip()

                        # Parse the value
                        parsed_value = _parse_pro_value(value)
                        data[current_section][key] = parsed_value

            except Exception as e:
                print(f"Warning: Error parsing line {line_num}: '{line}' - {e}")
                continue

    return data

def _parse_pro_value(value):
    """
    Parse a single value from a PRO file, handling IDL-specific syntax.
    """
    value = value.strip()

    # Handle special IDL values
    if value == '!VALUES.F_INFINITY':
        return float('inf')

    # Handle boolean values (IDL style)
    if value.lower() in ['0b', 'false']:
        return False
    elif value.lower() in ['1b', 'true']:
        return True

    # Handle quoted strings
    if (value.startswith("'") and value.endswith("'")) or \
       (value.startswith('"') and value.endswith('"')):
        return value[1:-1]

    # Handle arrays [val1, val2, ...]
    if value.startswith('[') and value.endswith(']'):
        return _parse_pro_array(value)

    # Handle replicate function: replicate(val, n)
    replicate_match = re.match(r'replicate\(([^,]+),\s*(\d+)\)', value, re.IGNORECASE)
    if replicate_match:
        val = _parse_pro_value(replicate_match.group(1))
        num = int(replicate_match.group(2))
        return [val] * num

    # *** MIGLIORAMENTO: Pattern per numeri più flessibile ***
    # Handle integers (with optional L suffix)
    if re.match(r'^-?\d+[lL]?$', value):
        return int(value.rstrip('lL'))

    # Handle floats (inclusi quelli che finiscono con .)
    if re.match(r'^-?\d*\.?\d*([eE][+-]?\d+)?[dD]?$', value) and any(c.isdigit() for c in value):
        try:
            return float(value.rstrip('dD'))
        except ValueError:
            pass

    # Handle scientific notation
    if re.match(r'^-?\d+[eE][+-]?\d+$', value):
        return float(value)

    # Handle mathematical expressions (e.g., 38.5/480)
    if '/' in value and re.match(r'^[\d\.\+\-\*/\(\)\s]+$', value):
        try:
            return eval(value)
        except:
            pass

    # Handle 'auto' and other special string values without quotes
    if value.lower() in ['auto']:
        return value.lower()

    # Default: return as string
    return value

def _parse_pro_array(array_str):
    """
    Parse a PRO array string like [val1, val2, val3].
    
    Args:
        array_str (str): Array string including brackets
        
    Returns:
        list: Parsed array
    """
    # Remove brackets
    content = array_str[1:-1].strip()
    if not content:
        return []

    # Split by comma
    elements = []
    parts = content.split(',')

    for part in parts:
        part = part.strip()
        if part:
            # Handle replicate within arrays
            if 'replicate(' in part.lower():
                replicate_match = re.match(r'replicate\(([^,]+),\s*(\d+)\)', part, re.IGNORECASE)
                if replicate_match:
                    val = _parse_pro_value(replicate_match.group(1))
                    num = int(replicate_match.group(2))
                    elements.extend([val] * num)
                    continue

            # Parse individual element
            elements.append(_parse_pro_value(part))

    return elements

def parse_params_file(file_path):
    """
    Parse a parameters file (YAML or .pro) and return its contents as a dictionary.

    Args:
        file_path (str): Path to the parameters file.

    Returns:
        dict: Parsed data as a dictionary.
    """
    if file_path.endswith('.yml') or file_path.endswith('.yaml'):
        with open(file_path, 'r') as file:
            return yaml.safe_load(file)
    elif file_path.endswith('.pro'):
        return parse_pro_file(file_path)
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

def prepare_interaction_matrix_params(params, wfs_type=None,
                                      wfs_index=None, dm_index=None):
    """
    Prepares parameters for synim.interaction_matrix from a SPECULA YAML configuration file.
    
    Args:
        params (dict): Dictionary with configuration parameters.
        wfs_type (str, optional): Type of WFS ('sh', 'pyr') or source type ('lgs', 'ngs', 'ref')
        wfs_index (int, optional): Index of the WFS to use (1-based)
        dm_index (int, optional): Index of the DM to use (1-based)
        
    Returns:
        dict: Parameters ready to be passed to synim.interaction_matrix
    """

    # Prepare the CalibManager
    main_params = params['main']
    cm = CalibManager(main_params['root_dir'])

    # Extract general parameters
    pixel_pupil = main_params['pixel_pupil']
    pixel_pitch = main_params['pixel_pitch']
    pup_diam_m = pixel_pupil * pixel_pitch

    # Determine if this is a simple or complex configuration
    simple_config = is_simple_config(params)

    # Extract all WFS and DM configurations
    wfs_list = extract_wfs_list(params)
    dm_list = extract_dm_list(params)

    # Load pupilstop and create pupil mask
    pup_mask = None
    if 'pupilstop' in params:
        pupilstop_params = params['pupilstop']
        pup_mask = load_pupilstop(cm, pupilstop_params, pixel_pupil, pixel_pitch, verbose=True)

    # If no pupilstop defined, create a default circular pupil
    if pup_mask is None:
        simul_params = SimulParams(pixel_pupil=pixel_pupil, pixel_pitch=pixel_pitch)
        pupilstop = Pupilstop(
            simul_params,
            mask_diam=1.0,
            obs_diam=0.0,
            target_device_idx=-1,
            precision=0
        )
        pup_mask = pupilstop.A

    # Find the appropriate DM based on dm_index
    selected_dm = None

    if dm_index is not None:
        # Try to find DM with specified index
        for dm in dm_list:
            if dm['index'] == str(dm_index):
                selected_dm = dm
                print(f"DM -- Using specified DM: {dm['name']}")
                break
    elif dm_index == 1 and "dm" in params:
        # Fallback: use 'dm' if 'dm1' not found and index==1
        dm_info = build_component_filename_part(
            params["dm"],
            'dm',
            n_modes_override=n_modes
        )
        filename_parts.append(dm_info)

    # If no DM found with specified index or no index specified, use first DM
    if selected_dm is None and dm_list:
        selected_dm = dm_list[0]
        print(f"DM -- Using first available DM: {selected_dm['name']}")

    if selected_dm is None:
        raise ValueError("No DM configuration found in the YAML file.")

    dm_key = selected_dm['name']
    dm_params = selected_dm['config']

    # Extract DM parameters
    dm_height = dm_params.get('height', 0)
    dm_rotation = dm_params.get('rotation', 0.0)

    # Load influence functions
    dm_array, dm_mask = load_influence_functions(cm, dm_params, pixel_pupil, verbose=True)

    if 'nmodes' in dm_params:
        nmodes = dm_params['nmodes']
        if dm_array.shape[2] > nmodes:
            print(f"     Trimming DM array to first {nmodes} modes")
            dm_array = dm_array[:, :, :nmodes]

    # WFS selection logic
    selected_wfs = None
    source_type = None

    print("WFS -- Looking for WFS parameters...")

    # Simple SCAO configuration case
    if simple_config:
        print("     Simple SCAO configuration detected")
        if len(wfs_list) > 0:
            selected_wfs = wfs_list[0]
            source_type = 'ngs'  # Default for simple configs
            print(f"     Using WFS: {selected_wfs['name']} of type {selected_wfs['type']}")
        else:
            raise ValueError("No WFS configuration found in the YAML file.")
    else:
        # Complex MCAO configuration
        print("     Complex MCAO configuration detected")

        # Case 1: wfs_type specifies the sensor type ('sh', 'pyr')
        if wfs_type in ['sh', 'pyr']:
            print(f"     Looking for WFS of type: {wfs_type}")
            matching_wfs = [wfs for wfs in wfs_list if wfs['type'] == wfs_type]

            if wfs_index is not None:
                # Try to find specific index
                for wfs in matching_wfs:
                    if wfs['index'] == str(wfs_index):
                        selected_wfs = wfs
                        print(f"     Found WFS with specified index: {wfs['name']}")
                        break

            # If no specific index found, use the first one
            if selected_wfs is None and matching_wfs:
                selected_wfs = matching_wfs[0]
                print(f"     Using first WFS of type {wfs_type}: {selected_wfs['name']}")

        # Case 2: wfs_type specifies the source type ('lgs', 'ngs', 'ref')
        elif wfs_type in ['lgs', 'ngs', 'ref']:
            source_type = wfs_type
            print(f"     Looking for WFS associated with {wfs_type} source")

            # Pattern for WFS names corresponding to the source type
            pattern = f"sh_{source_type}"
            matching_wfs = [wfs for wfs in wfs_list if pattern in wfs['name']]

            if wfs_index is not None:
                # Try to find specific index within the source type
                target_name = f"{pattern}{wfs_index}"
                for wfs in matching_wfs:
                    if wfs['name'] == target_name:
                        selected_wfs = wfs
                        print(f"     Found WFS with specified index: {wfs['name']}")
                        break

            # If no specific index found, use the first one
            if selected_wfs is None and matching_wfs:
                selected_wfs = matching_wfs[0]
                print(f"     Using first WFS for {wfs_type}: {selected_wfs['name']}")

        # Case 3: Only wfs_index is specified (no wfs_type)
        elif wfs_index is not None:
            print(f"     Looking for WFS with index: {wfs_index}")
            for wfs in wfs_list:
                if wfs['index'] == str(wfs_index):
                    selected_wfs = wfs
                    print(f"     Found WFS with specified index: {wfs['name']}")
                    break

        # Case 4: No specific criteria, use the first available WFS
        if selected_wfs is None and wfs_list:
            selected_wfs = wfs_list[0]
            print(f"     Using first available WFS: {selected_wfs['name']}")

    # If no WFS found, raise error
    if selected_wfs is None:
        raise ValueError("No matching WFS configuration found in the YAML file.")

    wfs_key = selected_wfs['name']
    wfs_params = selected_wfs['config']
    wfs_type_detected = selected_wfs['type']

    # Determine source type from WFS name if not already set
    if source_type is None:
        source_type = determine_source_type(wfs_key)

    print(f"     Source type determined: {source_type}")

    # Extract WFS parameters
    wfs_rotation = wfs_params.get('rotation', 0.0)
    wfs_translation = wfs_params.get('translation', [0.0, 0.0])
    wfs_magnification = wfs_params.get('magnification', 1.0)
    wfs_fov_arcsec = wfs_fov_from_config(wfs_params)

    if wfs_type_detected == 'sh':
        # Shack-Hartmann specific parameters
        wfs_nsubaps = wfs_params.get('subap_on_diameter', 0)
    else:
        # Pyramid specific parameters
        wfs_nsubaps = wfs_params.get('pup_diam', 0)

    # Load SubapData for valid subapertures
    idx_valid_sa = find_subapdata(cm, wfs_params, wfs_key, params, verbose=True)

    # Guide star parameters
    if source_type == 'lgs':
        # LGS is at finite height
        # Try to get height from source or use typical LGS height
        gs_height = None

        # Check if there's a specific source for this WFS and try to get height
        source_match = re.search(r'(lgs\d+)', wfs_key)
        if source_match:
            source_key = f'source_{source_match.group(1)}'
            if source_key in params:
                gs_height = params[source_key].get('height', None)

        # If still no height, use default
        if gs_height is None:
            gs_height = 90000.0  # Default LGS height in meters
    else:
        # NGS and REF are at infinite distance
        gs_height = float('inf')

    # Get source polar coordinates
    gs_pol_coo = extract_source_coordinates(params, wfs_key)

    # Return the prepared parameters
    return {
        'pup_diam_m': pup_diam_m,
        'pup_mask': pup_mask,
        'dm_array': dm_array,
        'dm_mask': dm_mask,
        'dm_height': dm_height,
        'dm_rotation': dm_rotation,
        'wfs_key': wfs_key,
        'wfs_type': wfs_type_detected,
        'wfs_nsubaps': wfs_nsubaps,
        'wfs_rotation': wfs_rotation,
        'wfs_translation': wfs_translation,
        'wfs_magnification': wfs_magnification,
        'wfs_fov_arcsec': wfs_fov_arcsec,
        'gs_pol_coo': gs_pol_coo,
        'gs_height': gs_height,
        'idx_valid_sa': idx_valid_sa,
        'dm_key': dm_key,
        'source_type': source_type
    }


def extract_wfs_list(config):
    """Extract all WFS configurations from config"""
    wfs_list = []

    # Find Pyramid WFSs
    for key in config:
        if key == 'pyramid' or (
            isinstance(key, str) and (key.startswith('pyramid') or key == 'pyr')
        ):
            wfs_list.append({
                'name': key,
                'type': 'pyr',
                'index': re.findall(r'\d+', key)[0] if re.findall(r'\d+', key) else '1',
                'config': config[key]
            })

    # Find Shack-Hartmann WFSs
    for key in config:
        if key == 'sh' or key.startswith('sh_') or key.startswith('sh'):
            wfs_list.append({
                'name': key,
                'type': 'sh', 
                'index': re.findall(r'\d+', key)[0] if re.findall(r'\d+', key) else '1',
                'config': config[key]
            })

    return wfs_list

def extract_dm_list(config):
    """Extract all DM configurations from config"""
    dm_list = []

    for key in config:
        if key == 'dm' or (isinstance(key, str) and key.startswith('dm')):
            dm_list.append({
                'name': key,
                'index': re.findall(r'\d+', key)[0] if re.findall(r'\d+', key) else '1',
                'config': config[key]
            })

    # If no DMs found, create a default one
    if len(dm_list) == 0 and 'dm' in config:
        dm_list.append({
            'name': 'dm',
            'index': '1',
            'config': config['dm']
        })

    return dm_list


def extract_layer_list(config):
    """Extract all layer configurations from config"""
    layer_list = []

    for key in config:
        if isinstance(key, str) and key.startswith('layer'):
            layer_list.append({
                'name': key,
                'index': re.findall(r'\d+', key)[0] if re.findall(r'\d+', key) else '1',
                'config': config[key]
            })

    return layer_list

def extract_source_info(config, wfs_name):
    """Extract source information related to a specific WFS"""
    source_info = {}

    # Check both formats: direct reference or connection through prop
    if 'inputs' in config[wfs_name] and 'in_ef' in config[wfs_name]['inputs']:
        source_ref = config[wfs_name]['inputs']['in_ef']

        # Format might be 'prop.out_source_lgs1_ef' or direct 'out_source_lgs1_ef'
        match = re.search(r'out_(\w+)_ef', source_ref)
        if match:
            source_name = match.group(1)

            # Find the source in the config
            if source_name in config:
                source = config[source_name]
                source_info['type'] = 'lgs' if 'lgs' in source_name else 'ngs'
                source_info['name'] = source_name
                if 'polar_coordinate' in source or 'polar_coordinates' in source:
                    source_info['pol_coords'] = source.get(
                        'polar_coordinate', source.get('polar_coordinates')
                    )
                if 'height' in source:
                    source_info['height'] = source['height']
                if 'wavelengthInNm' in source:
                    source_info['wavelength'] = source['wavelengthInNm']

    # Direct reference within source objects (for simple YAML)
    if not source_info and 'on_axis_source' in config:
        source_info['type'] = 'ngs'
        source_info['name'] = 'on_axis_source'
        if config['on_axis_source'].source.get('polar_coordinate', source.get('polar_coordinates')):
            source_info['pol_coords'] = config['on_axis_source'].get(
                'polar_coordinate', config['on_axis_source'].get('polar_coordinates')
            )
        else:
            source_info['pol_coords'] = [0, 0]
        source_info['wavelength'] = config['on_axis_source'].get('wavelengthInNm', 750)

    return source_info


def extract_opt_list(params):
    """
    Extract list of optical sources (source_optX) from parameters.
    
    Args:
        params (dict): Configuration dictionary
        
    Returns:
        list: List of dictionaries with optical source information:
            - 'name': Source name (e.g., 'source_opt1')
            - 'index': Source index (1-based, as integer)
            - 'config': Configuration dictionary for this source
    """
    # Try YAML format first
    proj_params = extract_projection_params(params)
    if proj_params is not None:
        # Convert projection format to opt_list format
        opt_list = []
        for i, src in enumerate(proj_params['opt_sources']):
            opt_list.append({
                'name': f'opt{i+1}',
                'index': i+1,
                'config': src  # Already has 'polar_coordinates' and 'weight'
            })
        return opt_list

    # Fallback: IDL- style source_optX format
    opt_list = []
    pattern = re.compile(r'^source_opt(\d+)$')
    for key in params.keys():
        match = pattern.match(key)
        if match:
            opt_index = int(match.group(1))
            opt_list.append({
                'name': key,
                'index': opt_index,
                'config': params[key]
            })

    # Sort by index
    opt_list.sort(key=lambda x: x['index'])
    return opt_list


def validate_opt_sources(params, verbose=False):
    """
    Validate optical source configurations.
    
    Args:
        params (dict): Configuration dictionary
        verbose (bool): Whether to print validation messages
        
    Returns:
        bool: True if all sources are valid
        
    Raises:
        ValueError: If validation fails
    """
    opt_list = extract_opt_list(params)

    if len(opt_list) == 0:
        if verbose:
            print("Warning: No optical sources (source_optX) found in configuration")
        return True

    for opt in opt_list:
        opt_name = opt['name']
        opt_config = opt['config']

        # Check required fields
        if 'polar_coordinate' in opt_config:
            pc_string = 'polar_coordinate'
        elif 'polar_coordinates' in opt_config:
            pc_string = 'polar_coordinates'
        else:
            raise ValueError(f"{opt_name} missing required field:"
                             f" polar_coordinate or polar_coordinates")

        # Validate polar_coordinate format
        pol_coo = opt_config[pc_string]
        if not isinstance(pol_coo, (list, tuple)) or len(pol_coo) != 2:
            raise ValueError(
                f"{opt_name}.{pc_string} must be [distance, angle], "
                f"got {pol_coo}"
            )

        # Set defaults for optional fields
        if 'height' not in opt_config:
            opt_config['height'] = float('inf')
            if verbose:
                print(f"{opt_name}: height not specified, using infinity (NGS)")

        if 'weight' not in opt_config:
            opt_config['weight'] = 1.0
            if verbose:
                print(f"{opt_name}: weight not specified, using 1.0")

    if verbose:
        print(f"✓ Validated {len(opt_list)} optical sources")
        for opt in opt_list:
            pol_coo = opt['config'][pc_string]
            height = opt['config']['height']
            weight = opt['config']['weight']
            h_str = f"{height:.0f}m" if not np.isinf(height) else "∞ (NGS)"
            print(f"  {opt['name']}: [{pol_coo[0]:.1f}\", {pol_coo[1]:.0f}°] "
                  f"h={h_str}, w={weight:.2f}")

    return True


def generate_im_filename(params_file, wfs_type=None,
                         wfs_index=None, dm_index=None,
                         layer_index=None, n_modes=None,
                         timestamp=False, verbose=False):
    """
    Generate the interaction matrix filename for a given WFS-DM/Layer combination.
    
    Args:
        params_file (str or dict): Path to configuration file or dictionary
        wfs_type (str, optional): Type of WFS source ('lgs', 'ngs', 'ref')
        wfs_index (int, optional): Index of the WFS (1-based)
        dm_index (int, optional): Index of the DM (1-based)
        layer_index (int, optional): Index of the Layer (1-based)
        n_modes (int, optional): Number of modes used (if limited)
        timestamp (bool): Whether to include a timestamp in the filename
        verbose (bool): Whether to print detailed information
        
    Returns:
        str: Generated filename
    """
    # Check that only one of dm_index or layer_index is specified
    if dm_index is not None and layer_index is not None:
        raise ValueError("Cannot specify both dm_index and layer_index")

    if dm_index is None and layer_index is None:
        raise ValueError("Must specify either dm_index or layer_index")

    # Load configuration if needed
    if isinstance(params_file, str):
        params = parse_params_file(params_file)
    else:
        params = params_file

    # Get main parameters
    main_params = params.get('main', {})

    # Find the WFS configuration
    wfs_key = None
    if wfs_type and wfs_index:
        wfs_key = f"sh_{wfs_type}{wfs_index}"
    elif wfs_type:
        wfs_key = f"sh_{wfs_type}1"  # Default to first WFS of this type

    if wfs_key and wfs_key not in params:
        # Try without index
        wfs_key = f"sh_{wfs_type}"
        if wfs_key not in params:
            if verbose:
                print(f"Warning: WFS key {wfs_key} not found in configuration")
            wfs_key = None

    # Build filename components
    filename_parts = ["IM_syn"]

    # Add source information
    if wfs_key:
        # Extract source configuration
        source_config = extract_source_config(params, wfs_key)
        source_info = build_source_filename_part(source_config)
        filename_parts.append(source_info)

        # Add WFS information
        wfs_params = params[wfs_key]
        wfs_info = build_wfs_filename_part(wfs_params, wfs_type)
        filename_parts.append(wfs_info)

    # Add DM or Layer information
    if dm_index is not None:
        dm_key = f"dm{dm_index}"
        if dm_key in params:
            dm_info = build_component_filename_part(
                params[dm_key],
                'dm',
                n_modes_override=n_modes
            )
            filename_parts.append(dm_info)
        elif dm_index == 1 and "dm" in params:
            # Fallback: use 'dm' if 'dm1' not found and index==1
            dm_info = build_component_filename_part(
                params["dm"],
                'dm',
                n_modes_override=n_modes
            )
            filename_parts.append(dm_info)
    else:
        layer_key = f"layer{layer_index}"
        if layer_key in params:
            layer_info = build_component_filename_part(
                params[layer_key],
                'layer',
                n_modes_override=n_modes
            )
            filename_parts.append(layer_info)

    # Add timestamp if requested
    if timestamp:
        from datetime import datetime
        timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename_parts.append(timestamp_str)

    # Join parts and add extension
    filename = "_".join(filename_parts) + ".fits"

    return filename


def extract_source_config(params, wfs_key):
    """
    Extract source configuration from WFS key.
    
    Args:
        params (dict): Full configuration
        wfs_key (str): WFS key (e.g., 'sh_lgs1')
        
    Returns:
        dict: Source configuration
    """
    # Extract source type and index from WFS key
    import re
    match = re.search(r'(lgs|ngs|ref)(\d*)', wfs_key)
    if match:
        source_type = match.group(1)
        source_index = match.group(2) if match.group(2) else '1'
        source_key = f"source_{source_type}{source_index}"

        if source_key in params:
            return params[source_key]

    # Default configuration
    return {
        'polar_coordinates': [0.0, 0.0],
        'height': float('inf') if 'ngs' in wfs_key or 'ref' in wfs_key else 90000.0
    }


def extract_projection_params(params):
    proj_params = params.get('projection', {})
    if not proj_params:
        return None

    reg_factor = proj_params.get('reg_factor', 1e-6)
    ifunc_inv_tag = proj_params['ifunc_inverse_tag']

    opt_sources = proj_params['opt_sources']
    polar_coordinates = opt_sources['polar_coordinates']  # list of [r, theta]
    weights = opt_sources['weights']

    if len(polar_coordinates) != len(weights):
        raise ValueError(f"Mismatch: {len(polar_coordinates)} coordinates vs"
                         f" {len(weights)} weights")

    # Build the projection parameters dictionary
    source_list = []
    for (r, theta), weight in zip(polar_coordinates, weights):
        source_list.append({
            'polar_coordinates': [r, theta],
            'weight': weight
        })

    return {
        'reg_factor': reg_factor,
        'ifunc_inverse_tag': ifunc_inv_tag,
        'opt_sources': source_list
    }


def generate_pm_filename(config_file, opt_index=None,
                        dm_index=None, layer_index=None,
                        timestamp=False, verbose=False):
    """
    Generate a specific projection matrix filename based on optical source and DM/layer indices.

    Args:
        config_file (str or dict): Path to YAML/PRO configuration file or config dictionary
        opt_index (int, optional): Index of the optical source to use (1-based)
        dm_index (int, optional): Index of the DM to use (1-based)
        layer_index (int, optional): Index of the layer to use (1-based)
        timestamp (bool, optional): Whether to include timestamp in the filename
        verbose (bool, optional): Whether to print verbose output

    Returns:
        str: Filename for the projection matrix with the specified parameters
    """
    # Check if both dm_index and layer_index are provided
    if dm_index is not None and layer_index is not None:
        raise ValueError("Cannot specify both dm_index and layer_index at the same time")

    if dm_index is None and layer_index is None:
        raise ValueError("Must specify either dm_index or layer_index")

    # Load configuration
    if isinstance(config_file, str):
        config = parse_params_file(config_file)
    else:
        config = config_file

    # Convert indices to strings for comparison
    opt_index_str = str(opt_index) if opt_index is not None else None
    dm_index_str = str(dm_index) if dm_index is not None else None
    layer_index_str = str(layer_index) if layer_index is not None else None

    # Find all optical sources
    opt_sources = extract_opt_list(config)

    # Extract DM and layer configurations
    dm_list = extract_dm_list(config) if dm_index is not None else []
    layer_list = extract_layer_list(config) if layer_index is not None else []

    # Filter optical sources based on opt_index
    if opt_index_str:
        filtered_sources = [src for src in opt_sources if src['index'] == int(opt_index_str)]
    else:
        filtered_sources = opt_sources

    # Filter DM or layer list based on provided index
    if dm_index_str:
        filtered_components = [dm for dm in dm_list if dm['index'] == dm_index_str]
        component_type = "dm"
    elif layer_index_str:
        filtered_components = [layer for layer in layer_list if layer['index'] == layer_index_str]
        component_type = "layer"
    else:
        # Default to first DM if no specific component requested
        filtered_components = dm_list
        component_type = "dm"

    # If we couldn't find matching source or component, return None
    if not filtered_sources or not filtered_components:
        if not filtered_sources:
            print("Warning: No matching optical source found with the specified parameters")
        if not filtered_components:
            print("Warning: No matching DM or layer found with the specified parameters")
        return None

    # Select the first source and component from the filtered lists
    selected_source = filtered_sources[0]
    selected_component = filtered_components[0]

    if verbose:
        print(f"Selected optical source: {selected_source['name']}"
              f" (index: {selected_source['index']})")
        print(f"Selected {component_type}: {selected_component['name']}"
              f" (index: {selected_component['index']})")

    # Extract source information
    source_config = selected_source['config']
    source_coords = source_config.get('polar_coordinates', None)
    if source_coords is None:
        source_coords = source_config.get('polar_coordinate', [0.0, 0.0])
    source_height = source_config.get('height', float('inf'))

    # Extract component configuration
    component_config = selected_component['config']
    component_height = component_config.get('height', 0)

    # *** Extract nmodes and start_mode ***
    component_nmodes = component_config.get('nmodes', None)
    component_start_mode = component_config.get('start_mode', 0)

    # Generate filename parts
    base_name = "PM_syn"
    parts = [base_name]

    # Source coordinates
    if source_coords and (source_coords[0] != 0.0 or source_coords[1] != 0.0):
        dist, angle = source_coords
        parts.append(f"pd{dist:.1f}a{angle:.0f}")
    else:
        parts.append("pd0.0a0")

    # Source height (only include if not infinite)
    if not np.isinf(source_height):
        parts.append(f"h{source_height:.0f}")

    # *** Component part - always use "dm" prefix, include modes info ***
    # Component height
    parts.append(f"dmH{component_height:.1f}")

    # *** Add modes information ***
    if component_nmodes is not None:
        if component_start_mode > 0:
            # Format: mn{nmodes}s{start_mode}
            # Actual modes used = nmodes - start_mode
            actual_modes = component_nmodes - component_start_mode
            parts.append(f"mn{actual_modes}s{component_start_mode}")
        else:
            # Format: mn{nmodes}
            parts.append(f"mn{component_nmodes}")

    # Add ifunc/m2c tags
    ifunc_tag = get_tag_or_object(component_config, 'ifunc')
    if ifunc_tag is not None:
        parts.append(f"ifunc_{ifunc_tag}")

    m2c_tag = get_tag_or_object(component_config, 'm2c')
    if m2c_tag is not None:
        parts.append(f"m2c_{m2c_tag}")

    # Add timestamp if requested
    if timestamp:
        import datetime
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        parts.append(ts)

    # Join all parts with underscores and add extension
    filename = "_".join(parts) + ".fits"
    return filename


def generate_pm_filenames(config_file, timestamp=False):
    """
    Generate projection matrix filenames for all optical source and DM combinations.
    
    Args:
        config_file (str or dict): Path to YAML/PRO file or config dictionary
        timestamp (bool, optional): Whether to include timestamp in filenames
        
    Returns:
        list: List of projection matrix filenames
    """
    # Load configuration
    if isinstance(config_file, str):
        config = parse_params_file(config_file)
    else:
        config = config_file

    # Find all optical sources
    opt_sources = extract_opt_list(config)

    # Extract all DM configurations
    dm_list = extract_dm_list(config)

    # Output list of filenames
    filenames = []

    # Generate filenames for all optical source and DM combinations
    for source in opt_sources:
        source_config = source['config']
        source_coords = source_config.get('polar_coordinate',
                                          source_config.get('polar_coordinates',
                                                            [0.0, 0.0]))
        source_height = source_config.get('height', float('inf'))

        for dm in dm_list:
            dm_config = dm['config']
            dm_height = dm_config.get('height', 0)

            # Generate filename parts
            base_name = "PM_syn"
            parts = [base_name]

            # Specific source identifier
            parts.append(f"opt{source['index']}")

            # Source coordinates
            if source_coords:
                dist, angle = source_coords
                parts.append(f"pd{dist:.1f}a{angle:.0f}")

            # Source height (only include if not infinite)
            if not np.isinf(source_height):
                parts.append(f"h{source_height:.0f}")

            # DM info
            parts.append(f"dm{dm['index']}")
            parts.append(f"dmH{dm_height}")

            # Add DM-specific parts
            parts.extend(build_component_filename_part(dm_config, 'dm'))

            # Add timestamp if requested
            if timestamp:
                ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                parts.append(ts)

            # Join all parts with underscores and add extension
            filename = "_".join(parts) + ".fits"
            filenames.append(filename)

    return filenames


def generate_cov_filename(component_config, pup_diam_m, r0, L0):
    """
    Generate a unique filename for a covariance matrix, handling _tag/_object.
    Args:
        component_config (dict): DM or layer config
        component_type (str): 'dm' or 'layer'
        pup_diam_m (float): pupil diameter in meters
        r0 (float): Fried parameter
        L0 (float): Outer scale
    Returns:
        str: filename
    """
    # Prefer m2c_tag, then ifunc_tag, then fallback
    base_tag = get_tag_or_object(component_config, 'm2c')
    if base_tag is None:
        base_tag = get_tag_or_object(component_config, 'ifunc')

    diam_str = f"{pup_diam_m:.1f}".strip()
    r0_str = f"{r0:.3f}".strip()
    L0_str = f"{L0:.1f}".strip()

    filename = f"covariance_{base_tag}_{diam_str}diam_{r0_str}r0_{L0_str}L0.fits"
    return filename, base_tag


def compute_layer_weights_from_turbulence(params,
                                          component_indices,
                                          component_type='layer',
                                          verbose=False):
    """
    Compute layer weights from atmospheric turbulence profile.
    
    Associates each layer/DM with the nearest turbulence height and assigns
    the corresponding CN² contribution as weight.
    
    Args:
        params (dict): Configuration dictionary with 'atmo' section
        component_indices (list): List of component indices to compute weights for
        component_type (str): Type of component ('layer' or 'dm')
        verbose (bool): Whether to print detailed information
        
    Returns:
        numpy.ndarray: Normalized weights for each component (sum to 1)
        
    Raises:
        ValueError: If 'atmo' section or required fields are missing
    """
    if 'atmo' not in params:
        raise ValueError("'atmo' section not found in configuration. Cannot compute layer weights.")

    atmo = params['atmo']

    cn2_key = 'cn2'
    if 'Cn2' in atmo:
        cn2_key = 'Cn2'
    elif 'CN2' in atmo:
        cn2_key = 'CN2'

    # Get turbulence profile
    if 'heights' not in atmo or cn2_key not in atmo:
        raise ValueError("'atmo.heights' and 'atmo.cn2' must be specified"
                         " for automatic weight computation")

    turb_heights = np.array(atmo['heights'])
    turb_cn2 = np.array(atmo[cn2_key])

    if len(turb_heights) != len(turb_cn2):
        raise ValueError(f"Mismatch: {len(turb_heights)} heights vs {len(turb_cn2)} CN² values")

    # Normalize CN² to sum to 1
    turb_cn2_norm = turb_cn2 / np.sum(turb_cn2)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Computing Layer Weights from Atmospheric Profile")
        print(f"{'='*60}")
        print(f"  Turbulence heights: {turb_heights}")
        print(f"  CN² values: {turb_cn2}")
        print(f"  Normalized CN²: {turb_cn2_norm}")

    # Get layer heights
    layer_heights = []
    for comp_idx in component_indices:
        comp_key = f'{component_type}{comp_idx}'
        if comp_key not in params:
            raise ValueError(f"Component {comp_key} not found in configuration")

        layer_height = params[comp_key].get('height', 0.0)
        layer_heights.append(layer_height)

    layer_heights = np.array(layer_heights)

    if verbose:
        print(f"\n  Layer heights: {layer_heights}")

    # Associate each turbulence height with its nearest layer
    # For each turbulence height, find which layer it's closest to
    weights = np.zeros(len(layer_heights))

    for j, turb_h in enumerate(turb_heights):
        # Find the layer closest to this turbulence height
        idx_nearest_layer = np.argmin(np.abs(layer_heights - turb_h))

        # Add this turbulence contribution to the nearest layer
        weights[idx_nearest_layer] += turb_cn2_norm[j]

    # Normalize weights to sum to 1 (should already be ~1, but ensure it)
    weights = weights / np.sum(weights)

    if verbose:
        print(f"\n  Final normalized weights: {weights}")
        print(f"{'='*60}\n")

    return weights


def compute_mmse_reconstructor(interaction_matrix, C_atm,
                              noise_variance=None, C_noise=None,
                              cinverse=False, use_inverse=False,
                              xp=np, dtype=np.float32,
                              verbose=False):
    """
    Compute the Minimum Mean Square Error (MMSE) reconstructor.
    
    Args:
        interaction_matrix (numpy.ndarray): Interaction matrix A relating modes to slopes
        C_atm (numpy.ndarray): Covariance matrix of atmospheric modes (Cx)
        noise_variance (list, optional): List of noise variances per WFS. 
                                        Used to build C_noise if C_noise is None.
        C_noise (numpy.ndarray, optional): Covariance matrix of measurement noise (Cz).
                                         If None, it's built from noise_variance.
        cinverse (bool, optional): If True, C_atm and C_noise are already inverted.
        use_inverse (bool, optional): If True, use standard inverse; otherwise, use pseudo-inverse.
        xp (module, optional): Numerical library to use (e.g., numpy or cupy).
        dtype (data-type, optional): Data type for computations.
        verbose (bool, optional): Whether to print detailed information during computation.
        
    Returns:
        numpy.ndarray: MMSE reconstructor matrix
    """
    if verbose:
        print("Starting MMSE reconstructor computation")

    # Setup matrices
    A = to_xp(xp, interaction_matrix, dtype=dtype)

    # Handle noise covariance matrix
    if C_noise is None and noise_variance is not None:
        n_slopes_total = A.shape[1]
        n_wfs = len(noise_variance)
        n_slopes_per_wfs = n_slopes_total // n_wfs

        if verbose:
            print(f"Building noise covariance matrix for {n_wfs}"
                  f" WFSs with {n_slopes_per_wfs} slopes each")

        C_noise = xp.zeros((n_slopes_total, n_slopes_total), dtype=dtype)
        for i in range(n_wfs):
            # Set the diagonal elements for this WFS
            start_idx = i * n_slopes_per_wfs
            end_idx = (i + 1) * n_slopes_per_wfs
            C_noise[start_idx:end_idx, start_idx:end_idx] = \
                noise_variance[i] * xp.eye(n_slopes_per_wfs, dtype=dtype)
    else:
        C_noise = to_xp(xp, C_noise, dtype=dtype) if C_noise is not None else None

    if C_atm is not None:
        C_atm = to_xp(xp, C_atm, dtype=dtype)
    else:
        raise ValueError("C_atm (atmospheric covariance matrix) must be provided")

    # Check dimensions
    if A.shape[1] != C_atm.shape[0]:
        raise ValueError(f"A ({A.shape}) and C_atm ({C_atm.shape}) must have compatible dimensions")

    if C_noise is not None and A.shape[0] != C_noise.shape[0]:
        raise ValueError(f"A ({A.shape}) and C_noise"
                         f" ({C_noise.shape}) must have compatible dimensions")

    # Compute inverses if needed
    if not cinverse:
        # Check if matrices are diagonal
        if C_noise is not None:
            is_diag_noise = xp.all(xp.abs(xp.diag(xp.diag(C_noise)) - C_noise) < 1e-10)

            if is_diag_noise:
                if verbose:
                    print("C_noise is diagonal, using optimized inversion")
                C_noise_inv = xp.diag(1.0 / xp.diag(C_noise))
            else:
                if verbose:
                    print("Inverting C_noise matrix")
                try:
                    C_noise_inv = xp.linalg.inv(C_noise)
                except xp.linalg.LinAlgError:
                    if verbose:
                        print("Warning: C_noise inversion failed, using pseudo-inverse")
                    C_noise_inv = xp.linalg.pinv(C_noise)
        else:
            # Default: identity matrix (no noise)
            if verbose:
                print("No C_noise provided, using identity matrix")
            C_noise_inv = xp.eye(A.shape[1], dtype=dtype)

        is_diag_atm = xp.all(xp.abs(xp.diag(xp.diag(C_atm)) - C_atm) < 1e-10)

        if is_diag_atm:
            if verbose:
                print("C_atm is diagonal, using optimized inversion")
            C_atm_inv = xp.diag(1.0 / xp.diag(C_atm))
        else:
            if verbose:
                print("Inverting C_atm matrix")
            if use_inverse:
                C_atm_inv = xp.linalg.inv(C_atm)
            else:
                if verbose:
                    print("Warning: Using pseudo-inverse")
                C_atm_inv = xp.linalg.pinv(C_atm)
    else:
        # Matrices are already inverted
        C_atm_inv = C_atm
        C_noise_inv = C_noise if C_noise is not None else xp.eye(A.shape[1], dtype=dtype)
    # Compute H = A' Cz^(-1) A + Cx^(-1)
    if verbose:
        print("Computing H = A' Cz^(-1) A + Cx^(-1)")

    # Check if C_noise_inv is scalar
    if isinstance(C_noise_inv, (int, float)) \
        or (hasattr(C_noise_inv, 'size') and C_noise_inv.size == 1):
        H = C_noise_inv * xp.dot(A.T, A) + C_atm_inv
    else:
        H = xp.dot(A.T, xp.dot(C_noise_inv, A)) + C_atm_inv

    # *** DIAGNOSTICS ***
    if verbose:
        print(f"\nH matrix:")
        print(f"  Shape: {H.shape}")
        print(f"  dtype: {H.dtype}")
        print(f"  Is finite: {xp.all(xp.isfinite(H))}")
        print(f"  Memory: {H.nbytes / 1024**3:.2f} GB")

    # Check for NaN/Inf
    if not xp.all(xp.isfinite(H)):
        raise ValueError("H contains NaN or Inf values!")

    # Compute H^(-1)
    if verbose:
        print("Inverting H")

    try:
        if use_inverse:
            H_inv = xp.linalg.inv(H)
        else:
            if verbose:
                print(f"Using pseudo-inverse")
            H_inv = xp.linalg.pinv(H)
    except Exception as e:
        print(f"ERROR during H inversion: {e}")
        print(f"H shape: {H.shape}, dtype: {H.dtype}")
        raise

    # Compute W = H^(-1) A' Cz^(-1)
    if verbose:
        print("Computing W = H^(-1) A' Cz^(-1)")

    # Check if C_noise_inv is scalar
    if isinstance(C_noise_inv, (int, float)) \
        or (hasattr(C_noise_inv, 'size') and C_noise_inv.size == 1):
        W_mmse = C_noise_inv * xp.dot(H_inv, A.T)
    else:
        W_mmse = xp.dot(H_inv, xp.dot(A.T, C_noise_inv))

    if verbose:
        print("MMSE reconstruction matrix computed")
        print(f"Matrix shape: {W_mmse.shape}")

    return cpuArray(W_mmse)


__all__ = [
    'parse_params_file',
    'is_simple_config',
    'extract_wfs_list',
    'extract_dm_list',
    'extract_layer_list',
    'extract_opt_list',
    'validate_opt_sources',
    'load_pupilstop',
    'load_influence_functions',
    'find_subapdata',
    'extract_projection_params',
    'extract_source_coordinates',
    'extract_source_info',
    'determine_source_type',
    'wfs_fov_from_config',
    'build_source_filename_part',
    'build_pupil_filename_part',
    'build_wfs_filename_part',
    'build_component_filename_part',
    'generate_im_filename',
    'generate_pm_filename',
    'generate_pm_filenames',
    'generate_cov_filename',
    'compute_layer_weights_from_turbulence',
    'compute_mmse_reconstructor',
]
