from synim import xp, cpuArray, to_xp, float_dtype

import matplotlib.pyplot as plt
from synim.utils import (
    apply_mask,
    rotshiftzoom_array,
    shiftzoom_from_source_dm_params,
    has_transformations,
    dm2d_to_3d,
    dm3d_to_2d
)

def transpose_base_array_for_specula(base_inv_array, pup_mask_original, verbose=False):
    """
    Transpose base_inv_array for SPECULA convention.
    
    CRITICAL: Requires the ORIGINAL (non-transposed) pupil mask to correctly
    extract and re-map pixel indices after transposition.
    
    Handles both 2D and 3D base arrays:
    - 3D arrays: simple transpose (height, width, nmodes) → (width, height, nmodes)
    - 2D arrays: reconstructs 3D using ORIGINAL mask, transposes, then re-extracts
    
    Parameters:
    - base_inv_array: 2D or 3D numpy array
        - 2D can be: (nmodes, npixels_total) [IFunc] or (npixels_total, nmodes) [IFuncInv]
        - 3D: (height, width, nmodes)
    - pup_mask_original: numpy 2D array, ORIGINAL (non-transposed) pupil mask
    - verbose: bool, print debug information
    
    Returns:
    - transposed_array: numpy array with swapped X-Y coordinates
    """

    # *** MODIFIED: Convert inputs to xp arrays with correct dtype ***
    base_inv_array = to_xp(xp, base_inv_array, dtype=float_dtype)
    pup_mask_original = to_xp(xp, pup_mask_original, dtype=float_dtype)

    if base_inv_array.ndim == 3:
        # ============================================================
        # SIMPLE CASE: 3D array - just transpose
        # ============================================================
        transposed = xp.transpose(base_inv_array, (1, 0, 2))

        if verbose:
            print(f"  Transposed 3D base: {base_inv_array.shape} → {transposed.shape}")

        return transposed

    elif base_inv_array.ndim == 2:
        # ============================================================
        # COMPLEX CASE: 2D array - use dm2d_to_3d and dm3d_to_2d
        # ============================================================

        n_rows, n_cols = base_inv_array.shape

        # Get original mask info
        pup_pixels_total = pup_mask_original.shape[0] * pup_mask_original.shape[1]
        n_valid_pixels = int(xp.sum(pup_mask_original > 0.5))

        # Determine format
        if n_cols == pup_pixels_total:
            # IFunc format: (nmodes, npixels_total)
            format_name = "IFunc (full pupil)"
            base_2d = base_inv_array

        elif n_rows == pup_pixels_total:
            # IFuncInv format: (npixels_total, nmodes)
            format_name = "IFuncInv (full pupil)"
            base_2d = base_inv_array.T  # Convert to IFunc format

        elif n_cols == n_valid_pixels:
            # IFunc format with only valid pixels: (nmodes, npixels_valid)
            format_name = "IFunc (valid pixels only)"
            base_2d = base_inv_array

        elif n_rows == n_valid_pixels:
            # IFuncInv format with only valid pixels: (npixels_valid, nmodes)
            format_name = "IFuncInv (valid pixels only)"
            base_2d = base_inv_array.T

        else:
            # Cannot determine format
            if verbose:
                print(f"  ERROR: Cannot determine format of 2D base")
                print(f"    Shape: {base_inv_array.shape}")
                print(f"    Total pixels: {pup_pixels_total}")
                print(f"    Valid pixels: {n_valid_pixels}")
            raise ValueError(f"Cannot determine format of 2D base array"
                           f" with shape {base_inv_array.shape}")

        if verbose:
            print(f"  Detected format: {format_name}")
            print(f"  Shape: {base_inv_array.shape}")

        # *** USE dm2d_to_3d TO RECONSTRUCT ***
        # base_2d is now always in IFunc format: (nmodes, npixels)
        base_3d_orig = dm2d_to_3d(base_2d, pup_mask_original, normalize=False)

        if verbose:
            print(f"  Reconstructed 3D with dm2d_to_3d: {base_3d_orig.shape}")

        # Transpose the 3D
        base_3d_transposed = xp.transpose(base_3d_orig, (1, 0, 2))

        if verbose:
            print(f"  Transposed 3D: {base_3d_transposed.shape}")

        if n_cols == n_valid_pixels or n_rows == n_valid_pixels:
            # Get transposed mask
            pup_mask_transposed = xp.transpose(pup_mask_original)

            # *** MODIFIED: USE dm3d_to_2d TO EXTRACT ***
            base_2d_transposed = dm3d_to_2d(base_3d_transposed, pup_mask_transposed)

            # normalize by number of valid pixels
            base_2d_transposed /= n_valid_pixels

            if verbose:
                print(f"  Re-extracted to 2D with dm3d_to_2d: {base_2d_transposed.shape}")

            # Return in same format as input
            if n_rows == n_valid_pixels:
                # Was IFuncInv, return transposed
                return base_2d_transposed.T
            else:
                # Was IFunc, return as-is
                return base_2d_transposed
        else:
            # Full pupil - return 3D transposed
            if verbose:
                print(f"  Returning full 3D transposed: {base_3d_transposed.shape}")
            return base_3d_transposed

    else:
        raise ValueError(f"base_inv_array must be 2D or 3D, got {base_inv_array.ndim}D")


def projection_matrix(pup_diam_m, pup_mask,
                      dm_array, dm_mask,
                      base_inv_array, dm_height,
                      dm_rotation, base_rotation,
                      base_translation, base_magnification,
                      gs_pol_coo, gs_height,
                      verbose=False, specula_convention=True,
                      specula_convention_inv=False):
    """
    Computes a projection matrix for DM modes onto a desired basis.
    Uses intelligent workflow selection like interaction_matrix.

    Parameters:
    - pup_diam_m: float, size in m of the side of the pupil
    - pup_mask: numpy 2D array, pupil mask (n_pup x n_pup)
    - dm_array: numpy 3D array, Deformable Mirror 2D shapes (n x n x n_dm_modes)
    - dm_mask: numpy 2D array, DM mask (n x n)
    - base_inv_array: numpy 2D or 3D array, inverted basis for projection
                      Can be:
                      - 2D: (nmodes, npixels_valid) - IFunc format
                      - 2D: (npixels_valid, nmodes) - IFuncInv format  
                      - 3D: (npix, npix, nmodes) - full 3D array
    - dm_height: float, conjugation altitude of the Deformable Mirror
    - dm_rotation: float, rotation in deg of the Deformable Mirror with respect to the pupil
    - base_rotation: float, rotation of the basis in deg
    - base_translation: tuple, translation of the basis (x, y) in pixels
    - base_magnification: tuple, magnification of the basis (x, y)
    - gs_pol_coo: tuple, polar coordinates of the guide star radius in arcsec and angle in deg
    - gs_height: float, altitude of the guide star
    - verbose: bool, optional, display verbose output
    - specula_convention: bool, optional, use SPECULA convention (transpose arrays)
    - specula_convention_inv: bool, optional, use SPECULA convention for base_inv_array

    Returns:
    - projection: numpy 2D array, projection matrix (n_dm_modes, n_base_modes)
    
    Workflow Selection:
    - SEPARATED: Used when EITHER DM OR Base has transformations (not both)
                 Applies transformations in 2 steps (more flexible)
    - COMBINED: Used when BOTH DM AND Base have transformations
                Applies transformations in 1 step (avoids double interpolation)
    """

    # ================================================================
    # STEP 1: Detect which transformations are present
    # ================================================================
    has_dm_transform = has_transformations(dm_rotation, (0, 0), (1, 1)) or \
                       gs_pol_coo != (0, 0) or dm_height != 0
    has_base_transform = has_transformations(base_rotation, base_translation, base_magnification)

    # Choose workflow: COMBINED only if BOTH have transformations
    use_combined = has_dm_transform and has_base_transform

    if verbose:
        print(f"\n{'='*60}")
        print(f"Projection Matrix Computation")
        print(f"{'='*60}")
        print(f"DM transformations: {has_dm_transform}")
        print(f"  - Height: {dm_height} m")
        print(f"  - Rotation: {dm_rotation}°")
        print(f"  - GS position: {gs_pol_coo}")
        print(f"Base transformations: {has_base_transform}")
        print(f"  - Rotation: {base_rotation}°")
        print(f"  - Translation: {base_translation}")
        print(f"  - Magnification: {base_magnification}")
        print(f"Using {'COMBINED' if use_combined else 'SEPARATED'} workflow")
        print(f"{'='*60}\n")

    # ================================================================
    # STEP 2: Apply SPECULA convention + Convert to xp with dtype
    # ================================================================
    if specula_convention:
        # *** Convert inputs to xp FIRST, then save original mask ***
        dm_array = to_xp(xp, dm_array, dtype=float_dtype)
        dm_mask = to_xp(xp, dm_mask, dtype=float_dtype)
        pup_mask = to_xp(xp, pup_mask, dtype=float_dtype)
        base_inv_array = to_xp(xp, base_inv_array, dtype=float_dtype)

        # Save ORIGINAL (non-transposed) mask for transpose_base_array_for_specula
        pup_mask_original = pup_mask.copy()

        # Now transpose
        dm_array = xp.transpose(dm_array, (1, 0, 2))
        dm_mask = xp.transpose(dm_mask)
        pup_mask = xp.transpose(pup_mask)

        if specula_convention_inv:
            # *** PASS xp ARRAY (not CPU!) ***
            base_inv_array = transpose_base_array_for_specula(
                base_inv_array,
                pup_mask_original,  # Already xp array with correct dtype
                verbose=False
            )
    else:
        # *** Still convert to xp even without SPECULA convention ***
        dm_array = to_xp(xp, dm_array, dtype=float_dtype)
        dm_mask = to_xp(xp, dm_mask, dtype=float_dtype)
        pup_mask = to_xp(xp, pup_mask, dtype=float_dtype)
        base_inv_array = to_xp(xp, base_inv_array, dtype=float_dtype)

    # ================================================================
    # STEP 3: Setup basic parameters
    # ================================================================
    pup_diam_pix = pup_mask.shape[0]
    pixel_pitch = pup_diam_m / pup_diam_pix

    if dm_mask.shape[0] != dm_array.shape[0]:
        raise ValueError('DM and mask arrays must have the same dimensions.')

    # Calculate DM transformations based on guide star geometry
    dm_translation, dm_magnification = shiftzoom_from_source_dm_params(
        gs_pol_coo, gs_height, dm_height, pixel_pitch
    )
    output_size = (pup_diam_pix, pup_diam_pix)

    # ================================================================
    # STEP 4: Apply transformations (COMBINED or SEPARATED)
    # ================================================================
    if use_combined:
        # ============================================================
        # COMBINED WORKFLOW: Apply DM + Base transformations together
        # ============================================================
        # This uses a SINGLE interpolation step, avoiding cumulative errors

        if verbose:
            print(f'[COMBINED] Applying DM+Base transformations in one step:')
            print(f'  DM translation: {dm_translation} pixels')
            print(f'  DM rotation: {dm_rotation}°')
            print(f'  DM magnification: {dm_magnification}')
            print(f'  Base translation: {base_translation} pixels')
            print(f'  Base rotation: {base_rotation}°')
            print(f'  Base magnification: {base_magnification}')

        # Transform DM array with BOTH DM and Base transformations
        trans_dm_array = rotshiftzoom_array(
            dm_array,
            dm_translation=dm_translation,      # From guide star geometry
            dm_rotation=dm_rotation,            # DM rotation
            dm_magnification=dm_magnification,  # From guide star geometry
            wfs_translation=base_translation,   # Base translation
            wfs_rotation=base_rotation,         # Base rotation
            wfs_magnification=base_magnification, # Base magnification
            output_size=output_size
        )

        # Transform DM mask (ONLY DM transformations, not base)
        trans_dm_mask = rotshiftzoom_array(
            dm_mask,
            dm_translation=dm_translation,
            dm_rotation=dm_rotation,
            dm_magnification=dm_magnification,
            wfs_translation=(0, 0),  # No base transformation for mask
            wfs_rotation=0,
            wfs_magnification=(1, 1),
            output_size=output_size
        )
        trans_dm_mask[trans_dm_mask < 0.5] = 0

        # Transform pupil mask (ONLY Base transformations, not DM)
        trans_pup_mask = rotshiftzoom_array(
            pup_mask,
            dm_translation=(0, 0),  # No DM transformation for pupil
            dm_rotation=0,
            dm_magnification=(1, 1),
            wfs_translation=base_translation,   # Base translation
            wfs_rotation=base_rotation,         # Base rotation
            wfs_magnification=base_magnification, # Base magnification
            output_size=output_size
        )
        trans_pup_mask[trans_pup_mask < 0.5] = 0

    else:
        # ============================================================
        # SEPARATED WORKFLOW: Apply DM and Base transformations separately
        # ============================================================
        # This uses TWO interpolation steps but is more flexible

        if verbose:
            print(f'[SEPARATED] Applying transformations in two steps:')
            print(f'  Step 1 - DM transformations:')
            print(f'    Translation: {dm_translation} pixels')
            print(f'    Rotation: {dm_rotation}°')
            print(f'    Magnification: {dm_magnification}')

        # Transform DM array (ONLY DM transformations)
        trans_dm_array = rotshiftzoom_array(
            dm_array,
            dm_translation=dm_translation,
            dm_rotation=dm_rotation,
            dm_magnification=dm_magnification,
            wfs_translation=(0, 0),  # No base transformation yet
            wfs_rotation=0,
            wfs_magnification=(1, 1),
            output_size=output_size
        )

        # Transform DM mask (ONLY DM transformations)
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

        # Transform pupil mask (ONLY Base transformations)
        trans_pup_mask = rotshiftzoom_array(
            pup_mask,
            dm_translation=(0, 0),  # No DM transformation
            dm_rotation=0,
            dm_magnification=(1, 1),
            wfs_translation=base_translation,   # Base translation
            wfs_rotation=base_rotation,         # Base rotation
            wfs_magnification=base_magnification, # Base magnification
            output_size=output_size
        )
        trans_pup_mask[trans_pup_mask < 0.5] = 0

        if has_base_transform and verbose:
            print(f'  Step 2 - Base transformations:')
            print(f'    Translation: {base_translation} pixels')
            print(f'    Rotation: {base_rotation}°')
            print(f'    Magnification: {base_magnification}')

    # ================================================================
    # STEP 5: Validate transformed arrays
    # ================================================================
    if xp.max(trans_dm_mask) <= 0:
        raise ValueError('Transformed DM mask is empty.')
    if xp.max(trans_pup_mask) <= 0:
        raise ValueError('Transformed pupil mask is empty.')

    # Apply DM mask to DM array
    trans_dm_array = apply_mask(trans_dm_array, trans_dm_mask, in_place=True)

    if verbose:
        print(f'  ✓ DM array transformed: {trans_dm_array.shape}')
        print(f'  ✓ DM mask valid pixels: {xp.sum(trans_dm_mask > 0.5)}')
        print(f'  ✓ Pupil mask valid pixels: {xp.sum(trans_pup_mask > 0.5)}')

    # ================================================================
    # STEP 6: Find valid pixels (intersection of DM and pupil)
    # ================================================================
    valid_mask = trans_pup_mask.copy() #trans_dm_mask * trans_pup_mask
    idx_valid = xp.where(valid_mask > 0.5)  # Returns tuple: (row_indices, col_indices)
    n_valid_pixels = len(idx_valid[0])

    if verbose:
        print(f'  ✓ Valid pixels (intersection): {n_valid_pixels}')

    # ================================================================
    # STEP 7: Extract valid pixel values from DM array
    # ================================================================
    dm_valid_values = trans_dm_array[idx_valid[0], idx_valid[1], :]
    # Result shape: (n_valid_pixels, n_modes)

    if verbose:
        print(f'  ✓ DM valid values extracted: {dm_valid_values.shape}')

    # ================================================================
    # STEP 8: Extract valid pixel values from base_inv_array
    # ================================================================
    # Handle different input formats:

    if base_inv_array.ndim == 2:
        # --------------------------------------------------------
        # 2D FORMAT: Could be IFunc or IFuncInv
        # --------------------------------------------------------
        n_rows, n_cols = base_inv_array.shape

        if n_cols == n_valid_pixels:
            # IFunc format: (nmodes, npixels_valid)
            # Each ROW is a mode, need to transpose
            base_valid_values = base_inv_array.T  # → (npixels_valid, nmodes)
            base_format = "IFunc (nmodes, npixels)"

        elif n_rows == n_valid_pixels:
            # IFuncInv format: (npixels_valid, nmodes)
            # Already in correct format!
            base_valid_values = base_inv_array
            base_format = "IFuncInv (npixels, nmodes)"

        else:
            raise ValueError(
                f"Base 2D shape {base_inv_array.shape} doesn't match "
                f"{n_valid_pixels} valid pixels. "
                f"Expected either ({n_valid_pixels}, nmodes) or (nmodes, {n_valid_pixels})"
            )

        if verbose:
            print(f'  Base format detected: {base_format}')
            print(f'  Base shape: {base_inv_array.shape} → {base_valid_values.shape}')

    elif base_inv_array.ndim == 3:
        # --------------------------------------------------------
        # 3D FORMAT: (height, width, n_modes)
        # --------------------------------------------------------
        # Extract using the SAME indices as DM
        base_valid_values = base_inv_array[idx_valid[0], idx_valid[1], :]
        # Result shape: (n_valid_pixels, n_modes_base)

        if verbose:
            print(f'  Base 3D: {base_inv_array.shape} → {base_valid_values.shape}')

    else:
        raise ValueError(f"base_inv_array must be 2D or 3D, got {base_inv_array.ndim}D")

    # ================================================================
    # STEP 9: Compute projection matrix
    # ================================================================
    # Matrix multiplication:
    # dm_valid_values:   (n_valid_pixels, n_dm_modes)
    # base_valid_values: (n_valid_pixels, n_base_modes)
    #
    # We want: projection = DM^T × Base
    # Result: (n_dm_modes, n_base_modes)

    # *** MODIFIED: Use xp.dot instead of np.dot ***
    projection = xp.dot(dm_valid_values.T, base_valid_values)

    if verbose:
        print(f'\n  ✓ PROJECTION COMPUTED: {projection.shape}')
        print(f'    DM modes: {projection.shape[0]}')
        print(f'    Base modes: {projection.shape[1]}')
        print(f'    Valid pixels used: {n_valid_pixels}')
        print(f'{"="*60}\n')

    # ================================================================
    # STEP 10: Optional display
    # ================================================================
    plot_debug = False
    if plot_debug:
        # *** MODIFIED: Convert to CPU for plotting if needed ***
        valid_mask_cpu = cpuArray(valid_mask)
        trans_dm_array_cpu = cpuArray(trans_dm_array)
        projection_cpu = cpuArray(projection)

        plt.figure(figsize=(8, 6))
        plt.imshow(valid_mask_cpu, cmap='gray')
        plt.title(f'Valid Pixels Mask ({n_valid_pixels} pixels)')
        plt.colorbar()

        # Display a couple of DM modes
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.imshow(trans_dm_array_cpu[:, :, 0], cmap='seismic')
        plt.title('Transformed DM Mode 0')
        plt.colorbar()
        plt.subplot(1, 2, 2)
        plt.imshow(trans_dm_array_cpu[:, :, 1], cmap='seismic')
        plt.title('Transformed DM Mode 1')
        plt.colorbar()

        # Display a couple of Base modes
        plt.figure(figsize=(12, 5))
        if base_inv_array.ndim == 3:
            base_inv_array_cpu = cpuArray(base_inv_array)
            plt.subplot(1, 2, 1)
            plt.imshow(base_inv_array_cpu[:, :, 0], cmap='seismic')
            plt.title('Base Mode 0')
            plt.colorbar()
            plt.subplot(1, 2, 2)
            plt.imshow(base_inv_array_cpu[:, :, 1], cmap='seismic')
            plt.title('Base Mode 1')
            plt.colorbar()
        else:
            plt.text(0.1, 0.5, 'Base modes not in 3D format', fontsize=12)
            plt.axis('off')

        # Display projection coefficients
        plt.figure(figsize=(10, 6))
        x = xp.arange(projection.shape[0])+1
        for i in range(min(5, projection.shape[1])):
            plt.plot(x, projection_cpu[:, i], label=f'Basis mode {i}', marker='o', markersize=3)
        plt.xscale('log')
        plt.legend()
        plt.title('Projection Coefficients')
        plt.xlabel('DM mode index')
        plt.ylabel('Coefficient')
        plt.grid(True)

        # Display projection matrix
        plt.figure(figsize=(10, 8))
        plt.imshow(projection_cpu, cmap='seismic', origin='lower', aspect='auto')
        plt.title(f'Projection Matrix ({projection.shape[0]} × {projection.shape[1]})')
        plt.xlabel('Basis mode index')
        plt.ylabel('DM mode index')
        plt.colorbar()
        plt.grid(True, alpha=0.3)
        plt.show()

    return projection
