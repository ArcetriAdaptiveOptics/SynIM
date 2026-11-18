"""Utility functions for SynIM simulations, including mask extrapolation"""

import numpy as np
import matplotlib.pyplot as plt
import synim as _synim
_synim._require_init(__name__)
from synim import (
    xp, cpuArray, to_xp, float_dtype, affine_transform, rotate, shift, zoom
)
from scipy.ndimage import binary_dilation
import scipy.sparse

# Labels for the extrapolation directions
directions_labels = ['Down (y+1)', 'Up (y-1)', 'Right (x+1)', 'Left (x-1)']

def calculate_extrapolation_indices_coeffs(mask, debug=False, debug_pixels=None):
    """
    Calculates indices and coefficients for extrapolating edge pixels of a mask.

    Parameters:
        mask (ndarray): Binary mask (True/1 inside, False/0 outside).
        debug (bool): If True, displays debug information and plots.
        debug_pixels (list): List of [y, x] coordinates for detailed debug output.

    Returns:
        tuple: (edge_pixels, reference_indices, coefficients)
            - edge_pixels: Linear indices of the edge pixels to extrapolate.
            - reference_indices: Array of reference pixel indices for extrapolation.
            - coefficients: Coefficients for linear extrapolation.
    """

    # *** MODIFIED: Ensure mask is CPU numpy ***
    mask = cpuArray(mask)

    # Convert the mask to boolean
    binary_mask = mask.astype(bool)

    # Identify edge pixels (outside but adjacent to the mask) using binary dilation
    dilated_mask = binary_dilation(binary_mask)
    edge_pixels = np.where(dilated_mask & ~binary_mask)
    edge_pixels_linear = np.ravel_multi_index(edge_pixels, mask.shape)

    if debug:
        print(f"Found {len(edge_pixels[0])} edge pixels to extrapolate.")

        # Plot the original mask and the edge pixels
        plt.figure(figsize=(10, 4))
        plt.subplot(121)
        plt.imshow(binary_mask, cmap='gray', interpolation='nearest')
        plt.title('Original Mask')

        plt.subplot(122)
        edge_mask = np.zeros_like(binary_mask)
        edge_mask[edge_pixels] = 1
        plt.imshow(binary_mask, cmap='gray', alpha=0.5, interpolation='nearest')
        plt.imshow(edge_mask, cmap='hot', alpha=0.5, interpolation='nearest')
        plt.title('Edge Pixels to Extrapolate (red)')
        plt.tight_layout()
        plt.show()

    # Preallocate arrays for reference indices and coefficients
    reference_indices = np.full((len(edge_pixels[0]), 8), -1, dtype=np.int32)
    coefficients = np.zeros((len(edge_pixels[0]), 8), dtype=np.float32)

    # Directions for extrapolation (y+1, y-1, x+1, x-1)
    directions = [
        (1, 0),  # y+1 (down)
        (-1, 0), # y-1 (up)
        (0, 1),  # x+1 (right)
        (0, -1)  # x-1 (left)
    ]

    # Iterate over each edge pixel
    for i, (y, x) in enumerate(zip(*edge_pixels)):
        # Check if this pixel is in the debug list
        is_debug_pixel = debug_pixels is not None and any(
            p[0] == y and p[1] == x for p in debug_pixels)
        if is_debug_pixel:
            print(f"\n[DEBUG] Detailed analysis for pixel [{y},{x}]:")

        # Linear indices of the reference points at distance 1 and 2 in each
        # direction, -1 where they are outside the image or the mask
        refs = []
        for dir_idx, (dy, dx) in enumerate(directions):
            y1, x1 = y + dy, x + dx
            y2, x2 = y + 2*dy, x + 2*dx
            valid_ref1 = (0 <= y1 < mask.shape[0] and
                          0 <= x1 < mask.shape[1] and
                          binary_mask[y1, x1])
            valid_ref2 = (valid_ref1 and
                          0 <= y2 < mask.shape[0] and
                          0 <= x2 < mask.shape[1] and
                          binary_mask[y2, x2])
            refs.append((y1 * mask.shape[1] + x1 if valid_ref1 else -1,
                         y2 * mask.shape[1] + x2 if valid_ref2 else -1))
            if is_debug_pixel:
                print(f"  Direction {directions_labels[dir_idx]}: "
                      f"Ref1 [{y1},{x1}] valid: {valid_ref1}, Ref2 [{y2},{x2}] valid: {valid_ref2}")

        # Linear extrapolation 2*P1 - P2, averaged over the directions with both
        # reference points; the value of the first valid neighbour only if there
        # is none. Averaging the two kinds would bias the result (by half a step
        # for a linear phase), so they are never mixed.
        linear = [d for d, (_, ref2) in enumerate(refs) if ref2 >= 0]
        if linear:
            for d in linear:
                reference_indices[i, 2*d:2*d + 2] = refs[d]
                coefficients[i, 2*d] = 2.0 / len(linear)
                coefficients[i, 2*d + 1] = -1.0 / len(linear)
        else:
            d = next(d for d, (ref1, _) in enumerate(refs) if ref1 >= 0)
            reference_indices[i, 2*d] = refs[d][0]
            coefficients[i, 2*d] = 1.0

        if is_debug_pixel:
            print("  Coefficients:", coefficients[i])

    if debug:
        print(f"Average valid directions per pixel:"
              f" {np.sum(coefficients != 0) / (len(edge_pixels[0]) * 2):.2f}")

        # Display coefficient matrix for the first 10 pixels
        if len(edge_pixels[0]) >= 10:
            print("\nCoefficients for the first 10 pixels:")
            for i in range(min(10, len(edge_pixels[0]))):
                print(f"Pixel {i} ({edge_pixels[0][i]}, {edge_pixels[1][i]}): {coefficients[i]}")
                print(f"Indices: {reference_indices[i]}")

    return edge_pixels_linear, reference_indices, coefficients


def apply_extrapolation(data, edge_pixels, reference_indices, coefficients, in_place=False):
    """
    Applies linear extrapolation to edge pixels using precalculated indices and coefficients.

    Parameters:
        data (ndarray): Input array to extrapolate.
        edge_pixels (ndarray): Linear indices of edge pixels to extrapolate.
        reference_indices (ndarray): Indices of reference pixels.
        coefficients (ndarray): Coefficients for linear extrapolation.
        in_place (bool): If True, modifies the input data array directly.

    Returns:
        ndarray: Array with extrapolated pixels.
    """
    # Decide whether to work on a copy or the original
    if in_place:
        result = data
    else:
        result = data.copy()

    # Handle 2D vs 3D arrays
    if data.ndim == 2:
        result_reshaped = result[..., np.newaxis]
        n_slices = 1
    else:
        result_reshaped = result
        n_slices = result.shape[2]

    flat_result = result_reshaped.reshape(-1, n_slices)
    flat_data = data.reshape(-1, n_slices) if data.ndim == 3 else data.reshape(-1, 1)

    # The explicit pairwise sum below assumes 8 neighbours per edge pixel
    if coefficients.shape[1] != 8 or reference_indices.shape[1] != 8:
        raise ValueError(f"Expected 8 reference pixels per edge pixel, got"
                         f" {coefficients.shape[1]} coefficients and"
                         f" {reference_indices.shape[1]} indices")

    # Vectorized extrapolation for all slices at once
    valid_ref_mask = reference_indices >= 0
    safe_ref_indices = xp.where(valid_ref_mask, reference_indices, 0)
    masked_coeffs = xp.where(valid_ref_mask, coefficients, 0.0)

    def contribution(j):
        # (n_edge, n_slices): coefficient j times reference pixel j
        return masked_coeffs[:, j:j + 1] * flat_data[safe_ref_indices[:, j], :]

    # Sum of the 8 contributions in the same order as numpy.sum over an axis
    # of length 8 (pairwise summation), so that the result is bit-identical
    # to the previous per-slice implementation on CPU
    extrap_values = ((contribution(0) + contribution(1)) + (contribution(2) + contribution(3))) \
        + ((contribution(4) + contribution(5)) + (contribution(6) + contribution(7)))
    flat_result[edge_pixels, :] = extrap_values

    # If in_place, result already points to data, so no need to return differently
    # If 2D input, squeeze back
    if data.ndim == 2:
        if in_place:
            # result is already data, just return it
            return result
        else:
            # result is a copy, squeeze and return
            return result_reshaped[..., 0]

    return result


def shiftzoom_from_source_dm_params(source_pol_coo, source_height, dm_height, pixel_pitch):
    """
    Compute the shift and zoom parameters for a DM based on the source coordinates and heights.
    
    Parameters:
    - source_pol_coo: tuple, (radius, angle) in polar coordinates
    - source_height: float, height of the source
    - dm_height: float, height of the DM
    - pixel_pitch: float, pixel pitch in meters

    Returns:
    - shift: tuple, (x_shift, y_shift) in pixels
    - zoom: tuple, (x_zoom, y_zoom) magnification factors
    """

    arcsec2rad = xp.pi/180/3600

    if xp.isinf(source_height):
        mag_factor = 1.0
    else:
        mag_factor = source_height/(source_height-dm_height)
    source_rec_coo_asec = polar_to_xy(source_pol_coo[0],source_pol_coo[1]*xp.pi/180)
    source_rec_coo_m = source_rec_coo_asec*dm_height*arcsec2rad
    # change sign to get the shift in the right direction considering
    # the convention applied in rotshiftzoom_array
    source_rec_coo_pix = -1 * source_rec_coo_m / pixel_pitch

    shift = tuple(source_rec_coo_pix)
    zoom = (mag_factor, mag_factor)

    return shift, zoom


def rotshiftzoom_matrix(input_shape, output_size, dm_translation=(0.0, 0.0),
                        dm_rotation=0.0, dm_magnification=(1.0, 1.0),
                        wfs_translation=(0.0, 0.0), wfs_rotation=0.0,
                        wfs_magnification=(1.0, 1.0),
                        wfs_anamorphosis_45=1.0):
    """
    2x2 matrix and 2-element offset of the affine transformation applied by
    rotshiftzoom_array (inverse mapping, as in affine_transform: the input
    coordinates of an output pixel o are matrix @ o + offset) to an array
    whose first two dimensions are input_shape, with an output of size
    output_size. The parameters are the same as in rotshiftzoom_array.
    """

    # Parameter handling: conversion of single values to tuples
    try:
        if not hasattr(dm_translation, '__len__') or len(dm_translation) != 2:
            dm_translation = (float(dm_translation), float(dm_translation))
    except (TypeError, ValueError):
        dm_translation = (0.0, 0.0)

    try:
        if not hasattr(wfs_translation, '__len__') or len(wfs_translation) != 2:
            wfs_translation = (float(wfs_translation), float(wfs_translation))
    except (TypeError, ValueError):
        wfs_translation = (0.0, 0.0)

    try:
        if not hasattr(dm_magnification, '__len__'):
            # If it is a single value, we create a tuple with two identical elements
            dm_magnification = (float(dm_magnification), float(dm_magnification))
        elif len(dm_magnification) != 2:
            # if it is a sequence but not of length 2
            dm_magnification = (float(dm_magnification[0]), float(dm_magnification[0]))
    except (TypeError, ValueError):
        dm_magnification = (1.0, 1.0)

    try:
        if not hasattr(wfs_magnification, '__len__'):
            # If it is a single value, we create a tuple with two identical elements
            wfs_magnification = (float(wfs_magnification), float(wfs_magnification))
        elif len(wfs_magnification) != 2:
            # if it is a sequence but not of length 2
            wfs_magnification = (float(wfs_magnification[0]), float(wfs_magnification[0]))
    except (TypeError, ValueError):
        wfs_magnification = (1.0, 1.0)

    # Center of the input array: pixel i is centred at i, so the centre of
    # N pixels is (N - 1) / 2 (the same convention as SPECULA AtmoPropagation)
    center = (xp.array(input_shape[:2]) - 1) / 2.0
    # Convert rotations to radians
    # Note: Inverting the sign of rotation to match the first function's direction
    dm_rot_rad = xp.deg2rad(-dm_rotation)  # Negative sign to reverse direction
    wfs_rot_rad = xp.deg2rad(-wfs_rotation)  # Negative sign to reverse direction
    # Create the transformation matrices
    # For DM transformation
    dm_scale_matrix = xp.array(
        [[1.0/dm_magnification[0], 0], [0, 1.0/dm_magnification[1]]]
    )
    dm_rot_matrix = xp.array(
        [[xp.cos(dm_rot_rad), -xp.sin(dm_rot_rad)], [xp.sin(dm_rot_rad), xp.cos(dm_rot_rad)]]
    )
    dm_matrix = xp.dot(dm_rot_matrix, dm_scale_matrix)

    # For WFS transformation
    wfs_scale_matrix = xp.array(
        [[1.0/wfs_magnification[0], 0], [0, 1.0/wfs_magnification[1]]]
    )
    wfs_rot_matrix = xp.array(
        [[xp.cos(wfs_rot_rad), -xp.sin(wfs_rot_rad)],
         [xp.sin(wfs_rot_rad), xp.cos(wfs_rot_rad)]]
    )

    # OPTICAL ORDER: The Shack-Hartmann sensor rotation is the FINAL step 
    # of the physical optical path. Therefore, in Inverse Mapping, we MUST un-rotate FIRST.
    # xp.dot(A, B) applies B first. So wfs_rot_matrix MUST be on the right!
    if wfs_anamorphosis_45 != 1.0:
        k = wfs_anamorphosis_45
        anam_45_matrix = xp.array([
            [(1.0 + k) / 2.0, (1.0 - k) / 2.0],
            [(1.0 - k) / 2.0, (1.0 + k) / 2.0]
        ])
        wfs_matrix = xp.dot(anam_45_matrix, xp.dot(wfs_scale_matrix, wfs_rot_matrix))
    else:
        wfs_matrix = xp.dot(wfs_scale_matrix, wfs_rot_matrix)

    # Combine transformations
    # INVERSE MAPPING ALGEBRA: To match a 2-step process (DM then WFS),
    # the inverse matrix must be M_DM_inv * M_WFS_inv.
    combined_matrix = xp.dot(dm_matrix, wfs_matrix)

    # Calculate offset
    output_center = (xp.array(output_size) - 1) / 2.0

    # The mathematical proof of the 2-step equivalence requires:
    # 1. WFS translation is multiplied by the FULL combined matrix.
    # 2. DM translation is multiplied ONLY by the DM rotation matrix.
    rotated_dm_translation = xp.dot(dm_rot_matrix, xp.array(dm_translation)[:2])
    scaled_wfs_translation = xp.dot(combined_matrix, xp.array(wfs_translation)[:2])

    offset = center - xp.dot(combined_matrix, output_center) \
        - rotated_dm_translation - scaled_wfs_translation

    return combined_matrix, offset


def rotshiftzoom_array(input_array, dm_translation=(0.0, 0.0),
                       dm_rotation=0.0, dm_magnification=(1.0, 1.0),
                       wfs_translation=(0.0, 0.0), wfs_rotation=0.0,
                       wfs_magnification=(1.0, 1.0),
                       wfs_anamorphosis_45=1.0,
                       output_size=None):
    """
    This function applies magnification, rotation, shift and resize of a
    2D or 3D numpy/cupy array using affine transformation.
    Rotation is applied in the same direction as the first function.

    Parameters:
    - input_array: numpy/cupy array, input data to be transformed
    - dm_translation: tuple, translation for DM (x, y)
    - dm_rotation: float, rotation angle for DM in degrees
    - dm_magnification: tuple, magnification factors for DM (x, y)
    - wfs_translation: tuple, translation for WFS (x, y)
    - wfs_rotation: float, rotation angle for WFS in degrees
    - wfs_magnification: tuple, magnification factors for WFS (x, y)
    - wfs_anamorphosis_45 (float): Diagonal anamorphosis factor.
      Values > 1 stretch along +45° diagonal and compress along -45° diagonal.
      Implemented as a shear transformation.
    - output_size: tuple, desired output size (height, width)

    Returns:
    - output: numpy/cupy array, transformed data
    """

    if xp.isnan(input_array).any():
        input_array = xp.nan_to_num(input_array, copy=True, nan=0.0, posinf=None, neginf=None)

    # Check if array is 2D or 3D
    is_3d = len(input_array.shape) == 3

    # resize
    if output_size is None:
        output_size = input_array.shape[:2]  # Only take the first two dimensions

    combined_matrix, offset = rotshiftzoom_matrix(
        input_array.shape[:2], output_size,
        dm_translation=dm_translation, dm_rotation=dm_rotation,
        dm_magnification=dm_magnification, wfs_translation=wfs_translation,
        wfs_rotation=wfs_rotation, wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45)

    # For 3D arrays, extend the transformation matrix to 3x3
    if is_3d:
        combined_matrix_3d = xp.eye(3)
        combined_matrix_3d[:2, :2] = combined_matrix
        combined_matrix = combined_matrix_3d
        offset_3d = xp.zeros(3, dtype=offset.dtype)
        offset_3d[:2] = offset
        offset = offset_3d

    # Apply transformation (scipy requires numpy)
    output = affine_transform(
        input_array,
        combined_matrix,
        offset=offset,
        output_shape=output_size if not is_3d else output_size + (input_array.shape[2],),
        order=1
    )

    return output


def affine_interpolation_matrix(input_shape, output_shape, matrix, offset, transposed=False):
    """
    Sparse (n_out, n_in) matrix of the bilinear interpolation done by
    affine_transform(input, matrix, offset, output_shape, order=1) (constant
    mode, cval=0) on a 2D array of shape input_shape: output.ravel() is
    the matrix times input.ravel(). Applied to the (n_pixels, n_modes)
    reshape of a cube, it transforms all the modes at once.

    As in scipy.ndimage, an output pixel whose input coordinates are outside
    [0, n - 1] on either axis is 0.

    transposed: if True, the input is the transpose of an array of shape
        input_shape[::-1] and the column indices refer to that array, so that
        the matrix applies to the original (not transposed) array.

    Returns a float64 scipy.sparse.csr_matrix (computed on CPU).
    """
    matrix = np.asarray(cpuArray(matrix), dtype=np.float64)
    offset = np.asarray(cpuArray(offset), dtype=np.float64)
    n0, n1 = input_shape[:2]
    if n0 < 2 or n1 < 2:
        raise ValueError(f'The input must have at least 2x2 pixels, got {input_shape}')

    out0, out1 = np.meshgrid(np.arange(output_shape[0], dtype=np.float64),
                             np.arange(output_shape[1], dtype=np.float64), indexing='ij')
    c0 = (matrix[0, 0] * out0 + matrix[0, 1] * out1 + offset[0]).ravel()
    c1 = (matrix[1, 0] * out0 + matrix[1, 1] * out1 + offset[1]).ravel()
    inside = (c0 >= 0) & (c0 <= n0 - 1) & (c1 >= 0) & (c1 <= n1 - 1)
    c0 = c0[inside]
    c1 = c1[inside]

    # Lower neighbour; at the last pixel (c = n - 1) the lower neighbour is
    # n - 2 with fraction 1, so that the upper neighbour is inside the array
    i0 = np.minimum(np.floor(c0).astype(np.intp), n0 - 2)
    i1 = np.minimum(np.floor(c1).astype(np.intp), n1 - 2)
    f0 = c0 - i0
    f1 = c1 - i1

    def column(a, b):
        return b * n0 + a if transposed else a * n1 + b

    # CSR arrays built directly: 4 neighbours for each output pixel inside the input
    cols = np.stack((column(i0, i1), column(i0, i1 + 1),
                     column(i0 + 1, i1), column(i0 + 1, i1 + 1)), axis=1)
    values = np.stack(((1 - f0) * (1 - f1), (1 - f0) * f1, f0 * (1 - f1), f0 * f1), axis=1)
    indptr = np.concatenate(([0], np.cumsum(4 * inside)))
    operator = scipy.sparse.csr_matrix(
        (values.ravel(), cols.ravel(), indptr),
        shape=(output_shape[0] * output_shape[1], n0 * n1))
    operator.eliminate_zeros()
    return operator


def dm3d_to_2d(dm_array, mask,
               xp_local=xp, float_dtype_local=float_dtype, block_size=256):
    """
    Convert a 3D DM influence function (n, n, n_modes) to a 2D array
    (n_modes, n_valid_pixels) using a mask.

    Each mode is normalized to unit RMS and then its mean (piston) is removed.

    The modes are processed in blocks of block_size: the temporary arrays
    have the size of one block, not of the whole dm_array.
    """

    # *** MODIFIED: Convert inputs to xp with correct dtype ***
    dm_array = to_xp(xp_local, dm_array, dtype=float_dtype_local)
    mask = to_xp(xp_local, mask, dtype=float_dtype_local)

    # Check if the mask is 2D
    if mask.ndim != 2:
        raise ValueError("The mask must be a 2D array.")
    # Check if the dm_array is 3D
    if dm_array.ndim != 3:
        raise ValueError("The dm_array must be a 3D array.")
    nmodes = dm_array.shape[2]
    # Valid pixels as rows of the (n * n, n_modes) view of the 3D array:
    # each one is a contiguous row, faster to gather than with 2D indices
    flat_idx = xp_local.flatnonzero(mask > 0)
    dm_rows = dm_array.reshape(-1, nmodes)
    dm_array_2d = xp_local.empty((nmodes, len(flat_idx)), dtype=float_dtype_local)
    block_size = max(1, int(block_size))
    pix_block = block_size * 64
    for start in range(0, len(flat_idx), pix_block):
        stop = min(start + pix_block, len(flat_idx))
        dm_array_2d[:, start:stop] = dm_rows[flat_idx[start:stop]].T
    normalize_modes(dm_array_2d, xp_local=xp_local, block_size=block_size)

    return dm_array_2d


def normalize_modes(modes_2d, xp_local=xp, block_size=256):
    """
    Normalize in place each mode (row) of a 2D array (n_modes, n_pixels) to
    unit RMS, then remove its mean (piston), as dm2d_to_3d and dm3d_to_2d do.

    The modes are processed in blocks of block_size rows. Returns modes_2d.
    """
    block_size = max(1, int(block_size))
    for start in range(0, modes_2d.shape[0], block_size):
        block = modes_2d[start:start + block_size]
        block /= xp_local.sqrt(xp_local.mean(block**2, axis=1, keepdims=True))
        block -= xp_local.mean(block, axis=1, keepdims=True)
    return modes_2d


def dm2d_to_3d(dm_array, mask, normalize=True,
               xp_local=xp, float_dtype_local=float_dtype, block_size=256):
    """
    Convert a 2D DM influence function (n_modes, n_valid_pixels) to a 3D
    array (n, n, n_modes) using a mask. The input array is not modified.

    If normalize is True, each mode is normalized to unit RMS and then its
    mean (piston) is removed.

    The modes are processed in blocks of block_size: the temporary arrays
    have the size of one block, not of the whole dm_array.
    """

    # *** Convert inputs to xp with correct dtype ***
    dm_array = to_xp(xp_local, dm_array, dtype=float_dtype_local)
    mask = to_xp(xp_local, mask, dtype=float_dtype_local)

    # Check if the mask is 2D
    if mask.ndim != 2:
        raise ValueError("The mask must be a 2D array.")
    # Check if the dm_array is 2D
    if dm_array.ndim != 2:
        raise ValueError("The dm_array must be a 2D array.")
    nmodes = dm_array.shape[0]
    idx = xp_local.where(mask > 0)
    dm_array_3d = xp_local.zeros((mask.shape[0], mask.shape[1], nmodes),
                                 dtype=float_dtype_local)
    block_size = max(1, int(block_size))
    for start in range(0, nmodes, block_size):
        stop = min(start + block_size, nmodes)
        # Copy of a block of modes: the caller's array is not modified
        block = dm_array[start:stop].copy()
        if normalize:
            # normalize by the RMS, then remove the mean (piston)
            block /= xp_local.sqrt(xp_local.mean(block**2, axis=1, keepdims=True))
            block -= xp_local.mean(block, axis=1, keepdims=True)
        dm_array_3d[idx[0], idx[1], start:stop] = block.T

    return dm_array_3d


def _to_xp_floating(array):
    """Convert to xp; keep floating point dtypes, convert others to float_dtype."""
    array = to_xp(xp, array)
    if not xp.issubdtype(array.dtype, xp.floating):
        array = array.astype(float_dtype)
    return array


def apply_mask(array, mask, norm=False, fill_value=None, in_place=False):
    """Apply a 2D or 3D mask to a 2D or 3D array.
    
    Parameters:
        array: Input array
        mask: Mask array
        norm: If True, normalize by mask
        fill_value: Value to fill masked regions
        in_place: If True, modify array in-place (only works if no dtype conversion needed)

    A floating point array keeps its dtype; other types (bool, int) are
    converted to float_dtype.
    """
    # *** Convert inputs to xp ***
    array = _to_xp_floating(array)
    mask = to_xp(xp, mask, dtype=float_dtype)

    # Broadcast mask for 3D arrays
    if array.ndim == 3 and mask.ndim == 2:
        norm_mask = mask[:, :, xp.newaxis]
    else:
        norm_mask = mask
    if norm:
        safe_mask = xp.where(norm_mask == 0, 1, norm_mask)
        norm_mask = 1.0 / safe_mask

    if fill_value is not None:
        result = xp.where(norm_mask, array, fill_value)
    else:
        result = array * norm_mask
        if norm and fill_value is None:
            # Set to 0 where mask was zero (to avoid inf)
            result = xp.where(
                mask if mask.ndim == array.ndim else mask[:, :, xp.newaxis],
                result, 0
            )

    if in_place:
        array[:] = result
        return array
    else:
        return result


def has_transformations(rotation, translation, magnification):
    """
    Helper function to check if there are any non-trivial transformations.
    
    Returns:
    - bool: True if there are transformations, False otherwise
    """
    # Check rotation
    has_rotation = rotation != 0.0

    # Check translation
    if hasattr(translation, '__len__'):
        has_translation = not all(t == 0.0 for t in translation)
    else:
        has_translation = translation != 0.0

    # Check magnification
    if hasattr(magnification, '__len__'):
        has_magnification = not all(m == 1.0 for m in magnification)
    else:
        has_magnification = magnification != 1.0

    return has_rotation or has_translation or has_magnification


def _nanmean_blocks_inplace(blocks):
    """
    numpy.nanmean(blocks, axis=(1, 3)) with the same operations as numpy
    (NaN replaced by 0, sum, count of the non-NaN values, division), but with
    the NaN values replaced in place instead of in a copy of the array.
    blocks is modified. Blocks without valid values give NaN, as in numpy.
    """
    nan_mask = np.isnan(blocks)
    np.copyto(blocks, 0, where=nan_mask)
    total = np.sum(blocks, axis=(1, 3))
    count = np.sum(~nan_mask, axis=(1, 3), dtype=np.intp)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.divide(total, count, out=total, casting='unsafe')


def rebin_matrix(m, M, dtype=None):
    """
    (M, m) matrix of the overlap between m pixels and M bins along one axis.

    Element [s, j] is the fraction of pixel j (interval [j, j + 1]) inside
    bin s (interval [s * m / M, (s + 1) * m / M]): each column sums to 1,
    each row to m / M. The boundaries are computed with integer arithmetic,
    so the overlaps are exact (0 or 1 when M divides m).
    """
    return to_xp(xp, _rebin_matrix_np(m, M), dtype=float_dtype if dtype is None else dtype)


def _rebin_matrix_np(m, M):
    """rebin_matrix(m, M) as a float64 numpy array."""
    s = np.arange(M)[:, np.newaxis]
    j = np.arange(m)[np.newaxis, :]
    low = np.maximum(s * m, j * M)
    high = np.minimum((s + 1) * m, (j + 1) * M)
    return np.clip(high - low, 0, None) / M


def rebin_fractional_sum(array, row_matrix, col_matrix):
    """
    Weighted sum on the first two axes of a 2D or 3D array:
    row_matrix @ array @ col_matrix.T for each slice (the weights are
    typically built with rebin_matrix). Returns a (M, N[, n_modes]) array.
    """
    # (M, n[, k]): the temporary array is M / m times the input
    rows = xp.tensordot(row_matrix, array, axes=(1, 0))
    # (N, M[, k]) -> (M, N[, k])
    out = xp.tensordot(col_matrix, rows, axes=(1, 1))
    return xp.ascontiguousarray(xp.swapaxes(out, 0, 1))


def _rebin_fractional(array, M, N, method, overwrite_input):
    """
    Compression of a floating point array by a non-integer factor: each
    pixel contributes to a bin with the fraction of its area inside the bin
    (sum), divided by the bin area (average) or by the area of the non-NaN
    pixels (nanmean, NaN where no valid pixel contributes).
    """
    m, n = array.shape[0:2]
    row_matrix = rebin_matrix(m, M, dtype=array.dtype)
    col_matrix = rebin_matrix(n, N, dtype=array.dtype)
    if method == 'sum':
        return rebin_fractional_sum(array, row_matrix, col_matrix)
    if method == 'average':
        out = rebin_fractional_sum(array, row_matrix, col_matrix)
        out /= (m / M) * (n / N)
        return out
    # nanmean
    invalid = xp.isnan(array)
    if not overwrite_input:
        array = array.copy()
    xp.copyto(array, 0, where=invalid)
    total = rebin_fractional_sum(array, row_matrix, col_matrix)
    weight = rebin_fractional_sum((~invalid).astype(array.dtype), row_matrix, col_matrix)
    del invalid
    out = total / xp.where(weight > 0, weight, 1)
    out[weight == 0] = xp.nan
    return out


def rebin(array, new_shape, method='average', overwrite_input=False):
    """
    Resize array to new dimensions.

    A floating point array keeps its dtype in both the expansion and the
    compression case. For other types (bool, int) the expansion returns the
    original dtype and the compression returns float_dtype.

    Parameters:
        array: 2D or 3D array (the first two axes are resized)
        new_shape: (M, N) output size of the first two axes
        method: 'sum', 'average' or 'nanmean' (compression case)
        overwrite_input: if True, the 'nanmean' compression may modify array
            (the NaN values are set to 0) instead of working on a copy.
            Same result, lower memory use.

    Compression by a non-integer factor (e.g. 480 pixels to 68 bins): each
    pixel contributes to a bin with the fraction of its area inside the bin
    (see rebin_matrix); 'average' divides by the bin area and 'nanmean' by
    the area of the non-NaN pixels.
    """

    # *** MODIFIED: Convert input to xp ***
    orig_dtype = array.dtype
    array = _to_xp_floating(array)

    if array.ndim == 1:
        array = array.reshape(array.shape[0], 1)

    shape = array.shape
    m, n = shape[0:2]
    M, N = new_shape

    if M > m or N > n:
        # Expansion case
        if M % m != 0 or N % n != 0:
            raise ValueError("New shape must be multiples of the input dimensions.")
        if array.ndim == 3:
            # *** MODIFIED: Use xp.tile ***
            rebinned_array = xp.tile(array, (M//m, N//n, 1))
        else:
            rebinned_array = xp.tile(array, (M//m, N//n))
        if orig_dtype != rebinned_array.dtype:
            rebinned_array = rebinned_array.astype(orig_dtype)
        return rebinned_array

    # Compression case
    if M == 0 or N == 0:
        raise ValueError("New shape dimensions must be greater than 0.")
    if method not in ('sum', 'average', 'nanmean'):
        raise ValueError(f"Unsupported method: {method}."
                         f" Use 'sum', 'average', or 'nanmean'.")

    # Non-integer ratio: exact area weighting (see rebin_matrix)
    if m % M != 0 or n % N != 0:
        return _rebin_fractional(array, M, N, method, overwrite_input)

    # (M, m//M, N, n//N[, n_modes]) view: the blocks are reduced on axes 1, 3
    blocks = array[:M*(m//M), :N*(n//N)].reshape((M, m//M, N, n//N) + array.shape[2:])

    if method == 'sum':
        rebinned_array = xp.sum(blocks, axis=(1, 3))
    elif method == 'average':
        rebinned_array = xp.mean(blocks, axis=(1, 3))
    elif xp.__name__ == 'cupy':
        # cupy.nanmean does not copy the array
        rebinned_array = xp.nanmean(blocks, axis=(1, 3))
    elif overwrite_input:
        rebinned_array = _nanmean_blocks_inplace(blocks)
    else:
        with xp.errstate(invalid='ignore'):
            rebinned_array = xp.nanmean(blocks, axis=(1, 3))

    return rebinned_array


def polar_to_xy(r,theta):
    # conversion polar to rectangular coordinates
    # theta is in rad
    return xp.array(( r * xp.cos(theta),r * xp.sin(theta) ))


def make_orto_modes(array):
    # return an othogonal 2D array

    size_array = xp.shape(array)

    if len(size_array) != 2:
        raise ValueError('Error in input data, the input array must have two dimensions.')

    if size_array[1] > size_array[0]:
        Q, R = xp.linalg.qr(array.T)
        Q = Q.T
    else:
        Q, R = xp.linalg.qr(array)

    return Q


def detect_spider_arms_azimuthal(
    pup_mask,
    n_bins=None,
    min_radial_fill=0.35,
    radial_range=None,
):
    """
    Detect radial spider arms from a 2D pupil mask using a 1D azimuthal profile
    in the mid-annulus (avoiding inner obstruction and outer boundary notches).

    Args:
        pup_mask (ndarray): 2D pupil mask (1 inside pupil, 0 on spider/outside).
        n_bins (int, optional): Number of azimuthal bins over ``[0, 360)``
            degrees. If ``None`` (default), chosen adaptively so each bin spans
            ~1 pixel of arc at the mid-annulus radius (detecting even 1-px thin
            spiders).
        min_radial_fill (float): Minimum radial obscuration fraction in the
            mid-annulus to classify an azimuthal bin as part of a spider arm.
        radial_range (tuple of float or None): Fractional radial interval
            ``(f_min, f_max)`` within ``[r_in, r_out]`` used to build the
            azimuthal profile, where ``r = r_in + f * (r_out - r_in)``. If
            ``None``, the full-transmission radial plateau is estimated
            automatically from the 1D radial profile of ``pup_mask``.

    Returns:
        tuple: ``(arms_geometry, center_yx, (r_in, r_out))``
            - ``arms_geometry``: List of dicts, sorted by ``azimuth_deg``, with:
                - ``'arm_index'`` (int): 1-based arm index.
                - ``'azimuth_deg'`` (float): Arm angle in degrees ``[0, 360)`` in
                  ``pup_mask`` coordinates ``atan2(row - cy, col - cx)``.
                - ``'u_yx'`` (tuple): Unit radial direction ``(uy, ux)`` in
                  ``(row, col)`` coordinates.
                - ``'n_yx'`` (tuple): Unit normal ``(ny, nx) = (-ux, uy)`` in
                  ``(row, col)`` coordinates.
                - ``'center_yx'`` (tuple): Arm centroid offset ``(y0, x0)``
                  relative to ``center_yx`` in pixels.
                - ``'half_width_px'`` (float): Estimated half-width in pixels.
                - ``'width_px'`` (float): Estimated full width in pixels.
            - ``center_yx``: ``(cy, cx)`` array center in pixels.
            - ``(r_in, r_out)``: Inner and outer pupil radii in pixels.
    """
    mask_bin = np.asarray(cpuArray(pup_mask)) > 0.5
    h, w = mask_bin.shape
    cy, cx = (h - 1) / 2.0, (w - 1) / 2.0
    if not np.any(mask_bin) or np.all(mask_bin):
        return [], (cy, cx), (0.0, 0.0)

    yy, xx = np.mgrid[:h, :w]
    dy, dx = yy - cy, xx - cx
    r = np.hypot(dx, dy)
    r_valid = r[mask_bin]
    r_in, r_out = float(r_valid.min()), float(r_valid.max())
    radial_span = r_out - r_in
    if radial_span < 4.0:
        return [], (cy, cx), (r_in, r_out)

    if radial_range is None:
        # Automatically find the radial plateau free of inner/outer boundary notches
        n_rbins = max(16, int(round(radial_span / 2.0)))
        rbins = np.linspace(r_in, r_out, n_rbins + 1)
        in_ring = (r >= r_in) & (r <= r_out)
        r_tot, _ = np.histogram(r[in_ring], bins=rbins)
        r_val, _ = np.histogram(r[mask_bin], bins=rbins)
        r_fill = r_val / np.maximum(1, r_tot)
        plateau = np.where(r_fill >= 0.85 * float(np.max(r_fill)))[0]
        if plateau.size >= 3:
            r_low = float(rbins[plateau[0] + 1])
            r_high = float(rbins[plateau[-1]])
        else:
            r_low = r_in + 0.15 * radial_span
            r_high = r_out - 0.15 * radial_span
        r_low = max(r_low, 0.08 * r_out)
    else:
        f_min, f_max = float(radial_range[0]), float(radial_range[1])
        if not (0.0 <= f_min < f_max <= 1.0):
            raise ValueError(
                f"radial_range must satisfy 0.0 <= f_min < f_max <= 1.0, got {radial_range}"
            )
        r_low = max(r_in + f_min * radial_span, f_min * r_out)
        r_high = r_in + f_max * radial_span
    if r_high <= r_low:
        return [], (cy, cx), (r_in, r_out)

    mid_ann = (r >= r_low) & (r <= r_high)
    spider_mid = mid_ann & (~mask_bin)
    if not np.any(spider_mid):
        return [], (cy, cx), (r_in, r_out)

    if n_bins is None:
        r_mid_ann = 0.5 * (r_low + r_high)
        n_bins = int(np.clip(round(2.0 * np.pi * r_mid_ann), 360, 2880))
    else:
        n_bins = int(n_bins)

    angles_all = np.degrees(np.arctan2(dy[mid_ann], dx[mid_ann])) % 360.0
    angles_sp = np.degrees(np.arctan2(dy[spider_mid], dx[spider_mid])) % 360.0

    bins = np.linspace(0.0, 360.0, n_bins + 1)
    tot_hist, _ = np.histogram(angles_all, bins=bins)
    sp_hist, _ = np.histogram(angles_sp, bins=bins)
    fill = sp_hist / np.maximum(1, tot_hist)

    # Bridge 1-bin discretization hops without diluting thin 1-px spider peaks:
    # sum of adjacent bin pair (for arms split across a bin boundary)
    fill_pair = fill + np.maximum(np.roll(fill, 1), np.roll(fill, -1))
    active = (fill >= float(min_radial_fill)) | (fill_pair >= 1.2 * float(min_radial_fill))
    # Close single-bin holes inside an arm peak
    active = active | (np.roll(active, 1) & np.roll(active, -1))
    if not np.any(active) or np.all(active):
        return [], (cy, cx), (r_in, r_out)

    first_inactive = int(np.argmin(active))
    active_rot = np.roll(active, -first_inactive)
    padded = np.concatenate([[False], active_rot, [False]])
    starts = np.where(padded[1:] & ~padded[:-1])[0]
    ends = np.where(~padded[1:] & padded[:-1])[0]

    bin_step = 360.0 / float(n_bins)
    margin_deg = max(1.0, 2.0 * bin_step)
    sp_y, sp_x = dy[spider_mid], dx[spider_mid]

    arms = []
    for s, e in zip(starts, ends):
        ang_start = ((s + first_inactive) * bin_step - margin_deg) % 360.0
        ang_end = ((e + first_inactive) * bin_step + margin_deg) % 360.0
        in_sector = ((angles_sp - ang_start) % 360.0) <= ((ang_end - ang_start) % 360.0)
        n_pix_sector = int(np.count_nonzero(in_sector))
        if n_pix_sector < 3:
            continue

        yk, xk = sp_y[in_sector], sp_x[in_sector]
        rk = np.hypot(xk, yk)
        if float(rk.max() - rk.min()) < 0.5 * (r_high - r_low):
            continue

        r_mid = 0.5 * (rk.min() + rk.max())
        if not np.any(rk >= r_mid) or not np.any(rk < r_mid):
            continue

        # Translation-invariant radial unit vector from inner half to outer half
        uy = float(yk[rk >= r_mid].mean() - yk[rk < r_mid].mean())
        ux = float(xk[rk >= r_mid].mean() - xk[rk < r_mid].mean())
        norm_u = np.hypot(ux, uy)
        if norm_u <= 0:
            continue
        uy, ux = uy / norm_u, ux / norm_u
        ny, nx = -ux, uy
        y0, x0 = float(yk.mean()), float(xk.mean())

        width_px = max(0.5, float(n_pix_sector) / (r_high - r_low))
        half_width = 0.5 * width_px
        ang_deg = float(np.degrees(np.arctan2(uy, ux)) % 360.0)

        arms.append({
            'azimuth_deg': ang_deg,
            'u_yx': (uy, ux),
            'n_yx': (ny, nx),
            'center_yx': (y0, x0),
            'half_width_px': half_width,
            'width_px': width_px,
        })

    arms.sort(key=lambda a: a['azimuth_deg'])
    for idx, arm in enumerate(arms, start=1):
        arm['arm_index'] = idx
    return arms, (cy, cx), (r_in, r_out)


def _resolve_spider_widths(spider_width_px, n_arms, default=3.0):
    """Expand scalar/None spider width or validate per-arm width sequence."""
    if spider_width_px is None:
        return [float(default)] * n_arms
    if hasattr(spider_width_px, '__len__'):
        widths = [float(w) for w in spider_width_px]
        if len(widths) != n_arms:
            raise ValueError(
                f"spider_width_px length ({len(widths)}) must match number of arms ({n_arms})"
            )
        return widths
    return [float(spider_width_px)] * n_arms


def compute_spider_subaperture_weights_and_normals(
    pup_mask,
    wfs_nsubaps,
    wfs_rotation=0.0,
    wfs_translation=(0.0, 0.0),
    wfs_magnification=(1.0, 1.0),
    wfs_anamorphosis_90=1.0,
    wfs_anamorphosis_45=1.0,
    idx_valid_sa=None,
    dilation_px=0.0,
    smooth_sigma_px=0.0,
    spider_angles_deg=None,
    spider_width_px=None,
    n_bins=None,
    min_radial_fill=0.35,
    radial_range=None,
    arms_geometry=None,
    normalize_per_arm=False,
    specula_convention=True,
):
    """
    Compute per-subaperture weights and unit normal vectors for each spider arm
    in the WFS reference frame using the 1D azimuthal/angular spider model and
    the analytical inverse affine matrix ``rotshiftzoom_matrix``.

    The unit normal ``normal = (n_0, n_1)`` for each arm is returned in the exact
    slope block ordering of SynIM (where block 0 is the first half of the WFS
    slopes and block 1 is the second half):
    - With ``specula_convention=True``, ``pup_mask`` is transposed before the WFS
      transformation and SynIM concatenates ``(signal_y, signal_x)`` where
      ``signal_y`` is the derivative along ``axis=0`` (row) and ``signal_x`` is
      the derivative along ``axis=1`` (col) of the transformed array. Therefore
      ``n_0 = n_row`` and ``n_1 = n_col``.
    - With ``specula_convention=False``, ``n_0 = n_col`` and ``n_1 = n_row``.

    Args:
        pup_mask (ndarray): 2D pupil mask (n_pup, n_pup).
        wfs_nsubaps (int): Number of subapertures across the pupil diameter.
        wfs_rotation (float): WFS rotation angle in degrees.
        wfs_translation (tuple): WFS translation ``(dx, dy)`` in pixels.
        wfs_magnification (tuple or float): WFS magnification ``(mx, my)`` or scalar.
        wfs_anamorphosis_90 (float): WFS 90-deg anamorphosis factor (applied if
            ``wfs_magnification`` is scalar).
        wfs_anamorphosis_45 (float): WFS 45-deg anamorphosis factor.
        idx_valid_sa (ndarray, optional): Valid subaperture indices.
        dilation_px (float): Flat plateau half-width in pupil pixels around each
            spider arm before Gaussian roll-off (0 = strictly the spider bar).
        smooth_sigma_px (float): Gaussian roll-off standard deviation in pupil
            pixels beyond ``dilation_px`` (0 = sharp step on the spider bar).
        spider_angles_deg (sequence of float, optional): Explicit radial arm
            angles in degrees ``[0, 360)`` in the ``pup_mask`` frame
            (``atan2(row - cy, col - cx)``). If ``None``, angles and widths are
            detected automatically from the 1D azimuthal profile via
            ``detect_spider_arms_azimuthal``.
        spider_width_px (float or sequence of float, optional): Explicit spider
            arm width(s) in pixels (scalar or one per arm).
        n_bins (int, optional): Number of azimuthal bins for automatic detection.
        min_radial_fill (float): Minimum radial fill for automatic detection.
        radial_range (tuple of float or None): Fractional radial interval
            ``(f_min, f_max)`` within ``[r_in, r_out]`` for automatic detection
            (or ``None`` for automatic radial plateau estimation).
        arms_geometry (tuple, optional): Pre-computed output
            ``(arms, center_yx, (r_in, r_out))`` from
            ``detect_spider_arms_azimuthal(pup_mask, ...)`` to avoid repeating
            detection across multiple WFSs.
        normalize_per_arm (bool): If True, normalize each arm's subaperture
            weights so that the maximum weight along the arm is 1.0. Default is
            False so that coarse WFSs (e.g. 1x1, 8x8) are not inflated to 1.0.
        specula_convention (bool): Whether to follow SPECULA axis conventions.

    Returns:
        list of dict: One dictionary per spider arm with keys:
            - ``'arm_index'`` (int): 1-based arm index.
            - ``'normal'`` (np.ndarray): Unit normal ``(n_0, n_1)`` in WFS slope
              ordering.
            - ``'angle_deg'`` (float): Line orientation angle in degrees in the
              WFS slope frame (orthogonal to ``normal``).
            - ``'azimuth_deg'`` (float): Radial arm azimuth in the ``pup_mask``
              frame.
            - ``'width_px'`` (float): Arm full width in pupil pixels.
            - ``'weights'`` (np.ndarray): 1D array of weights for valid
              subapertures, shape ``(n_valid_sa,)``.
            - ``'weights_2d'`` (np.ndarray): 2D array of weights on the
              ``(wfs_nsubaps, wfs_nsubaps)`` grid.
    """
    from synim.geometry import _valid_sa_indices

    mask_bin = np.asarray(cpuArray(pup_mask)) > 0.5
    if specula_convention:
        mask_in = np.transpose(mask_bin)
        wfs_translation_local = (-1.0 * float(wfs_translation[1]),
                                 -1.0 * float(wfs_translation[0]))
    else:
        mask_in = mask_bin
        wfs_translation_local = (float(wfs_translation[0]), float(wfs_translation[1]))

    # Detect (or construct) arms in the non-transposed pup_mask frame so
    # azimuth_deg and arm_index are always consistent across conventions, then
    # map (row, col) -> (col, row) when specula_convention=True.
    if spider_angles_deg is not None:
        h_orig, w_orig = mask_bin.shape
        cy_orig, cx_orig = (h_orig - 1) / 2.0, (w_orig - 1) / 2.0
        yy, xx = np.mgrid[:h_orig, :w_orig]
        r_valid = np.hypot(xx[mask_bin] - cx_orig, yy[mask_bin] - cy_orig)
        r_in = float(r_valid.min()) if r_valid.size else 0.0
        r_out = float(r_valid.max()) if r_valid.size else min(h_orig, w_orig) / 2.0

        angles_list = [float(a) for a in spider_angles_deg]
        widths_list = _resolve_spider_widths(spider_width_px, len(angles_list), default=3.0)

        arms_orig = []
        for idx_a, (ang, w_px) in enumerate(zip(angles_list, widths_list), start=1):
            rad = np.deg2rad(ang)
            uy_orig, ux_orig = float(np.sin(rad)), float(np.cos(rad))
            arms_orig.append({
                'arm_index': idx_a,
                'azimuth_deg': ang % 360.0,
                'u_yx': (uy_orig, ux_orig),
                'n_yx': (-ux_orig, uy_orig),
                'center_yx': (0.0, 0.0),
                'half_width_px': 0.5 * w_px,
                'width_px': w_px,
            })
    else:
        if arms_geometry is not None:
            arms_orig_raw, (cy_orig, cx_orig), (r_in, r_out) = arms_geometry
            arms_orig = [dict(a) for a in arms_orig_raw]
        else:
            arms_orig, (cy_orig, cx_orig), (r_in, r_out) = detect_spider_arms_azimuthal(
                mask_bin,
                n_bins=n_bins,
                min_radial_fill=min_radial_fill,
                radial_range=radial_range,
            )
        if spider_width_px is not None:
            widths_list = _resolve_spider_widths(spider_width_px, len(arms_orig))
            for arm, w_px in zip(arms_orig, widths_list):
                arm['width_px'] = w_px
                arm['half_width_px'] = 0.5 * w_px

    if not arms_orig:
        return []

    cy, cx = (cx_orig, cy_orig) if specula_convention else (cy_orig, cx_orig)
    arms_geom = []
    for arm in arms_orig:
        uy_o, ux_o = arm['u_yx']
        y0_o, x0_o = arm['center_yx']
        if specula_convention:
            uy, ux = ux_o, uy_o
            y0, x0 = x0_o, y0_o
        else:
            uy, ux = uy_o, ux_o
            y0, x0 = y0_o, x0_o
        arms_geom.append({
            'arm_index': arm['arm_index'],
            'azimuth_deg': arm['azimuth_deg'],
            'u_yx': (uy, ux),
            'n_yx': (-ux, uy),
            'center_yx': (y0, x0),
            'half_width_px': arm['half_width_px'],
            'width_px': arm['width_px'],
        })

    if not hasattr(wfs_magnification, '__len__'):
        wfs_mag_tuple = (
            float(wfs_magnification),
            float(wfs_magnification) * float(wfs_anamorphosis_90),
        )
    else:
        wfs_mag_tuple = (float(wfs_magnification[0]), float(wfs_magnification[1]))

    h, w = mask_in.shape
    matrix, offset = rotshiftzoom_matrix(
        (h, w),
        (h, w),
        dm_translation=(0.0, 0.0),
        dm_rotation=0.0,
        dm_magnification=(1.0, 1.0),
        wfs_translation=wfs_translation_local,
        wfs_rotation=float(wfs_rotation),
        wfs_magnification=wfs_mag_tuple,
        wfs_anamorphosis_45=float(wfs_anamorphosis_45),
    )
    A = np.asarray(cpuArray(matrix), dtype=np.float64)
    b_off = np.asarray(cpuArray(offset), dtype=np.float64)
    center_in = np.array([cy, cx], dtype=np.float64)
    dc = b_off - center_in

    # Input coordinates relative to pupil center: dx_in(r, c) = A @ [r, c] + dc
    rows = np.arange(h, dtype=np.float64)[:, None]
    cols = np.arange(w, dtype=np.float64)[None, :]
    in_y = A[0, 0] * rows + A[0, 1] * cols + dc[0]
    in_x = A[1, 0] * rows + A[1, 1] * cols + dc[1]
    r_in_sq = in_y * in_y + in_x * in_x
    in_annulus = (r_in_sq >= (r_in * r_in)) & (r_in_sq <= (r_out * r_out))

    idx = _valid_sa_indices(idx_valid_sa, wfs_nsubaps, specula_convention)

    arms_info = []
    for arm in arms_geom:
        u_in = np.array(arm['u_yx'], dtype=np.float64)
        n_in = np.array(arm['n_yx'], dtype=np.float64)
        p0_in = np.array(arm['center_yx'], dtype=np.float64)
        hw = float(arm['half_width_px'])

        # Analytical gradient of transverse coordinate in WFS (row, col)
        g_rc = A.T @ n_in
        g_norm = float(np.linalg.norm(g_rc))
        if g_norm <= 0:
            continue
        n_rc = g_rc / g_norm
        t_rc = A.T @ u_in

        if specula_convention:
            normal = np.array([n_rc[0], n_rc[1]], dtype=np.float64)
        else:
            normal = np.array([n_rc[1], n_rc[0]], dtype=np.float64)

        # Tangent strictly orthogonal to normal in WFS slope coordinates
        tangent = np.array([-normal[1], normal[0]], dtype=np.float64)
        angle_deg = float(np.degrees(np.arctan2(tangent[1], tangent[0])) % 180.0)

        # 1D outer sums for longitudinal and transverse coordinates (in pupil px)
        s_par = (t_rc[0] * rows + t_rc[1] * cols) + float(u_in @ dc)
        d_perp = np.abs(
            (g_rc[0] * rows + g_rc[1] * cols) + float(n_in @ (dc - p0_in))
        )
        # Small 1e-3 px tolerance avoids float tie-breaking at half-integer pixel centers
        d_strip = np.maximum(0.0, d_perp - (hw + 1e-3))
        valid_half_plane = in_annulus & (s_par >= 0.0)

        if smooth_sigma_px > 0:
            arm_profile = np.where(
                valid_half_plane,
                np.exp(
                    -0.5 * (np.maximum(0.0, d_strip - float(dilation_px)) / float(smooth_sigma_px)) ** 2
                ),
                0.0,
            )
            arm_profile[arm_profile < 1e-4] = 0.0
        else:
            arm_profile = (valid_half_plane & (d_strip <= float(dilation_px))).astype(np.float64)

        w_2d = cpuArray(rebin(arm_profile, (wfs_nsubaps, wfs_nsubaps), method='average'))
        w_flat = w_2d.flatten()
        w_sa = w_flat if idx is None else w_flat[idx]
        w_sa = np.asarray(w_sa, dtype=np.float64)

        if normalize_per_arm:
            max_w = float(np.max(w_sa)) if w_sa.size > 0 else 0.0
            if max_w > 0:
                w_sa = w_sa / max_w
                w_2d = w_2d / max_w

        arms_info.append({
            'arm_index': arm['arm_index'],
            'normal': normal,
            'angle_deg': angle_deg,
            'azimuth_deg': arm['azimuth_deg'],
            'width_px': arm['width_px'],
            'weights': w_sa,
            'weights_2d': w_2d,
        })

    return arms_info


def update_wfs_noise_cov_inv_for_spider(
    C_noise_inv_wfs,
    arms_info,
    alpha=10.0,
    sigma_noise2_base=None,
    weight_threshold=1e-3,
    in_place=False,
    check_block_diagonal=True,
):
    """
    Update a single WFS inverse noise covariance matrix with a directional
    spider/LWE penalty using the vectorized 2x2 Sherman-Morrison formula.

    Models an additional slope variance ``sigma_{j,k}^2 * n_k @ n_k.T`` along the
    unit normal ``n_k`` to spider arm ``k`` for each subaperture ``j``, where
    ``sigma_{j,k}^2 = alpha * w_{j,k} * sigma_noise2_base``.

    Preserves exact bit-for-bit symmetry and the 2x2 block-diagonal sparsity
    pattern (at most 2 non-zeros per row) required by ``_sparse_rows`` in
    ``compute_mmse_reconstructor``. Only subapertures touched by at least one
    arm (``w > weight_threshold``) are modified.

    Args:
        C_noise_inv_wfs (ndarray): Inverse noise covariance for one WFS, shape
            ``(2 * n_sub, 2 * n_sub)``.
        arms_info (list of dict): Output of
            ``compute_spider_subaperture_weights_and_normals``.
        alpha (float): Dimensionless penalty factor relative to
            ``sigma_noise2_base``. If ``alpha <= 0``, returns a copy (or the
            input if ``in_place=True``) unchanged.
        sigma_noise2_base (float, optional): Base noise variance in slope^2
            units. If ``None``, inferred automatically as ``1.0 / max(lambda_max)``
            across all 2x2 subaperture blocks.
        weight_threshold (float): Minimum subaperture weight to apply update.
        in_place (bool): If True, modify ``C_noise_inv_wfs`` in place.
        check_block_diagonal (bool): If True, verify that ``C_noise_inv_wfs``
            only has non-zeros on the 2x2 subaperture blocks ``(j, j + n_sub)``.

    Returns:
        ndarray: Updated inverse noise covariance matrix for the WFS.
    """
    out = C_noise_inv_wfs if in_place else C_noise_inv_wfs.copy()
    if alpha <= 0 or not arms_info:
        return out

    n_slopes = out.shape[0]
    if n_slopes % 2 != 0 or out.shape[1] != n_slopes:
        raise ValueError(f"Expected square (2*n_sub, 2*n_sub) matrix, got {out.shape}")

    n_sub = n_slopes // 2
    out_dtype = out.dtype

    idx_0 = np.arange(n_sub)
    idx_1 = idx_0 + n_sub

    a = np.asarray(cpuArray(out[idx_0, idx_0]), dtype=np.float64).copy()
    d = np.asarray(cpuArray(out[idx_1, idx_1]), dtype=np.float64).copy()
    b01 = np.asarray(cpuArray(out[idx_0, idx_1]), dtype=np.float64)
    b10 = np.asarray(cpuArray(out[idx_1, idx_0]), dtype=np.float64)

    if check_block_diagonal:
        # Check that no non-zero entries exist outside the 2x2 blocks (j, j + n_sub)
        module = np if isinstance(out, np.ndarray) else _synim.xp
        total_nnz = int(module.count_nonzero(out))
        block_nnz = (
            int(np.count_nonzero(a)) + int(np.count_nonzero(d)) +
            int(np.count_nonzero(b01)) + int(np.count_nonzero(b10))
        )
        if total_nnz != block_nnz:
            raise ValueError(
                "C_noise_inv_wfs has non-zero entries outside the 2x2 subaperture "
                "blocks (j, j + n_sub); Sherman-Morrison 2x2 update requires a "
                "2x2 block-diagonal precision matrix."
            )

    b = 0.5 * (b01 + b10)

    if sigma_noise2_base is None:
        if n_sub == 0:
            return out
        # Maximum eigenvalue of each 2x2 block M_j = [[a_j, b_j], [b_j, d_j]]
        lambda_max = 0.5 * ((a + d) + np.hypot(a - d, 2.0 * b))
        max_prec = float(np.max(lambda_max))
        if max_prec <= 0:
            return out
        sigma_noise2_base = 1.0 / max_prec
    else:
        sigma_noise2_base = float(sigma_noise2_base)

    modified = np.zeros(n_sub, dtype=bool)
    for arm in arms_info:
        n0, n1 = float(arm['normal'][0]), float(arm['normal'][1])
        w = np.asarray(arm['weights'], dtype=np.float64)
        if w.shape[0] != n_sub:
            raise ValueError(
                f"Arm weights length ({w.shape[0]}) does not match n_sub ({n_sub})"
            )

        active = w > weight_threshold
        if not np.any(active):
            continue

        modified |= active
        sigma2_add = float(alpha) * w[active] * sigma_noise2_base
        a_act = a[active]
        b_act = b[active]
        d_act = d[active]

        v0 = a_act * n0 + b_act * n1
        v1 = b_act * n0 + d_act * n1
        den = 1.0 + sigma2_add * (n0 * v0 + n1 * v1)
        scale = sigma2_add / den

        a[active] = a_act - scale * (v0 * v0)
        d[active] = d_act - scale * (v1 * v1)
        b[active] = b_act - scale * (v0 * v1)

    if not np.any(modified):
        return out

    m_idx0 = idx_0[modified]
    m_idx1 = idx_1[modified]
    a_cast = a[modified].astype(out_dtype)
    d_cast = d[modified].astype(out_dtype)
    b_cast = b[modified].astype(out_dtype)

    if isinstance(out, np.ndarray):
        out[m_idx0, m_idx0] = a_cast
        out[m_idx1, m_idx1] = d_cast
        out[m_idx0, m_idx1] = b_cast
        out[m_idx1, m_idx0] = b_cast
    else:
        m_idx0_dev = to_xp(_synim.xp, m_idx0)
        m_idx1_dev = to_xp(_synim.xp, m_idx1)
        out[m_idx0_dev, m_idx0_dev] = to_xp(_synim.xp, a_cast, dtype=out_dtype)
        out[m_idx1_dev, m_idx1_dev] = to_xp(_synim.xp, d_cast, dtype=out_dtype)
        b_dev = to_xp(_synim.xp, b_cast, dtype=out_dtype)
        out[m_idx0_dev, m_idx1_dev] = b_dev
        out[m_idx1_dev, m_idx0_dev] = b_dev

    return out


__all__ = [
    # Mask creation and manipulation
    'make_orto_modes',
    'apply_mask',

    # DM array conversions
    'dm3d_to_2d',
    'normalize_modes',
    'dm2d_to_3d',

    # Array transformations
    'rebin',
    'rebin_matrix',
    'rebin_fractional_sum',
    'rotshiftzoom_array',
    'shiftzoom_from_source_dm_params',
    'has_transformations',

    # Extrapolation functions
    'apply_extrapolation',
    'calculate_extrapolation_indices_coeffs',

    # Spider / LWE noise covariance functions
    'detect_spider_arms_azimuthal',
    'compute_spider_subaperture_weights_and_normals',
    'update_wfs_noise_cov_inv_for_spider',
]