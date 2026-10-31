import math
import warnings

import numpy as np
import scipy.sparse

import synim as _synim
_synim._require_init(__name__)
from synim import xp, cpuArray, to_xp, float_dtype
import matplotlib.pyplot as plt
from synim.utils import (
    apply_mask,
    rebin,
    rebin_matrix,
    rebin_fractional_sum,
    _rebin_matrix_np,
    affine_interpolation_matrix,
    rotshiftzoom_array,
    rotshiftzoom_matrix,
    shiftzoom_from_source_dm_params,
    apply_extrapolation,
    calculate_extrapolation_indices_coeffs
)


def _without_nan(array):
    """array with the NaN values set to 0 (a copy only if there are NaN values)."""
    if xp.isnan(array).any():
        return xp.nan_to_num(array, nan=0.0)
    return array


def _warn_if_nan(im, name='interaction matrix'):
    """Warn if an interaction matrix (numpy or backend array) contains NaN values."""
    module = np if isinstance(im, np.ndarray) else xp
    nan_values = module.isnan(im)
    n_nan = int(module.count_nonzero(nan_values))
    if n_nan:
        n_rows = int(module.count_nonzero(nan_values.any(axis=1)))
        warnings.warn(f'The {name} contains {n_nan} NaN values in {n_rows} of'
                      f' {im.shape[0]} slopes: check the pupil and DM masks and the'
                      f' valid subapertures (idx_valid_sa).', RuntimeWarning, stacklevel=3)


def _repair_interpolated_phase(data, mask, threshold=0.999999, in_place=False):
    """
    Vital repair for edge artifacts. 
    Restores pixels that collapsed towards 0 due to bilinear interpolation 
    by using a strict threshold logic to select only uncorrupted pixels.

    If in_place is True, data is modified and returned (no copy of the array).
    """
    if mask is not None:
        if mask.max() < threshold:
            raise ValueError(f'Mask max value is {mask.max()},'
                             f' expected binary mask with values 0 and 1.')

        # 1. Strict binarization: isolate only the pure "core" that escaped interpolation
        strict_mask = xp.where(mask >= threshold, 1, 0)

        # 2. Zero out the softened edges corrupted by interpolation
        if in_place:
            data_repaired = data
            data_repaired[strict_mask == 0] = 0
        else:
            data_repaired = apply_mask(data, strict_mask, fill_value=0)

        # 3. Calculate indices and extrapolate using ONLY the strict mask
        edge_pixels, reference_indices, coefficients = calculate_extrapolation_indices_coeffs(
            cpuArray(strict_mask), debug=False, debug_pixels=None
        )

        edge_pixels = to_xp(xp, edge_pixels, dtype=xp.int32)
        reference_indices = to_xp(xp, reference_indices, dtype=xp.int32)
        coefficients = to_xp(xp, coefficients, dtype=float_dtype)

        # 4. Apply extrapolation, overwriting the corrupted edges with restored values
        data_repaired = apply_extrapolation(
            data_repaired, edge_pixels, reference_indices, coefficients, in_place=True
        )

        return data_repaired, strict_mask
    return data, mask


def compute_telsum_with_extrapolation(data, mask=None, wfs_nsubaps=None, verbose=False,
                                      in_place=False, pup_mask=None):
    """
    Computes the raw telescoping sum (average phase difference per pixel) for each subaperture.
    Acts as the telescoping sum equivalent of computing continuous derivatives.

    The difference between two adjacent pixels is used when both are valid
    for the phase (mask >= 0.999999, the DM mask), and is 0 otherwise.
    The differences are averaged over the pixel pairs (baselines) inside
    each subaperture with both pixels in the pupil (pup_mask > 0), as the
    derivatives are averaged over the pupil pixels. Without pup_mask the
    average is over the pairs valid for the phase. When the number of
    pixels is not a multiple of wfs_nsubaps, see _telsum_fractional.

    If in_place is True, the edge repair may modify data instead of a copy.
    """
    if wfs_nsubaps is None:
        raise ValueError("wfs_nsubaps must be provided to compute telescoping sum.")

    pup_diam_pix = data.shape[0]
    if pup_diam_pix < 2 * wfs_nsubaps:
        raise ValueError(f'The telescoping sum needs at least 2 pixels per subaperture:'
                         f' {pup_diam_pix} pixels, {wfs_nsubaps} subapertures.')
    if pup_diam_pix % wfs_nsubaps != 0 or data.shape[1] % wfs_nsubaps != 0:
        if verbose:
            print(f"  * Telescoping Sum: fractional subapertures"
                  f" ({pup_diam_pix / wfs_nsubaps:.3f} pixels per subaperture)")
        return _telsum_fractional(data, mask, wfs_nsubaps, pup_mask=pup_mask)
    N = pup_diam_pix // wfs_nsubaps

    # Repair interpolated edges using strict logic
    if mask is not None:
        repaired_data, strict_mask = _repair_interpolated_phase(data, mask, in_place=in_place)
    else:
        repaired_data = data
        strict_mask = xp.ones((data.shape[0], data.shape[1]), dtype=bool)
    valid_mask = strict_mask > 0
    # Pixels over which the differences are averaged
    weight_mask = valid_mask if pup_mask is None else pup_mask > 0

    # Setup 4D arrays for the Telescoping Sum
    is_3d = repaired_data.ndim == 3
    shape_4d = (wfs_nsubaps, N, wfs_nsubaps, N) + ((1,) if is_3d else ())
    p = repaired_data.reshape((wfs_nsubaps, N, wfs_nsubaps, N) + repaired_data.shape[2:])
    m = valid_mask.reshape(shape_4d)
    w = weight_mask.reshape(shape_4d)

    # Fast 4D Telescoping Sum. The weights (number of valid baselines) are
    # converted to the data dtype, so that the result keeps the data precision.
    # The differences are zeroed in place where not valid (same values as
    # xp.where(valid, d, 0.0)) and x and y are computed one after the other,
    # so that only one array of differences is allocated at a time.
    dx = p[:, :, :, 1:] - p[:, :, :, :-1]
    weight_dx = w[:, :, :, 1:] & w[:, :, :, :-1]
    xp.copyto(dx, 0.0, where=~(m[:, :, :, 1:] & m[:, :, :, :-1] & weight_dx))
    sum_dx = xp.sum(dx, axis=(1, 3))
    weight_dx = xp.sum(weight_dx, axis=(1, 3)).astype(sum_dx.dtype)
    del dx

    dy = p[:, 1:, :, :] - p[:, :-1, :, :]
    weight_dy = w[:, 1:, :, :] & w[:, :-1, :, :]
    xp.copyto(dy, 0.0, where=~(m[:, 1:, :, :] & m[:, :-1, :, :] & weight_dy))
    sum_dy = xp.sum(dy, axis=(1, 3))
    weight_dy = xp.sum(weight_dy, axis=(1, 3)).astype(sum_dy.dtype)
    del dy

    # Normalize to get Delta Phi per valid baseline
    raw_telsum_x = xp.where(weight_dx > 0, (sum_dx / xp.where(weight_dx > 0, weight_dx, 1.0)), 0.0)
    raw_telsum_y = xp.where(weight_dy > 0, (sum_dy / xp.where(weight_dy > 0, weight_dy, 1.0)), 0.0)

    return raw_telsum_x, raw_telsum_y


def _telsum_fractional(data, mask, wfs_nsubaps, threshold=0.999999, pup_mask=None):
    """
    Telescoping sum when the number of pixels is not a multiple of
    wfs_nsubaps. A baseline (pair of adjacent valid pixels j, j + 1) belongs
    to a subaperture with weight min(f_j, f_j+1), where f is the fraction of
    the pixel area inside the subaperture along the difference axis (the
    pixel pairs split by a subaperture boundary contribute to both
    subapertures with the fraction shared by the two pixels), and with
    weight f along the other axis. With an integer ratio these weights
    select the baselines inside each subaperture, as in the block
    computation. No interpolation of the data is needed. The differences
    and the pixels over which they are averaged are selected as in
    compute_telsum_with_extrapolation (mask and pup_mask).
    """
    dtype = data.dtype
    if mask is not None:
        if mask.max() < threshold:
            raise ValueError(f'Mask max value is {mask.max()},'
                             f' expected binary mask with values 0 and 1.')
        valid = (mask >= threshold).astype(dtype)
    else:
        valid = xp.ones(data.shape[:2], dtype=dtype)
    weight_mask = valid if pup_mask is None else (pup_mask > 0).astype(dtype)
    is_3d = data.ndim == 3

    def weights(n):
        fraction = rebin_matrix(n, wfs_nsubaps, dtype=dtype)
        return fraction, xp.minimum(fraction[:, 1:], fraction[:, :-1])

    rows, rows_baselines = weights(data.shape[0])
    cols, cols_baselines = weights(data.shape[1])

    def mean_difference(axis):
        if axis == 1:
            diff = data[:, 1:] - data[:, :-1]
            used = weight_mask[:, 1:] * weight_mask[:, :-1]
            zero = (valid[:, 1:] * valid[:, :-1] * used) == 0
            row_matrix, col_matrix = rows, cols_baselines
        else:
            diff = data[1:] - data[:-1]
            used = weight_mask[1:] * weight_mask[:-1]
            zero = (valid[1:] * valid[:-1] * used) == 0
            row_matrix, col_matrix = rows_baselines, cols
        # differences zeroed where not valid (in place, NaN safe)
        xp.copyto(diff, 0.0, where=zero[:, :, xp.newaxis] if is_3d else zero)
        del zero
        total = rebin_fractional_sum(diff, row_matrix, col_matrix)
        del diff
        weight = rebin_fractional_sum(used, row_matrix, col_matrix)
        if is_3d:
            weight = weight[:, :, xp.newaxis]
        return xp.where(weight > 0, total / xp.where(weight > 0, weight, 1), 0.0)

    return mean_difference(1), mean_difference(0)


def _valid_sa_indices(idx_valid_sa, wfs_nsubaps, specula_convention):
    """
    Indices (numpy int32) of the valid subapertures in the row-major list of
    all the subapertures, or None to select all of them. idx_valid_sa holds
    either (row, column) pairs (transposed in the SPECULA convention) or
    linear indices.
    """
    if idx_valid_sa is None:
        return None
    idx_valid_sa = cpuArray(idx_valid_sa)
    if specula_convention and len(idx_valid_sa.shape) > 1 and idx_valid_sa.shape[1] == 2:
        sa_2d = np.zeros((wfs_nsubaps, wfs_nsubaps))
        sa_2d[idx_valid_sa[:, 0], idx_valid_sa[:, 1]] = 1
        idx_temp = np.where(np.transpose(sa_2d) > 0)
        idx_valid_sa_new = np.zeros_like(idx_valid_sa)
        idx_valid_sa_new[:, 0] = idx_temp[0]
        idx_valid_sa_new[:, 1] = idx_temp[1]
    else:
        idx_valid_sa_new = idx_valid_sa

    if len(idx_valid_sa_new.shape) > 1 and idx_valid_sa_new.shape[1] == 2:
        return (idx_valid_sa_new[:, 0] * wfs_nsubaps + idx_valid_sa_new[:, 1]).astype(np.int32)
    return idx_valid_sa_new.astype(np.int32)


def _compute_slopes_from_telsum(raw_telsum_x, raw_telsum_y, pup_mask, dm_mask,
                               wfs_nsubaps, wfs_fov_arcsec, pup_diam_m,
                               idx_valid_sa, verbose, specula_convention):
    """
    Formats the raw telescoping sum into the final 1D slopes array.
    Signature is fully symmetrical to _compute_slopes_from_derivatives.
    """
    is_3d = raw_telsum_x.ndim == 3

    # 1. Pixels per subaperture (the telescoping sum is a difference per pixel)
    pup_diam_pix = pup_mask.shape[0]
    N = pup_diam_pix / wfs_nsubaps

    # 2. Rebin masks to compute valid subapertures (identical to derivative method)
    pup_mask = _without_nan(pup_mask)
    dm_mask = _without_nan(dm_mask)

    pup_mask_sa = rebin(pup_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    dm_mask_sa = rebin(dm_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    combined_mask_sa = (dm_mask_sa > 0.0) & (pup_mask_sa > 0.0)

    # 3. Apply mask
    telsum_x = xp.where(combined_mask_sa[:, :, xp.newaxis] if is_3d \
              else combined_mask_sa, raw_telsum_x, 0.0)
    telsum_y = xp.where(combined_mask_sa[:, :, xp.newaxis] if is_3d \
              else combined_mask_sa, raw_telsum_y, 0.0)

    # 4. Reshape to 2D matrix
    wfs_signal_x_2d = telsum_x.reshape((-1, telsum_x.shape[2] if is_3d else 1))
    wfs_signal_y_2d = telsum_y.reshape((-1, telsum_y.shape[2] if is_3d else 1))

    # 5. Select valid subapertures using provided indices
    idx = _valid_sa_indices(idx_valid_sa, wfs_nsubaps, specula_convention)
    if idx is not None:
        wfs_signal_x_2d = wfs_signal_x_2d[to_xp(xp, idx), :]
        wfs_signal_y_2d = wfs_signal_y_2d[to_xp(xp, idx), :]

    # 6. Concatenate X and Y arrays
    if specula_convention:
        im = xp.concatenate((wfs_signal_y_2d, wfs_signal_x_2d))
    else:
        im = xp.concatenate((wfs_signal_x_2d, wfs_signal_y_2d))

    # 7. Convert to slope units and apply N scaling
    coeff = 1e-9 / (pup_diam_m / wfs_nsubaps) * 206265
    coeff *= 1 / (0.5 * wfs_fov_arcsec)
    coeff *= N

    im = im * coeff

    if verbose:
        print(f'  G-tilt slopes formatted, shape: {im.shape}')

    return im


def _gradient(data, axis):
    """
    Same result as xp.gradient(data, axis=axis, edge_order=1) with unit
    spacing (central differences inside, one-sided differences at the edges,
    same floating point operations), computed directly in the output array
    without temporary arrays of the size of data.
    """
    def sl(start, stop):
        index = [slice(None)] * data.ndim
        index[axis] = slice(start, stop)
        return tuple(index)

    out = xp.empty_like(data)
    xp.subtract(data[sl(2, None)], data[sl(None, -2)], out=out[sl(1, -1)])
    out[sl(1, -1)] /= 2.0
    xp.subtract(data[sl(1, 2)], data[sl(0, 1)], out=out[sl(0, 1)])
    xp.subtract(data[sl(-1, None)], data[sl(-2, -1)], out=out[sl(-1, None)])
    return out


def compute_derivatives_with_extrapolation(data, mask=None, in_place=False):
    """
    Compute x and y derivatives (as numpy.gradient) on a 2D or 3D numpy array
    if mask is present does an extrapolation to avoid issue at the edges
    
    Parameters:
    - data: numpy 3D array
    - mask: optional, numpy 2D array, mask
    - in_place: bool, if True the edge repair modifies data instead of a copy

    Returns:
    - dx: numpy 3D array, x derivative (NaN outside the mask)
    - dy: numpy 3D array, y derivative (NaN outside the mask)
    """

    if mask is not None:
        data, mask = _repair_interpolated_phase(data, mask, in_place=in_place)

    # Compute x derivative
    dx = _gradient(data, axis=1)

    # Compute y derivative
    dy = _gradient(data, axis=0)

    if mask is not None:
        # Gracefully handle both 2D and 3D arrays
        is_3d = dx.ndim == 3
        idx = xp.ravel(xp.array(xp.where(mask.flatten() == 0)))

        dx_2d = dx.reshape((-1, dx.shape[2] if is_3d else 1))
        dx_2d[idx, :] = xp.nan

        dy_2d = dy.reshape((-1, dy.shape[2] if is_3d else 1))
        dy_2d[idx, :] = xp.nan

        dx = dx_2d.reshape(dx.shape)
        dy = dy_2d.reshape(dy.shape)

    return dx, dy


def _compute_slopes_from_derivatives(derivatives_x, derivatives_y, pup_mask, dm_mask,
                                     wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa,
                                     verbose, specula_convention):
    """
    Common function to compute slopes from derivatives.

    The slope of a subaperture is the mean of the derivatives over the pupil
    pixels (pup_mask > 0) inside it. When the number of pixels is not a
    multiple of wfs_nsubaps, a pixel shared by two subapertures contributes
    to each with the fraction of its area inside it (see rebin).

    derivatives_x and derivatives_y are modified in place (NaN values and
    values outside the pupil set to 0), to avoid copies of the arrays. For a
    given pupil mask these operations are idempotent: calling the function
    again on the same derivatives gives the same slopes. pup_mask and
    dm_mask are not modified.
    """

    # Clean up masks (a copy only if they contain NaN)
    pup_mask = _without_nan(pup_mask)
    dm_mask = _without_nan(dm_mask)

    # Rebin masks to WFS resolution
    pup_mask_sa = rebin(pup_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    pup_mask_sa = pup_mask_sa / xp.max(pup_mask_sa) if xp.max(pup_mask_sa) > 0 else pup_mask_sa

    dm_mask_sa = rebin(dm_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    if xp.max(dm_mask_sa) <= 0:
        raise ValueError('DM mask is empty after rebinning.')
    dm_mask_sa = dm_mask_sa / xp.max(dm_mask_sa)

    # Clean derivatives (in place): NaN values (outside the DM mask) and
    # values outside the pupil set to 0
    outside_pupil = pup_mask == 0
    for derivatives in (derivatives_x, derivatives_y):
        xp.copyto(derivatives, 0.0, where=xp.isnan(derivatives))
        derivatives[outside_pupil] = 0.0

    # Mean over the pupil pixels: sum of the derivatives divided by the
    # number (area) of pupil pixels of each subaperture. With an integer
    # ratio this is the same as a nanmean with NaN outside the pupil.
    weight = rebin((~outside_pupil).astype(derivatives_x.dtype),
                   (wfs_nsubaps, wfs_nsubaps), method='sum')
    weight = xp.where(weight > 0, weight, 1)
    if derivatives_x.ndim == 3:
        weight = weight[:, :, xp.newaxis]
    scale_factor = derivatives_x.shape[0] / wfs_nsubaps

    wfs_signal_x = rebin(derivatives_x, (wfs_nsubaps, wfs_nsubaps), method='sum') / weight \
        * scale_factor
    wfs_signal_y = rebin(derivatives_y, (wfs_nsubaps, wfs_nsubaps), method='sum') / weight \
        * scale_factor

    # Combined mask
    combined_mask_sa = (dm_mask_sa > 0.0) & (pup_mask_sa > 0.0)

    # Apply mask
    wfs_signal_x = apply_mask(wfs_signal_x, combined_mask_sa, fill_value=0)
    wfs_signal_y = apply_mask(wfs_signal_y, combined_mask_sa, fill_value=0)

    # Check if data is 3D (modes) or 2D (single phase screen)
    is_3d = wfs_signal_x.ndim == 3

    # Reshape gracefully handling both 2D and 3D arrays
    wfs_signal_x_2d = wfs_signal_x.reshape((-1, wfs_signal_x.shape[2] if is_3d else 1))
    wfs_signal_y_2d = wfs_signal_y.reshape((-1, wfs_signal_y.shape[2] if is_3d else 1))

    # Select valid subapertures
    idx = _valid_sa_indices(idx_valid_sa, wfs_nsubaps, specula_convention)
    if idx is not None:
        wfs_signal_x_2d = wfs_signal_x_2d[to_xp(xp, idx), :]
        wfs_signal_y_2d = wfs_signal_y_2d[to_xp(xp, idx), :]

    # Concatenate
    if specula_convention:
        im = xp.concatenate((wfs_signal_y_2d, wfs_signal_x_2d))
    else:
        im = xp.concatenate((wfs_signal_x_2d, wfs_signal_y_2d))

    # Convert to slope units
    pup_diam_pix = pup_mask.shape[0]
    pixel_pitch = pup_diam_m / pup_diam_pix
    coeff = 1e-9 / (pup_diam_m / wfs_nsubaps) * 206265
    coeff *= 1 / (0.5 * wfs_fov_arcsec)
    im = im * coeff

    if verbose:
        print(f'  Slopes computed, shape: {im.shape}')

    return im


def _check_dm_mask_shape(dm_array, dm_mask, specula_convention):
    """Check that the DM array and mask have the same size along the first axis
    of the transformed arrays (the second one in the SPECULA convention)."""
    axis = 1 if specula_convention else 0
    if dm_mask.shape[axis] != dm_array.shape[axis]:
        raise ValueError('DM and mask arrays must have the same dimensions.')


def _transform_mask(mask, transform, threshold=0.5):
    """
    Mask transformed as rotshiftzoom_array(mask, **transform), with the
    values below threshold set to 0. The interpolation is computed on CPU in
    float64 (affine_interpolation_matrix), so that the masks, and the
    operators that depend on them through their thresholds, are the same on
    CPU and GPU (cupyx affine_transform uses single precision coordinates
    for single precision arrays). Returns a float_dtype array on the target
    device.
    """
    transform = dict(transform)
    output_size = transform.pop('output_size')
    mask = np.nan_to_num(cpuArray(mask).astype(np.float64), nan=0.0)
    matrix, offset = rotshiftzoom_matrix(mask.shape, output_size, **transform)
    trans_mask = (affine_interpolation_matrix(mask.shape, output_size, matrix, offset)
                  @ mask.ravel()).reshape(output_size)
    trans_mask[trans_mask < threshold] = 0
    return to_xp(xp, trans_mask, dtype=float_dtype)


def _combined_transformation(pup_diam_m, pup_mask, dm_mask, dm_height, dm_rotation,
                             gs_pol_coo, gs_height, wfs_rotation, wfs_translation,
                             wfs_mag_global, wfs_anamorphosis_90, wfs_anamorphosis_45,
                             specula_convention, verbose=False):
    """
    Parameters of the combined DM + WFS transformation and transformed masks.

    Returns:
    - transform: dict of rotshiftzoom_array keyword arguments that apply the
      combined transformation to the DM array (transposed in the SPECULA
      convention)
    - trans_dm_mask: transformed DM mask (values below 0.5 set to 0)
    - trans_pup_mask: pupil mask with the WFS transformations only (values
      below 0.5 set to 0)
    """
    # *** Compute WFS magnification including anamorphosis at 90° ***
    wfs_magnification = (wfs_mag_global, wfs_mag_global * wfs_anamorphosis_90)

    # The masks are transformed on CPU (see _transform_mask)
    dm_mask = np.asarray(cpuArray(dm_mask), dtype=float_dtype)
    pup_mask = np.asarray(cpuArray(pup_mask), dtype=float_dtype)

    if specula_convention:
        dm_mask = np.transpose(dm_mask)
        pup_mask = np.transpose(pup_mask)
        wfs_translation_local = (-1*wfs_translation[1], -1*wfs_translation[0])
    else:
        wfs_translation_local = wfs_translation

    pup_diam_pix = pup_mask.shape[0]
    pixel_pitch = pup_diam_m / pup_diam_pix

    dm_translation, dm_magnification = shiftzoom_from_source_dm_params(
        source_pol_coo=gs_pol_coo,
        source_height=gs_height,
        dm_height=dm_height,
        pixel_pitch=pixel_pitch
    )
    output_size = (pup_diam_pix, pup_diam_pix)

    if verbose:
        print(f'Combined DM+WFS transformations:')
        print(f'  DM translation: {dm_translation} pixels')
        print(f'  DM rotation: {dm_rotation} deg')
        print(f'  DM magnification: {dm_magnification}')
        print(f'  WFS translation: {wfs_translation} pixels')
        print(f'  WFS rotation: {wfs_rotation} deg')
        print(f'  WFS magnification: {wfs_magnification}')

    transform = dict(
        dm_translation=dm_translation,
        dm_rotation=dm_rotation,
        dm_magnification=dm_magnification,
        wfs_translation=wfs_translation_local,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45,
        output_size=output_size
    )

    # DM mask (same transformation as the DM array)
    trans_dm_mask = _transform_mask(dm_mask, transform)

    # Pupil mask (only WFS transformations)
    trans_pup_mask = _transform_mask(pup_mask, dict(
        transform, dm_translation=(0, 0), dm_rotation=0, dm_magnification=(1, 1)))

    if xp.max(trans_dm_mask) <= 0:
        raise ValueError('Transformed DM mask is empty.')
    if xp.max(trans_pup_mask) <= 0:
        raise ValueError('Transformed pupil mask is empty.')

    return transform, trans_dm_mask, trans_pup_mask


def apply_dm_transformations_combined(pup_diam_m, pup_mask, dm_array, dm_mask,
                                      dm_height, dm_rotation,
                                      gs_pol_coo, gs_height,
                                      wfs_rotation, wfs_translation,
                                      wfs_nsubaps=None,
                                      wfs_mag_global=1.0,
                                      wfs_anamorphosis_90=1.0,
                                      wfs_anamorphosis_45=1.0,
                                      slope_method='derivatives',
                                      specula_convention=True,
                                      verbose=False):
    """
    Apply DM and WFS transformations COMBINED (single interpolation step).
    This avoids cumulative interpolation errors when both DM and WFS have rotations.

    The returned trans_dm_array is the transformed DM array after masking and
    edge repair (the repair is done in place to avoid a copy of the cube).
    """

    dm_array = to_xp(xp, dm_array, dtype=float_dtype)
    _check_dm_mask_shape(dm_array, dm_mask, specula_convention)
    if specula_convention:
        dm_array = xp.transpose(dm_array, (1, 0, 2))

    transform, trans_dm_mask, trans_pup_mask = _combined_transformation(
        pup_diam_m, pup_mask, dm_mask, dm_height, dm_rotation, gs_pol_coo, gs_height,
        wfs_rotation, wfs_translation, wfs_mag_global, wfs_anamorphosis_90,
        wfs_anamorphosis_45, specula_convention, verbose=verbose)

    # Apply ALL transformations in one step
    trans_dm_array = rotshiftzoom_array(dm_array, **transform)

    # trans_dm_array is a new array: mask it and repair its edges in place
    # (same values as apply_mask, without a copy of the cube)
    if trans_dm_array.ndim == 3:
        trans_dm_array *= trans_dm_mask[:, :, xp.newaxis]
    else:
        trans_dm_array *= trans_dm_mask

    # Compute derivatives or telescoping sum on already-transformed array
    if slope_method == 'derivatives':
        derivatives_x, derivatives_y = compute_derivatives_with_extrapolation(
            trans_dm_array, mask=trans_dm_mask, in_place=True
        )
    elif slope_method == 'telsum':
        derivatives_x, derivatives_y = compute_telsum_with_extrapolation(
            trans_dm_array, mask=trans_dm_mask, wfs_nsubaps=wfs_nsubaps, verbose=verbose,
            in_place=True, pup_mask=trans_pup_mask
        )
    else:
        raise ValueError(f"Unknown slope_method: {slope_method}")

    if verbose:
        print(f'  Combined transformation applied, shape: {trans_dm_array.shape}')
        print(f'  Slopes/Derivatives computed')

    return trans_dm_array, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y


def apply_wfs_transformations_combined(derivatives_x, derivatives_y, trans_pup_mask, dm_mask,
                                       wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa=None,
                                       slope_method='derivatives', specula_convention=True,
                                       verbose=False, in_place=False):
    """
    Compute slopes from pre-transformed derivatives (for combined workflow).
    No additional transformations needed.

    in_place: with slope_method='derivatives', if True derivatives_x and
        derivatives_y are modified (no copy of the arrays): calling the
        function again on the same derivatives gives correct results only
        with the same pupil mask (trans_pup_mask). If False (default) the
        function works on copies and the inputs are not modified.
    """

    # Derivatives are already transformed - route to the appropriate formatter
    if slope_method == 'derivatives':
        if not in_place:
            derivatives_x = derivatives_x.copy()
            derivatives_y = derivatives_y.copy()
        return _compute_slopes_from_derivatives(
            derivatives_x, derivatives_y, trans_pup_mask, dm_mask,
            wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa,
            verbose, specula_convention
        )
    elif slope_method == 'telsum':
        return _compute_slopes_from_telsum(
            derivatives_x, derivatives_y, trans_pup_mask, dm_mask,
            wfs_nsubaps, wfs_fov_arcsec, pup_diam_m,
            idx_valid_sa, verbose, specula_convention
        )
    else:
        raise ValueError(f"Unknown slope_method: {slope_method}")


def _gradient_operator(n):
    """(n, n) sparse matrix of _gradient along an axis of length n."""
    gradient = scipy.sparse.diags([-0.5, 0.5], [-1, 1], shape=(n, n), format='lil')
    gradient[0, :2] = [-1, 1]
    gradient[n - 1, n - 2:] = [-1, 1]
    return gradient.tocsr()


def _repair_operator(strict_mask):
    """
    Edge repair factor R (n_pixels, n_pixels) of the interaction matrix
    operator: same result as _repair_interpolated_phase. The pixels of
    strict_mask are kept, the edge pixels around it are extrapolated from up
    to 8 of its pixels (calculate_extrapolation_indices_coeffs), all the
    others are 0.
    """
    strict_pixels = np.flatnonzero(strict_mask)
    edge_pixels, reference_indices, coefficients = calculate_extrapolation_indices_coeffs(
        strict_mask.astype(np.int64), debug=False, debug_pixels=None)
    used = (reference_indices >= 0) & (coefficients != 0)
    edge_rows = np.broadcast_to(np.asarray(edge_pixels)[:, np.newaxis], used.shape)
    rows = np.concatenate((strict_pixels, edge_rows[used]))
    cols = np.concatenate((strict_pixels, reference_indices[used]))
    values = np.concatenate((np.ones(len(strict_pixels)),
                             coefficients[used].astype(np.float64)))
    n = strict_mask.size
    return scipy.sparse.csr_matrix((values, (rows, cols)), shape=(n, n))


def _slope_operators(trans_dm_mask, trans_pup_mask, wfs_nsubaps, slope_method,
                     threshold=0.999999):
    """
    Factors M_sa . G . R . D_mask of the interaction matrix operator (see
    interaction_matrix): sparse (wfs_nsubaps**2, n_pixels) float64 operators
    from the transformed DM array to the x and y slopes of all the
    subapertures (row-major), before the selection of the valid
    subapertures and the conversion to slope units (_wfs_operator).
    Same slopes as the explicit computation: masking, edge repair and
    derivatives (_compute_slopes_from_derivatives), or differences of valid
    adjacent pixels (_telsum_fractional, _compute_slopes_from_telsum), and
    mean over each subaperture with the pixel fractions of rebin_matrix.
    """
    mask = np.nan_to_num(cpuArray(trans_dm_mask).astype(np.float64), nan=0.0)
    pup_mask = np.nan_to_num(cpuArray(trans_pup_mask).astype(np.float64), nan=0.0)
    if mask.max() < threshold:
        raise ValueError(f'Mask max value is {mask.max()},'
                         f' expected binary mask with values 0 and 1.')
    n0, n1 = mask.shape
    strict_mask = mask >= threshold
    in_pupil = (pup_mask > 0).astype(np.float64)

    # Pixel fraction inside each subaperture along each axis (see rebin)
    rows = _rebin_matrix_np(n0, wfs_nsubaps)
    cols = _rebin_matrix_np(n1, wfs_nsubaps)
    dm_mask_sa = rows @ mask @ cols.T
    combined_mask_sa = (dm_mask_sa > 0) & (rows @ pup_mask @ cols.T > 0)

    def diag(values):
        return scipy.sparse.diags(values.ravel().astype(np.float64), format='csr')

    def kron(a, b):
        return scipy.sparse.kron(scipy.sparse.csr_matrix(a), scipy.sparse.csr_matrix(b),
                                 format='csr')

    def adjacent_differences(n):
        return scipy.sparse.diags([-1.0, 1.0], [0, 1], shape=(n - 1, n))

    masking = diag(mask)
    if slope_method == 'derivatives':
        if dm_mask_sa.max() <= 0:
            raise ValueError('DM mask is empty after rebinning.')
        # see _compute_slopes_from_derivatives
        repair = _repair_operator(strict_mask) @ masking
        weight = rows @ in_pupil @ cols.T
        weight = np.where(weight > 0, weight, 1)
        scale = np.where(combined_mask_sa, (n0 / wfs_nsubaps) / weight, 0)
        mean = diag(scale) @ kron(rows, cols) @ diag(strict_mask & (in_pupil > 0))
        gradient_x = kron(scipy.sparse.identity(n0), _gradient_operator(n1))
        gradient_y = kron(_gradient_operator(n0), scipy.sparse.identity(n1))
        return mean @ gradient_x @ repair, mean @ gradient_y @ repair

    if slope_method == 'telsum':
        if n0 < 2 * wfs_nsubaps:
            raise ValueError(f'The telescoping sum needs at least 2 pixels per subaperture:'
                             f' {n0} pixels, {wfs_nsubaps} subapertures.')
        # see _telsum_fractional and _compute_slopes_from_telsum
        def operator(axis):
            if axis == 1:
                used = in_pupil[:, 1:] * in_pupil[:, :-1]
                valid = strict_mask[:, 1:] & strict_mask[:, :-1]
                row_matrix, col_matrix = rows, np.minimum(cols[:, 1:], cols[:, :-1])
                difference = kron(scipy.sparse.identity(n0), adjacent_differences(n1))
            else:
                used = in_pupil[1:] * in_pupil[:-1]
                valid = strict_mask[1:] & strict_mask[:-1]
                row_matrix, col_matrix = np.minimum(rows[:, 1:], rows[:, :-1]), cols
                difference = kron(adjacent_differences(n0), scipy.sparse.identity(n1))
            weight = row_matrix @ used @ col_matrix.T
            scale = np.where(combined_mask_sa & (weight > 0),
                             (n0 / wfs_nsubaps) / np.where(weight > 0, weight, 1), 0)
            mean = diag(scale) @ kron(row_matrix, col_matrix) @ diag(valid * used)
            return mean @ difference @ masking

        return operator(1), operator(0)

    raise ValueError(f"Unknown slope_method: {slope_method}")


def _wfs_operator(slopes_x, slopes_y, wfs_nsubaps, wfs_fov_arcsec, pup_diam_m,
                  idx_valid_sa, specula_convention):
    """
    Factors C . P_valid of the interaction matrix operator: selection of the
    valid subapertures, x/y order of the convention and conversion to slope
    units, applied to the slope operators of _slope_operators.
    """
    idx = _valid_sa_indices(idx_valid_sa, wfs_nsubaps, specula_convention)
    if idx is not None:
        slopes_x = slopes_x[idx]
        slopes_y = slopes_y[idx]
    pair = (slopes_y, slopes_x) if specula_convention else (slopes_x, slopes_y)
    coeff = 1e-9 / (pup_diam_m / wfs_nsubaps) * 206265
    coeff *= 1 / (0.5 * wfs_fov_arcsec)
    return scipy.sparse.vstack(pair, format='csr') * coeff


def _transformation_operator(transform, dm_shape, specula_convention):
    """
    Factor T of the interaction matrix operator: sparse bilinear interpolation
    matrix of the combined transformation (_combined_transformation), applied
    to the DM array as given (the SPECULA transposition is folded into the
    column indices, see affine_interpolation_matrix).
    """
    transform = dict(transform)
    output_size = transform.pop('output_size')
    input_shape = (dm_shape[1], dm_shape[0]) if specula_convention else tuple(dm_shape[:2])
    matrix, offset = rotshiftzoom_matrix(input_shape, output_size, **transform)
    return affine_interpolation_matrix(input_shape, output_size, matrix, offset,
                                       transposed=specula_convention)


# CSR matrix times C-contiguous dense matrix on GPU, one thread per output
# element: deterministic (cuSPARSE SpMM is not) and without the F-contiguous
# copy of the DM array that cuSPARSE needs
_CSR_DENSE_KERNEL_CODE = """
    const long long row = i / n_cols;
    const long long col = i % n_cols;
    T acc = 0;
    for (I k = indptr[row]; k < indptr[row + 1]; k++) {
        acc += data[k] * dense[(long long)indices[k] * n_cols + col];
    }
    out = acc;
"""
_csr_dense_kernel = None


def _apply_operator(operator, dm_array):
    """
    IM = A . Phi: operator A times dm_array seen as Phi, a
    (n_pixels, n_modes) matrix (a reshape, not a copy), on the target device.
    On GPU the product uses _CSR_DENSE_KERNEL_CODE (deterministic, reads the
    C-ordered cube directly). Non-finite values of dm_array are replaced as
    in rotshiftzoom_array (nan_to_num, on a copy).
    """
    n_modes = dm_array.shape[2] if dm_array.ndim == 3 else 1
    dm_2d = xp.ascontiguousarray(dm_array).reshape(-1, n_modes)
    # The sum is non-finite if any value is (NaN or inf); the full check
    # avoids a copy if the sum only overflows
    if not xp.isfinite(xp.sum(dm_2d)) and not xp.isfinite(dm_2d).all():
        dm_2d = xp.nan_to_num(dm_2d, copy=True, nan=0.0, posinf=None, neginf=None)
    operator = scipy.sparse.csr_matrix(operator, dtype=float_dtype)
    operator.sum_duplicates()  # same summation order on CPU and GPU
    operator.eliminate_zeros()  # e.g. cancelled terms of the telescoping sums
    if xp is np:
        return np.ascontiguousarray(operator @ dm_2d)

    global _csr_dense_kernel
    if _csr_dense_kernel is None:
        _csr_dense_kernel = xp.ElementwiseKernel(
            'raw T data, raw I indptr, raw I indices, raw T dense, int64 n_cols',
            'T out', _CSR_DENSE_KERNEL_CODE, 'synim_csr_dense_product')
    out = xp.empty((operator.shape[0], n_modes), dtype=dm_2d.dtype)
    _csr_dense_kernel(xp.asarray(operator.data), xp.asarray(operator.indptr),
                      xp.asarray(operator.indices), dm_2d, n_modes, out)
    return out


def interaction_matrix(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                       wfs_nsubaps, wfs_fov_arcsec, gs_pol_coo, gs_height,
                       wfs_rotation, wfs_translation, wfs_mag_global,
                       wfs_anamorphosis_90=1.0, wfs_anamorphosis_45=1.0,
                       idx_valid_sa=None, slope_method='derivatives',
                       specula_convention=True, verbose=False,
                       display=False):
    """
    Computes an interaction matrix.

    All the geometric transformations (DM shift, rotation and magnification
    from the guide star direction and height, WFS shift, rotation,
    magnification and anamorphosis) are combined into a single affine
    transformation applied to the DM modes (the phase), with one
    interpolation step. The slopes (derivatives or telescoping sum) are then
    computed on the transformed phase.

    Method (see "Computation Workflows" in docs/general.rst). The masks do
    not depend on the DM mode, and once they are fixed every step is linear
    in the phase phi of a mode (raveled to n_pix pixels), so the whole chain
    is one sparse matrix A:

        s = C . P_valid . M_sa . G . R . D_mask . T . phi = A . phi

    T        bilinear interpolation of the combined transformation
             (_transformation_operator)
    D_mask   multiplication by the transformed DM mask
    R        edge repair (_repair_operator; derivatives only)
    G        derivatives, or differences of adjacent pixels (telescoping sum)
    M_sa     mean over each, possibly fractional, subaperture
             (D_mask, G and M_sa: _slope_operators)
    P_valid, C  valid subapertures, x/y order, slope units (_wfs_operator)

    With the DM cube seen as Phi (n_pix x n_modes, a reshape), the
    interaction matrix is the single product IM = A . Phi (_apply_operator).
    A is built once on CPU in float64 from the masks and rounded once to
    float_dtype; it does not depend on n_modes and is very sparse (the mean
    of the derivatives over a subaperture telescopes to its border pixels),
    so no transformed copy of the cube or of its derivatives is allocated.
    The explicit step-by-step computation (apply_dm_transformations_combined
    + apply_wfs_transformations_combined) gives the same results up to
    rounding, and is kept as the reference and for debugging.

    Note: WFS transformations must not be applied to derivative maps. The
    gradient of the transformed phase includes the Jacobian of the
    transformation (it mixes the x and y components for rotations and 45 deg
    anamorphosis and scales them for magnification and 90 deg anamorphosis),
    which a spatial resampling of the derivative maps does not account for.
    """

    if verbose:
        print(f"\n{'='*60}")
        print(f"Interaction Matrix Computation")
        print(f"{'='*60}\n")

    dm_array = to_xp(xp, dm_array, dtype=float_dtype)
    _check_dm_mask_shape(dm_array, dm_mask, specula_convention)

    transform, trans_dm_mask, trans_pup_mask = _combined_transformation(
        pup_diam_m, pup_mask, dm_mask, dm_height, dm_rotation, gs_pol_coo, gs_height,
        wfs_rotation, wfs_translation, wfs_mag_global, wfs_anamorphosis_90,
        wfs_anamorphosis_45, specula_convention, verbose=verbose)

    # The whole computation (transformation, masking, edge repair, slopes)
    # is a sparse linear operator applied to the DM array
    slopes_x, slopes_y = _slope_operators(trans_dm_mask, trans_pup_mask, wfs_nsubaps,
                                          slope_method)
    transformation = _transformation_operator(transform, dm_array.shape, specula_convention)
    operator = _wfs_operator(slopes_x, slopes_y, wfs_nsubaps, wfs_fov_arcsec, pup_diam_m,
                             idx_valid_sa, specula_convention) @ transformation
    im = _apply_operator(operator, dm_array)

    if verbose:
        print(f'  Interaction matrix computed with a sparse operator'
              f' ({operator.nnz} non-zero elements), shape: {im.shape}')

    if display:
        idx_plot = [2, 5]
        pup_mask_cpu = cpuArray(trans_pup_mask)
        trans_dm_mask_cpu = cpuArray(trans_dm_mask)
        # Transformed and masked DM shapes of the two modes
        dm_modes = cpuArray(dm_array.reshape(transformation.shape[1], -1)[:, idx_plot])
        trans_dm_array_cpu = (transformation @ dm_modes).reshape(
            trans_dm_mask_cpu.shape + (len(idx_plot),)) * trans_dm_mask_cpu[:, :, np.newaxis]
        fig, axs = plt.subplots(2, 2)
        im3 = axs[0, 0].imshow(pup_mask_cpu, cmap='seismic')
        axs[0, 1].imshow(trans_dm_mask_cpu, cmap='seismic')
        axs[1, 0].imshow(trans_dm_array_cpu[:, :, 0], cmap='seismic')
        axs[1, 1].imshow(trans_dm_array_cpu[:, :, 1], cmap='seismic')
        fig.suptitle(f'Mask, DM mask, DM shapes (modes {idx_plot[0]} and {idx_plot[1]})')
        fig.colorbar(im3, ax=axs.ravel().tolist(), fraction=0.02)
        plt.show()

    _warn_if_nan(im)
    return im


def _wfs_magnification_params(wfs_config):
    """
    Return (wfs_mag_global, wfs_anamorphosis_90) from the 'magnification'
    entry of a WFS configuration: a (x, y) pair or a single value.
    """
    magnification = wfs_config.get('magnification', (1.0, 1.0))
    if isinstance(magnification, (tuple, list)) and len(magnification) == 2:
        wfs_mag_global = math.sqrt(magnification[0] * magnification[1])
        wfs_anamorphosis_90 = magnification[1] / magnification[0] \
            if magnification[0] != 0 else 1.0
    else:
        wfs_mag_global = float(magnification)
        wfs_anamorphosis_90 = 1.0
    return wfs_mag_global, wfs_anamorphosis_90


def _as_key(value):
    """Hashable, type-independent version of a scalar or a sequence of numbers."""
    if hasattr(value, '__len__'):
        return tuple(float(v) for v in value)
    return float(value)


def _dm_transformation_key(wfs_params, slope_method):
    """
    Parameters that determine the transformed DM array and the derivatives:
    WFS that share them can share the computation.
    """
    return (
        _as_key(wfs_params['gs_pol_coo']), _as_key(wfs_params['gs_height']),
        _as_key(wfs_params['rotation']), _as_key(wfs_params['translation']),
        _as_key(wfs_params['mag_global']), _as_key(wfs_params['anamorphosis_90']),
        _as_key(wfs_params['anamorphosis_45']), int(wfs_params['nsubaps']), slope_method,
    )


def interaction_matrices_multi_wfs(pup_diam_m, pup_mask,
                                   dm_array, dm_mask,
                                   dm_height, dm_rotation,
                                   wfs_configs, gs_pol_coo=None,
                                   gs_height=None,
                                   slope_method='derivatives',
                                   specula_convention=True,
                                   im_on_cpu=False,
                                   minimize_memory=False,
                                   verbose=False):
    """
    Computes interaction matrices for multiple WFS configurations.
    
    Each WFS can have its own guide star position (gs_pol_coo) and height (gs_height).
    The DM array is converted to the target device once; each WFS is
    computed as in interaction_matrix (single combined transformation, as a
    sparse operator applied to the DM array).
    WFS with the same guide star, WFS transformations, number of subapertures
    and slope method share the transformation and the slope operators,
    which are computed once for the whole group.
    
    Parameters:
    - pup_diam_m: float, pupil diameter in meters
    - pup_mask: numpy 2D array, pupil mask (n_pup x n_pup)
    - dm_array: numpy 3D array, DM modes (n x n x n_dm_modes)
    - dm_mask: numpy 2D array, DM mask (n x n)
    - dm_height: float, DM conjugation altitude
    - dm_rotation: float, DM rotation in degrees
    - wfs_configs: list of dict, each containing WFS parameters
    - gs_pol_coo: tuple or None (DEPRECATED)
    - gs_height: float or None (DEPRECATED)
    - slope_method: str, 'derivatives' or 'telsum'
    - specula_convention: bool, optional
    - im_on_cpu: bool, optional, force output interaction matrices on CPU
    - minimize_memory: bool, optional, kept for compatibility: no transformed
      copy of the DM array is allocated
    - verbose: bool, optional
    
    Returns:
    - im_dict: dict, interaction matrices keyed by WFS name or index
    - derivatives_info: dict with metadata about the computation:
        'workflow': 'combined'
        'groups': list of lists of WFS names sharing the DM transformation
        'n_dm_transformations': number of DM transformations computed
    """

    if verbose:
        print(f"\n{'='*60}")
        print(f"Computing interaction matrices for {len(wfs_configs)} WFS")
        print(f"{'='*60}")

    # Check if using deprecated global gs_pol_coo/gs_height
    use_global_gs = gs_pol_coo is not None and gs_height is not None

    if use_global_gs and verbose:
        print("WARNING: Using global gs_pol_coo and gs_height for all WFS (deprecated)")
        print("         Consider specifying gs_pol_coo and gs_height in each wfs_config")

    # Extract gs_pol_coo and gs_height for each WFS (checked before computing)
    wfs_gs_info = []
    for i, wfs_config in enumerate(wfs_configs):
        if use_global_gs:
            wfs_gs_pol_coo = gs_pol_coo
            wfs_gs_height = gs_height
        else:
            if 'gs_pol_coo' not in wfs_config:
                raise ValueError(f"WFS {i}: 'gs_pol_coo' must be"
                                 f" specified in wfs_config when gs_pol_coo=None")
            if 'gs_height' not in wfs_config:
                raise ValueError(f"WFS {i}: 'gs_height' must be"
                                 f" specified in wfs_config when gs_height=None")

            wfs_gs_pol_coo = wfs_config['gs_pol_coo']
            wfs_gs_height = wfs_config['gs_height']

        wfs_gs_info.append((wfs_gs_pol_coo, wfs_gs_height))

    # Parameters of each WFS
    wfs_params = []
    for i, wfs_config in enumerate(wfs_configs):
        wfs_mag_global, wfs_anamorphosis_90 = _wfs_magnification_params(wfs_config)
        wfs_params.append(dict(
            name=wfs_config.get('name', f'wfs_{i}'),
            nsubaps=wfs_config['nsubaps'],
            rotation=wfs_config.get('rotation', 0.0),
            translation=wfs_config.get('translation', (0.0, 0.0)),
            mag_global=wfs_mag_global,
            anamorphosis_90=wfs_anamorphosis_90,
            anamorphosis_45=wfs_config.get('anamorphosis_45', 1.0),
            fov_arcsec=wfs_config['fov_arcsec'],
            idx_valid_sa=wfs_config.get('idx_valid_sa', None),
            gs_pol_coo=wfs_gs_info[i][0],
            gs_height=wfs_gs_info[i][1],
        ))

    # WFS with the same guide star, WFS transformations, number of
    # subapertures and slope method share the transformed DM array and the
    # derivatives: they are computed once per group (in the order of first
    # appearance), then the slopes are computed for each WFS of the group.
    groups = {}
    for i, params in enumerate(wfs_params):
        groups.setdefault(_dm_transformation_key(params, slope_method), []).append(i)

    # Convert the DM array once, instead of once per WFS
    dm_array = to_xp(xp, dm_array, dtype=float_dtype)
    _check_dm_mask_shape(dm_array, dm_mask, specula_convention)

    operators = [None] * len(wfs_params)
    for i_group, members in enumerate(groups.values()):
        first = wfs_params[members[0]]

        if verbose:
            names = ', '.join(str(wfs_params[i]['name']) for i in members)
            print(f"  [group {i_group+1}/{len(groups)}] {names}:")
            print(f"    Subapertures: {first['nsubaps']}x{first['nsubaps']}")
            print(f"    GS: {first['gs_pol_coo']}, height: {first['gs_height']} m")

        transform, trans_dm_mask, trans_pup_mask = _combined_transformation(
            pup_diam_m, pup_mask, dm_mask, dm_height, dm_rotation,
            first['gs_pol_coo'], first['gs_height'], first['rotation'],
            first['translation'], first['mag_global'], first['anamorphosis_90'],
            first['anamorphosis_45'], specula_convention)
        slopes_x, slopes_y = _slope_operators(trans_dm_mask, trans_pup_mask,
                                              first['nsubaps'], slope_method)
        transformation = _transformation_operator(transform, dm_array.shape,
                                                  specula_convention)

        # Only the subaperture selection and the slope units differ in a group
        for i in members:
            params = wfs_params[i]
            operators[i] = _wfs_operator(
                slopes_x, slopes_y, params['nsubaps'], params['fov_arcsec'], pup_diam_m,
                params['idx_valid_sa'], specula_convention) @ transformation

    # All the interaction matrices with a single product, reading the DM
    # array once. Each row of the product is computed independently (scipy
    # on CPU, _CSR_DENSE_KERNEL_CODE on GPU), so the results are identical
    # to interaction_matrix.
    im_all = _apply_operator(scipy.sparse.vstack(operators, format='csr'), dm_array)
    bounds = np.cumsum([0] + [operator.shape[0] for operator in operators])
    im_list = []
    for params, start, stop in zip(wfs_params, bounds[:-1], bounds[1:]):
        im = im_all[start:stop]
        im_list.append(cpuArray(im) if im_on_cpu else im)
        if verbose:
            print(f"    {params['name']}: IM shape {im.shape}")

    im_dict = {params['name']: im for params, im in zip(wfs_params, im_list)}
    for name, im in im_dict.items():
        _warn_if_nan(im, name=f'interaction matrix of {name}')

    derivatives_info = {
        'workflow': 'combined',
        # WFS names of each group sharing the DM transformation and slope operators
        'groups': [[wfs_params[i]['name'] for i in members] for members in groups.values()],
        'n_dm_transformations': len(groups),
    }

    if verbose:
        print(f"\n{'='*60}")
        print(f"Completed {len(im_dict)} interaction matrices"
              f" ({len(groups)} DM transformations)")
        print(f"{'='*60}\n")

    return im_dict, derivatives_info


def compute_subaperture_illumination(pup_mask, wfs_nsubaps, wfs_rotation=0.0,
                                    wfs_translation=(0.0, 0.0),
                                    wfs_magnification=(1.0, 1.0),
                                    idx_valid_sa=None,
                                    specula_convention=True,
                                    verbose=False):
    """
    Compute the relative illumination of valid subapertures.
    
    This is useful for weighting the noise covariance matrix based on the 
    actual flux received by each subaperture (edge subapertures receive less light).
    
    Parameters:
    - pup_mask: numpy 2D array, pupil mask
    - wfs_nsubaps: int, number of subapertures along diameter
    - wfs_rotation: float, WFS rotation in degrees
    - wfs_translation: tuple, WFS translation (x, y) in pixels
    - wfs_magnification: tuple, WFS magnification (x, y)
    - idx_valid_sa: array, indices of valid subapertures (2D format preferred)
    - specula_convention: bool, whether to use SPECULA convention
    - verbose: bool, whether to print information
    
    Returns:
    - illumination: 1D array, relative illumination of each valid subaperture (normalized to max=1)
    """

    # *** Convert to target device ***
    pup_mask = to_xp(xp, pup_mask, dtype=float_dtype)

    if specula_convention:
        pup_mask = xp.transpose(pup_mask)

    output_size = pup_mask.shape

    # Apply WFS transformations to pupil mask using the EXACT output size (W, W)
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

    # We must keep the fractionally interpolated pixels (anti-aliased edges) 
    # to correctly integrate the flux per subaperture (e.g. across spider arms).

    if xp.max(trans_pup_mask) <= 0:
        raise ValueError('Transformed pupil mask is empty.')

    # Rebin to WFS resolution - use 'sum' to get total flux per subaperture
    pup_mask_sa = rebin(trans_pup_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')

    # Normalize to theoretical maximum (fully illuminated subaperture)
    max_illumination = xp.max(pup_mask_sa)
    if max_illumination > 0:
        pup_mask_sa = pup_mask_sa / max_illumination

    # Flatten to 1D
    illumination_2d = pup_mask_sa.flatten()

    # Select only valid subapertures
    idx = _valid_sa_indices(idx_valid_sa, wfs_nsubaps, specula_convention)
    illumination = illumination_2d if idx is None else illumination_2d[to_xp(xp, idx)]

    # Convert to CPU for return
    illumination = cpuArray(illumination)

    if verbose:
        print(f"Subaperture illumination statistics:")
        print(f"  Min: {illumination.min():.3f}")
        print(f"  Max: {illumination.max():.3f}")
        print(f"  Mean: {illumination.mean():.3f}")
        print(f"  Std: {illumination.std():.3f}")

    return illumination
