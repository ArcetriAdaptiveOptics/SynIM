"""
Reference implementation of the interaction matrix: the explicit,
step-by-step computation on the DM cube.

For each mode the phase is transformed with the combined DM + WFS
transformation (one interpolation step), multiplied by the transformed DM
mask, repaired at the edges, differentiated (derivatives or telescoping
sum) and averaged over each subaperture; then the valid subapertures are
selected and converted to slope units. Each step returns arrays that can be
inspected.

This is the implementation where the physics is defined: every change to
the computation is made here first, checked with the physical tests, then
ported to the optimized implementation (synim.sparse_operator), and the
equivalence test (test/test_sparse_operator.py) must pass. See "Reference
and optimized implementations" in docs/general.rst.

The geometry shared with the optimized implementation (transformation
parameters, transformed masks, valid subapertures, slope units) is in
synim.geometry.
"""
import synim as _synim
_synim._require_init(__name__)
from synim import xp, cpuArray, to_xp, float_dtype
from synim.utils import (
    apply_mask,
    rebin,
    rebin_matrix,
    rebin_fractional_sum,
    rotshiftzoom_array,
    apply_extrapolation,
    calculate_extrapolation_indices_coeffs
)
from synim.geometry import (
    _without_nan,
    _valid_sa_indices,
    _slope_units,
    _check_dm_mask_shape,
    _combined_transformation,
)


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
    coeff = _slope_units(pup_diam_m, wfs_nsubaps, wfs_fov_arcsec)
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
    im = im * _slope_units(pup_diam_m, wfs_nsubaps, wfs_fov_arcsec)

    if verbose:
        print(f'  Slopes computed, shape: {im.shape}')

    return im


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
    trans_dm_mask = to_xp(xp, trans_dm_mask)
    trans_pup_mask = to_xp(xp, trans_pup_mask)

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



def interaction_matrix(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                       wfs_nsubaps, wfs_fov_arcsec, gs_pol_coo, gs_height,
                       wfs_rotation, wfs_translation, wfs_mag_global,
                       wfs_anamorphosis_90=1.0, wfs_anamorphosis_45=1.0,
                       idx_valid_sa=None, slope_method='derivatives',
                       specula_convention=True, verbose=False):
    """
    Interaction matrix computed step by step on the DM cube
    (apply_dm_transformations_combined, then
    apply_wfs_transformations_combined). Same parameters and result as
    synim.synim.interaction_matrix with method='reference'.
    """
    _, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
        apply_dm_transformations_combined(
            pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
            gs_pol_coo, gs_height, wfs_rotation, wfs_translation,
            wfs_nsubaps=wfs_nsubaps, wfs_mag_global=wfs_mag_global,
            wfs_anamorphosis_90=wfs_anamorphosis_90,
            wfs_anamorphosis_45=wfs_anamorphosis_45, slope_method=slope_method,
            specula_convention=specula_convention, verbose=verbose)
    im = apply_wfs_transformations_combined(
        derivatives_x, derivatives_y, trans_pup_mask, trans_dm_mask, wfs_nsubaps,
        wfs_fov_arcsec, pup_diam_m, idx_valid_sa=idx_valid_sa, slope_method=slope_method,
        specula_convention=specula_convention, verbose=verbose, in_place=True)

    if verbose:
        print(f'  Interaction matrix computed step by step (reference), shape: {im.shape}')
    return im


def interaction_matrices(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                         wfs_params, groups, slope_method, specula_convention,
                         verbose=False):
    """
    Interaction matrices of several WFS (see
    synim.synim.interaction_matrices_multi_wfs, method='reference').
    wfs_params is the list of the WFS parameter dicts and groups the lists
    of indices of the WFS that share the transformed DM array and the
    derivatives, both built by interaction_matrices_multi_wfs. The
    derivatives are computed once per group, then the slopes of each WFS of
    the group. Returns the list of the interaction matrices, in the order
    of wfs_params.
    """
    im_list = [None] * len(wfs_params)
    for members in groups:
        first = wfs_params[members[0]]
        _, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
            apply_dm_transformations_combined(
                pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                first['gs_pol_coo'], first['gs_height'], first['rotation'],
                first['translation'], wfs_nsubaps=first['nsubaps'],
                wfs_mag_global=first['mag_global'],
                wfs_anamorphosis_90=first['anamorphosis_90'],
                wfs_anamorphosis_45=first['anamorphosis_45'], slope_method=slope_method,
                specula_convention=specula_convention)
        # The in-place cleaning of the derivatives is idempotent for the
        # pupil mask of the group (see apply_wfs_transformations_combined)
        for i in members:
            params = wfs_params[i]
            im_list[i] = apply_wfs_transformations_combined(
                derivatives_x, derivatives_y, trans_pup_mask, trans_dm_mask,
                params['nsubaps'], params['fov_arcsec'], pup_diam_m,
                idx_valid_sa=params['idx_valid_sa'], slope_method=slope_method,
                specula_convention=specula_convention, in_place=True)
        del derivatives_x, derivatives_y
    return im_list
