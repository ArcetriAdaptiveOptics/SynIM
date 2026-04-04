from synim import xp, cpuArray, to_xp, float_dtype
import matplotlib.pyplot as plt
from synim.utils import (
    apply_mask,
    rebin,
    rotshiftzoom_array,
    shiftzoom_from_source_dm_params,
    apply_extrapolation,
    calculate_extrapolation_indices_coeffs,
    has_transformations
)


def _repair_interpolated_phase(data, mask, threshold=0.999999):
    """
    Vital repair for edge artifacts. 
    Restores pixels that collapsed towards 0 due to bilinear interpolation 
    by using a strict threshold logic to select only uncorrupted pixels.
    """
    if mask is not None:
        if mask.max() < threshold:
            raise ValueError(f'Mask max value is {mask.max()},'
                             f' expected binary mask with values 0 and 1.')

        # 1. Strict binarization: isolate only the pure "core" that escaped interpolation
        strict_mask = xp.where(mask >= threshold, 1, 0)

        # 2. Zero out the softened edges corrupted by interpolation
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


def compute_gtilt_with_extrapolation(data, mask=None, wfs_nsubaps=None, verbose=False):
    """
    Computes the raw G-tilt (average phase difference per pixel) for each subaperture.
    Acts as the G-tilt equivalent of computing continuous derivatives.
    """
    if wfs_nsubaps is None:
        raise ValueError("wfs_nsubaps must be provided to compute G-tilt.")

    pup_diam_pix = data.shape[0]
    N = max(int(xp.ceil(pup_diam_pix / wfs_nsubaps)), 2)
    W = N * wfs_nsubaps

    # ON-THE-FLY INTERPOLATION
    if pup_diam_pix != W:
        if verbose:
            print(f"  * G-Tilt: Interpolating grid from {pup_diam_pix} to {W} (N={N})")
        mag = W / pup_diam_pix
        data = rotshiftzoom_array(data, dm_magnification=(mag, mag), output_size=(W, W))
        if mask is not None:
            mask = rotshiftzoom_array(mask, dm_magnification=(mag, mag), output_size=(W, W))
            mask[mask < 0.5] = 0

    # Repair interpolated edges using strict logic
    if mask is not None:
        repaired_data, strict_mask = _repair_interpolated_phase(data, mask)
    else:
        repaired_data = data
        strict_mask = xp.ones((data.shape[0], data.shape[1]), dtype=float_dtype)

    # Setup 4D arrays for the Telescoping Sum
    is_3d = repaired_data.ndim == 3
    if is_3d:
        n_modes = repaired_data.shape[2]
        p = repaired_data.reshape(wfs_nsubaps, N, wfs_nsubaps, N, n_modes)
        m = strict_mask.reshape(wfs_nsubaps, N, wfs_nsubaps, N, 1)
    else:
        p = repaired_data.reshape(wfs_nsubaps, N, wfs_nsubaps, N)
        m = strict_mask.reshape(wfs_nsubaps, N, wfs_nsubaps, N)

    # Fast 4D Telescoping Sum
    dx = p[:, :, :, 1:] - p[:, :, :, :-1]
    valid_dx = m[:, :, :, 1:] & m[:, :, :, :-1]
    sum_dx = xp.sum(xp.where(valid_dx, dx, 0.0), axis=(1, 3))
    weight_dx = xp.sum(valid_dx, axis=(1, 3))

    dy = p[:, 1:, :, :] - p[:, :-1, :, :]
    valid_dy = m[:, 1:, :, :] & m[:, :-1, :, :]
    sum_dy = xp.sum(xp.where(valid_dy, dy, 0.0), axis=(1, 3))
    weight_dy = xp.sum(valid_dy, axis=(1, 3))

    # Normalize to get Delta Phi per valid baseline
    raw_gtilt_x = xp.where(weight_dx > 0, (sum_dx / xp.where(weight_dx > 0, weight_dx, 1.0)), 0.0)
    raw_gtilt_y = xp.where(weight_dy > 0, (sum_dy / xp.where(weight_dy > 0, weight_dy, 1.0)), 0.0)

    return raw_gtilt_x, raw_gtilt_y


def _compute_slopes_from_gtilt(raw_gtilt_x, raw_gtilt_y, pup_mask, dm_mask,
                               wfs_nsubaps, wfs_fov_arcsec, pup_diam_m,
                               idx_valid_sa, verbose, specula_convention):
    """
    Formats the raw G-tilt into the final 1D slopes array.
    Signature is fully symmetrical to _compute_slopes_from_derivatives.
    """
    is_3d = raw_gtilt_x.ndim == 3

    # 1. Deduce N internally
    pup_diam_pix = pup_mask.shape[0]
    N = max(int(xp.ceil(pup_diam_pix / wfs_nsubaps)), 2)

    # 2. Rebin masks to compute valid subapertures (identical to derivative method)
    if xp.isnan(pup_mask).any():
        xp.nan_to_num(pup_mask, copy=False, nan=0.0)
    if xp.isnan(dm_mask).any():
        xp.nan_to_num(dm_mask, copy=False, nan=0.0)

    pup_mask_sa = rebin(pup_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    dm_mask_sa = rebin(dm_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    combined_mask_sa = (dm_mask_sa > 0.0) & (pup_mask_sa > 0.0)

    # 3. Apply mask
    gtilt_x = xp.where(combined_mask_sa[:, :, xp.newaxis] if is_3d \
              else combined_mask_sa, raw_gtilt_x, 0.0)
    gtilt_y = xp.where(combined_mask_sa[:, :, xp.newaxis] if is_3d \
              else combined_mask_sa, raw_gtilt_y, 0.0)

    # 4. Reshape to 2D matrix
    wfs_signal_x_2d = gtilt_x.reshape((-1, gtilt_x.shape[2] if is_3d else 1))
    wfs_signal_y_2d = gtilt_y.reshape((-1, gtilt_y.shape[2] if is_3d else 1))

    # 5. Select valid subapertures using provided indices
    if idx_valid_sa is not None:
        if specula_convention and len(idx_valid_sa.shape) > 1 and idx_valid_sa.shape[1] == 2:
            sa_2d = xp.zeros((wfs_nsubaps, wfs_nsubaps), dtype=float_dtype)
            sa_2d[idx_valid_sa[:, 0], idx_valid_sa[:, 1]] = 1
            sa_2d = xp.transpose(sa_2d)
            idx_temp = xp.where(sa_2d > 0)
            idx_valid_sa_new = xp.zeros_like(idx_valid_sa)
            idx_valid_sa_new[:, 0] = idx_temp[0]
            idx_valid_sa_new[:, 1] = idx_temp[1]
        else:
            idx_valid_sa_new = idx_valid_sa

        if len(idx_valid_sa_new.shape) > 1 and idx_valid_sa_new.shape[1] == 2:
            width = wfs_nsubaps
            linear_indices = idx_valid_sa_new[:, 0] * width + idx_valid_sa_new[:, 1]
            wfs_signal_x_2d = wfs_signal_x_2d[linear_indices.astype(xp.int32), :]
            wfs_signal_y_2d = wfs_signal_y_2d[linear_indices.astype(xp.int32), :]
        else:
            wfs_signal_x_2d = wfs_signal_x_2d[idx_valid_sa_new.astype(xp.int32), :]
            wfs_signal_y_2d = wfs_signal_y_2d[idx_valid_sa_new.astype(xp.int32), :]

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
        print(f'  ✓ G-tilt slopes formatted, shape: {im.shape}')

    return im


def compute_derivatives_with_extrapolation(data, mask=None):
    """
    Compute x and y derivatives using numpy.gradient on a 2D or 3D numpy array
    if mask is present does an extrapolation to avoid issue at the edges
    
    Parameters:
    - data: numpy 3D array
    - mask: optional, numpy 2D array, mask

    Returns:
    - dx: numpy 3D array, x derivative
    - dy: numpy 3D array, y derivative
    """

    if mask is not None:
        data, mask = _repair_interpolated_phase(data, mask)

    # Compute x derivative
    dx = xp.gradient(data, axis=(1), edge_order=1)

    # Compute y derivative
    dy = xp.gradient(data, axis=(0), edge_order=1)

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
    Used by both separated and combined workflows.
    """

    # Clean up masks
    if xp.isnan(pup_mask).any():
        xp.nan_to_num(pup_mask, copy=False, nan=0.0)
    if xp.isnan(dm_mask).any():
        xp.nan_to_num(dm_mask, copy=False, nan=0.0)

    # Rebin masks to WFS resolution
    pup_mask_sa = rebin(pup_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    pup_mask_sa = pup_mask_sa / xp.max(pup_mask_sa) if xp.max(pup_mask_sa) > 0 else pup_mask_sa

    dm_mask_sa = rebin(dm_mask, (wfs_nsubaps, wfs_nsubaps), method='sum')
    if xp.max(dm_mask_sa) <= 0:
        raise ValueError('DM mask is empty after rebinning.')
    dm_mask_sa = dm_mask_sa / xp.max(dm_mask_sa)

    # Clean derivatives
    if xp.isnan(derivatives_x).any():
        xp.nan_to_num(derivatives_x, copy=False, nan=0.0)
    if xp.isnan(derivatives_y).any():
        xp.nan_to_num(derivatives_y, copy=False, nan=0.0)

    # Apply pupil mask
    trans_der_x = apply_mask(derivatives_x, pup_mask, fill_value=xp.nan)
    trans_der_y = apply_mask(derivatives_y, pup_mask, fill_value=xp.nan)

    # Rebin derivatives
    # Since we use 'nanmean', the average is already correctly normalized by the valid area.
    scale_factor = trans_der_x.shape[0] / wfs_nsubaps

    wfs_signal_x = rebin(trans_der_x, (wfs_nsubaps, wfs_nsubaps), method='nanmean') * scale_factor
    wfs_signal_y = rebin(trans_der_y, (wfs_nsubaps, wfs_nsubaps), method='nanmean') * scale_factor

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
    if idx_valid_sa is not None:
        if specula_convention and len(idx_valid_sa.shape) > 1 and idx_valid_sa.shape[1] == 2:
            # *** sa_2d should use float_dtype (it's a mask with 0/1 values) ***
            sa_2d = xp.zeros((wfs_nsubaps, wfs_nsubaps), dtype=float_dtype)
            sa_2d[idx_valid_sa[:, 0], idx_valid_sa[:, 1]] = 1
            sa_2d = xp.transpose(sa_2d)
            idx_temp = xp.where(sa_2d > 0)
            # *** But idx_valid_sa_new should keep integer type (indices!) ***
            idx_valid_sa_new = xp.zeros_like(idx_valid_sa)  # Keep original dtype (int)
            idx_valid_sa_new[:, 0] = idx_temp[0]
            idx_valid_sa_new[:, 1] = idx_temp[1]
        else:
            idx_valid_sa_new = idx_valid_sa

        if len(idx_valid_sa_new.shape) > 1 and idx_valid_sa_new.shape[1] == 2:
            width = wfs_nsubaps
            linear_indices = idx_valid_sa_new[:, 0] * width + idx_valid_sa_new[:, 1]
            # *** Ensure indices are integers ***
            wfs_signal_x_2d = wfs_signal_x_2d[linear_indices.astype(xp.int32), :]
            wfs_signal_y_2d = wfs_signal_y_2d[linear_indices.astype(xp.int32), :]
        else:
            # *** Ensure indices are integers ***
            wfs_signal_x_2d = wfs_signal_x_2d[idx_valid_sa_new.astype(xp.int32), :]
            wfs_signal_y_2d = wfs_signal_y_2d[idx_valid_sa_new.astype(xp.int32), :]

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
        print(f'  ✓ Slopes computed, shape: {im.shape}')

    return im


def apply_dm_transformations_separated(pup_diam_m, pup_mask, dm_array, dm_mask,
                                       dm_height, dm_rotation,
                                       gs_pol_coo, gs_height,
                                       wfs_nsubaps=None,
                                       slope_method='derivatives',
                                       specula_convention=True,
                                       verbose=False):
    """
    Apply ONLY DM transformations (for separated workflow).
    Returns slopes/derivatives that need WFS transformations applied separately.
    """

    # *** Convert inputs to target device with correct dtype ***
    dm_array = to_xp(xp, dm_array, dtype=float_dtype)
    dm_mask = to_xp(xp, dm_mask, dtype=float_dtype)
    pup_mask = to_xp(xp, pup_mask, dtype=float_dtype)

    if specula_convention:
        dm_array = xp.transpose(dm_array, (1, 0, 2))
        dm_mask = xp.transpose(dm_mask)
        pup_mask = xp.transpose(pup_mask)

    pup_diam_pix = pup_mask.shape[0]
    pixel_pitch = pup_diam_m / pup_diam_pix

    if dm_mask.shape[0] != dm_array.shape[0]:
        raise ValueError('DM and mask arrays must have the same dimensions.')

    dm_translation, dm_magnification = shiftzoom_from_source_dm_params(
        gs_pol_coo, gs_height, dm_height, pixel_pitch
    )
    output_size = (pup_diam_pix, pup_diam_pix)

    if verbose:
        print(f'DM transformations (separated):')
        print(f'  Translation: {dm_translation} pixels')
        print(f'  Rotation: {dm_rotation} deg')
        print(f'  Magnification: {dm_magnification}')

    # Apply ONLY DM transformations
    trans_dm_array = rotshiftzoom_array(
        dm_array,
        dm_translation=dm_translation,
        dm_rotation=dm_rotation,
        dm_magnification=dm_magnification,
        wfs_translation=(0, 0),
        wfs_rotation=0,
        wfs_magnification=(1, 1),
        output_size=output_size
    )

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

    if xp.max(trans_dm_mask) <= 0:
        raise ValueError('Transformed DM mask is empty.')

    trans_dm_array = apply_mask(trans_dm_array, trans_dm_mask)

    # Compute derivatives or gtilt on DM-transformed array
    if slope_method == 'derivatives':
        derivatives_x, derivatives_y = compute_derivatives_with_extrapolation(
            trans_dm_array, mask=trans_dm_mask
        )
    elif slope_method == 'gtilt':
        derivatives_x, derivatives_y = compute_gtilt_with_extrapolation(
            trans_dm_array, mask=trans_dm_mask, wfs_nsubaps=wfs_nsubaps, verbose=verbose
        )
    else:
        raise ValueError(f"Unknown slope_method: {slope_method}")

    if verbose:
        print(f'  ✓ DM array transformed, shape: {trans_dm_array.shape}')
        print(f'  ✓ Slopes/Derivatives computed')

    # Return the ORIGINAL pupil mask (not transformed), so WFS transformations can be applied later
    return trans_dm_array, trans_dm_mask, pup_mask, derivatives_x, derivatives_y


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
    """

    # *** Compute WFS magnification including anamorphosis at 90° ***
    wfs_magnification = (wfs_mag_global, wfs_mag_global * wfs_anamorphosis_90)

    # *** Convert inputs to target device with correct dtype ***
    dm_array = to_xp(xp, dm_array, dtype=float_dtype)
    dm_mask = to_xp(xp, dm_mask, dtype=float_dtype)
    pup_mask = to_xp(xp, pup_mask, dtype=float_dtype)

    if specula_convention:
        dm_array = xp.transpose(dm_array, (1, 0, 2))
        dm_mask = xp.transpose(dm_mask)
        pup_mask = xp.transpose(pup_mask)
        wfs_translation_local = (-1*wfs_translation[1], -1*wfs_translation[0])
    else:
        wfs_translation_local = wfs_translation

    pup_diam_pix = pup_mask.shape[0]
    pixel_pitch = pup_diam_m / pup_diam_pix

    if dm_mask.shape[0] != dm_array.shape[0]:
        raise ValueError('DM and mask arrays must have the same dimensions.')

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

    # Apply ALL transformations in one step
    trans_dm_array = rotshiftzoom_array(
        dm_array,
        dm_translation=dm_translation,
        dm_rotation=dm_rotation,
        dm_magnification=dm_magnification,
        wfs_translation=wfs_translation_local,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45,
        output_size=output_size
    )

    # DM mask (only DM transformations)
    trans_dm_mask = rotshiftzoom_array(
        dm_mask,
        dm_translation=dm_translation,
        dm_rotation=dm_rotation,
        dm_magnification=dm_magnification,
        wfs_translation=wfs_translation_local,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45,
        output_size=output_size
    )
    trans_dm_mask[trans_dm_mask < 0.5] = 0

    # Pupil mask (only WFS transformations)
    trans_pup_mask = rotshiftzoom_array(
        pup_mask,
        dm_translation=(0, 0),
        dm_rotation=0,
        dm_magnification=(1, 1),
        wfs_translation=wfs_translation_local,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45,
        output_size=output_size
    )
    trans_pup_mask[trans_pup_mask < 0.5] = 0

    if xp.max(trans_dm_mask) <= 0:
        raise ValueError('Transformed DM mask is empty.')
    if xp.max(trans_pup_mask) <= 0:
        raise ValueError('Transformed pupil mask is empty.')

    trans_dm_array = apply_mask(trans_dm_array, trans_dm_mask)

    # Compute derivatives or gtilt on already-transformed array
    if slope_method == 'derivatives':
        derivatives_x, derivatives_y = compute_derivatives_with_extrapolation(
            trans_dm_array, mask=trans_dm_mask
        )
    elif slope_method == 'gtilt':
        derivatives_x, derivatives_y = compute_gtilt_with_extrapolation(
            trans_dm_array, mask=trans_dm_mask, wfs_nsubaps=wfs_nsubaps, verbose=verbose
        )
    else:
        raise ValueError(f"Unknown slope_method: {slope_method}")

    if verbose:
        print(f'  ✓ Combined transformation applied, shape: {trans_dm_array.shape}')
        print(f'  ✓ Slopes/Derivatives computed')

    return trans_dm_array, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y


def apply_wfs_transformations_separated(derivatives_x, derivatives_y,
                                        pup_mask, dm_mask,
                                        wfs_nsubaps, wfs_fov_arcsec,
                                        pup_diam_m, wfs_rotation,
                                        wfs_translation, wfs_mag_global,
                                        wfs_anamorphosis_90=1.0,
                                        wfs_anamorphosis_45=1.0,
                                        idx_valid_sa=None,
                                        slope_method='derivatives',
                                        specula_convention=True,
                                        verbose=False):
    """
    Apply WFS transformations to derivatives (for separated workflow).
    """

    output_size = pup_mask.shape

    # *** Compute WFS magnification including anamorphosis at 90° ***
    wfs_magnification = (wfs_mag_global, wfs_mag_global * wfs_anamorphosis_90)

    if specula_convention:
        wfs_translation_local = (-1*wfs_translation[1], -1*wfs_translation[0])
    else:
        wfs_translation_local = wfs_translation

    # Transform pupil mask
    trans_pup_mask = rotshiftzoom_array(
        pup_mask,
        dm_translation=(0, 0),
        dm_rotation=0,
        dm_magnification=(1, 1),
        wfs_translation=wfs_translation_local,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45,
        output_size=output_size
    )
    trans_pup_mask[trans_pup_mask < 0.5] = 0

    if xp.max(trans_pup_mask) <= 0:
        raise ValueError('Transformed pupil mask is empty.')

    if verbose:
        print(f'WFS transformations (separated):')
        print(f'  Translation: {wfs_translation} pixels')
        print(f'  Rotation: {wfs_rotation} deg')
        print(f'  Magnification: {wfs_magnification}')

    # Transform derivatives
    trans_der_x = rotshiftzoom_array(
        derivatives_x,
        dm_translation=(0, 0),
        dm_rotation=0,
        dm_magnification=(1, 1),
        wfs_translation=wfs_translation_local,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45,
        output_size=output_size
    )

    trans_der_y = rotshiftzoom_array(
        derivatives_y,
        dm_translation=(0, 0),
        dm_rotation=0,
        dm_magnification=(1, 1),
        wfs_translation=wfs_translation_local,
        wfs_rotation=wfs_rotation,
        wfs_magnification=wfs_magnification,
        wfs_anamorphosis_45=wfs_anamorphosis_45,
        output_size=output_size
    )

    # Route to the appropriate formatter
    if slope_method == 'derivatives':
        return _compute_slopes_from_derivatives(
            trans_der_x, trans_der_y, trans_pup_mask, dm_mask,
            wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa,
            verbose, specula_convention
        )
    elif slope_method == 'gtilt':
        return _compute_slopes_from_gtilt(
            trans_der_x, trans_der_y, trans_pup_mask, dm_mask,
            wfs_nsubaps, wfs_fov_arcsec, pup_diam_m,
            idx_valid_sa, verbose, specula_convention
        )
    else:
        raise ValueError(f"Unknown slope_method: {slope_method}")


def apply_wfs_transformations_combined(derivatives_x, derivatives_y, trans_pup_mask, dm_mask,
                                       wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa=None,
                                       slope_method='derivatives', specula_convention=True,
                                       verbose=False):
    """
    Compute slopes from pre-transformed derivatives (for combined workflow).
    No additional transformations needed.
    """

    # Derivatives are already transformed - route to the appropriate formatter
    if slope_method == 'derivatives':
        return _compute_slopes_from_derivatives(
            derivatives_x, derivatives_y, trans_pup_mask, dm_mask,
            wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa,
            verbose, specula_convention
        )
    elif slope_method == 'gtilt':
        return _compute_slopes_from_gtilt(
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
                       specula_convention=True, verbose=False,
                       display=False):
    """
    Computes interaction matrix using intelligent workflow selection.
    
    Automatically chooses between:
    - Separated workflow: When transformations are only in DM OR WFS (2 interpolation steps)
    - Combined workflow: When both DM and WFS have transformations (1 interpolation step)
    """

    # *** Convert idx_valid_sa if provided ***
    if idx_valid_sa is not None:
        idx_valid_sa = to_xp(xp, idx_valid_sa)

    # Detect which transformations are present
    has_dm_transform = has_transformations(dm_rotation, (0, 0), (1, 1)) or \
                       gs_pol_coo != (0, 0) or dm_height != 0
    has_wfs_transform = has_transformations(wfs_rotation, wfs_translation, wfs_mag_global)

    # Choose workflow
    use_combined = has_dm_transform and has_wfs_transform and slope_method == 'derivatives'

    if verbose:
        print(f"\n{'='*60}")
        print(f"Interaction Matrix Computation")
        print(f"{'='*60}")
        print(f"DM transformations: {has_dm_transform}")
        print(f"WFS transformations: {has_wfs_transform}")
        print(f"Using {'COMBINED' if use_combined else 'SEPARATED'} workflow")
        if wfs_anamorphosis_45 != 1.0:
            print(f"- Note: WFS anamorphosis at 45° is set to {wfs_anamorphosis_45}, "
                  f"which will be applied in both workflows.")
        print(f"{'='*60}\n")

    if use_combined:
        # Combined workflow: single interpolation
        trans_dm_array, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
            apply_dm_transformations_combined(
                pup_diam_m=pup_diam_m,
                pup_mask=pup_mask,
                dm_array=dm_array,
                dm_mask=dm_mask,
                gs_pol_coo=gs_pol_coo,
                gs_height=gs_height,
                dm_height=dm_height,
                dm_rotation=dm_rotation,
                wfs_rotation=wfs_rotation,
                wfs_translation=wfs_translation,
                wfs_nsubaps=wfs_nsubaps,
                slope_method=slope_method,
                wfs_mag_global=wfs_mag_global,
                wfs_anamorphosis_90=wfs_anamorphosis_90,
                wfs_anamorphosis_45=wfs_anamorphosis_45,
                verbose=verbose,
                specula_convention=specula_convention
            )

        im = apply_wfs_transformations_combined(
            derivatives_x, derivatives_y, trans_pup_mask, trans_dm_mask,
            wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa=idx_valid_sa,
            slope_method=slope_method, verbose=verbose, specula_convention=specula_convention
        )
    else:
        # Separated workflow: two interpolation steps (more flexible)
        trans_dm_array, trans_dm_mask, pup_mask_conv, derivatives_x, derivatives_y = \
            apply_dm_transformations_separated(
                pup_diam_m, pup_mask, dm_array, dm_mask,
                dm_height, dm_rotation,
                gs_pol_coo, gs_height,
                wfs_nsubaps=wfs_nsubaps, slope_method=slope_method,
                verbose=verbose, specula_convention=specula_convention
            )

        im = apply_wfs_transformations_separated(
            derivatives_x=derivatives_x,
            derivatives_y=derivatives_y,
            pup_mask=pup_mask_conv,
            dm_mask=trans_dm_mask,
            wfs_nsubaps=wfs_nsubaps,
            wfs_fov_arcsec=wfs_fov_arcsec,
            pup_diam_m=pup_diam_m,
            wfs_rotation=wfs_rotation,
            wfs_translation=wfs_translation,
            wfs_mag_global=wfs_mag_global,
            wfs_anamorphosis_90=wfs_anamorphosis_90,
            wfs_anamorphosis_45=wfs_anamorphosis_45,
            idx_valid_sa=idx_valid_sa,
            slope_method=slope_method,
            verbose=verbose, specula_convention=specula_convention
        )

    if display:
        idx_plot = [2, 5]
        pup_mask_cpu = cpuArray(pup_mask_conv)
        trans_dm_mask_cpu = cpuArray(trans_dm_mask)
        trans_dm_array_cpu = cpuArray(trans_dm_array)
        fig, axs = plt.subplots(2, 2)
        im3 = axs[0, 0].imshow(pup_mask_cpu, cmap='seismic')
        axs[0, 1].imshow(trans_dm_mask_cpu, cmap='seismic')
        axs[1, 0].imshow(trans_dm_array_cpu[:, :, idx_plot[0]], cmap='seismic')
        axs[1, 1].imshow(trans_dm_array_cpu[:, :, idx_plot[1]], cmap='seismic')
        fig.suptitle(f'Mask, DM mask, DM shapes (modes {idx_plot[0]} and {idx_plot[1]})')
        fig.colorbar(im3, ax=axs.ravel().tolist(), fraction=0.02)
        plt.show()

    return im


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
    - slope_method: str, 'derivatives' or 'gtilt'
    - specula_convention: bool, optional
    - im_on_cpu: bool, optional, force output interaction matrices on CPU
    - minimize_memory: bool, optional, delete intermediate variables to save memory
    - verbose: bool, optional
    
    Returns:
    - im_dict: dict, interaction matrices keyed by WFS name or index
    - derivatives_info: dict with metadata about the computation
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

    # Extract gs_pol_coo and gs_height for each WFS
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
                                 f"specified in wfs_config when gs_height=None")

            wfs_gs_pol_coo = wfs_config['gs_pol_coo']
            wfs_gs_height = wfs_config['gs_height']

        wfs_gs_info.append((wfs_gs_pol_coo, wfs_gs_height))

    # Check if all WFS see DM from the same direction
    all_gs_same = all(gs_info == wfs_gs_info[0] for gs_info in wfs_gs_info)

    # Detect WFS transformations
    wfs_transforms = []
    has_wfs_transform_list = []

    for config in wfs_configs:
        wfs_rot = config.get('rotation', 0.0)
        wfs_trans = config.get('translation', (0.0, 0.0))
        wfs_mag = config.get('magnification', (1.0, 1.0))
        wfs_transforms.append((wfs_rot, wfs_trans, wfs_mag))
        has_wfs_transform_list.append(
            has_transformations(wfs_rot, wfs_trans, wfs_mag)
        )

    # Check if all WFS transforms are identical
    all_wfs_same = all(t == wfs_transforms[0] for t in wfs_transforms)

    # Check if ANY WFS has transformations
    any_wfs_transform = any(has_wfs_transform_list)

    # Detect DM transformations (relative to first GS position)
    gs_pol_coo_ref, gs_height_ref = wfs_gs_info[0]
    has_dm_transform = has_transformations(dm_rotation, (0, 0), (1, 1)) or \
                       gs_pol_coo_ref != (0, 0) or gs_height_ref != 0 or dm_height != 0

    # *** Improved workflow decision logic ***
    # SEPARATED workflow can be used when:
    # 1. All WFS see DM from same direction AND have same transforms
    # 2. AND we can compute derivatives once and reuse them
    #    This is possible when:
    #    - ONLY DM transforms (no WFS transforms)
    #    - OR ONLY WFS transforms (no DM transforms)
    #    - NOT BOTH (that would require COMBINED)

    # XOR logic: exactly one of them has transforms, not both
    can_separate_transforms = (has_dm_transform and not any_wfs_transform) or \
                             (not has_dm_transform and any_wfs_transform)

    use_separated = all_gs_same and all_wfs_same and can_separate_transforms \
                and slope_method == 'derivatives'

    if verbose:
        print(f"All WFS see DM from same direction: {all_gs_same}")
        print(f"All WFS have same transforms: {all_wfs_same}")
        print(f"DM transformations present: {has_dm_transform}")
        print(f"WFS transformations present: {any_wfs_transform}")
        print(f"Can separate transforms (XOR): {can_separate_transforms}")

        if not all_gs_same:
            print(f"  → Different GS positions require per-WFS computation")
        if not all_wfs_same:
            print(f"  → Different WFS transforms require per-WFS computation")
        if not can_separate_transforms:
            print(f"  → Both DM and WFS transforms present")

        print(f"Using {'SEPARATED' if use_separated else 'COMBINED (per-WFS)'} workflow")
        print(f"{'='*60}\n")

    im_dict = {}
    derivatives_info = {
        'workflow': 'separated' if use_separated else 'combined',
        'all_gs_same': all_gs_same,
        'all_wfs_same': all_wfs_same,
        'has_dm_transform': has_dm_transform,
        'any_wfs_transform': any_wfs_transform,
        'can_separate': can_separate_transforms
    }

    if use_separated:
        # SEPARATED WORKFLOW: Compute transformations once, then reuse
        if verbose:
            print("[SEPARATED WORKFLOW]")
            print("  Step 1/2: Computing DM transformations ONCE...")

        # Use first WFS's gs_pol_coo and gs_height (they're all the same)
        gs_pol_coo_ref, gs_height_ref = wfs_gs_info[0]
        # Needs the wfs_nsubaps of the first WFS if it's identical
        wfs_nsubaps_ref = wfs_configs[0]['nsubaps']

        trans_dm_array, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
            apply_dm_transformations_separated(
                pup_diam_m, pup_mask, dm_array, dm_mask,
                dm_height, dm_rotation,
                gs_pol_coo_ref, gs_height_ref,
                wfs_nsubaps=wfs_nsubaps_ref, slope_method=slope_method,
                verbose=verbose,
                specula_convention=specula_convention
            )

        if verbose:
            print(f"  ✓ DM transformed: {trans_dm_array.shape}")
            print(f"  ✓ Slopes/Derivatives computed: {derivatives_x.shape}")
            if any_wfs_transform:
                print(f"\n[SEPARATED] Step 2/2: Applying WFS transformations to each WFS...")
            else:
                print(f"\n[SEPARATED] Step 2/2: Computing slopes for each WFS...")

        # Store derivatives for potential reuse
        derivatives_info['derivatives_x'] = derivatives_x
        derivatives_info['derivatives_y'] = derivatives_y
        derivatives_info['dm_array'] = trans_dm_array
        derivatives_info['dm_mask'] = trans_dm_mask

        # Apply WFS transformations for each WFS
        for i, wfs_config in enumerate(wfs_configs):
            wfs_name = wfs_config.get('name', f'wfs_{i}')
            wfs_nsubaps = wfs_config['nsubaps']
            wfs_rotation = wfs_config.get('rotation', 0.0)
            wfs_translation = wfs_config.get('translation', (0.0, 0.0))
            wfs_magnification = wfs_config.get('magnification', (1.0, 1.0))
            if isinstance(wfs_magnification, (tuple, list)) and len(wfs_magnification) == 2:
                wfs_mag_global = xp.sqrt(wfs_magnification[0] * wfs_magnification[1])
                wfs_anamorphosis_90 = wfs_magnification[1] / wfs_magnification[0] \
                    if wfs_magnification[0] != 0 else 1.0
            else:
                wfs_mag_global = float(wfs_magnification)
                wfs_anamorphosis_90 = 1.0
            wfs_anamorphosis_45 = wfs_config.get('anamorphosis_45', 1.0)
            wfs_fov_arcsec = wfs_config['fov_arcsec']
            idx_valid_sa = wfs_config.get('idx_valid_sa', None)

            if verbose:
                print(f"\n  [{i+1}/{len(wfs_configs)}] {wfs_name}:")
                print(f"    Subapertures: {wfs_nsubaps}x{wfs_nsubaps}")
                print(f"    FOV: {wfs_fov_arcsec}''")

            im = apply_wfs_transformations_separated(
                derivatives_x=derivatives_x,
                derivatives_y=derivatives_y,
                pup_mask=trans_pup_mask,
                dm_mask=trans_dm_mask,
                wfs_nsubaps=wfs_nsubaps,
                wfs_fov_arcsec=wfs_fov_arcsec,
                pup_diam_m=pup_diam_m,
                wfs_rotation=wfs_rotation,
                wfs_translation=wfs_translation,
                wfs_mag_global=wfs_mag_global,
                wfs_anamorphosis_90=wfs_anamorphosis_90,
                wfs_anamorphosis_45=wfs_anamorphosis_45,
                idx_valid_sa=idx_valid_sa,
                slope_method=slope_method,
                verbose=False,  # Suppress inner verbose
                specula_convention=specula_convention
            )

            im_dict[wfs_name] = im

            if verbose:
                print(f"    ✓ IM shape: {im.shape}")

    else:
        # COMBINED WORKFLOW: Each WFS computed independently
        # But within each WFS, use SEPARATED if XOR condition is met
        if verbose:
            print("[COMBINED WORKFLOW - Per-WFS Computation]")
            if not all_gs_same:
                print("  Reason: Different guide star positions")
            if not all_wfs_same:
                print("  Reason: Different WFS transformations")
            print()

        for i, wfs_config in enumerate(wfs_configs):
            wfs_name = wfs_config.get('name', f'wfs_{i}')
            wfs_nsubaps = wfs_config['nsubaps']
            wfs_rotation = wfs_config.get('rotation', 0.0)
            wfs_translation = wfs_config.get('translation', (0.0, 0.0))
            wfs_magnification = wfs_config.get('magnification', (1.0, 1.0))
            if isinstance(wfs_magnification, (tuple, list)) and len(wfs_magnification) == 2:
                wfs_mag_global = xp.sqrt(wfs_magnification[0] * wfs_magnification[1])
                wfs_anamorphosis_90 = wfs_magnification[1] / wfs_magnification[0] \
                    if wfs_magnification[0] != 0 else 1.0
            else:
                wfs_mag_global = float(wfs_magnification)
                wfs_anamorphosis_90 = 1.0
            wfs_anamorphosis_45 = wfs_config.get('anamorphosis_45', 1.0)
            wfs_fov_arcsec = wfs_config['fov_arcsec']
            idx_valid_sa = wfs_config.get('idx_valid_sa', None)

            # Get WFS-specific gs_pol_coo and gs_height
            gs_pol_coo_wfs, gs_height_wfs = wfs_gs_info[i]

            # *** For this WFS, determine if we can use SEPARATED workflow ***
            has_dm_transform_wfs = has_transformations(dm_rotation, (0, 0), (1, 1)) or \
                                   gs_pol_coo_wfs != (0, 0) or gs_height_wfs != 0 or dm_height != 0
            has_wfs_transform_wfs = has_transformations(
                wfs_rotation, wfs_translation, wfs_magnification
            )

            # XOR: only one type of transform
            use_separated_this_wfs = ((has_dm_transform_wfs and not has_wfs_transform_wfs) or \
                                    (not has_dm_transform_wfs and has_wfs_transform_wfs)) \
                                    and slope_method == 'derivatives'  # Can only reuse derivatives if using derivatives method

            if verbose:
                print(f"  [{i+1}/{len(wfs_configs)}] {wfs_name}:")
                print(f"    Subapertures: {wfs_nsubaps}x{wfs_nsubaps}, FOV:"
                      f" {wfs_fov_arcsec}''")
                print(f"    GS: {gs_pol_coo_wfs}, height: {gs_height_wfs} m")
                print(f"    DM transforms: {has_dm_transform_wfs}, WFS transforms:"
                      f" {has_wfs_transform_wfs}")
                print(f"    → Using {'SEPARATED' if use_separated_this_wfs else 'COMBINED'}"
                      f" for this WFS")

            # Call the appropriate workflow function for THIS WFS
            if use_separated_this_wfs:
                # SEPARATED: Two interpolation steps
                trans_dm_array, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
                    apply_dm_transformations_separated(
                        pup_diam_m=pup_diam_m,
                        pup_mask=pup_mask,
                        dm_array=dm_array,
                        dm_mask=dm_mask,
                        dm_height=dm_height,
                        dm_rotation=dm_rotation,
                        gs_pol_coo=gs_pol_coo_wfs,
                        gs_height=gs_height_wfs,
                        wfs_nsubaps=wfs_nsubaps, slope_method=slope_method,
                        verbose=False,
                        specula_convention=specula_convention
                    )

                im = apply_wfs_transformations_separated(
                    derivatives_x=derivatives_x,
                    derivatives_y=derivatives_y,
                    pup_mask=trans_pup_mask,
                    dm_mask=trans_dm_mask,
                    wfs_nsubaps=wfs_nsubaps,
                    wfs_fov_arcsec=wfs_fov_arcsec,
                    pup_diam_m=pup_diam_m,
                    wfs_rotation=wfs_rotation,
                    wfs_translation=wfs_translation,
                    wfs_mag_global=wfs_mag_global,
                    wfs_anamorphosis_90=wfs_anamorphosis_90,
                    wfs_anamorphosis_45=wfs_anamorphosis_45,
                    idx_valid_sa=idx_valid_sa,
                    slope_method=slope_method,
                    verbose=False,
                    specula_convention=specula_convention
                )
            else:
                # COMBINED: Single interpolation step
                trans_dm_array, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
                    apply_dm_transformations_combined(
                        pup_diam_m=pup_diam_m,
                        pup_mask=pup_mask,
                        dm_array=dm_array,
                        dm_mask=dm_mask,
                        dm_height=dm_height,
                        dm_rotation=dm_rotation,
                        gs_pol_coo=gs_pol_coo_wfs,
                        gs_height=gs_height_wfs,
                        wfs_rotation=wfs_rotation,
                        wfs_translation=wfs_translation,
                        wfs_nsubaps=wfs_nsubaps, slope_method=slope_method,
                        wfs_mag_global=wfs_mag_global,
                        wfs_anamorphosis_90=wfs_anamorphosis_90,
                        wfs_anamorphosis_45=wfs_anamorphosis_45,
                        verbose=False,
                        specula_convention=specula_convention
                    )

                im = apply_wfs_transformations_combined(
                    derivatives_x, derivatives_y, trans_pup_mask, trans_dm_mask,
                    wfs_nsubaps, wfs_fov_arcsec, pup_diam_m, idx_valid_sa=idx_valid_sa,
                    slope_method=slope_method, verbose=False,
                    specula_convention=specula_convention
                )

            if minimize_memory:
                # free as much memory as possible
                del trans_dm_array, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y

            if im_on_cpu:
                #  move im to CPU
                im_dict[wfs_name] = cpuArray(im)
            else:
                im_dict[wfs_name] = im

            if verbose:
                print(f"    ✓ IM shape: {im.shape}")

    display = False  # Disable display for multi-WFS case
                     # Please note that it is not compatible with minimize_memory=True
                     # as the variables would be deleted before plotting
    if display:
        idx_plot = [0, 2, 5]
        trans_dm_array_cpu = cpuArray(trans_dm_array)
        derivatives_x_cpu = cpuArray(apply_mask(derivatives_x, trans_pup_mask, fill_value=xp.nan))
        derivatives_y_cpu = cpuArray(apply_mask(derivatives_y, trans_pup_mask, fill_value=xp.nan))
        for idx in idx_plot:
            # 3 plots: DM shape, derivative x, derivative y
            fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            im3 = axs[0].imshow(trans_dm_array_cpu[:, :, idx], cmap='seismic')
            axs[1].imshow(derivatives_x_cpu[:, :, idx], cmap='seismic')
            axs[2].imshow(derivatives_y_cpu[:, :, idx], cmap='seismic')
            fig.suptitle(f'DM shapes (modes {idx_plot[0]} and {idx_plot[1]})')
            fig.colorbar(im3, ax=axs.ravel().tolist(), fraction=0.02)
        plt.show()

    if verbose:
        print(f"\n{'='*60}")
        print(f"Completed {len(im_dict)} interaction matrices")
        print(f"Overall workflow: {derivatives_info['workflow'].upper()}")
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
    - idx_valid_sa: array, indices of valid subapertures
    - specula_convention: bool, whether to use SPECULA convention
    - verbose: bool, whether to print information
    
    Returns:
    - illumination: 1D array, relative illumination of each valid subaperture (normalized to max=1)
    """

    # *** Convert to target device ***
    pup_mask = to_xp(xp, pup_mask, dtype=float_dtype)
    if idx_valid_sa is not None:
        idx_valid_sa = to_xp(xp, idx_valid_sa)

    if specula_convention:
        pup_mask = xp.transpose(pup_mask)

    output_size = pup_mask.shape

    # Apply WFS transformations to pupil mask
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
    trans_pup_mask[trans_pup_mask < 0.5] = 0

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
    if idx_valid_sa is not None:
        if specula_convention and len(idx_valid_sa.shape) > 1 and idx_valid_sa.shape[1] == 2:
            # Convert SPECULA format indices
            sa_2d = xp.zeros((wfs_nsubaps, wfs_nsubaps), dtype=float_dtype)
            sa_2d[idx_valid_sa[:, 0], idx_valid_sa[:, 1]] = 1
            sa_2d = xp.transpose(sa_2d)
            idx_temp = xp.where(sa_2d > 0)
            idx_valid_sa_new = xp.zeros_like(idx_valid_sa)
            idx_valid_sa_new[:, 0] = idx_temp[0]
            idx_valid_sa_new[:, 1] = idx_temp[1]
        else:
            idx_valid_sa_new = idx_valid_sa

        if len(idx_valid_sa_new.shape) > 1 and idx_valid_sa_new.shape[1] == 2:
            width = wfs_nsubaps
            linear_indices = idx_valid_sa_new[:, 0] * width + idx_valid_sa_new[:, 1]
            illumination = illumination_2d[linear_indices.astype(xp.int32)]
        else:
            illumination = illumination_2d[idx_valid_sa_new.astype(xp.int32)]
    else:
        # Use all subapertures
        illumination = illumination_2d

    # Convert to CPU for return
    illumination = cpuArray(illumination)

    if verbose:
        print(f"Subaperture illumination statistics:")
        print(f"  Min: {xp.min(illumination):.3f}")
        print(f"  Max: {xp.max(illumination):.3f}")
        print(f"  Mean: {xp.mean(illumination):.3f}")
        print(f"  Std: {xp.std(illumination):.3f}")

    return illumination
