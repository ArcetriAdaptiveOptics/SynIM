"""
Geometry shared by the two implementations of the interaction matrix
(synim.reference and synim.sparse_operator): parameters of the combined DM +
WFS transformation, transformed masks, selection of the valid subapertures
and conversion to slope units.

Both implementations use these functions, so an error here is not detected
by the test that checks that the two implementations agree
(test/test_sparse_operator.py): it is covered by the physical tests of the
reference implementation (see "Reference and optimized implementations" in
docs/general.rst).
"""
import warnings

import numpy as np
import scipy.ndimage

import synim as _synim
_synim._require_init(__name__)
from synim import xp, cpuArray, float_dtype
from synim.utils import rotshiftzoom_matrix, shiftzoom_from_source_dm_params


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


def _slope_units(pup_diam_m, wfs_nsubaps, wfs_fov_arcsec):
    """
    Slope units: factor from the phase difference across one subaperture
    (in nm: the mean difference per pixel times the number of pixels per
    subaperture) to the tilt in arcsec, normalized to half the WFS field of
    view.
    """
    coeff = 1e-9 / (pup_diam_m / wfs_nsubaps) * 206265
    coeff *= 1 / (0.5 * wfs_fov_arcsec)
    return coeff


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
    float64 (scipy.ndimage.affine_transform), so that the masks, and the
    operators that depend on them through their thresholds, are the same on
    CPU and GPU (cupyx affine_transform uses single precision coordinates
    for single precision arrays). Returns a float_dtype numpy array: the
    operators are built on CPU, only the explicit path needs the masks on
    the target device.
    """
    transform = dict(transform)
    output_size = transform.pop('output_size')
    mask = np.nan_to_num(cpuArray(mask).astype(np.float64), nan=0.0)
    matrix, offset = rotshiftzoom_matrix(mask.shape, output_size, **transform)
    trans_mask = scipy.ndimage.affine_transform(
        mask, cpuArray(matrix), offset=cpuArray(offset), output_shape=output_size, order=1)
    trans_mask[trans_mask < threshold] = 0
    return trans_mask.astype(float_dtype)


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
    The masks are numpy arrays (see _transform_mask).
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

    if np.max(trans_dm_mask) <= 0:
        raise ValueError('Transformed DM mask is empty.')
    if np.max(trans_pup_mask) <= 0:
        raise ValueError('Transformed pupil mask is empty.')

    return transform, trans_dm_mask, trans_pup_mask

