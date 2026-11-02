"""
Optimized implementation of the interaction matrix: the whole computation
as one sparse linear operator applied to the DM cube.

The masks do not depend on the DM mode, and once they are fixed every step
of the reference implementation (synim.reference) is linear in the phase
phi of a mode (raveled to n_pix pixels), so the whole chain is one sparse
matrix A:

    s = C . P_valid . M_sa . G . R . D_mask . T . phi = A . phi

T        bilinear interpolation of the combined transformation
         (_transformation_operator)
D_mask   multiplication by the transformed DM mask
R        edge repair (_repair_operator; derivatives only)
G        derivatives, or differences of adjacent pixels (telescoping sum)
M_sa     mean over each, possibly fractional, subaperture
         (D_mask, G and M_sa: _slope_operators)
P_valid, C  valid subapertures, x/y order, slope units (_wfs_operator)

With the DM cube seen as Phi (n_pix x n_modes, a reshape), the interaction
matrix is the single product IM = A . Phi (_apply_operator). See
"Computation Workflows" in docs/general.rst.

This module is a port of synim.reference, which defines the physics:
changes are made there first and then ported here, and the equivalence
test (test/test_sparse_operator.py) must pass. See "Reference and optimized
implementations" in docs/general.rst. The geometry shared with the
reference implementation is in synim.geometry.
"""
import functools

import numpy as np
import scipy.sparse

import synim as _synim
_synim._require_init(__name__)
from synim import xp, cpuArray, to_xp, float_dtype
from synim.utils import (
    _rebin_matrix_np,
    affine_interpolation_matrix,
    rotshiftzoom_matrix,
    calculate_extrapolation_indices_coeffs
)
from synim.geometry import (
    _valid_sa_indices,
    _slope_units,
    _check_dm_mask_shape,
    _combined_transformation,
)


@functools.lru_cache(maxsize=8)
def _difference_operators(n0, n1, slope_method):
    """
    Sparse x and y operators on an n0 x n1 grid: derivatives as _gradient
    (central differences, one-sided at the edges) or, for the telescoping
    sum, differences of adjacent pixels. They depend only on the grid size,
    so they are cached: all the WFS of a calibration share them. The
    returned matrices are shared by all the callers: do not modify them in
    place (e.g. sort_indices, eliminate_zeros, data *= ...).
    """
    def one_axis(n):
        if slope_method == 'telsum':
            return scipy.sparse.diags([-1.0, 1.0], [0, 1], shape=(n - 1, n))
        gradient = scipy.sparse.diags([-0.5, 0.5], [-1, 1], shape=(n, n), format='lil')
        gradient[0, :2] = [-1, 1]
        gradient[n - 1, n - 2:] = [-1, 1]
        return gradient

    return (scipy.sparse.kron(scipy.sparse.identity(n0), one_axis(n1), format='csr'),
            scipy.sparse.kron(one_axis(n0), scipy.sparse.identity(n1), format='csr'))


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
    the module docstring): sparse (wfs_nsubaps**2, n_pixels) float64 operators
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
        gradient_x, gradient_y = _difference_operators(n0, n1, slope_method)
        return mean @ gradient_x @ repair, mean @ gradient_y @ repair

    if slope_method == 'telsum':
        if n0 < 2 * wfs_nsubaps:
            raise ValueError(f'The telescoping sum needs at least 2 pixels per subaperture:'
                             f' {n0} pixels, {wfs_nsubaps} subapertures.')
        # see _telsum_fractional and _compute_slopes_from_telsum
        difference_x, difference_y = _difference_operators(n0, n1, slope_method)

        def operator(axis):
            if axis == 1:
                used = in_pupil[:, 1:] * in_pupil[:, :-1]
                valid = strict_mask[:, 1:] & strict_mask[:, :-1]
                row_matrix, col_matrix = rows, np.minimum(cols[:, 1:], cols[:, :-1])
                difference = difference_x
            else:
                used = in_pupil[1:] * in_pupil[:-1]
                valid = strict_mask[1:] & strict_mask[:-1]
                row_matrix, col_matrix = np.minimum(rows[:, 1:], rows[:, :-1]), cols
                difference = difference_y
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
    coeff = _slope_units(pup_diam_m, wfs_nsubaps, wfs_fov_arcsec)
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
                       specula_convention=True, verbose=False):
    """
    Interaction matrix computed as a sparse operator applied to the DM
    array. Same parameters and result as synim.synim.interaction_matrix
    with method='operator'.
    """
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
    return im


def interaction_matrices(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                         wfs_params, groups, slope_method, specula_convention,
                         verbose=False):
    """
    Interaction matrices of several WFS (see
    synim.synim.interaction_matrices_multi_wfs, method='operator').
    wfs_params is the list of the WFS parameter dicts and groups the lists
    of indices of the WFS that share the transformation and the slope
    operators, both built by interaction_matrices_multi_wfs. The operators
    are built once per group, and the operators of all the WFS are applied
    to the DM array with a single product. Returns the list of the
    interaction matrices, in the order of wfs_params.
    """
    operators = [None] * len(wfs_params)
    for members in groups:
        first = wfs_params[members[0]]
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
    return [im_all[start:stop] for start, stop in zip(bounds[:-1], bounds[1:])]
