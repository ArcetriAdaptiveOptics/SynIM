"""
Interaction matrices of Shack-Hartmann WFS, synthesized from the DM
influence functions (public API).

The interaction matrix has two implementations, selected with the method
argument of interaction_matrix and interaction_matrices_multi_wfs:

- 'operator' (default): the whole computation as one sparse linear operator
  applied to the DM cube (synim.sparse_operator). Fast, and it allocates no
  transformed copy of the cube.
- 'reference': the explicit, step-by-step computation on the DM cube
  (synim.reference). It defines the physics, and its intermediate arrays
  (transformed phase, derivatives) can be inspected for debugging.

The geometry shared by the two implementations is in synim.geometry. Every
change to the computation is made first in the reference implementation,
then ported to the optimized one; test/test_sparse_operator.py checks that
the two agree, and the physical tests check the reference. See "Reference
and optimized implementations" in docs/general.rst.
"""
import math

import numpy as np

import synim as _synim
_synim._require_init(__name__)
from synim import xp, cpuArray, to_xp, float_dtype
import matplotlib.pyplot as plt
from synim.utils import rebin, rotshiftzoom_array
from synim import reference, sparse_operator

# The functions below were defined in this module before the reference and
# optimized implementations were separated: they are imported here so that
# synim.synim.<name> keeps working.
from synim.geometry import (  # noqa: F401
    _without_nan,
    _warn_if_nan,
    _valid_sa_indices,
    _check_dm_mask_shape,
    _transform_mask,
    _combined_transformation,
)
from synim.reference import (  # noqa: F401
    _repair_interpolated_phase,
    compute_telsum_with_extrapolation,
    _telsum_fractional,
    _compute_slopes_from_telsum,
    _gradient,
    compute_derivatives_with_extrapolation,
    _compute_slopes_from_derivatives,
    apply_dm_transformations_combined,
    apply_wfs_transformations_combined,
)
from synim.sparse_operator import (  # noqa: F401
    _difference_operators,
    _repair_operator,
    _slope_operators,
    _wfs_operator,
    _transformation_operator,
    _CSR_DENSE_KERNEL_CODE,
    _apply_operator,
)

# Implementations of the interaction matrix, selected by the method argument
METHODS = {
    'operator': sparse_operator,
    'reference': reference,
}


def _implementation(method):
    """Module implementing the interaction matrix with the given method."""
    try:
        return METHODS[method]
    except (KeyError, TypeError):
        raise ValueError(f"Unknown method: {method!r}, expected one of"
                         f" {', '.join(repr(m) for m in METHODS)}") from None


def _display_transformed_dm(pup_diam_m, pup_mask, dm_array, dm_mask, transformation_args,
                            specula_convention, idx_plot=(2, 5)):
    """
    Plot the transformed pupil and DM masks and two transformed and masked
    DM modes (display option of interaction_matrix).
    """
    idx_plot = list(idx_plot)
    transform, trans_dm_mask, trans_pup_mask = _combined_transformation(
        pup_diam_m, pup_mask, dm_mask, *transformation_args, specula_convention)
    modes = dm_array[:, :, idx_plot]
    if specula_convention:
        modes = xp.transpose(modes, (1, 0, 2))
    trans_modes = cpuArray(rotshiftzoom_array(modes, **transform)) \
        * trans_dm_mask[:, :, np.newaxis]
    fig, axs = plt.subplots(2, 2)
    im3 = axs[0, 0].imshow(trans_pup_mask, cmap='seismic')
    axs[0, 1].imshow(trans_dm_mask, cmap='seismic')
    axs[1, 0].imshow(trans_modes[:, :, 0], cmap='seismic')
    axs[1, 1].imshow(trans_modes[:, :, 1], cmap='seismic')
    fig.suptitle(f'Mask, DM mask, DM shapes (modes {idx_plot[0]} and {idx_plot[1]})')
    fig.colorbar(im3, ax=axs.ravel().tolist(), fraction=0.02)
    plt.show()


def interaction_matrix(pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
                       wfs_nsubaps, wfs_fov_arcsec, gs_pol_coo, gs_height,
                       wfs_rotation, wfs_translation, wfs_mag_global,
                       wfs_anamorphosis_90=1.0, wfs_anamorphosis_45=1.0,
                       idx_valid_sa=None, slope_method='derivatives',
                       specula_convention=True, verbose=False,
                       display=False, method='operator'):
    """
    Computes an interaction matrix.

    All the geometric transformations (DM shift, rotation and magnification
    from the guide star direction and height, WFS shift, rotation,
    magnification and anamorphosis) are combined into a single affine
    transformation applied to the DM modes (the phase), with one
    interpolation step. The slopes (derivatives or telescoping sum) are then
    computed on the transformed phase.

    method selects the implementation (see "Computation Workflows" and
    "Reference and optimized implementations" in docs/general.rst):

    - 'operator' (default): the masks do not depend on the DM mode, and
      once they are fixed every step is linear in the phase, so the whole
      chain is one sparse matrix A applied to the DM cube, IM = A . Phi
      (synim.sparse_operator). No transformed copy of the cube or of its
      derivatives is allocated.
    - 'reference': the explicit step-by-step computation on the
      transformed cube (apply_dm_transformations_combined +
      apply_wfs_transformations_combined, synim.reference). It needs memory
      for the transformed cube and its derivatives; use it for validation
      and debugging.

    The two methods give the same results up to rounding.

    Note: WFS transformations must not be applied to derivative maps. The
    gradient of the transformed phase includes the Jacobian of the
    transformation (it mixes the x and y components for rotations and 45 deg
    anamorphosis and scales them for magnification and 90 deg anamorphosis),
    which a spatial resampling of the derivative maps does not account for.
    """
    implementation = _implementation(method)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Interaction Matrix Computation")
        print(f"{'='*60}\n")

    dm_array = to_xp(xp, dm_array, dtype=float_dtype)
    _check_dm_mask_shape(dm_array, dm_mask, specula_convention)

    im = implementation.interaction_matrix(
        pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
        wfs_nsubaps, wfs_fov_arcsec, gs_pol_coo, gs_height,
        wfs_rotation, wfs_translation, wfs_mag_global,
        wfs_anamorphosis_90=wfs_anamorphosis_90, wfs_anamorphosis_45=wfs_anamorphosis_45,
        idx_valid_sa=idx_valid_sa, slope_method=slope_method,
        specula_convention=specula_convention, verbose=verbose)

    if display:
        _display_transformed_dm(
            pup_diam_m, pup_mask, dm_array, dm_mask,
            (dm_height, dm_rotation, gs_pol_coo, gs_height, wfs_rotation, wfs_translation,
             wfs_mag_global, wfs_anamorphosis_90, wfs_anamorphosis_45),
            specula_convention)

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
                                   verbose=False,
                                   method='operator'):
    """
    Computes interaction matrices for multiple WFS configurations.
    
    Each WFS can have its own guide star position (gs_pol_coo) and height (gs_height).
    The DM array is converted to the target device once; each WFS is
    computed as in interaction_matrix (single combined transformation).
    WFS with the same guide star, WFS transformations, number of subapertures
    and slope method share the transformation and the slope operators
    (method='operator') or the transformed DM array and its derivatives
    (method='reference'), which are computed once for the whole group.
    With method='operator' the operators of all the WFS are applied to the
    DM array with a single product.
    
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
    - minimize_memory: bool, optional, kept for compatibility and ignored
      (method='operator' allocates no transformed copy of the DM array)
    - verbose: bool, optional
    - method: str, optional, 'operator' (default) or 'reference', the
      implementation (see interaction_matrix)

    Returns:
    - im_dict: dict, interaction matrices keyed by WFS name or index
    - derivatives_info: dict with metadata about the computation:
        'workflow': 'combined'
        'groups': list of lists of WFS names sharing the DM transformation
        'n_dm_transformations': number of DM transformations computed
    """

    implementation = _implementation(method)

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

    if verbose:
        for i_group, members in enumerate(groups.values()):
            first = wfs_params[members[0]]
            names = ', '.join(str(wfs_params[i]['name']) for i in members)
            print(f"  [group {i_group+1}/{len(groups)}] {names}:")
            print(f"    Subapertures: {first['nsubaps']}x{first['nsubaps']}")
            print(f"    GS: {first['gs_pol_coo']}, height: {first['gs_height']} m")

    ims = implementation.interaction_matrices(
        pup_diam_m, pup_mask, dm_array, dm_mask, dm_height, dm_rotation,
        wfs_params, list(groups.values()), slope_method, specula_convention,
        verbose=verbose)

    im_list = []
    for params, im in zip(wfs_params, ims):
        im_list.append(cpuArray(im) if im_on_cpu else im)
        if verbose:
            print(f"    {params['name']}: IM shape {im.shape}")

    im_dict = {params['name']: im for params, im in zip(wfs_params, im_list)}
    for name, im in im_dict.items():
        _warn_if_nan(im, name=f'interaction matrix of {name}')

    derivatives_info = {
        'workflow': 'combined',
        # WFS names of each group sharing the DM transformation
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
