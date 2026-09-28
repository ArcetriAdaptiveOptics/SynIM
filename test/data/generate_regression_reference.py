"""
Generate test/data/regression_reference.npz, used by test/test_regression_reference.py.

The reference file was generated on commit 12ea3de (before the removal of
the separated workflow), on CPU, in single precision:

    python test/data/generate_regression_reference.py

Interaction matrices are always computed with the COMBINED workflow (DM and
WFS transformations applied to the phase in a single interpolation, then
derivatives / telescoping sum). On 12ea3de some geometries were routed to the
separated workflow, which applied the WFS transformation to the derivative
maps and is wrong for magnification and anamorphosis (the Jacobian of the
transformation is missing). For those cases the reference is computed by
calling the combined functions of 12ea3de directly. For the geometries where
the separated workflow was correct (no WFS transformation, pure WFS shift)
the script also reports the difference between the two workflows.

Regenerating the file on a later commit uses the public interaction_matrix,
which only has the combined workflow.

The telescoping sum entries (im_telsum_*) were regenerated when the
telescoping sum was changed to average over the pupil pixels, as the
derivatives, instead of over the DM mask: they differ from 12ea3de in the
subapertures where the pupil and the DM mask differ (edge and central
obstruction). All the other entries are unchanged.
"""
import os
import sys

import numpy as np

import synim
# Initialize on CPU when run as a script; when imported by the tests SynIM is
# already initialized (possibly on GPU, see test/__init__.py)
if synim.xp is None:
    synim.init(device_idx=-1, precision=1)
import synim.synim as synim_core  # noqa: E402
import synim.synpm as synpm  # noqa: E402
import synim.utils as synim_utils  # noqa: E402
import synim.params_utils as params_utils  # noqa: E402

OUTPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'regression_reference.npz')

N_PIX = 96
N_SUBAPS = 12
N_MODES = 30
PUP_DIAM_M = 8.0
FOV_ARCSEC = 4.0

# name: interaction_matrix geometry parameters
GEOMETRIES = {
    'none': dict(),
    'dm_only': dict(dm_height=5000.0, dm_rotation=2.0, gs_pol_coo=(30.0, 45.0)),
    'wfs_shift': dict(wfs_translation=(1.3, -0.7)),
    'wfs_magnification': dict(wfs_mag_global=1.05),
    'wfs_anamorphosis_90': dict(wfs_anamorphosis_90=1.05),
    'wfs_anamorphosis_45': dict(wfs_anamorphosis_45=1.05),
    'dm_and_wfs_anamorphosis': dict(dm_height=5000.0, dm_rotation=1.0,
                                    gs_pol_coo=(30.0, 45.0), wfs_anamorphosis_45=0.97),
    'wfs_rotation': dict(wfs_rotation=4.0),
    'dm_and_wfs': dict(dm_height=8000.0, dm_rotation=1.0, gs_pol_coo=(20.0, 10.0),
                       gs_height=90e3, wfs_rotation=3.0, wfs_translation=(0.5, 0.2),
                       wfs_mag_global=0.98),
}
DEFAULTS = dict(dm_height=0.0, dm_rotation=0.0, gs_pol_coo=(0.0, 0.0), gs_height=np.inf,
                wfs_rotation=0.0, wfs_translation=(0.0, 0.0), wfs_mag_global=1.0,
                wfs_anamorphosis_90=1.0, wfs_anamorphosis_45=1.0)


def make_inputs():
    n = N_PIX
    yy, xx = np.mgrid[:n, :n]
    r = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5)
    pup_mask = ((r < n / 2) & (r > n / 8)).astype(np.float32)
    dm_mask = (r < n / 2 + 4).astype(np.float32)
    rng = np.random.default_rng(0)
    x = (xx - n / 2) / (n / 2)
    y = (yy - n / 2) / (n / 2)
    modes = [rng.normal() * x ** (k % 4) * y ** (k // 4 % 4)
             + rng.normal() * np.sin((k + 1) * x) * np.cos((k % 5) * y)
             for k in range(N_MODES)]
    dm_array = np.stack(modes, axis=2).astype(np.float32)
    return pup_mask, dm_mask, dm_array


def inputs_checksum(pup_mask, dm_mask, dm_array):
    dm64 = dm_array.astype(np.float64)
    return np.array([pup_mask.sum(), dm_mask.sum(), dm64.sum(),
                     np.abs(dm64).sum(), dm64[::7, ::5, ::3].sum()])


def combined_im(pup_mask, dm_mask, dm_array, geometry, slope_method, idx_valid_sa=None):
    """Interaction matrix with the combined workflow."""
    p = dict(DEFAULTS, **geometry)
    if not hasattr(synim_core, 'apply_dm_transformations_combined'):
        return synim_core.interaction_matrix(
            PUP_DIAM_M, pup_mask, dm_array, dm_mask, wfs_nsubaps=N_SUBAPS,
            wfs_fov_arcsec=FOV_ARCSEC, idx_valid_sa=idx_valid_sa,
            slope_method=slope_method, **p)
    if idx_valid_sa is not None:
        idx_valid_sa = synim.to_xp(synim.xp, idx_valid_sa)
    _, trans_dm_mask, trans_pup_mask, der_x, der_y = \
        synim_core.apply_dm_transformations_combined(
            PUP_DIAM_M, pup_mask, dm_array, dm_mask,
            p['dm_height'], p['dm_rotation'], p['gs_pol_coo'], p['gs_height'],
            p['wfs_rotation'], p['wfs_translation'], wfs_nsubaps=N_SUBAPS,
            wfs_mag_global=p['wfs_mag_global'],
            wfs_anamorphosis_90=p['wfs_anamorphosis_90'],
            wfs_anamorphosis_45=p['wfs_anamorphosis_45'],
            slope_method=slope_method)
    return synim_core.apply_wfs_transformations_combined(
        der_x, der_y, trans_pup_mask, trans_dm_mask, N_SUBAPS, FOV_ARCSEC,
        PUP_DIAM_M, idx_valid_sa=idx_valid_sa, slope_method=slope_method)


def main():
    pup_mask, dm_mask, dm_array = make_inputs()
    # The inputs are regenerated by the test with make_inputs(); only a
    # checksum is stored, to detect a change of the inputs.
    ref = {'inputs_checksum': inputs_checksum(pup_mask, dm_mask, dm_array)}

    illumination = synim_core.compute_subaperture_illumination(pup_mask, N_SUBAPS)
    ref['illumination'] = synim.cpuArray(illumination)
    # Valid subapertures (2D indices, SPECULA convention)
    valid = np.argwhere(ref['illumination'].reshape(N_SUBAPS, N_SUBAPS).T > 0.5)
    ref['idx_valid_sa'] = valid.astype(np.int32)

    for slope_method in ('derivatives', 'telsum'):
        for name, geometry in GEOMETRIES.items():
            im = combined_im(pup_mask, dm_mask, dm_array, geometry, slope_method)
            # Stored in float32 (on 12ea3de telsum returned float64 by mistake)
            ref[f'im_{slope_method}_{name}'] = synim.cpuArray(im).astype(np.float32)
            # Report the difference with the public function where it works
            try:
                im_public = synim.cpuArray(synim_core.interaction_matrix(
                    PUP_DIAM_M, pup_mask, dm_array, dm_mask, wfs_nsubaps=N_SUBAPS,
                    wfs_fov_arcsec=FOV_ARCSEC, slope_method=slope_method,
                    **dict(DEFAULTS, **geometry)))
                im_ref = ref[f'im_{slope_method}_{name}']
                diff = np.nanmax(np.abs(im_public - im_ref)) / np.nanmax(np.abs(im_ref))
                print(f'{slope_method:11s} {name:24s} interaction_matrix vs combined:'
                      f' {diff:.2e}')
            except ValueError as err:
                print(f'{slope_method:11s} {name:24s} interaction_matrix fails: {err}')

    # Valid subaperture selection (SPECULA convention)
    ref['im_derivatives_dm_and_wfs_idx_valid'] = synim.cpuArray(combined_im(
        pup_mask, dm_mask, dm_array, GEOMETRIES['dm_and_wfs'], 'derivatives',
        idx_valid_sa=ref['idx_valid_sa']))

    # Multi-WFS: different guide stars and WFS transformations
    configs = [
        dict(name='a', nsubaps=N_SUBAPS, fov_arcsec=FOV_ARCSEC,
             gs_pol_coo=(30.0, 0.0), gs_height=90e3),
        dict(name='b', nsubaps=N_SUBAPS, fov_arcsec=FOV_ARCSEC,
             gs_pol_coo=(30.0, 120.0), gs_height=90e3, rotation=10.0,
             translation=(0.3, -0.1)),
    ]
    im_dict, _ = synim_core.interaction_matrices_multi_wfs(
        PUP_DIAM_M, pup_mask, dm_array, dm_mask, 8000.0, 1.0, configs)
    for key, value in im_dict.items():
        ref[f'multi_{key}'] = synim.cpuArray(value)

    # Projection matrix
    idx = np.where(pup_mask > 0.5)
    base = dm_array[idx[0], idx[1], :10].astype(np.float64)
    base_inv = np.linalg.pinv(base).astype(np.float32)
    ref['projection'] = synim.cpuArray(synpm.projection_matrix(
        PUP_DIAM_M, pup_mask, dm_array, dm_mask, base_inv, 5000.0, 1.0,
        0.0, (0, 0), (1, 1), (20, 10), np.inf))

    # 2D <-> 3D influence function conversion
    dm_2d = synim_utils.dm3d_to_2d(dm_array[:, :, :3].copy(), dm_mask, xp_local=np,
                                   float_dtype_local=np.float32)
    ref['dm3d_to_2d'] = synim.cpuArray(dm_2d)
    dm_3d = synim.cpuArray(synim_utils.dm2d_to_3d(
        dm_2d.copy(), dm_mask, xp_local=np, float_dtype_local=np.float32))
    ref['dm2d_to_3d_on_mask'] = dm_3d[dm_mask > 0]   # zeros outside the mask
    ref['dm2d_to_3d_outside_mask_abs_sum'] = np.array(np.abs(dm_3d[dm_mask == 0]).sum())

    # Reconstructors (double precision inputs)
    im = ref['im_derivatives_none'].astype(np.float64)
    c_atm = np.diag(1.0 / np.arange(1, N_MODES + 1) ** 1.5)
    c_noise = np.eye(im.shape[0]) * 1e-2
    ref['reconstructor_mmse'] = params_utils.compute_mmse_reconstructor(
        im, c_atm, C_noise=c_noise, xp=np, dtype=np.float64)
    ref['reconstructor_pinv'] = params_utils.compute_pseudoinverse_reconstructor(
        im, xp=np, dtype=np.float64)

    np.savez_compressed(OUTPUT, **ref)
    print(f'\nSaved {len(ref)} arrays to {OUTPUT}'
          f' ({os.path.getsize(OUTPUT) / 1024:.0f} KB)')


if __name__ == '__main__':
    sys.exit(main())
