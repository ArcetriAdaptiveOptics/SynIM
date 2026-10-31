"""
interaction_matrix and interaction_matrices_multi_wfs compute the
interaction matrix as a sparse linear operator applied to the DM array
(see "Computation Workflows" in docs/general.rst). The explicit computation
on the transformed DM array (apply_dm_transformations_combined followed by
apply_wfs_transformations_combined) is kept as the independent reference:
the two must give the same results, up to rounding errors, for both slope
methods, both conventions, integer and non-integer numbers of pixels per
subaperture, with and without a selection of the valid subapertures, and
for each geometric term separately (DM rotation, DM shift and magnification
from the guide star, WFS rotation, shift, magnification and both
anamorphoses) and combined.
"""
import itertools
import unittest

import numpy as np

import synim
import synim.synim as synim_core
import synim.utils as synim_utils

N_PIX = 96
N_MODES = 6
PUP_DIAM_M = 8.0
FOV_ARCSEC = 3.0
# Integer (96 / 12) and non-integer (96 / 13) pixels per subaperture
N_SUBAPS = (12, 13)
# Tolerance relative to max(|IM|). On GPU in single precision the explicit
# computation is less accurate: cupyx affine_transform computes the
# interpolation coordinates of the DM array in single precision.
if synim.float_dtype != np.float32:
    RTOL = 1e-10
elif synim.xp is np:
    RTOL = 1e-5
else:
    RTOL = 1e-4

NO_TRANSFORMATION = dict(
    dm_height=0.0, dm_rotation=0.0, gs_pol_coo=(0.0, 0.0), gs_height=np.inf,
    wfs_rotation=0.0, wfs_translation=(0.0, 0.0), wfs_mag_global=1.0,
    wfs_anamorphosis_90=1.0, wfs_anamorphosis_45=1.0)

# One geometric term at a time, then all of them
GEOMETRIES = {
    'none': {},
    'dm_rotation': dict(dm_rotation=7.0),
    'dm_shift': dict(dm_height=8000.0, gs_pol_coo=(30.0, 40.0)),
    'lgs_on_axis': dict(dm_height=8000.0, gs_height=90e3),
    'lgs_off_axis': dict(dm_height=8000.0, gs_pol_coo=(30.0, 40.0), gs_height=90e3),
    'wfs_rotation': dict(wfs_rotation=12.0),
    'wfs_shift': dict(wfs_translation=(0.7, -1.3)),
    'wfs_magnification': dict(wfs_mag_global=1.02),
    'wfs_anamorphosis_90': dict(wfs_anamorphosis_90=1.02),
    'wfs_anamorphosis_45': dict(wfs_anamorphosis_45=1.02),
    'all': dict(dm_height=8000.0, dm_rotation=-5.0, gs_pol_coo=(25.0, 200.0),
                gs_height=90e3, wfs_rotation=-9.0, wfs_translation=(-0.4, 0.9),
                wfs_mag_global=0.99, wfs_anamorphosis_90=1.01,
                wfs_anamorphosis_45=0.985),
}


class TestSparseOperator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        n = N_PIX
        yy, xx = np.mgrid[:n, :n]
        r = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5)
        cls.pup_mask = ((r < n / 2 - 3) & (r > n / 8)).astype(np.float32)
        cls.dm_mask = (r < n / 2 - 1).astype(np.float32)
        cls.dm_array = (np.random.default_rng(0).standard_normal((n, n, N_MODES))
                        * cls.dm_mask[:, :, np.newaxis]).astype(np.float32)
        rng = np.random.default_rng(1)
        cls.idx_valid_sa = {n_subaps: np.argwhere(rng.random((n_subaps, n_subaps)) > 0.3)
                            for n_subaps in N_SUBAPS}

    def _options(self):
        """Slope method, convention, number of subapertures, valid subapertures."""
        for slope_method, specula_convention, n_subaps, select in itertools.product(
                ('derivatives', 'telsum'), (True, False), N_SUBAPS, (False, True)):
            yield dict(slope_method=slope_method, specula_convention=specula_convention,
                       wfs_nsubaps=n_subaps,
                       idx_valid_sa=self.idx_valid_sa[n_subaps] if select else None)

    def _explicit(self, geometry, options):
        """Interaction matrix computed explicitly on the transformed DM array."""
        options = dict(options)
        idx_valid_sa = options.pop('idx_valid_sa')
        if idx_valid_sa is not None:
            idx_valid_sa = synim.xp.asarray(idx_valid_sa)
        _, trans_dm_mask, trans_pup_mask, derivatives_x, derivatives_y = \
            synim_core.apply_dm_transformations_combined(
                PUP_DIAM_M, self.pup_mask, self.dm_array, self.dm_mask,
                **geometry, **options)
        return synim.cpuArray(synim_core.apply_wfs_transformations_combined(
            derivatives_x, derivatives_y, trans_pup_mask, trans_dm_mask,
            options['wfs_nsubaps'], FOV_ARCSEC, PUP_DIAM_M, idx_valid_sa=idx_valid_sa,
            slope_method=options['slope_method'],
            specula_convention=options['specula_convention'], in_place=True))

    def assert_same(self, im, expected):
        im = synim.cpuArray(im)
        self.assertEqual(im.shape, expected.shape)
        np.testing.assert_allclose(im, expected, rtol=0, atol=RTOL * np.abs(expected).max())

    def test_interaction_matrix(self):
        for name, terms in GEOMETRIES.items():
            geometry = dict(NO_TRANSFORMATION, **terms)
            for options in self._options():
                with self.subTest(geometry=name, **{k: v is not None if k == 'idx_valid_sa'
                                                    else v for k, v in options.items()}):
                    im = synim_core.interaction_matrix(
                        PUP_DIAM_M, self.pup_mask, self.dm_array, self.dm_mask,
                        wfs_fov_arcsec=FOV_ARCSEC, **geometry, **options)
                    self.assert_same(im, self._explicit(geometry, options))

    def test_interaction_matrices_multi_wfs(self):
        """All the geometries with the same DM parameters in one call, one WFS each."""
        dm_params = {}
        for name, terms in GEOMETRIES.items():
            geometry = dict(NO_TRANSFORMATION, **terms)
            dm_params.setdefault((geometry['dm_height'], geometry['dm_rotation']),
                                 []).append((name, geometry))

        for (dm_height, dm_rotation), geometries in dm_params.items():
            for options in self._options():
                configs = []
                for name, geometry in geometries:
                    mag, anam_90 = geometry['wfs_mag_global'], geometry['wfs_anamorphosis_90']
                    configs.append(dict(
                        name=name, nsubaps=options['wfs_nsubaps'], fov_arcsec=FOV_ARCSEC,
                        gs_pol_coo=geometry['gs_pol_coo'], gs_height=geometry['gs_height'],
                        rotation=geometry['wfs_rotation'],
                        translation=geometry['wfs_translation'],
                        magnification=(mag / np.sqrt(anam_90), mag * np.sqrt(anam_90)),
                        anamorphosis_45=geometry['wfs_anamorphosis_45'],
                        idx_valid_sa=options['idx_valid_sa']))
                im_dict, _ = synim_core.interaction_matrices_multi_wfs(
                    PUP_DIAM_M, self.pup_mask, self.dm_array, self.dm_mask,
                    dm_height, dm_rotation, configs, slope_method=options['slope_method'],
                    specula_convention=options['specula_convention'])

                for config, (name, geometry) in zip(configs, geometries):
                    with self.subTest(geometry=name, **{
                            k: v is not None if k == 'idx_valid_sa' else v
                            for k, v in options.items()}):
                        # same magnification parameters as the multi-WFS function
                        mag, anam_90 = synim_core._wfs_magnification_params(config)
                        expected = self._explicit(
                            dict(geometry, wfs_mag_global=mag, wfs_anamorphosis_90=anam_90),
                            options)
                        self.assert_same(im_dict[name], expected)

    def test_interpolation_matrix(self):
        """affine_interpolation_matrix gives the same result as rotshiftzoom_array."""
        rng = np.random.default_rng(2)
        cube = rng.standard_normal((70, 64, 3))
        transform = dict(dm_translation=(2.0, -1.5), dm_rotation=12.0,
                         dm_magnification=(1.05, 0.98), wfs_translation=(0.4, 0.7),
                         wfs_rotation=-33.0, wfs_magnification=(1.01, 0.99),
                         wfs_anamorphosis_45=1.02)
        for transposed in (False, True):
            with self.subTest(transposed=transposed):
                source = np.transpose(cube, (1, 0, 2)) if transposed else cube
                expected = synim.cpuArray(synim_utils.rotshiftzoom_array(
                    synim.xp.asarray(source, dtype=synim.float_dtype),
                    output_size=(60, 60), **transform))
                matrix, offset = synim_utils.rotshiftzoom_matrix(
                    source.shape[:2], (60, 60), **transform)
                operator = synim_utils.affine_interpolation_matrix(
                    source.shape[:2], (60, 60), matrix, offset, transposed=transposed)
                result = (operator @ cube.reshape(-1, 3)).reshape(60, 60, 3)
                np.testing.assert_allclose(result, expected, rtol=0,
                                           atol=RTOL * np.abs(expected).max())

    def test_non_finite_values(self):
        """NaN and inf values of the DM array are replaced with nan_to_num, as
        rotshiftzoom_array does in the explicit computation."""
        geometry = dict(NO_TRANSFORMATION, **GEOMETRIES['all'])
        for value in (np.nan, np.inf):
            with self.subTest(value=value):
                dm_array = self.dm_array.copy()
                # a pixel inside the pupil, used by the operator
                dm_array[N_PIX // 2, N_PIX // 2 + 30, 1] = value
                im, expected = (synim.cpuArray(synim_core.interaction_matrix(
                    PUP_DIAM_M, self.pup_mask, array, self.dm_mask, wfs_nsubaps=12,
                    wfs_fov_arcsec=FOV_ARCSEC, **geometry))
                    for array in (dm_array, np.nan_to_num(dm_array)))
                np.testing.assert_array_equal(im, expected)
                self.assertFalse(np.isnan(im).any())

if __name__ == '__main__':
    unittest.main()
