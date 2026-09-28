"""
interaction_matrices_multi_wfs: WFS with the same guide star, WFS
transformations, number of subapertures and slope method share the DM
transformation and the derivatives. The results must be identical to the
interaction matrices computed one WFS at a time.
"""
import unittest

import numpy as np

import synim
import synim.synim as synim_core

N_PIX = 96
N_SUBAPS = 12


class TestMultiWfsGroups(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        n = N_PIX
        yy, xx = np.mgrid[:n, :n]
        r = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5)
        cls.pup_mask = ((r < n / 2) & (r > n / 8)).astype(np.float32)
        cls.dm_mask = (r < n / 2 + 4).astype(np.float32)
        cls.dm_array = np.random.default_rng(0).standard_normal((n, n, 8)).astype(np.float32)
        illumination = synim_core.compute_subaperture_illumination(cls.pup_mask, N_SUBAPS)
        cls.idx_valid_sa = np.argwhere(
            illumination.reshape(N_SUBAPS, N_SUBAPS).T > 0.5).astype(np.int32)
        common = dict(nsubaps=N_SUBAPS, gs_pol_coo=(30.0, 0.0), gs_height=90e3,
                      rotation=10.0, translation=(0.3, -0.1))
        cls.configs = [
            dict(common, name='a', fov_arcsec=4.0),
            # same DM transformation as 'a', given with other sequence types
            dict(common, name='b', fov_arcsec=4.0, gs_pol_coo=[30, 0],
                 translation=np.array([0.3, -0.1])),
            dict(common, name='c', fov_arcsec=4.0, gs_pol_coo=(30.0, 120.0)),
            # FoV and valid subapertures do not change the DM transformation
            dict(common, name='d', fov_arcsec=6.0, idx_valid_sa=cls.idx_valid_sa),
            dict(common, name='e', fov_arcsec=4.0, magnification=(1.02, 1.0)),
        ]

    def _single(self, config, slope_method):
        mag = config.get('magnification', (1.0, 1.0))
        return synim.cpuArray(synim_core.interaction_matrix(
            8.0, self.pup_mask, self.dm_array, self.dm_mask, 8000.0, 1.0,
            config['nsubaps'], config['fov_arcsec'], tuple(config['gs_pol_coo']),
            config['gs_height'], config['rotation'], tuple(config['translation']),
            float(np.sqrt(mag[0] * mag[1])), wfs_anamorphosis_90=mag[1] / mag[0],
            idx_valid_sa=config.get('idx_valid_sa'), slope_method=slope_method))

    def test_groups_and_results(self):
        for slope_method in ('derivatives', 'telsum'):
            with self.subTest(slope_method=slope_method):
                im_dict, info = synim_core.interaction_matrices_multi_wfs(
                    8.0, self.pup_mask, self.dm_array, self.dm_mask, 8000.0, 1.0,
                    self.configs, slope_method=slope_method)
                self.assertEqual(info['workflow'], 'combined')
                self.assertEqual(info['groups'], [['a', 'b', 'd'], ['c'], ['e']])
                self.assertEqual(info['n_dm_transformations'], 3)
                self.assertEqual(list(im_dict), ['a', 'b', 'c', 'd', 'e'])
                for config in self.configs:
                    np.testing.assert_array_equal(
                        synim.cpuArray(im_dict[config['name']]),
                        self._single(config, slope_method),
                        err_msg=f"WFS {config['name']}")

    def test_im_on_cpu(self):
        im_dict, _ = synim_core.interaction_matrices_multi_wfs(
            8.0, self.pup_mask, self.dm_array, self.dm_mask, 8000.0, 1.0,
            self.configs[:2], im_on_cpu=True)
        for im in im_dict.values():
            self.assertIsInstance(im, np.ndarray)


if __name__ == '__main__':
    unittest.main()
