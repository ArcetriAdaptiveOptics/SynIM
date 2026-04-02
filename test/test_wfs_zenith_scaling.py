import unittest
from unittest.mock import patch
import numpy as np

from synim.params_manager import ParamsManager
from synim.params_utils import generate_im_filename


class TestWfsZenithScaling(unittest.TestCase):
    def setUp(self):
        self.mock_config = {
            'main': {
                'root_dir': '/tmp',
                'pixel_pupil': 64,
                'pixel_pitch': 0.01,
                'zenithAngleInDeg': 60.0,
            },
            'pupilstop': {
                'mask_diam': 1.0,
                'obs_diam': 0.0,
            },
            'dm1': {'class': 'DM', 'height': 0.0},
            'sh_lgs1': {'class': 'ShSlopec'},
            'source_lgs1': {
                'height': 90000.0,
                'polar_coordinates': [10.0, 0.0],
            },
        }

        self.pm = ParamsManager(self.mock_config, verbose=False)
        self.pm.wfs_list = [
            {
                'name': 'sh_lgs1',
                'type': 'sh',
                'index': '1',
                'config': {'subap_on_diameter': 8, 'fov': 4.0},
            },
        ]

    @patch('synim.params_manager.find_subapdata', return_value=None)
    @patch.object(ParamsManager, 'get_component_params')
    def test_prepare_params_scales_gs_height_for_interaction_matrix(
        self, mock_get_component_params, _mock_subap
    ):
        mock_get_component_params.return_value = {
            'dm_array': np.zeros((8, 8, 2), dtype=np.float32),
            'dm_mask': np.ones((8, 8), dtype=np.float32),
            'dm_height': 0.0,
            'dm_rotation': 0.0,
            'component_key': 'dm1',
        }

        params = self.pm.prepare_interaction_matrix_params(
            wfs_type='lgs', wfs_index=1, dm_index=1
        )

        self.assertAlmostEqual(params['gs_height'], 180000.0, places=6)

    @patch('synim.params_manager.find_subapdata', return_value=None)
    @patch.object(ParamsManager, 'get_component_params')
    @patch('synim.params_manager.generate_im_filename', return_value='im_test.fits')
    @patch('synim.params_manager.Intmat')
    @patch('synim.params_manager.synim.interaction_matrices_multi_wfs')
    def test_wfs_configs_scales_gs_height_before_multi_wfs_call(
        self,
        mock_multi,
        _mock_intmat,
        _mock_gen_name,
        mock_get_component_params,
        _mock_subap,
    ):
        mock_get_component_params.return_value = {
            'dm_array': np.zeros((8, 8, 2), dtype=np.float32),
            'dm_mask': np.ones((8, 8), dtype=np.float32),
            'dm_height': 0.0,
            'dm_rotation': 0.0,
            'component_key': 'dm1',
        }

        mock_multi.return_value = ({'sh_lgs1': np.zeros((4, 2), dtype=np.float32)},
                                   {'workflow': 'unified'})

        self.pm.compute_interaction_matrices(
            output_im_dir='/tmp',
            output_rec_dir='/tmp',
            wfs_type='lgs',
            overwrite=True,
            verbose=False,
            display=False,
        )

        called_wfs_configs = mock_multi.call_args.kwargs['wfs_configs']
        self.assertEqual(len(called_wfs_configs), 1)
        self.assertAlmostEqual(called_wfs_configs[0]['gs_height'], 180000.0, places=6)

    def test_im_filename_uses_zenith_scaled_height_integer_meters(self):
        filename = generate_im_filename(
            self.mock_config,
            wfs_type='lgs',
            wfs_index=1,
            dm_index=1,
        )

        self.assertIn('h180000', filename)
        self.assertNotIn('h180000.', filename)
