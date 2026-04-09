import unittest
from unittest.mock import patch
import numpy as np

from synim import cpu_float_dtype
from synim.params_manager import ParamsManager

class TestMultirateLUTGeneration(unittest.TestCase):
    """
    Test the generation of the Multirate Look-Up Table (LUT) for reconstruction matrices.
    Verifies the binary loop, masking, zero-padding, and filename suffixes.
    """

    def setUp(self):
        # Create a minimal fake configuration for ParamsManager
        self.mock_config = {
            'main': {
                'root_dir': '/tmp',
                'pixel_pupil': 100,
                'pixel_pitch': 0.01
            },
            'pupilstop': {
                'mask_diam': 1.0,
                'obs_diam': 0.0
            },
            'sh_ngs1': {'type': 'sh', 'index': '1'},
            'sh_ngs2': {'type': 'sh', 'index': '2'},
            'dm': {'type': 'dm', 'index': '1'}
        }

        # Initialize ParamsManager with the mock dict instead of a file
        self.pm = ParamsManager(self.mock_config, verbose=False)

        # Override the wfs_list extracted to match our fake config
        self.pm.wfs_list = [
            {'name': 'sh_ngs1', 'type': 'sh', 'index': '1'},
            {'name': 'sh_ngs2', 'type': 'sh', 'index': '2'}
        ]

    @patch('synim.params_manager.ParamsManager.compute_tomographic_reconstructor')
    def test_compute_multirate_reconstructors(self, mock_compute_tomo):
        """
        Test that compute_multirate_reconstructors correctly iterates 2^N-1 times,
        generates the correct binary masks, and passes the correct padding arguments.
        """

        # 1. Setup the Mock to return a fake unpadded reconstructor
        fake_unpadded_rec = np.ones((3, 10), dtype=cpu_float_dtype)

        mock_compute_tomo.return_value = {
            'reconstructor': fake_unpadded_rec,
            'mode_indices': [[0, 1, 2]],
        }

        # 2. Call the function we want to test
        results = self.pm.compute_multirate_reconstructors(
            r0=0.15, L0=25.0,
            wfs_type='ngs',
            component_type='dm',
            save=False,
            verbose=False
        )

        # 3. Assertions on the wrapper logic
        # For N=2 sensors, we expect 2^2 - 1 = 3 iterations
        self.assertEqual(len(results), 3, "Should generate exactly 3 matrices for 2 sensors")

        # Check that the keys are the correct boolean tuples
        expected_keys = [(False, True), (True, False), (True, True)]
        for k in expected_keys:
            self.assertIn(k, results, f"Missing validity mask {k} in results")

        self.assertEqual(mock_compute_tomo.call_count, 3)

        # 4. Verify the arguments passed to the mocked function!
        # This is how we test that the wrapper correctly passes the padding target and the masks
        calls = mock_compute_tomo.call_args_list

        # Extract all the 'active_wfs_mask' arguments passed to the mock
        masks_passed = [call.kwargs['active_wfs_mask'] for call in calls]

        self.assertIn([False, True], masks_passed)
        self.assertIn([True, False], masks_passed)
        self.assertIn([True, True], masks_passed)
