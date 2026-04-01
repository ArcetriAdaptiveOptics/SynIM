import unittest
import numpy as np

from synim.params_manager import ParamsManager


class TestReconstructorModeCut(unittest.TestCase):

    def test_split_mode_columns_respects_absolute_mode_index(self):
        # Component 1 starts from mode 2, component 2 from mode 0.
        mode_indices = [
            [2, 3, 4],
            [0, 1, 2, 3],
        ]

        keep_cols, drop_cols = ParamsManager._split_mode_columns_by_absolute_index(
            mode_indices,
            n_low_modes_cut=2,
        )

        # Global columns are: [2,3,4,0,1,2,3]
        # Drop only absolute modes < 2 -> columns for [0,1] in component 2.
        np.testing.assert_array_equal(drop_cols, np.array([3, 4]))
        np.testing.assert_array_equal(keep_cols, np.array([0, 1, 2, 5, 6]))

    def test_restore_reconstructor_with_zero_rows(self):
        rec_reduced = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
        keep_cols = np.array([1, 3], dtype=int)

        rec_full = ParamsManager._restore_reconstructor_with_zero_rows(
            rec_reduced=rec_reduced,
            keep_cols=keep_cols,
            n_full_modes=5,
        )

        self.assertEqual(rec_full.shape, (5, 2))
        np.testing.assert_array_equal(rec_full[1], rec_reduced[0])
        np.testing.assert_array_equal(rec_full[3], rec_reduced[1])
        np.testing.assert_array_equal(rec_full[0], np.zeros(2, dtype=np.float32))
        np.testing.assert_array_equal(rec_full[2], np.zeros(2, dtype=np.float32))
        np.testing.assert_array_equal(rec_full[4], np.zeros(2, dtype=np.float32))

    def test_low_mode_cut_requires_offline_filter(self):
        pm = ParamsManager.__new__(ParamsManager)
        pm.wfs_list = [{'name': 'ngs1', 'index': 1}]
        pm.params = {
            'slopec_ngs1': {
                'filtmat_tag': 'dummy_filter'
            }
        }
        pm.ngs_recon_params = {
            'sigma2_in_nm2': 0.0,
            'noise_elong_model': False,
            'naThicknessInM': 0.0,
            'tGparameter': 0.0,
            'n_low_modes_zero_filt': 3,
        }

        n_cut = pm._get_low_modes_cut_if_applicable(
            wfs_type='ngs',
            apply_filter=True,
            active_wfs_mask=None,
            verbose=False,
        )
        self.assertEqual(n_cut, 3)

        with self.assertRaises(ValueError):
            pm._get_low_modes_cut_if_applicable(
                wfs_type='ngs',
                apply_filter=False,
                active_wfs_mask=None,
                verbose=False,
            )
