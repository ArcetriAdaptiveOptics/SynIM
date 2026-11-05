"""
Tests of ParamsManager.compute_tomographic_projection_matrix (PMs read one
optical source at a time) against the explicit formula

    tpdm_pdm = sum_i w_i P_dm,i P_dm,i^T        tpdm_pl = sum_i w_i P_dm,i P_l,i^T
    P_opt    = (tpdm_pdm + reg I)^-1 tpdm_pl

computed independently in float64 from the same random PMs written to disk.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

from test import on_backend  # noqa: F401  (initializes SynIM before the imports below)

import specula
from specula.data_objects.intmat import Intmat

from synim.params_manager import ParamsManager
from synim.params_utils import generate_pm_filename

N_PUPIL_MODES = 30
REG = 1e-8
# single precision library vs float64 reference; tpdm_pdm is well conditioned
RTOL = 2e-3
ATOL = 2e-4


def make_config():
    """3 sources, 2 DMs (the 2nd one with start_mode>0), 3 layers (the 2nd with start_mode>0)."""
    return {
        'main': {'root_dir': '/tmp', 'pixel_pupil': 20, 'pixel_pitch': 0.05},
        'pupilstop': {'mask_diam': 1.0, 'obs_diam': 0.0},
        'projection': {
            'reg_factor': REG,
            'ifunc_inverse_tag': 'unused',
            'opt_sources': {
                'polar_coordinates': [[0.0, 0.0], [30.0, 0.0], [30.0, 120.0]],
                'weights': [1.0, 2.0, 5.0],
            },
        },
        'dm1': {'height': 0.0, 'nmodes': 8},
        'dm2': {'height': 6000.0, 'nmodes': 10, 'start_mode': 3},      # 7 modes in the PM
        'layer1': {'height': 100.0, 'nmodes': 12},
        'layer2': {'height': 7000.0, 'nmodes': 9, 'start_mode': 2},    # 7 modes in the PM
        'layer3': {'height': 12000.0, 'nmodes': 6},
    }


DM_AVAIL = [8, 7]
LAYER_AVAIL = [12, 7, 6]
LAYER_START = [0, 2, 0]
WEIGHTS = np.array([1.0, 2.0, 5.0])


class TestTomographicProjectionMatrix(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.config = make_config()
        self.pm = ParamsManager(self.config, verbose=False)
        rng = np.random.default_rng(1234)
        # pm_dm[i][j]: (DM_AVAIL[j], N_PUPIL_MODES), source i; same for layers
        self.pm_dm = [[rng.standard_normal((n, N_PUPIL_MODES)).astype(np.float32)
                       for n in DM_AVAIL] for _ in WEIGHTS]
        self.pm_layer = [[rng.standard_normal((n, N_PUPIL_MODES)).astype(np.float32)
                          for n in LAYER_AVAIL] for _ in WEIGHTS]
        for i in range(len(WEIGHTS)):
            for j in range(len(DM_AVAIL)):
                self._save(self.pm_dm[i][j], i + 1, dm_index=j + 1)
            for j in range(len(LAYER_AVAIL)):
                self._save(self.pm_layer[i][j], i + 1, layer_index=j + 1)

    def _path(self, opt_index, dm_index=None, layer_index=None):
        name = generate_pm_filename(self.config, opt_index=opt_index,
                                    dm_index=dm_index, layer_index=layer_index)
        return os.path.join(self.tmp, name)

    def _save(self, mat, opt_index, dm_index=None, layer_index=None):
        Intmat(mat, target_device_idx=-1, precision=1).save(
            self._path(opt_index, dm_index, layer_index), overwrite=True)

    def _set_reconstruction(self, n_modes_dm, n_modes_layer, n_low):
        self.pm.ref_n_modes_dm = n_modes_dm
        self.pm.ref_n_modes_layer = n_modes_layer
        self.pm.ref_recon_params = {'n_low_modes_zero_filt': n_low}

    def _expected(self, n_modes_dm, n_modes_layer, n_low, zero_layer=None):
        """Explicit formula in float64. zero_layer: layer index (0-based) replaced by zeros."""
        dm_cols = []
        off = 0
        for j, avail in enumerate(DM_AVAIL):
            keep = min(n_modes_dm[j], avail) if n_modes_dm else avail
            dm_cols += list(range(off, off + keep))
            off += avail
        layer_cols, target_cols = [], []
        off = tgt = 0
        for j, avail in enumerate(LAYER_AVAIL):
            keep = min(n_modes_layer[j], avail) if n_modes_layer else avail
            for k in range(keep):
                if LAYER_START[j] + k >= n_low:
                    layer_cols.append(off + k)
                    target_cols.append(tgt + k)
            off += avail
            tgt += keep
        n_target = tgt

        w = WEIGHTS / WEIGHTS.sum()
        pdm_all, pl_all = [], []
        for i in range(len(WEIGHTS)):
            pdm_all.append(np.concatenate(self.pm_dm[i], axis=0).astype(np.float64))
            blocks = [np.zeros_like(b) if j == zero_layer else b
                      for j, b in enumerate(self.pm_layer[i])]
            pl_all.append(np.concatenate(blocks, axis=0).astype(np.float64))
        tpdm_pdm = sum(w[i] * pdm_all[i][dm_cols] @ pdm_all[i][dm_cols].T for i in range(3))
        tpdm_pl = sum(w[i] * pdm_all[i][dm_cols] @ pl_all[i][layer_cols].T for i in range(3))
        reduced = np.linalg.solve(tpdm_pdm + REG * np.eye(len(dm_cols)), tpdm_pl)
        p_opt = np.zeros((len(dm_cols), n_target))
        p_opt[:, target_cols] = reduced
        return p_opt, target_cols

    def _run(self):
        return self.pm.compute_tomographic_projection_matrix(
            output_dir=self.tmp, save=False, wfs_type='ref', verbose=False)

    def test_no_truncation_no_low_modes_cut(self):
        self._set_reconstruction([], [], 0)
        p_opt, _, _, info = self._run()
        expected, _ = self._expected([], [], 0)
        self.assertEqual(p_opt.shape, (15, 25))
        np.testing.assert_allclose(p_opt, expected, rtol=RTOL, atol=ATOL)
        self.assertEqual(info['n_dm_modes_target'], 15)
        self.assertEqual(info['n_layer_modes_target'], 25)

    def test_truncation_and_low_modes_zero_columns(self):
        n_dm, n_layer, n_low = [5, 4], [10, 6, 3], 4
        self._set_reconstruction(n_dm, n_layer, n_low)
        p_opt, _, _, _ = self._run()
        expected, target_cols = self._expected(n_dm, n_layer, n_low)
        self.assertEqual(p_opt.shape, (9, 19))
        # layer1 abs modes 0..3 -> cols 0..3; layer2 abs modes 2,3 -> local 0,1 -> cols 10,11;
        # layer3 (3 modes kept, abs 0..2) -> cols 16..18
        zero_cols = sorted(set(range(19)) - set(target_cols))
        self.assertEqual(zero_cols, [0, 1, 2, 3, 10, 11, 16, 17, 18])
        np.testing.assert_array_equal(p_opt[:, zero_cols], 0.0)
        # the kept columns are not trivially zero
        self.assertTrue(np.all(np.abs(p_opt[:, target_cols]).max(axis=0) > 0))
        np.testing.assert_allclose(p_opt, expected, rtol=RTOL, atol=ATOL)

    def test_truncation_lists_shorter_than_components(self):
        # components beyond the list are ignored (DM2 and layer3 dropped)
        n_dm, n_layer = [6], [8, 5]
        self._set_reconstruction(n_dm, n_layer, 0)
        p_opt, _, _, _ = self._run()
        expected, _ = self._expected([6, 0], [8, 5, 0], 0)
        self.assertEqual(p_opt.shape, (6, 13))
        np.testing.assert_allclose(p_opt, expected, rtol=RTOL, atol=ATOL)

    def test_returned_full_pms_are_first_source_only(self):
        self._set_reconstruction([5, 4], [10, 6, 3], 4)
        _, pm_full_dm, pm_full_layer, _ = self._run()
        self.assertEqual(pm_full_dm.shape, (1, 15, N_PUPIL_MODES))
        self.assertEqual(pm_full_layer.shape, (1, 25, N_PUPIL_MODES))
        np.testing.assert_array_equal(pm_full_dm[0], np.concatenate(self.pm_dm[0], axis=0))
        np.testing.assert_array_equal(pm_full_layer[0], np.concatenate(self.pm_layer[0], axis=0))

    def test_missing_layer_file_is_zeros(self):
        self._set_reconstruction([], [], 0)
        # layer2 missing for every source
        for i in range(1, 4):
            os.remove(self._path(i, layer_index=2))
        p_opt, _, pm_full_layer, _ = self._run()
        expected, _ = self._expected([], [], 0, zero_layer=1)
        np.testing.assert_allclose(p_opt, expected, rtol=RTOL, atol=ATOL)
        # the zero block (12:19 of the layer modes) has no effect: columns are exactly zero
        np.testing.assert_array_equal(p_opt[:, 12:19], 0.0)
        np.testing.assert_array_equal(pm_full_layer[0, 12:19], 0.0)

    def test_missing_dm_file_raises(self):
        self._set_reconstruction([], [], 0)
        os.remove(self._path(2, dm_index=2))
        with self.assertRaises(FileNotFoundError):
            self._run()

    def test_wrong_number_of_rows_raises(self):
        self._set_reconstruction([], [], 0)
        # DM2 has 7 modes (nmodes - start_mode); write 10 rows for source 3
        rng = np.random.default_rng(0)
        self._save(rng.standard_normal((10, N_PUPIL_MODES)).astype(np.float32), 3, dm_index=2)
        with self.assertRaises(ValueError):
            self._run()

    def test_wrong_number_of_rows_layer_raises(self):
        self._set_reconstruction([], [], 0)
        rng = np.random.default_rng(0)
        self._save(rng.standard_normal((5, N_PUPIL_MODES)).astype(np.float32), 1, layer_index=3)
        with self.assertRaises(ValueError):
            self._run()


if __name__ == '__main__':
    unittest.main()
