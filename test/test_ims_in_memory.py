"""
Interaction matrices kept in memory: compute_interaction_matrices(in_memory=True)
writes nothing, and assemble_interaction_matrices(im_arrays=...) gives the
same result as the file-based workflow.
"""
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

import synim
import synim.params_manager as params_manager
from synim.params_manager import ParamsManager

HERE = os.path.dirname(os.path.abspath(__file__))
YAML_FILE = os.path.join(HERE, 'params_scao_sh_test.yml')
CALIB_DIR = os.path.join(HERE, 'calib')


def _all_subaps(*args, **kwargs):
    # valid subapertures of the 8x8 SH of the test configuration (the
    # SubapData file is generated only by test_sh_intmat)
    n = 8
    yy, xx = np.mgrid[:n, :n]
    keep = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5) < n / 2
    return np.column_stack(np.nonzero(keep)).astype(np.int64)


class TestInteractionMatricesInMemory(unittest.TestCase):

    def setUp(self):
        patcher = mock.patch.object(params_manager, 'find_subapdata', side_effect=_all_subaps)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def _pm(self):
        pm = ParamsManager(YAML_FILE, root_dir=CALIB_DIR, verbose=False)
        return pm, pm.wfs_list[0]['name']

    def test_same_as_files(self):
        pm, wfs = self._pm()
        files_dir = os.path.join(self.tmp.name, 'files')
        pm.compute_interaction_matrices(output_im_dir=files_dir, output_rec_dir=files_dir,
                                        wfs_type=wfs, component_type='dm', overwrite=True)
        expected = pm.assemble_interaction_matrices(wfs_type=wfs, output_im_dir=files_dir,
                                                    component_type='dm')

        pm, wfs = self._pm()
        unused_dir = os.path.join(self.tmp.name, 'unused')
        ims = pm.compute_interaction_matrices(output_im_dir=unused_dir, output_rec_dir=unused_dir,
                                              wfs_type=wfs, component_type='dm', in_memory=True)
        self.assertFalse(os.path.exists(unused_dir))
        self.assertEqual(sorted(ims), sorted(os.listdir(files_dir)))
        for im in ims.values():
            self.assertIsInstance(im, synim.xp.ndarray)

        result = pm.assemble_interaction_matrices(wfs_type=wfs, component_type='dm', im_arrays=ims)
        np.testing.assert_array_equal(result[0], expected[0])
        self.assertEqual(result[1:], expected[1:])

    def test_missing_matrix(self):
        pm, wfs = self._pm()
        with self.assertRaises(KeyError):
            pm.assemble_interaction_matrices(wfs_type=wfs, component_type='dm', im_arrays={})


if __name__ == '__main__':
    unittest.main()
