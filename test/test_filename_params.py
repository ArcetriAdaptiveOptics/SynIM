"""
ParamsManager builds the IM / PM filenames from the parameters it already
parsed, without reading the YAML file again, and the names are the same as
those built from the file.
"""
import os
import unittest
from unittest import mock

import synim.params_utils as params_utils
from synim.params_manager import ParamsManager
from synim.params_utils import generate_im_filename

HERE = os.path.dirname(os.path.abspath(__file__))
YAML_FILE = os.path.join(HERE, 'params_scao_sh_test.yml')
CALIB_DIR = os.path.join(HERE, 'calib')


class TestFilenamesFromParsedParams(unittest.TestCase):

    def setUp(self):
        self.pm = ParamsManager(YAML_FILE, root_dir=CALIB_DIR, verbose=False)

    def test_im_filename_does_not_parse_file(self):
        kwargs = dict(wfs_type='ngs', wfs_index=1, dm_index=1)
        expected = generate_im_filename(YAML_FILE, **kwargs)
        with mock.patch.object(params_utils, 'parse_params_file',
                               side_effect=AssertionError('YAML parsed again')):
            name = self.pm.generate_im_filename(**kwargs)
        self.assertEqual(name, expected)


if __name__ == '__main__':
    unittest.main()
