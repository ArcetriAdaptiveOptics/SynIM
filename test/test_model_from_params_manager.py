import os
import unittest
import numpy as np

from synim.params_utils import parse_params_file
from synim.registration.model import System

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


class TestSystemFromParamsManager(unittest.TestCase):
    """
    Read-only integration test: `System.from_params_manager` against the
    existing SynIM test yml configs, reusing `synim.params_utils`'s own
    extraction helpers (no modification to params_manager.py/params_utils.py).
    """

    def test_shift_test_config(self):
        config = parse_params_file(os.path.join(TEST_DIR, "params_scao_sh_shift_test.yml"))
        system = System.from_params_manager(config)

        self.assertEqual(list(system.wfss.keys()), ["sh"])
        self.assertEqual(list(system.dms.keys()), ["dm"])

        wfs = system.wfss["sh"]
        self.assertAlmostEqual(wfs.shift[0], 10.0)
        self.assertAlmostEqual(wfs.shift[1], 0.0)
        self.assertAlmostEqual(wfs.rotation, 0.0)
        self.assertAlmostEqual(wfs.magnification, 1.0)

        dm = system.dms["dm"]
        self.assertAlmostEqual(dm.height, 0.0)

        # On-axis NGS by default in this config.
        np.testing.assert_allclose(wfs.guide_star.position, (0.0, 0.0), atol=1e-8)
        self.assertTrue(np.isinf(wfs.guide_star.height))

        self.assertAlmostEqual(system.pixel_pitch, config["main"]["pixel_pitch"])

    def test_rot_test_config(self):
        config = parse_params_file(os.path.join(TEST_DIR, "params_scao_sh_rot_test.yml"))
        system = System.from_params_manager(config)

        wfs = system.wfss["sh"]
        self.assertAlmostEqual(wfs.rotation, 15.0)

    def test_dm_and_wfs_name_filtering(self):
        config = parse_params_file(os.path.join(TEST_DIR, "params_scao_sh_shift_test.yml"))
        system = System.from_params_manager(config, wfs_names=["sh"], dm_names=["dm"])
        self.assertEqual(set(system.wfss.keys()), {"sh"})
        self.assertEqual(set(system.dms.keys()), {"dm"})

    def test_local_params_are_computable_after_loading(self):
        config = parse_params_file(os.path.join(TEST_DIR, "params_scao_sh_shift_test.yml"))
        system = System.from_params_manager(config)
        for wfs_name, dm_name in system.pairs():
            shift, rotation, magnification, anam45, anam90 = system.local_params(wfs_name, dm_name)
            self.assertTrue(np.all(np.isfinite(shift)))


if __name__ == "__main__":
    unittest.main()
