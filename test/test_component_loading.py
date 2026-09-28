"""
Loading of DM / layer influence functions in ParamsManager:
- load_influence_functions builds the 3D array with the requested array
  module and dtype (on GPU the cube is built directly on the device)
- get_component_params(use_cache=False) neither reads nor writes the cache

The file I/O is replaced by in-memory objects.
"""
import unittest
from unittest import mock

import numpy as np

import synim.params_utils as params_utils
import synim.params_manager as params_manager
from synim.utils import dm2d_to_3d


def _mask(n=24, radius=10):
    yy, xx = np.mgrid[:n, :n]
    return (np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5) < radius).astype(np.float32)


class _FakeIFunc:
    def __init__(self, influence_function, mask):
        self.influence_function = influence_function
        self.mask_inf_func = mask


class _FakeCalibManager:
    def filename(self, kind, tag):
        return f'/nonexistent/{kind}/{tag}.fits'


class TestLoadInfluenceFunctions(unittest.TestCase):

    def setUp(self):
        self.mask = _mask()
        n_valid = int(self.mask.sum())
        self.ifunc_2d = np.random.default_rng(0).standard_normal((6, n_valid)).astype(np.float32)

    def _load(self, **kwargs):
        fake = _FakeIFunc(self.ifunc_2d.copy(), self.mask.copy())
        with mock.patch.object(params_utils.IFunc, 'restore', return_value=fake):
            return params_utils.load_influence_functions(
                _FakeCalibManager(), {'ifunc_tag': 'test'}, pixel_pupil=24, **kwargs)

    def test_default_builds_on_cpu(self):
        dm_array, dm_mask = self._load()
        expected = dm2d_to_3d(self.ifunc_2d, self.mask, xp_local=np,
                              float_dtype_local=params_utils.cpu_float_dtype)
        self.assertIsInstance(dm_array, np.ndarray)
        np.testing.assert_array_equal(dm_array, expected)
        np.testing.assert_array_equal(dm_mask, self.mask)

    def test_requested_module_and_dtype(self):
        calls = []
        original = params_utils.dm2d_to_3d

        def spy(*args, **kwargs):
            calls.append((kwargs.get('xp_local'), kwargs.get('float_dtype_local')))
            return original(*args, **kwargs)

        with mock.patch.object(params_utils, 'dm2d_to_3d', side_effect=spy):
            dm_array, _ = self._load(xp_local=np, float_dtype_local=np.float64)
        self.assertEqual(calls, [(np, np.float64)])
        self.assertEqual(dm_array.dtype, np.float64)


class TestGetComponentParamsCache(unittest.TestCase):

    def setUp(self):
        self.pm = params_manager.ParamsManager.__new__(params_manager.ParamsManager)
        self.pm.params = {'dm1': {'ifunc_tag': 'test', 'height': 0.0, 'rotation': 0.0}}
        self.pm.cm = _FakeCalibManager()
        self.pm.pixel_pupil = 24
        self.pm.verbose = False
        self.pm.dm_cache = {}
        self.loads = []
        mask = _mask()
        cube = np.ones((24, 24, 3), dtype=np.float32)

        def fake_load(cm, params, pixel_pupil, **kwargs):
            self.loads.append(kwargs)
            return cube.copy(), mask.copy()

        patcher = mock.patch.object(params_manager, 'load_influence_functions',
                                    side_effect=fake_load)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_cache_used_by_default(self):
        first = self.pm.get_component_params(1, xp_local=np, float_dtype_local=np.float32)
        second = self.pm.get_component_params(1, xp_local=np, float_dtype_local=np.float32)
        self.assertIs(first, second)
        self.assertEqual(len(self.loads), 1)
        self.assertEqual(len(self.pm.dm_cache), 1)

    def test_no_cache(self):
        first = self.pm.get_component_params(1, xp_local=np, float_dtype_local=np.float32,
                                             use_cache=False)
        second = self.pm.get_component_params(1, xp_local=np, float_dtype_local=np.float32,
                                              use_cache=False)
        self.assertIsNot(first, second)
        self.assertEqual(len(self.loads), 2)
        self.assertEqual(self.pm.dm_cache, {})
        np.testing.assert_array_equal(first['dm_array'], second['dm_array'])

    def test_no_cache_does_not_read_existing_entry(self):
        cached = self.pm.get_component_params(1, xp_local=np, float_dtype_local=np.float32)
        fresh = self.pm.get_component_params(1, xp_local=np, float_dtype_local=np.float32,
                                             use_cache=False)
        self.assertIsNot(cached, fresh)
        self.assertEqual(len(self.loads), 2)

    def test_array_module_passed_to_loading(self):
        self.pm.get_component_params(1, xp_local=np, float_dtype_local=np.float64,
                                     use_cache=False)
        self.assertIs(self.loads[-1]['xp_local'], np)
        self.assertIs(self.loads[-1]['float_dtype_local'], np.float64)


if __name__ == '__main__':
    unittest.main()
