"""
Tests for the interleaved-slopes guard in SynIM.

SynIM orders slopes as two blocks (all subapertures of one axis, then the
other), so ParamsManager.get_wfs_params must refuse WFSs whose slopec uses
interleave=True, and find_slopec_key must locate the right slopec section.
"""
import copy
import os
import unittest

import yaml

import specula
specula.init(device_idx=-1, precision=1)

from synim.params_manager import ParamsManager
from synim.params_utils import find_slopec_key

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_YML = os.path.join(TEST_DIR, 'params_scao_sh_test.yml')


class TestFindSlopecKey(unittest.TestCase):

    def test_sh_prefix(self):
        params = {'slopec_lgs1': {}, 'slopec': {}}
        self.assertEqual(find_slopec_key('sh_lgs1', params), 'slopec_lgs1')

    def test_pyramid_prefix(self):
        params = {'slopec_ngs1': {}}
        self.assertEqual(find_slopec_key('pyramid_ngs1', params), 'slopec_ngs1')

    def test_numeric_index(self):
        params = {'slopec1': {}, 'slopec2': {}}
        self.assertEqual(find_slopec_key('sh1', params), 'slopec1')
        self.assertEqual(find_slopec_key('sh2', params), 'slopec2')

    def test_fallback_to_generic_slopec(self):
        self.assertEqual(find_slopec_key('sh', {'slopec': {}}), 'slopec')
        # specific section missing -> generic one
        self.assertEqual(find_slopec_key('sh_lgs1', {'slopec': {}}), 'slopec')

    def test_none_when_nothing_matches(self):
        self.assertIsNone(find_slopec_key('sh_lgs1', {'slopec_ngs1': {}}))
        self.assertIsNone(find_slopec_key('sh', {}))

    def test_specific_key_wins_over_generic(self):
        params = {'slopec': {}, 'slopec_lgs1': {}}
        self.assertEqual(find_slopec_key('sh_lgs1', params), 'slopec_lgs1')


class TestInterleaveGuard(unittest.TestCase):

    def setUp(self):
        with open(BASE_YML) as f:
            self.params = yaml.safe_load(f)

    def _manager(self, params):
        return ParamsManager(copy.deepcopy(params),
                             root_dir=os.path.join(TEST_DIR, 'calib'),
                             require_subapdata_file=False)

    def _two_wfs_params(self, interleave_1, interleave_2):
        """Two NGS SH WFSs (sh_ngs1/slopec_ngs1, sh_ngs2/slopec_ngs2), each with its own slopec."""
        p = self.params
        p['sh_ngs1'] = p.pop('sh')
        p['sh_ngs2'] = copy.deepcopy(p['sh_ngs1'])
        p['slopec_ngs1'] = p.pop('slopec')
        p['slopec_ngs2'] = copy.deepcopy(p['slopec_ngs1'])
        p['slopec_ngs1']['interleave'] = interleave_1
        p['slopec_ngs2']['interleave'] = interleave_2
        return p

    def test_interleave_true_raises(self):
        self.params['slopec']['interleave'] = True
        with self.assertRaises(NotImplementedError) as ctx:
            self._manager(self.params).get_wfs_params(wfs_type='ngs')
        self.assertIn('slopec', str(ctx.exception))

    def test_interleave_false_ok(self):
        self.params['slopec']['interleave'] = False
        self._manager(self.params).get_wfs_params(wfs_type='ngs')

    def test_interleave_absent_ok(self):
        self.assertNotIn('interleave', self.params['slopec'])
        self._manager(self.params).get_wfs_params(wfs_type='ngs')

    def test_guard_is_per_wfs(self):
        p = self._two_wfs_params(interleave_1=True, interleave_2=False)
        mgr = self._manager(p)
        # WFS 2 is not interleaved: must work although WFS 1 is
        mgr.get_wfs_params(wfs_type='ngs', wfs_index=2)
        with self.assertRaises(NotImplementedError) as ctx:
            mgr.get_wfs_params(wfs_type='ngs', wfs_index=1)
        self.assertIn('slopec_ngs1', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
