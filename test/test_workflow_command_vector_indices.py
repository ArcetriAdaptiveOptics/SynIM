"""Tests for the MORFEO command vector indices of the calibration workflow."""
import copy
import importlib.util
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
WORKFLOW_PATH = ROOT / 'calibration_workflows' / 'run_calib_workflow_morfeo.py'
PARAMS_PATH = ROOT / 'calibration_workflows' / 'params_morfeo_simplified.yml'

_spec = importlib.util.spec_from_file_location('run_calib_workflow_morfeo', WORKFLOW_PATH)
workflow = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(workflow)

command_vector_indices = workflow.command_vector_indices
MORFEOCalibrationWorkflow = workflow.MORFEOCalibrationWorkflow


def make_config(nmodes=(1260, 1026, 1026), start_modes=(0, 5, 2)):
    """Minimal config with all the sections needed by the method."""
    config = {}
    for i, (n, k) in enumerate(zip(nmodes, start_modes), start=1):
        config[f'dm{i}'] = {'nmodes': n}
        if k is not None:
            config[f'dm{i}']['start_mode'] = k
        config[f'slicer{i}'] = {'slice_args': [-1, -1], 'class': 'BaseSlicer'}
    config['modal_combination'] = {'plate_scale_idx': -1}
    config['iir_lgs'] = {'int_gain': [0.0, 0.25, 0.0, 0.25], 'n_modes': [1, 2, 3, 4]}
    return config


def make_workflow(config):
    """Workflow instance without running the constructor."""
    wf = MORFEOCalibrationWorkflow.__new__(MORFEOCalibrationWorkflow)
    wf.config = config
    wf.logs = []
    wf._log = wf.logs.append
    return wf


class TestCommandVectorIndicesFunction(unittest.TestCase):

    def test_reference_values(self):
        """Values hard-coded in params_morfeo_simplified.yml."""
        idx = command_vector_indices([1260, 1026, 1026], [0, 5, 2])
        self.assertEqual(idx['slices'], [[0, 1260], [1260, 2281], [2281, 3305]])
        self.assertEqual(idx['plate_scale_idx'], 2281)
        self.assertEqual(idx['iir_n_modes'], [3, 2278, 3, 1021])

    def test_morfeo_4519(self):
        idx = command_vector_indices([4519, 1026, 1026], [0, 5, 2])
        self.assertEqual(idx['slices'], [[0, 4519], [4519, 5540], [5540, 6564]])
        self.assertEqual(idx['plate_scale_idx'], 5540)
        self.assertEqual(idx['iir_n_modes'], [3, 5537, 3, 1021])

    def test_morfeo_4654(self):
        idx = command_vector_indices([4654, 1026, 1026], [0, 5, 2])
        self.assertEqual(idx['slices'], [[0, 4654], [4654, 5675], [5675, 6699]])
        self.assertEqual(idx['plate_scale_idx'], 5675)
        self.assertEqual(idx['iir_n_modes'], [3, 5672, 3, 1021])

    def test_total_length_invariant(self):
        for nmodes, starts in [([1260, 1026, 1026], [0, 5, 2]),
                               ([4519, 1026, 1026], [0, 5, 2]),
                               ([4654, 1026, 1026], [0, 5, 2]),
                               ([10, 7, 8], [1, 3, 0]),
                               ([3, 5, 3], [0, 5, 0])]:
            idx = command_vector_indices(nmodes, starts)
            self.assertEqual(sum(idx['iir_n_modes']), idx['slices'][2][1])
            # Slices are contiguous and start at 0
            self.assertEqual(idx['slices'][0][0], 0)
            self.assertEqual(idx['slices'][0][1], idx['slices'][1][0])
            self.assertEqual(idx['slices'][1][1], idx['slices'][2][0])

    def test_boundary_minimum_ok(self):
        """N1 = 3, N2 = 0, N3 = 3 is the smallest valid case."""
        idx = command_vector_indices([3, 5, 3], [0, 5, 0])
        self.assertEqual(idx['slices'], [[0, 3], [3, 3], [3, 6]])
        self.assertEqual(idx['plate_scale_idx'], 3)
        self.assertEqual(idx['iir_n_modes'], [3, 0, 3, 0])

    def test_dm1_too_few_modes(self):
        with self.assertRaises(ValueError):
            command_vector_indices([2, 1026, 1026], [0, 5, 2])

    def test_dm3_too_few_modes_with_start_mode(self):
        with self.assertRaises(ValueError):
            command_vector_indices([1260, 1026, 4], [0, 5, 2])

    def test_dm2_negative_modes(self):
        with self.assertRaises(ValueError):
            command_vector_indices([1260, 4, 1026], [0, 5, 2])


class TestUpdateCommandVectorIndicesMethod(unittest.TestCase):

    def test_writes_values_and_keeps_other_keys(self):
        config = make_config()
        wf = make_workflow(config)
        wf._update_command_vector_indices()
        self.assertEqual(config['slicer1']['slice_args'], [0, 1260])
        self.assertEqual(config['slicer2']['slice_args'], [1260, 2281])
        self.assertEqual(config['slicer3']['slice_args'], [2281, 3305])
        self.assertEqual(config['modal_combination']['plate_scale_idx'], 2281)
        self.assertEqual(config['iir_lgs']['n_modes'], [3, 2278, 3, 1021])
        self.assertEqual(config['iir_lgs']['int_gain'], [0.0, 0.25, 0.0, 0.25])
        self.assertEqual(config['slicer1']['class'], 'BaseSlicer')
        self.assertEqual(config['dm2']['nmodes'], 1026)

    def test_morfeo_4519(self):
        config = make_config(nmodes=(4519, 1026, 1026))
        make_workflow(config)._update_command_vector_indices()
        self.assertEqual(config['slicer3']['slice_args'], [5540, 6564])
        self.assertEqual(config['modal_combination']['plate_scale_idx'], 5540)
        self.assertEqual(config['iir_lgs']['n_modes'], [3, 5537, 3, 1021])

    def test_missing_start_mode_is_zero(self):
        config = make_config(start_modes=(0, None, 2))
        self.assertNotIn('start_mode', config['dm2'])
        make_workflow(config)._update_command_vector_indices()
        # N2 = 1026 instead of 1021
        self.assertEqual(config['slicer2']['slice_args'], [1260, 2286])
        self.assertEqual(config['modal_combination']['plate_scale_idx'], 2286)
        self.assertEqual(config['iir_lgs']['n_modes'], [3, 1260 - 3 + 1026, 3, 1021])

    def test_none_start_mode_is_zero(self):
        config = make_config(start_modes=(None, 5, 2))
        config['dm1']['start_mode'] = None
        make_workflow(config)._update_command_vector_indices()
        self.assertEqual(config['slicer1']['slice_args'], [0, 1260])

    def test_missing_section_leaves_config_unchanged(self):
        for section in ['dm1', 'dm2', 'dm3', 'slicer1', 'slicer2', 'slicer3',
                        'modal_combination', 'iir_lgs']:
            with self.subTest(section=section):
                config = make_config()
                del config[section]
                before = copy.deepcopy(config)
                wf = make_workflow(config)
                wf._update_command_vector_indices()
                self.assertEqual(config, before)
                self.assertEqual(len(wf.logs), 1)
                self.assertIn(section, wf.logs[0])
                self.assertIn('not updated', wf.logs[0])

    def test_nmodes_none_leaves_config_unchanged(self):
        for dm in ['dm1', 'dm2', 'dm3']:
            with self.subTest(dm=dm):
                config = make_config()
                config[dm]['nmodes'] = None
                before = copy.deepcopy(config)
                wf = make_workflow(config)
                wf._update_command_vector_indices()
                self.assertEqual(config, before)
                self.assertEqual(len(wf.logs), 1)
                self.assertIn('nmodes not set', wf.logs[0])

    def test_iir_three_blocks_leaves_config_unchanged(self):
        config = make_config()
        config['iir_lgs']['n_modes'] = [3, 100, 3]
        before = copy.deepcopy(config)
        wf = make_workflow(config)
        wf._update_command_vector_indices()
        self.assertEqual(config, before)
        self.assertEqual(len(wf.logs), 1)
        self.assertIn('4 blocks', wf.logs[0])

    def test_invalid_nmodes_raises(self):
        config = make_config(nmodes=(2, 1026, 1026))
        with self.assertRaises(ValueError):
            make_workflow(config)._update_command_vector_indices()

    def test_params_morfeo_simplified_yaml_is_consistent(self):
        """The hard-coded YAML values equal the computed ones."""
        with open(PARAMS_PATH) as f:
            config = yaml.safe_load(f)
        for key, n in [('dm1', 1260), ('dm2', 1026), ('dm3', 1026)]:
            self.assertIsNone(config[key]['nmodes'])
            config[key]['nmodes'] = n
        before = copy.deepcopy(config)
        wf = make_workflow(config)
        wf._update_command_vector_indices()
        for key in ['slicer1', 'slicer2', 'slicer3']:
            self.assertEqual(config[key]['slice_args'], before[key]['slice_args'])
        self.assertEqual(config['modal_combination']['plate_scale_idx'],
                         before['modal_combination']['plate_scale_idx'])
        self.assertEqual(config['iir_lgs']['n_modes'], before['iir_lgs']['n_modes'])
        self.assertEqual(config, before)
        self.assertIn('Command vector indices:', wf.logs[0])


if __name__ == '__main__':
    unittest.main()
