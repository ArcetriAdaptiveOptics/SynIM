"""
compute_tomographic_reconstructor(mmse_on_device=True): the assembled IM and the
noise covariance are requested on the compute device, the MMSE runs there, and
the reconstructor (back on the host) matches the default CPU path.

The IMs, covariances and noise covariance are mocked with small random
matrices; compute_mmse_reconstructor is the real one.
"""
import unittest
from unittest import mock

import numpy as np

import synim
from synim import cpu_float_dtype
from synim.params_manager import ParamsManager

N_SLOPES, N_MODES = 60, 20


def _matrices():
    rng = np.random.default_rng(5)
    im = rng.standard_normal((N_SLOPES, N_MODES)).astype(cpu_float_dtype)
    b = rng.standard_normal((N_MODES, N_MODES))
    c_atm_inv = (b @ b.T / N_MODES + np.eye(N_MODES)).astype(cpu_float_dtype)
    b = rng.standard_normal((N_SLOPES, N_SLOPES)) * 0.1
    c_noise_inv = (b @ b.T / N_SLOPES + np.eye(N_SLOPES)).astype(cpu_float_dtype)
    return im, c_atm_inv, c_noise_inv


def _on(to_host, a):
    return a if to_host else synim.to_xp(synim.xp, a)


class TestTomographicReconstructorOnDevice(unittest.TestCase):

    def _run(self, n_low_modes_cut=0, **kwargs):
        im, c_atm_inv, c_noise_inv = _matrices()
        pm = ParamsManager.__new__(ParamsManager)
        pm.verbose = False
        pm.im_dir = pm.rec_dir = pm.cov_dir = None
        pm.lgs_recon_params = {'noise_elong_model': True}
        mocks = {
            'compute_interaction_matrices': mock.Mock(return_value={}),
            'assemble_interaction_matrices': mock.Mock(
                side_effect=lambda **kw: (_on(kw.get('to_host', True), im.copy()),
                                          [N_SLOPES], [list(range(N_MODES))], [1])),
            'compute_covariance_matrices': mock.Mock(
                return_value={'C_atm_blocks': None, 'component_indices': None}),
            'assemble_covariance_matrix': mock.Mock(return_value=c_atm_inv.copy()),
            'build_noise_covariance': mock.Mock(
                side_effect=lambda **kw: _on(kw.get('to_host', True), c_noise_inv.copy())),
            '_get_low_modes_cut_if_applicable': mock.Mock(return_value=n_low_modes_cut),
        }
        for name, m in mocks.items():
            setattr(pm, name, m)
        result = pm.compute_tomographic_reconstructor(r0=0.15, L0=25.0, wfs_type='lgs',
                                                      component_type='layer', **kwargs)
        return result, mocks

    def test_same_reconstructor_as_cpu_path(self):
        self._check(n_low_modes_cut=0)

    def test_low_mode_cut(self):
        # columns selected on the device, dropped modes restored as zero rows
        R = self._check(n_low_modes_cut=2)
        np.testing.assert_array_equal(R[:2], 0)
        self.assertTrue(np.all(np.any(R[2:] != 0, axis=1)))

    def _check(self, n_low_modes_cut):
        expected, mocks = self._run(n_low_modes_cut)
        self.assertTrue(mocks['assemble_interaction_matrices'].call_args.kwargs['to_host'])
        self.assertTrue(mocks['build_noise_covariance'].call_args.kwargs['to_host'])

        result, mocks = self._run(n_low_modes_cut, mmse_on_device=True)
        self.assertFalse(mocks['assemble_interaction_matrices'].call_args.kwargs['to_host'])
        self.assertFalse(mocks['build_noise_covariance'].call_args.kwargs['to_host'])
        self.assertIsInstance(result['im_full'], synim.xp.ndarray)
        self.assertIsInstance(result['C_noise'], synim.xp.ndarray)

        R = result['reconstructor']
        self.assertIsInstance(R, np.ndarray)
        R_ref = expected['reconstructor']
        self.assertEqual(R.shape, (N_MODES, N_SLOPES))
        np.testing.assert_allclose(R, R_ref, rtol=0, atol=1e-4 * np.abs(R_ref).max())
        return R


if __name__ == '__main__':
    unittest.main()
