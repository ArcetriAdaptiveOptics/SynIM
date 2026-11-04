"""
Check the slope ordering of the illumination weights in
ParamsManager._build_noise_covariance_impl.

The slopes of a WFS are two blocks (all subapertures of one axis, then all of
the other), so the weight of slope k and slope k + n must both be ill_k
(np.tile(ill, 2)), not interleaved (np.repeat(ill, 2)).

Checked, with a mocked illumination (1.0 everywhere except a few subapertures
at 0.3) and a mocked elongated precision matrix:
  - diagonal model: variance of slopes k and k + n is base / ill_k, others untouched;
  - elongated model: the precision matrix stays symmetric (its x-y couplings
    link i and i + n) and rows k, k + n are scaled by ill_k;
  - several WFSs with different numbers of valid subapertures.
"""
import unittest
from unittest import mock

import numpy as np

import specula
specula.init(device_idx=-1, precision=1)

import synim
synim.init(device_idx=-1, precision=1)

import synim.synim as synim_mod  # the module referenced as `synim` in params_manager
from synim.params_manager import ParamsManager

LOW = 0.3
BASE_SIGMA2_NM2 = 1.0


def _make_pm(n_valid_per_wfs, low_idx_per_wfs):
    """Minimal ParamsManager (no __init__) with mocked WFS parameters."""
    pm = ParamsManager.__new__(ParamsManager)
    pm.wfs_list = [{'name': f'lgs{i+1}', 'index': i + 1}
                   for i in range(len(n_valid_per_wfs))]
    pm.pup_mask = np.ones((4, 4))
    pm.params = {
        'main': {'pixel_pupil': 8, 'pixel_pitch': 0.1},
        'sh_lgs1': {'subap_on_diameter': 4, 'subap_npx': 6,
                    'sensor_pxscale': 0.5},
    }
    pm.lgs_recon_params = {'sigma2_in_nm2': BASE_SIGMA2_NM2}

    def get_wfs_params(wfs_type, idx, xp_local=np):
        n = n_valid_per_wfs[idx - 1]
        return {'idx_valid_sa': np.arange(n * 2).reshape(n, 2),
                'wfs_nsubaps': 4, 'wfs_rotation': 0.0,
                'wfs_translation': (0.0, 0.0), 'wfs_magnification': 1.0}

    pm.get_wfs_params = get_wfs_params

    ills = []
    for n, low in zip(n_valid_per_wfs, low_idx_per_wfs):
        ill = np.ones(n)
        ill[list(low)] = LOW
        ills.append(ill)
    pm._ills = ills
    return pm


def _run(pm, n_valid, use_elongation, elong_blocks=None):
    ill_iter = iter(pm._ills)
    n_total = 2 * sum(n_valid)
    patches = [mock.patch.object(
        synim_mod, 'compute_subaperture_illumination',
        side_effect=lambda **kw: next(ill_iter).copy())]
    if elong_blocks is not None:
        patches.append(mock.patch.object(
            pm, '_build_elongated_noise_covariance',
            create=True, return_value=_blockdiag(elong_blocks)))
    for p in patches:
        p.start()
    try:
        return pm._build_noise_covariance_impl(
            False, 'lgs', None, n_total, None, None, use_elongation,
            None, True)
    finally:
        mock.patch.stopall()


def _blockdiag(blocks):
    n = sum(b.shape[0] for b in blocks)
    out = np.zeros((n, n))
    o = 0
    for b in blocks:
        out[o:o + b.shape[0], o:o + b.shape[0]] = b
        o += b.shape[0]
    return out


def _elong_block(n, a=2.0, b=0.5, c=3.0):
    """Per-WFS precision with x-y coupling between slope i and i + n."""
    eye = np.eye(n)
    return np.block([[a * eye, b * eye], [b * eye, c * eye]])


class TestNoiseCovIllumination(unittest.TestCase):

    def test_diagonal_same_weight_for_both_axes(self):
        n = 5
        pm = _make_pm([n], [[1, 3]])
        C, diag = _run(pm, [n], use_elongation=False)
        w = diag['illumination_concat']
        expected = np.tile(pm._ills[0], 2)
        np.testing.assert_allclose(w, expected)

        prec = np.diag(C)
        base = prec[0]  # fully illuminated slope
        low = {1, 3}
        for k in range(n):
            f = LOW if k in low else 1.0
            self.assertAlmostEqual(prec[k] / base, f)
            self.assertAlmostEqual(prec[k + n] / base, f)

    def test_elongated_precision_symmetric_and_row_scaled(self):
        n = 5
        pm = _make_pm([n], [[1, 3]])
        P = _elong_block(n)
        C, _ = _run(pm, [n], use_elongation=True, elong_blocks=[P])
        np.testing.assert_allclose(C, C.T, atol=1e-12)
        w = np.tile(pm._ills[0], 2)
        np.testing.assert_allclose(C, w[:, None] * P)

    def test_multi_wfs_different_sizes(self):
        n_valid = [4, 6]
        pm = _make_pm(n_valid, [[0], [2, 5]])
        blocks = [_elong_block(n) for n in n_valid]
        C, diag = _run(pm, n_valid, use_elongation=True, elong_blocks=blocks)
        expected_w = np.concatenate([np.tile(i, 2) for i in pm._ills])
        np.testing.assert_allclose(diag['illumination_concat'], expected_w)
        np.testing.assert_allclose(C, C.T, atol=1e-12)
        np.testing.assert_allclose(C, expected_w[:, None] * _blockdiag(blocks))

        # diagonal model on the same geometry
        pm = _make_pm(n_valid, [[0], [2, 5]])
        Cd, _ = _run(pm, n_valid, use_elongation=False)
        prec = np.diag(Cd)
        base = prec.max()
        np.testing.assert_allclose(prec / base, expected_w)


if __name__ == '__main__':
    unittest.main()
