"""
Launcher in the telescope frame (``absolute_coordinates: True``) in
ParamsManager._build_elongated_noise_covariance.

With the flag, the launcher ``tel_position`` is in the telescope frame and
SynIM passes the WFS rotation (sh_lgsN ``rotation``, fallback
``rotAnglePhInDeg``, default 0) to SPECULA calc_noise_cov_elong, which rotates
the launcher into the WFS frame: p_wfs = R(+theta) p_tel. Without the flag the
call is unchanged (no ``wfs_rotation_in_deg`` keyword), so that older SPECULA
versions keep working.

Physical check: rotating the WFS by +90 deg maps a launcher at +x of the
telescope onto +y of the WFS, so the precision matrix must equal the one of a
launcher at +y given in the WFS frame (no flag). No SPECULA simulation is run.

The tests that call the real calc_noise_cov_elong with a rotation need a
SPECULA version with wfs_rotation_in_deg (SPECULA #892): they are skipped
with an older one. The others replace calc_noise_cov_elong and always run.
"""
import contextlib
import functools
import inspect
import unittest
from unittest import mock

import numpy as np

import synim.params_manager as pm_module
from synim.params_manager import ParamsManager

SPECULA_HAS_ROTATION = 'wfs_rotation_in_deg' in \
    inspect.signature(pm_module.calc_noise_cov_elong).parameters
needs_rotation = unittest.skipUnless(
    SPECULA_HAS_ROTATION, 'SPECULA calc_noise_cov_elong() without wfs_rotation_in_deg (#892)')

N_PIX = 48
N_SUBAPS = 8
DIAMETER_M = 8.0
OFFSET_M = DIAMETER_M / 2


def _valid_subaps():
    """Asymmetric round mask, valid subapertures as [column, row] pairs."""
    i, j = np.mgrid[:N_SUBAPS, :N_SUBAPS]
    mask = np.hypot(i - 3.5, j - 3.5) < 4.2
    mask[0, 2] = False
    return np.array([(jj, ii) for ii in range(N_SUBAPS) for jj in range(N_SUBAPS)
                     if mask[ii, jj]])


def _make_pm(launcher_xy_m, absolute=False, wfs_extra=None):
    """Minimal ParamsManager (no __init__) for one LGS WFS."""
    idx_valid_sa = _valid_subaps()
    sh = {'subap_on_diameter': N_SUBAPS, 'subap_npx': 8, 'sensor_pxscale': 0.5,
          'laser_launch_tel_ref': 'launcher1'}
    sh.update(wfs_extra or {})
    launcher = {'tel_position': [launcher_xy_m[0], launcher_xy_m[1], 0.0],
                'beacon_focus': 90000.0, 'spot_size': 1.0}
    if absolute:
        launcher['absolute_coordinates'] = True
    pm = ParamsManager.__new__(ParamsManager)
    pm.params = {
        'main': {'pixel_pupil': N_PIX, 'pixel_pitch': DIAMETER_M / N_PIX},
        'sh_lgs1': sh,
        'launcher1': launcher,
    }
    # tGparameter = 0: no truncated subapertures (divides by the elongation)
    pm.lgs_recon_params = {'sigma2_in_nm2': 1e6, 'naThicknessInM': 10000.0,
                           'tGparameter': 0.0}
    pm.get_wfs_params = lambda wfs_type, idx, xp_local=np: {
        'idx_valid_sa': idx_valid_sa, 'wfs_nsubaps': N_SUBAPS, 'wfs_fov_arcsec': 4.0}
    pm._count_wfs = lambda wfs_type: 1
    return pm


def _precision(launcher_xy_m, absolute=False, wfs_extra=None):
    pm = _make_pm(launcher_xy_m, absolute, wfs_extra)
    return np.asarray(pm._build_elongated_noise_covariance('lgs'), dtype=float)


def _recording_wrapper(real, calls):
    """Real function that records its kwargs; functools.wraps keeps the signature."""
    @functools.wraps(real)
    def wrapper(*args, **kwargs):
        calls.append(kwargs)
        return real(*args, **kwargs)
    return wrapper


def _recording_new_specula(calls):
    """
    calc_noise_cov_elong of a SPECULA with wfs_rotation_in_deg (in its
    signature), recording its kwargs: the rotation is not applied, the
    result is the one of the installed function without it.
    """
    real = pm_module.calc_noise_cov_elong

    def new_calc_noise_cov_elong(*args, wfs_rotation_in_deg=None, **kwargs):
        calls.append(dict(kwargs, wfs_rotation_in_deg=wfs_rotation_in_deg))
        return real(*args, **kwargs)
    return new_calc_noise_cov_elong


class TestNoiseCovAbsoluteLauncher(unittest.TestCase):

    @needs_rotation
    def test_rotation_90_equals_launcher_on_y(self):
        absolute = _precision((OFFSET_M, 0.0), True, {'rotAnglePhInDeg': 90.0})
        in_wfs_frame = _precision((0.0, OFFSET_M))
        np.testing.assert_allclose(absolute, in_wfs_frame, rtol=1e-5,
                                   atol=1e-6 * np.abs(in_wfs_frame).max())

    @needs_rotation
    def test_rotation_30_equals_launcher_in_wfs_frame(self):
        from specula.data_objects.laser_launch_telescope import launcher_in_wfs_frame
        pos, _ = launcher_in_wfs_frame([OFFSET_M, 0.0, 0.0], [0.0, 0.0], 30.0)
        absolute = _precision((OFFSET_M, 0.0), True, {'rotAnglePhInDeg': 30.0})
        in_wfs_frame = _precision((pos[0], pos[1]))
        np.testing.assert_allclose(absolute, in_wfs_frame, rtol=1e-5,
                                   atol=1e-6 * np.abs(in_wfs_frame).max())

    @needs_rotation
    def test_rotation_key_has_priority_over_rotAnglePhInDeg(self):
        absolute = _precision((OFFSET_M, 0.0), True,
                              {'rotation': 90.0, 'rotAnglePhInDeg': 0.0})
        in_wfs_frame = _precision((0.0, OFFSET_M))
        np.testing.assert_allclose(absolute, in_wfs_frame, rtol=1e-5,
                                   atol=1e-6 * np.abs(in_wfs_frame).max())

    @needs_rotation
    def test_not_vacuous(self):
        """The rotation changes the result: absolute != unrotated launcher."""
        absolute = _precision((OFFSET_M, 0.0), True, {'rotAnglePhInDeg': 90.0})
        unrotated = _precision((OFFSET_M, 0.0))
        scale = np.abs(unrotated).max()
        self.assertGreater(np.abs(absolute - unrotated).max(), 0.1 * scale)

    @needs_rotation
    def test_rotation_0_equals_no_flag(self):
        for extra in ({'rotAnglePhInDeg': 0.0}, {}):
            with self.subTest(extra=extra):
                np.testing.assert_allclose(
                    _precision((OFFSET_M, 0.0), True, extra),
                    _precision((OFFSET_M, 0.0)), rtol=1e-6, atol=0)

    def test_no_flag_does_not_pass_rotation_keyword(self):
        calls = []
        wrapper = _recording_wrapper(pm_module.calc_noise_cov_elong, calls)
        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.object(pm_module, 'calc_noise_cov_elong', wrapper))
            # a rotation in the WFS config must be ignored without the flag
            _precision((OFFSET_M, 0.0), False, {'rotAnglePhInDeg': 90.0})
        self.assertEqual(len(calls), 1)
        self.assertNotIn('wfs_rotation_in_deg', calls[0])

    def test_flag_passes_wfs_rotation(self):
        calls = []
        wrapper = _recording_new_specula(calls)
        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.object(pm_module, 'calc_noise_cov_elong', wrapper))
            _precision((OFFSET_M, 0.0), True, {'rotAnglePhInDeg': 30.0})
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]['wfs_rotation_in_deg'], 30.0)

    def test_flag_default_rotation_is_zero(self):
        calls = []
        wrapper = _recording_new_specula(calls)
        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.object(pm_module, 'calc_noise_cov_elong', wrapper))
            _precision((OFFSET_M, 0.0), True)
        self.assertEqual(calls[0]['wfs_rotation_in_deg'], 0.0)

    def test_flag_with_old_specula_raises(self):
        def old_calc_noise_cov_elong(diameter_in_m, zenith_angle_in_deg, na_thickness_in_m,
                                     launcher_coord_in_m, sub_aps_index, n_sub_aps,
                                     sub_aps_fov, sh_spot_fwhm, sigma_noise2, t_g_parameter,
                                     h_in_m=None, user_pofile_xy=None, theta=None,
                                     only_diag=False, eta_is_not_one=False, display=False):
            raise AssertionError('must not be called')

        pm = _make_pm((OFFSET_M, 0.0), True)
        with contextlib.ExitStack() as stack:
            stack.enter_context(
                mock.patch.object(pm_module, 'calc_noise_cov_elong', old_calc_noise_cov_elong))
            with self.assertRaises(ValueError):
                pm._build_elongated_noise_covariance('lgs')

    def test_no_flag_with_old_specula_works(self):
        real = pm_module.calc_noise_cov_elong

        def old_calc_noise_cov_elong(*args, **kwargs):
            self.assertNotIn('wfs_rotation_in_deg', kwargs)
            return real(*args, **kwargs)

        pm = _make_pm((OFFSET_M, 0.0), False)
        with contextlib.ExitStack() as stack:
            stack.enter_context(
                mock.patch.object(pm_module, 'calc_noise_cov_elong', old_calc_noise_cov_elong))
            precision = pm._build_elongated_noise_covariance('lgs')
        self.assertTrue(np.all(np.isfinite(precision)))


if __name__ == '__main__':
    unittest.main()
