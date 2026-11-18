"""
Physical and numerical tests for spider arm extraction and directional
noise covariance weighting (extract_spider_arms,
compute_spider_subaperture_weights_and_normals,
update_wfs_noise_cov_inv_for_spider, apply_spider_to_noise_covariance).

1. Arm extraction on an annular pupil with a 6-arm spider, on a segmented
   hexagonal-edge pupil (mimicking E-ELT M1 concavities), and on a spiderless
   annulus.
2. Physical alignment check: for a phase screen with piston steps across the
   petals, the slope vector (s_0, s_1) computed by synim.interaction_matrix on
   subapertures along each spider arm must be parallel to the arm's unit normal
   returned by compute_spider_subaperture_weights_and_normals, for both
   specula_convention=True and False, rotated/translated/magnified WFSs, and
   with/without idx_valid_sa selection.
3. Numerical check of the 2x2 Sherman-Morrison update: exact agreement with
   explicit 2x2 block inversion (C_j + sum_k sigma_{j,k}^2 n_k n_k^T)^(-1),
   preservation of tangential precision, bit-for-bit symmetry, and rejection
   of non-2x2-block-diagonal matrices.
4. End-to-end integration with apply_spider_to_noise_covariance, _sparse_rows
   (symmetric=True), and compute_mmse_reconstructor (petal step response
   decreases monotonically with alpha).
"""
import unittest

import numpy as np
from scipy.ndimage import label

import synim
if synim.xp is None:
    synim.init(device_idx=-1, precision=1)
import synim.synim as synim_core
from synim.params_utils import _sparse_rows, apply_spider_to_noise_covariance, compute_mmse_reconstructor
from synim.utils import (
    compute_spider_subaperture_weights_and_normals,
    detect_spider_arms_azimuthal,
    update_wfs_noise_cov_inv_for_spider,
)


def _make_six_arm_pupil_mask(
    n_pix=240, obsratio=0.283, spider_width_px=3.0, serrated_edges=False
):
    """Annular pupil mask with 3 spider lines (6 arms) at 0, 60, 120 deg."""
    yy, xx = np.mgrid[:n_pix, :n_pix]
    cy = (n_pix - 1) / 2.0
    cx = (n_pix - 1) / 2.0
    dy = yy - cy
    dx = xx - cx
    r_max = n_pix / 2.0 - 1.0
    r_min = obsratio * r_max
    r = np.hypot(dx, dy)
    theta = np.arctan2(dy, dx)

    if serrated_edges:
        # Simulate hexagonal M1 segment notches along inner and outer edges
        r_out_edge = r_max - 2.5 * np.abs(np.sin(24.0 * theta))
        r_in_edge = r_min + 2.0 * np.abs(np.sin(12.0 * theta))
        annulus = (r <= r_out_edge) & (r >= r_in_edge)
    else:
        annulus = (r <= r_max) & (r >= r_min)

    spider = np.zeros_like(annulus)
    for angle_deg in (0.0, 60.0, 120.0):
        ang = np.deg2rad(angle_deg)
        dist_line = np.abs(-np.sin(ang) * dx + np.cos(ang) * dy)
        spider |= dist_line <= (spider_width_px / 2.0)

    mask = (annulus & (~spider)).astype(np.float32)
    return mask, annulus


class _MockParamsManager:
    """Minimal mock of ParamsManager for testing apply_spider_to_noise_covariance."""

    def __init__(self, pup_mask, wfs_configs):
        self.pup_mask = pup_mask
        self.wfs_list = [
            {'name': f'sh_lgs{i + 1}', 'type': 'sh', 'index': str(i + 1), 'config': cfg}
            for i, cfg in enumerate(wfs_configs)
        ]
        self._wfs_configs = wfs_configs
        self.lgs_recon_params = {'noise_elong_model': False}

    def get_wfs_params(self, wfs_type, wfs_index, xp_local=np):
        return self._wfs_configs[wfs_index - 1]


class TestSpiderNoiseCovariance(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.n_pix = 240
        cls.n_subaps = 30
        cls.pup_diam_m = 38.5
        cls.fov_arcsec = 16.0
        cls.pup_mask, cls.annulus = _make_six_arm_pupil_mask(
            n_pix=cls.n_pix, obsratio=0.283, spider_width_px=3.2
        )
        cls.serrated_mask, cls.serrated_annulus = _make_six_arm_pupil_mask(
            n_pix=cls.n_pix, obsratio=0.283, spider_width_px=3.2, serrated_edges=True
        )
        labeled_petals, n_petals = label(cls.pup_mask > 0.5)
        assert n_petals == 6, f"Expected 6 petals, got {n_petals}"
        cls.petal_modes = np.stack(
            [(labeled_petals == p).astype(np.float32) for p in range(1, 7)],
            axis=-1,
        )
        cls.dm_mask = cls.annulus.astype(np.float32)

    def test_detect_spider_arms_on_serrated_and_spiderless_masks(self):
        for mask in (self.pup_mask, self.serrated_mask):
            arms, _, _ = detect_spider_arms_azimuthal(mask)
            self.assertEqual(len(arms), 6)
            for arm, exp_ang in zip(arms, [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]):
                diff_deg = abs((arm['azimuth_deg'] - exp_ang + 180.0) % 360.0 - 180.0)
                self.assertLess(diff_deg, 0.4)

        # Spiderless mask returns 0 arms
        arms_empty, _, _ = detect_spider_arms_azimuthal(self.annulus.astype(np.float32))
        self.assertEqual(len(arms_empty), 0)

    def test_detect_spider_arms_azimuthal_unequal_widths_and_4arms(self):
        """
        1D azimuthal detector accurately recovers arm angles and widths even
        when one arm is thicker (like E-ELT's 90-deg crane arm), the spider
        center is shifted by 0.5 px, and edges are serrated; also works for
        a 4-arm spider at arbitrary angles.
        """
        n = 320
        yy, xx = np.mgrid[:n, :n]
        # Shift spider origin by +0.5 px relative to (n-1)/2
        cy, cx = n / 2.0, n / 2.0
        dy, dx = yy - cy, xx - cx
        r = np.hypot(dx, dy)
        theta = np.arctan2(dy, dx)
        r_max = n / 2.0 - 2.0 - 2.0 * np.abs(np.sin(18.0 * theta))
        r_min = 0.283 * (n / 2.0) + 1.5 * np.abs(np.sin(12.0 * theta))
        annulus = (r <= r_max) & (r >= r_min)

        expected_angles = [30.0, 90.0, 150.0, 210.0, 270.0, 330.0]
        expected_widths = [3.0, 7.0, 3.0, 3.0, 3.0, 3.0]
        spider = np.zeros_like(annulus)
        for ang_deg, w_px in zip(expected_angles, expected_widths):
            ang = np.deg2rad(ang_deg)
            s_par = np.sin(ang) * dy + np.cos(ang) * dx
            d_perp = np.abs(-np.cos(ang) * dy + np.sin(ang) * dx)
            spider |= (s_par > 0.0) & (d_perp <= 0.5 * w_px)

        mask = (annulus & (~spider)).astype(np.float32)
        arms, _, _ = detect_spider_arms_azimuthal(mask)
        self.assertEqual(len(arms), 6)
        for arm, exp_ang, exp_w in zip(arms, expected_angles, expected_widths):
            self.assertAlmostEqual(arm['azimuth_deg'], exp_ang, delta=0.35)
            self.assertAlmostEqual(arm['width_px'], exp_w, delta=0.45)

        # 4-arm spider at [22.5, 112.5, 202.5, 292.5] deg
        spider_4 = np.zeros_like(annulus)
        angles_4 = [22.5, 112.5, 202.5, 292.5]
        for ang_deg in angles_4:
            ang = np.deg2rad(ang_deg)
            s_par = np.sin(ang) * dy + np.cos(ang) * dx
            d_perp = np.abs(-np.cos(ang) * dy + np.sin(ang) * dx)
            spider_4 |= (s_par > 0.0) & (d_perp <= 1.8)
        arms_4, _, _ = detect_spider_arms_azimuthal((annulus & (~spider_4)).astype(np.float32))
        self.assertEqual(len(arms_4), 4)
        for arm, exp_ang in zip(arms_4, angles_4):
            self.assertAlmostEqual(arm['azimuth_deg'], exp_ang, delta=0.35)

        # Automatic radial plateau (radial_range=None) and custom radial_range=(0.2, 0.8)
        for r_range in (None, (0.2, 0.8)):
            arms_rr, _, _ = detect_spider_arms_azimuthal(mask, radial_range=r_range)
            self.assertEqual(len(arms_rr), 6)
            for arm, exp_ang in zip(arms_rr, expected_angles):
                self.assertAlmostEqual(arm['azimuth_deg'], exp_ang, delta=0.35)

        # Invalid radial_range raises ValueError
        with self.assertRaises(ValueError):
            detect_spider_arms_azimuthal(mask, radial_range=(0.8, 0.2))

        # Thin spider (1.2 px wide on 480x480 grid) is reliably detected with adaptive n_bins
        thin_mask, _ = _make_six_arm_pupil_mask(
            n_pix=480, obsratio=0.283, spider_width_px=1.2, serrated_edges=True
        )
        thin_arms, _, _ = detect_spider_arms_azimuthal(thin_mask)
        self.assertEqual(len(thin_arms), 6)
        for arm, exp_ang in zip(thin_arms, [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]):
            diff_deg = abs((arm['azimuth_deg'] - exp_ang + 180.0) % 360.0 - 180.0)
            self.assertLess(diff_deg, 0.5)

    def test_explicit_spider_angles_match_auto_detection(self):
        """
        Explicitly passing spider_angles_deg and spider_width_px matches
        automatic azimuthal detection for both specula_convention=True and False.
        """
        angles_6 = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
        for spec_conv in (True, False):
            for rot_deg in (0.0, 14.2, -6.2):
                auto_arms = compute_spider_subaperture_weights_and_normals(
                    pup_mask=self.pup_mask,
                    wfs_nsubaps=self.n_subaps,
                    wfs_rotation=rot_deg,
                    dilation_px=1.0,
                    smooth_sigma_px=2.0,
                    specula_convention=spec_conv,
                )
                explicit_arms = compute_spider_subaperture_weights_and_normals(
                    pup_mask=self.pup_mask,
                    wfs_nsubaps=self.n_subaps,
                    wfs_rotation=rot_deg,
                    dilation_px=1.0,
                    smooth_sigma_px=2.0,
                    spider_angles_deg=angles_6,
                    spider_width_px=3.2,
                    specula_convention=spec_conv,
                )
                self.assertEqual(len(auto_arms), 6)
                self.assertEqual(len(explicit_arms), 6)
                for auto_a, exp_a in zip(auto_arms, explicit_arms):
                    # Arm order and azimuth_deg are consistent in both conventions
                    self.assertEqual(auto_a['arm_index'], exp_a['arm_index'])
                    diff_az = abs((auto_a['azimuth_deg'] - exp_a['azimuth_deg'] + 180.0) % 360.0 - 180.0)
                    self.assertLess(diff_az, 0.5)
                    dot_n = abs(float(np.dot(exp_a['normal'], auto_a['normal'])))
                    self.assertGreater(dot_n, 0.9995)
                    np.testing.assert_allclose(
                        exp_a['weights'], auto_a['weights'], rtol=0, atol=0.06
                    )

        # Sequence spider_width_px (one width per arm) and length mismatch check
        arms_seq = compute_spider_subaperture_weights_and_normals(
            pup_mask=self.pup_mask,
            wfs_nsubaps=self.n_subaps,
            spider_angles_deg=angles_6,
            spider_width_px=[3.0, 7.0, 3.0, 3.0, 3.0, 3.0],
        )
        self.assertEqual([a['width_px'] for a in arms_seq], [3.0, 7.0, 3.0, 3.0, 3.0, 3.0])
        with self.assertRaises(ValueError):
            compute_spider_subaperture_weights_and_normals(
                pup_mask=self.pup_mask,
                wfs_nsubaps=self.n_subaps,
                spider_angles_deg=angles_6,
                spider_width_px=[3.0, 7.0],
            )

    def test_zero_dilation_and_sigma_touches_only_spider_subapertures(self):
        """
        With dilation_px=0.0 and smooth_sigma_px=0.0, non-zero weights must
        appear strictly on subapertures that physically overlap the spider arms,
        while adjacent non-overlapping subapertures stay at 0.
        """
        from synim.geometry import _transform_mask

        spider_only = (self.annulus & (self.pup_mask <= 0.5)).astype(np.float64)

        for rot_deg in (0.0, 6.2, -6.2):
            arms_strict = compute_spider_subaperture_weights_and_normals(
                pup_mask=self.pup_mask,
                wfs_nsubaps=self.n_subaps,
                wfs_rotation=rot_deg,
                dilation_px=0.0,
                smooth_sigma_px=0.0,
                specula_convention=True,
            )
            arms_wide = compute_spider_subaperture_weights_and_normals(
                pup_mask=self.pup_mask,
                wfs_nsubaps=self.n_subaps,
                wfs_rotation=rot_deg,
                dilation_px=2.0,
                smooth_sigma_px=3.0,
                specula_convention=True,
            )

            transform_kwargs = dict(
                dm_translation=(0.0, 0.0),
                dm_rotation=0.0,
                dm_magnification=(1.0, 1.0),
                wfs_translation=(0.0, 0.0),
                wfs_rotation=rot_deg,
                wfs_magnification=(1.0, 1.0),
                wfs_anamorphosis_45=1.0,
                output_size=self.pup_mask.shape,
            )
            # Combined transformed spider mask in the WFS frame (transposed for specula_convention)
            trans_spider = _transform_mask(np.transpose(spider_only), transform_kwargs, threshold=0.25) > 0
            has_spider_sa = synim.cpuArray(
                synim.utils.rebin(trans_spider.astype(np.float64),
                                  (self.n_subaps, self.n_subaps), method='sum')
            ) > 0.0

            strict_nonzero = np.zeros((self.n_subaps, self.n_subaps), dtype=bool)
            wide_nonzero = np.zeros((self.n_subaps, self.n_subaps), dtype=bool)
            for a_s, a_w in zip(arms_strict, arms_wide):
                strict_nonzero |= (a_s['weights_2d'] > 0.0)
                wide_nonzero |= (a_w['weights_2d'] > 0.0)

            # Every weighted subaperture in strict mode must actually overlap the spider
            self.assertFalse(np.any(strict_nonzero & (~has_spider_sa)))
            # And wide mode must be a strict superset of strict mode
            self.assertTrue(np.all(wide_nonzero[strict_nonzero]))
            self.assertGreater(np.count_nonzero(wide_nonzero), np.count_nonzero(strict_nonzero))

    def test_coarse_wfs_not_inflated_by_default(self):
        """With normalize_per_arm=False, 1x1 WFS weights stay small (not inflated to 1.0)."""
        arms_1x1 = compute_spider_subaperture_weights_and_normals(
            pup_mask=self.pup_mask,
            wfs_nsubaps=1,
            dilation_px=1.0,
            smooth_sigma_px=1.5,
        )
        self.assertEqual(len(arms_1x1), 6)
        for arm in arms_1x1:
            self.assertLess(float(np.max(arm['weights'])), 0.1)

    def test_normals_match_step_phase_slopes(self):
        """
        For piston steps across the 6 petals, the slopes computed by
        interaction_matrix on subapertures crossed by spider arm k must point
        along the normal vector returned for arm k.
        """
        rng = np.random.default_rng(7)
        idx_valid_sa = np.argwhere(rng.random((self.n_subaps, self.n_subaps)) > 0.25)

        cases = [
            dict(spec_conv=True, rot_deg=0.0, trans=(0.0, 0.0), mag=1.0, anam90=1.0, valid_sa=None),
            dict(spec_conv=True, rot_deg=6.2, trans=(0.0, 0.0), mag=1.0, anam90=1.0, valid_sa=idx_valid_sa),
            dict(spec_conv=True, rot_deg=14.2, trans=(2.5, -1.8), mag=1.04, anam90=1.03, valid_sa=idx_valid_sa),
            dict(spec_conv=True, rot_deg=-6.2, trans=(-1.5, 2.0), mag=0.97, anam90=0.98, valid_sa=None),
            dict(spec_conv=False, rot_deg=6.2, trans=(2.0, -1.5), mag=1.03, anam90=1.02, valid_sa=idx_valid_sa),
            dict(spec_conv=False, rot_deg=-6.2, trans=(0.0, 0.0), mag=1.0, anam90=1.0, valid_sa=None),
        ]

        for case in cases:
            with self.subTest(**{k: v if k != 'valid_sa' else (v is not None) for k, v in case.items()}):
                im = synim.cpuArray(
                    synim_core.interaction_matrix(
                        pup_diam_m=self.pup_diam_m,
                        pup_mask=self.pup_mask,
                        dm_array=self.petal_modes,
                        dm_mask=self.dm_mask,
                        dm_height=0.0,
                        dm_rotation=0.0,
                        wfs_nsubaps=self.n_subaps,
                        wfs_fov_arcsec=self.fov_arcsec,
                        gs_pol_coo=(0.0, 0.0),
                        gs_height=np.inf,
                        wfs_rotation=case['rot_deg'],
                        wfs_translation=case['trans'],
                        wfs_mag_global=case['mag'],
                        wfs_anamorphosis_90=case['anam90'],
                        idx_valid_sa=case['valid_sa'],
                        slope_method='derivatives',
                        specula_convention=case['spec_conv'],
                    )
                )
                arms_info = compute_spider_subaperture_weights_and_normals(
                    pup_mask=self.pup_mask,
                    wfs_nsubaps=self.n_subaps,
                    wfs_rotation=case['rot_deg'],
                    wfs_translation=case['trans'],
                    wfs_magnification=case['mag'],
                    wfs_anamorphosis_90=case['anam90'],
                    idx_valid_sa=case['valid_sa'],
                    dilation_px=0.0,
                    smooth_sigma_px=0.0,
                    specula_convention=case['spec_conv'],
                )
                self.assertEqual(len(arms_info), 6)
                n_sa = im.shape[0] // 2

                for arm in arms_info:
                    w = arm['weights']
                    self.assertEqual(w.shape, (n_sa,))
                    n_vec = arm['normal']
                    t_vec = np.array([-n_vec[1], n_vec[0]])
                    np.testing.assert_allclose(np.linalg.norm(n_vec), 1.0, atol=1e-12)

                    idx_on_arm = np.where(w > 0.5 * np.max(w))[0]
                    self.assertGreater(len(idx_on_arm), 0)

                    s0 = im[idx_on_arm, :]
                    s1 = im[idx_on_arm + n_sa, :]
                    proj_n = np.linalg.norm(s0 * n_vec[0] + s1 * n_vec[1])
                    proj_t = np.linalg.norm(s0 * t_vec[0] + s1 * t_vec[1])
                    alignment = proj_n / np.hypot(proj_n, proj_t)
                    self.assertGreater(alignment, 0.985)

                # Opposite arms must be collinear (3 pairs of normals)
                normals = np.array([a['normal'] for a in arms_info])
                dots = np.abs(normals @ normals.T)
                for k in range(6):
                    self.assertEqual(np.count_nonzero(dots[k] > 0.995), 2)

    def test_sherman_morrison_matches_explicit_2x2_inverse(self):
        """
        update_wfs_noise_cov_inv_for_spider must match explicit 2x2 block
        inversion of (C_j + sum_k sigma_{j,k}^2 n_k n_k^T)^(-1), preserve
        bit-exact symmetry, and pass _sparse_rows(symmetric=True).
        """
        n_sa = 120
        rng = np.random.default_rng(123)

        a = rng.uniform(8.0, 20.0, size=n_sa)
        d = rng.uniform(8.0, 20.0, size=n_sa)
        b = rng.uniform(-4.0, 4.0, size=n_sa)

        C_inv = np.zeros((2 * n_sa, 2 * n_sa), dtype=np.float32)
        idx0 = np.arange(n_sa)
        idx1 = idx0 + n_sa
        C_inv[idx0, idx0] = a.astype(np.float32)
        C_inv[idx1, idx1] = d.astype(np.float32)
        C_inv[idx0, idx1] = b.astype(np.float32)
        C_inv[idx1, idx0] = b.astype(np.float32)

        arms_info = [
            {
                'arm_index': 1,
                'normal': np.array([0.5, -np.sqrt(0.75)]),
                'weights': rng.uniform(0.0, 1.0, size=n_sa),
            },
            {
                'arm_index': 2,
                'normal': np.array([0.8, 0.6]),
                'weights': rng.choice([0.0, 0.5, 1.0], size=n_sa),
            },
        ]
        alpha = 15.0
        sigma_base = 0.08

        C_inv_mod = update_wfs_noise_cov_inv_for_spider(
            C_inv,
            arms_info=arms_info,
            alpha=alpha,
            sigma_noise2_base=sigma_base,
            weight_threshold=1e-3,
            in_place=False,
        )

        self.assertTrue(np.array_equal(C_inv_mod, C_inv_mod.T))
        self.assertLessEqual(int(np.count_nonzero(C_inv_mod, axis=1).max()), 2)
        sparse_res = _sparse_rows(C_inv_mod, np, dtype=np.float32)
        self.assertIsNotNone(sparse_res)
        self.assertTrue(sparse_res[2])  # symmetric=True

        for j in range(n_sa):
            M_j = np.array(
                [[C_inv[j, j], C_inv[j, j + n_sa]],
                 [C_inv[j + n_sa, j], C_inv[j + n_sa, j + n_sa]]],
                dtype=np.float64,
            )
            C_j = np.linalg.inv(M_j)
            for arm in arms_info:
                w_j = arm['weights'][j]
                if w_j > 1e-3:
                    n_k = arm['normal'][:, None]
                    C_j = C_j + (alpha * w_j * sigma_base) * (n_k @ n_k.T)
            M_j_expected = np.linalg.inv(C_j).astype(np.float32)
            M_j_actual = np.array(
                [[C_inv_mod[j, j], C_inv_mod[j, j + n_sa]],
                 [C_inv_mod[j + n_sa, j], C_inv_mod[j + n_sa, j + n_sa]]],
                dtype=np.float32,
            )
            np.testing.assert_allclose(M_j_actual, M_j_expected, rtol=1e-5, atol=1e-6)

    def test_edge_cases_and_block_diagonal_check(self):
        n_sa = 4
        C_inv = np.eye(2 * n_sa, dtype=np.float32) * 5.0
        arms = [{'arm_index': 1, 'normal': np.array([1.0, 0.0]), 'weights': np.ones(n_sa)}]

        # alpha <= 0 returns untouched copy
        C_zero = update_wfs_noise_cov_inv_for_spider(C_inv, arms, alpha=0.0, in_place=False)
        np.testing.assert_array_equal(C_zero, C_inv)
        self.assertIsNot(C_zero, C_inv)

        # Wrong weight length raises ValueError
        with self.assertRaises(ValueError):
            update_wfs_noise_cov_inv_for_spider(
                C_inv, [{'arm_index': 1, 'normal': np.array([1.0, 0.0]), 'weights': np.ones(n_sa + 1)}]
            )

        # Off-block-diagonal entry raises ValueError when check_block_diagonal=True
        C_dense = C_inv.copy()
        C_dense[0, 1] = 0.25
        C_dense[1, 0] = 0.25
        with self.assertRaises(ValueError):
            update_wfs_noise_cov_inv_for_spider(C_dense, arms, alpha=10.0)

    def test_directional_selectivity_at_high_alpha(self):
        """
        For large alpha on a single arm, the precision along the arm normal n_k
        is suppressed while the precision along the tangent t_k is preserved
        for an isotropic 2x2 block.
        """
        n_sa = 2
        C_inv = np.eye(2 * n_sa, dtype=np.float64) * 10.0
        n_vec = np.array([np.cos(0.7), np.sin(0.7)])
        t_vec = np.array([-np.sin(0.7), np.cos(0.7)])
        arms_info = [
            {
                'arm_index': 1,
                'normal': n_vec,
                'weights': np.array([1.0, 0.0]),
            }
        ]
        C_mod = update_wfs_noise_cov_inv_for_spider(C_inv, arms_info, alpha=1e6)

        M0 = np.array([[C_mod[0, 0], C_mod[0, 2]], [C_mod[2, 0], C_mod[2, 2]]])
        self.assertLess(float(n_vec @ M0 @ n_vec), 1e-4)
        np.testing.assert_allclose(float(t_vec @ M0 @ t_vec), 10.0, rtol=1e-10)

        M1 = np.array([[C_mod[1, 1], C_mod[1, 3]], [C_mod[3, 1], C_mod[3, 3]]])
        np.testing.assert_allclose(M1, 10.0 * np.eye(2), rtol=0, atol=0)

    def test_end_to_end_apply_spider_and_mmse_petal_suppression(self):
        """
        End-to-end test with a multi-WFS mock ParamsManager, _sparse_rows, and
        compute_mmse_reconstructor: increasing alpha monotonically suppresses
        reconstructor response to petal piston steps while preserving smooth tilt.
        """
        yy_sa, xx_sa = np.mgrid[:self.n_subaps, :self.n_subaps]
        r_sa = np.hypot(xx_sa - (self.n_subaps - 1) / 2.0, yy_sa - (self.n_subaps - 1) / 2.0)
        idx_valid = np.argwhere((r_sa <= self.n_subaps / 2.0 - 1.0) & (r_sa >= 0.283 * self.n_subaps / 2.0))

        wfs_configs = [
            dict(
                wfs_nsubaps=self.n_subaps,
                wfs_rotation=rot,
                wfs_translation=(0.0, 0.0),
                wfs_magnification=(1.0, 1.0),
                wfs_anamorphosis_90=1.0,
                wfs_anamorphosis_45=1.0,
                idx_valid_sa=idx_valid,
            )
            for rot in (6.2, -6.2)
        ]
        pm_mock = _MockParamsManager(self.serrated_mask, wfs_configs)
        n_slopes_wfs = len(idx_valid) * 2
        n_slopes_total = 2 * n_slopes_wfs

        # Build block-diagonal 2x2 C_noise_inv
        C_base = np.eye(n_slopes_total, dtype=np.float32) * 100.0
        C_orig = C_base.copy()

        # Check active_wfs_mask validation
        with self.assertRaises(ValueError):
            apply_spider_to_noise_covariance(pm_mock, C_base, active_wfs_mask=[True])

        # Build IM for 6 petal step modes and 2 smooth global tilts separately
        yy, xx = np.mgrid[:self.n_pix, :self.n_pix].astype(np.float32)
        labeled_petals, _ = label(self.serrated_mask > 0.5)
        petals = np.stack([(labeled_petals == p).astype(np.float32) for p in range(1, 7)], axis=-1)
        tilts = np.stack([xx - self.n_pix / 2.0, yy - self.n_pix / 2.0], axis=-1)
        dm_mask = self.serrated_annulus.astype(np.float32)

        def _assemble_im(mode_cube):
            return np.concatenate([
                synim.cpuArray(
                    synim_core.interaction_matrix(
                        pup_diam_m=self.pup_diam_m,
                        pup_mask=self.serrated_mask,
                        dm_array=mode_cube,
                        dm_mask=dm_mask,
                        dm_height=0.0,
                        dm_rotation=0.0,
                        wfs_nsubaps=self.n_subaps,
                        wfs_fov_arcsec=self.fov_arcsec,
                        gs_pol_coo=(0.0, 0.0),
                        gs_height=np.inf,
                        wfs_rotation=cfg['wfs_rotation'],
                        wfs_translation=cfg['wfs_translation'],
                        wfs_mag_global=1.0,
                        idx_valid_sa=idx_valid,
                        specula_convention=True,
                    )
                )
                for cfg in wfs_configs
            ], axis=0)

        im_petals = _assemble_im(petals)
        im_tilts = _assemble_im(tilts)
        C_atm_petals_inv = np.eye(petals.shape[-1], dtype=np.float32) * 1e-2
        C_atm_tilts_inv = np.eye(tilts.shape[-1], dtype=np.float32) * 1e-4

        petal_norms = []
        tilt_norms = []
        for alpha in (0.0, 10.0, 100.0):
            C_mod, diag = apply_spider_to_noise_covariance(
                pm_mock,
                C_noise_inv=C_base,
                wfs_type='lgs',
                alpha=alpha,
                dilation_px=1.5,
                smooth_sigma_px=2.0,
                in_place=False,
                return_diagnostics=True,
            )
            # Input untouched when in_place=False
            self.assertTrue(np.array_equal(C_base, C_orig))
            self.assertTrue(np.array_equal(C_mod, C_mod.T))
            self.assertEqual(len(diag['wfs']), 2)
            self.assertEqual(len(diag['wfs'][0]['arms_info']), 6)

            sparse_res = _sparse_rows(C_mod, np, dtype=np.float32)
            self.assertIsNotNone(sparse_res)
            self.assertTrue(sparse_res[2])

            R_petals = compute_mmse_reconstructor(
                im_petals, C_atm_petals_inv, C_noise=C_mod, cinverse=True, xp=np, dtype=np.float32
            )
            R_tilts = compute_mmse_reconstructor(
                im_tilts, C_atm_tilts_inv, C_noise=C_mod, cinverse=True, xp=np, dtype=np.float32
            )
            petal_norms.append(float(np.linalg.norm(R_petals @ im_petals)))
            tilt_norms.append(float(np.linalg.norm(R_tilts @ im_tilts)))

        # Monotonic suppression of petal-step response as alpha increases
        self.assertGreater(petal_norms[0], petal_norms[1])
        self.assertGreater(petal_norms[1], petal_norms[2])
        # Smooth tilt modes remain well reconstructed
        self.assertGreater(tilt_norms[2] / tilt_norms[0], 0.99)

        # Display option saves diagnostic plot for all WFSs and for 'first' only
        import os
        import tempfile
        import matplotlib
        matplotlib.use('Agg')
        with tempfile.TemporaryDirectory() as tmpdir:
            for disp_mode, fname in ((True, 'all.png'), ('first', 'first.png')):
                out_png = os.path.join(tmpdir, fname)
                apply_spider_to_noise_covariance(
                    pm_mock,
                    C_noise_inv=C_base,
                    wfs_type='lgs',
                    alpha=10.0,
                    display=disp_mode,
                    display_save_path=out_png,
                )
                self.assertTrue(os.path.exists(out_png))
                self.assertGreater(os.path.getsize(out_png), 0)


if __name__ == '__main__':
    unittest.main()
