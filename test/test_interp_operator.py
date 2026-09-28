"""
Equivalence tests for the sparse bilinear operator used by rotshiftzoom_array
(interp='sparse') against the original affine_transform implementation.

The reference is a frozen copy of rotshiftzoom_array (test/reference_rotshiftzoom.py),
so these tests keep checking the original behaviour even if the library code
changes. Tests run on the active synim backend (CPU by default).
"""
import unittest
import numpy as np

import synim
import synim.synim as synim_core
import synim.synpm as synpm
from synim.utils import (
    rotshiftzoom_array,
    build_bilinear_operator,
    set_interp_method,
    get_interp_method,
)
from synim import cpuArray

try:
    from test.reference_rotshiftzoom import rotshiftzoom_array_reference
except ImportError:  # running from inside the test directory
    from reference_rotshiftzoom import rotshiftzoom_array_reference


# Geometries covering the transformations used by SynIM, including cases
# where input positions fall exactly on the grid border (identity, integer
# shifts, 90/180 degree rotations) and cases with parts of the output
# mapping outside the input (de-magnification, large shifts).
GEOMETRIES = {
    'identity': {},
    'integer_shift': dict(dm_translation=(3.0, -2.0)),
    'fractional_shift': dict(dm_translation=(1.37, -0.42)),
    'large_shift': dict(dm_translation=(15.5, 9.25)),
    'rotation_small': dict(dm_rotation=7.3),
    'rotation_45': dict(dm_rotation=45.0),
    'rotation_90': dict(dm_rotation=90.0),
    'rotation_180': dict(dm_rotation=180.0),
    'magnification': dict(dm_magnification=(1.1, 0.9)),
    'demagnification': dict(dm_magnification=(0.8, 0.8)),
    'wfs_all': dict(wfs_translation=(0.3, 1.2), wfs_rotation=20.0,
                    wfs_magnification=(1.02, 0.97)),
    'anamorphosis_45': dict(wfs_anamorphosis_45=1.05),
    'dm_and_wfs': dict(dm_translation=(-2.1, 0.7), dm_rotation=3.0,
                       dm_magnification=(1.03, 1.03), wfs_translation=(0.5, -0.25),
                       wfs_rotation=-31.0, wfs_magnification=(0.98, 1.01),
                       wfs_anamorphosis_45=0.97),
    'larger_output': dict(dm_rotation=12.0, output_size=(53, 53)),
    'smaller_output': dict(dm_rotation=-5.0, output_size=(33, 33)),
    'non_square_output': dict(dm_rotation=17.0, dm_translation=(0.4, 0.1),
                              output_size=(37, 49)),
}

# Tolerances relative to max(|input|). float32: the weights and the products
# are rounded to float32 (affine_transform accumulates in float64).
RTOL = {np.float32: 2e-6, np.float64: 1e-12}


def _random_cube(shape, dtype, seed=0):
    rng = np.random.default_rng(seed)
    return synim.xp.asarray(rng.standard_normal(shape).astype(dtype))


def _max_rel_diff(a, b, scale):
    return float(np.max(np.abs(cpuArray(a) - cpuArray(b)))) / scale


class TestSparseOperatorEquivalence(unittest.TestCase):
    """rotshiftzoom_array(interp='sparse') vs the frozen affine reference."""

    def tearDown(self):
        set_interp_method('affine')

    def _check(self, cube, geometry, dtype):
        ref = rotshiftzoom_array_reference(cube, **geometry)
        out = rotshiftzoom_array(cube, interp='sparse', **geometry)
        self.assertEqual(out.shape, ref.shape)
        self.assertEqual(out.dtype, ref.dtype)
        scale = float(np.max(np.abs(cpuArray(cube))))
        # Whole array, border included: a pixel classified differently at the
        # border would give an O(1) difference and fail this check.
        self.assertLessEqual(_max_rel_diff(out, ref, scale), RTOL[dtype])

    def test_geometries_square_input(self):
        for dtype in (np.float32, np.float64):
            cube = _random_cube((41, 41, 5), dtype)
            for name, geometry in GEOMETRIES.items():
                with self.subTest(geometry=name, dtype=dtype.__name__):
                    self._check(cube, geometry, dtype)

    def test_geometries_non_square_input(self):
        for dtype in (np.float32, np.float64):
            cube = _random_cube((40, 52, 4), dtype, seed=1)
            for name, geometry in GEOMETRIES.items():
                with self.subTest(geometry=name, dtype=dtype.__name__):
                    self._check(cube, geometry, dtype)

    def test_even_size_input(self):
        cube = _random_cube((64, 64, 3), np.float32, seed=2)
        for name, geometry in GEOMETRIES.items():
            with self.subTest(geometry=name):
                self._check(cube, geometry, np.float32)

    def test_single_slice(self):
        cube = _random_cube((41, 41, 1), np.float32, seed=3)
        self._check(cube, GEOMETRIES['dm_and_wfs'], np.float32)

    def test_transposed_input_view(self):
        # SPECULA convention: the pipeline passes a transposed view (1, 0, 2)
        base = _random_cube((41, 45, 4), np.float64, seed=4)
        cube = synim.xp.transpose(base, (1, 0, 2))
        self.assertFalse(cube.flags.c_contiguous)
        for name in ('identity', 'rotation_90', 'dm_and_wfs', 'non_square_output'):
            with self.subTest(geometry=name):
                self._check(cube, GEOMETRIES[name], np.float64)

    def test_generic_non_contiguous_input(self):
        base = _random_cube((41, 41, 8), np.float64, seed=5)
        cube = base[:, :, ::2]   # non-contiguous, not a transposed view
        self.assertFalse(cube.flags.c_contiguous)
        self._check(cube, GEOMETRIES['dm_and_wfs'], np.float64)

    def test_nan_input(self):
        cube = _random_cube((41, 41, 3), np.float64, seed=6)
        cube[5:9, 10:12, 1] = np.nan
        ref = rotshiftzoom_array_reference(cube, **GEOMETRIES['dm_and_wfs'])
        out = rotshiftzoom_array(cube, interp='sparse', **GEOMETRIES['dm_and_wfs'])
        self.assertFalse(np.isnan(cpuArray(out)).any())
        self.assertLessEqual(_max_rel_diff(out, ref, float(np.nanmax(np.abs(cube)))), 1e-12)

    def test_input_not_modified(self):
        cube = _random_cube((41, 41, 3), np.float32, seed=7)
        copy = cube.copy()
        rotshiftzoom_array(cube, interp='sparse', **GEOMETRIES['dm_and_wfs'])
        self.assertTrue(np.array_equal(cpuArray(cube), cpuArray(copy)))


class TestDefaultBehaviourUnchanged(unittest.TestCase):
    """The default path must stay bit-identical to the original implementation."""

    def tearDown(self):
        set_interp_method('affine')

    def test_default_is_affine(self):
        self.assertEqual(get_interp_method(), 'affine')

    def test_default_bit_identical_3d(self):
        cube = _random_cube((41, 41, 4), np.float32)
        for name, geometry in GEOMETRIES.items():
            with self.subTest(geometry=name):
                ref = rotshiftzoom_array_reference(cube, **geometry)
                out = rotshiftzoom_array(cube, **geometry)
                self.assertTrue(np.array_equal(cpuArray(out), cpuArray(ref)))

    def test_2d_always_affine(self):
        # 2D arrays (masks) keep using affine_transform, whatever the method
        img = _random_cube((41, 41), np.float32)
        set_interp_method('sparse')
        for name, geometry in GEOMETRIES.items():
            with self.subTest(geometry=name):
                ref = rotshiftzoom_array_reference(img, **geometry)
                out = rotshiftzoom_array(img, **geometry)
                self.assertTrue(np.array_equal(cpuArray(out), cpuArray(ref)))

    def test_set_interp_method(self):
        set_interp_method('sparse')
        self.assertEqual(get_interp_method(), 'sparse')
        cube = _random_cube((41, 41, 3), np.float64)
        ref = rotshiftzoom_array_reference(cube, **GEOMETRIES['dm_and_wfs'])
        out = rotshiftzoom_array(cube, **GEOMETRIES['dm_and_wfs'])  # uses default
        self.assertLessEqual(_max_rel_diff(out, ref, float(np.max(np.abs(cube)))), 1e-12)

    def test_invalid_method(self):
        with self.assertRaises(ValueError):
            set_interp_method('cubic')
        with self.assertRaises(ValueError):
            rotshiftzoom_array(_random_cube((8, 8, 2), np.float32), interp='cubic')


class TestOperatorStructure(unittest.TestCase):
    """Properties of the sparse operator itself."""

    def test_rows_and_weights(self):
        n = 41
        theta = np.deg2rad(13.0)
        matrix = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta), np.cos(theta)]]) / 1.07
        offset = np.array([n / 2, n / 2]) - matrix @ np.array([n / 2, n / 2]) + [0.3, -0.8]
        op = build_bilinear_operator((n, n), (n, n), matrix, offset, dtype=np.float64)
        self.assertEqual(op.shape, (n * n, n * n))

        nnz_per_row = np.diff(op.indptr)
        self.assertLessEqual(int(nnz_per_row.max()), 4)
        self.assertTrue((op.data > 0).all())

        # Rows of output pixels mapping inside the input sum to 1 (partition
        # of unity); rows mapping outside are empty (constant mode, cval=0)
        ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
        y = (offset[0] + ii * matrix[0, 0]) + jj * matrix[0, 1]
        x = (offset[1] + ii * matrix[1, 0]) + jj * matrix[1, 1]
        inside = ((y >= 0) & (y <= n - 1) & (x >= 0) & (x <= n - 1)).ravel()
        row_sums = np.asarray(op.sum(axis=1)).ravel()
        np.testing.assert_allclose(row_sums[inside], 1.0, rtol=0, atol=1e-12)
        self.assertTrue((nnz_per_row[~inside] == 0).all())

    def test_identity_operator(self):
        n = 17
        op = build_bilinear_operator((n, n), (n, n), np.eye(2), np.zeros(2),
                                     dtype=np.float64)
        self.assertEqual(op.nnz, n * n)
        np.testing.assert_array_equal(op.toarray(), np.eye(n * n))

    def test_border_rounding_matches_scipy(self):
        # Input positions that fall on the grid border up to rounding: whether
        # the pixel is interpolated or set to cval depends on the exact
        # floating point operation order. build_bilinear_operator must use
        # the same order as scipy's affine_transform. The offset is chosen so
        # that y is exactly 0 with a different summation order and a tiny
        # negative number with scipy's order (or vice versa).
        from scipy.ndimage import affine_transform as scipy_affine_transform
        n = 41
        rng = np.random.default_rng(3)
        img = rng.standard_normal((n, n))
        n_cases = 0
        for _ in range(2000):
            theta = rng.uniform(0, 2 * np.pi)
            matrix = np.array([[np.cos(theta), -np.sin(theta)],
                               [np.sin(theta), np.cos(theta)]]) * rng.uniform(0.8, 1.2)
            i, j = rng.integers(0, n, 2)
            offset = np.array([-(i * matrix[0, 0] + j * matrix[0, 1]),
                               rng.uniform(5, 30)])
            y_scipy_order = (offset[0] + i * matrix[0, 0]) + j * matrix[0, 1]
            if y_scipy_order >= 0:
                continue  # rounding does not matter for this case
            n_cases += 1
            ref = scipy_affine_transform(img, matrix, offset=offset, order=1)
            op = build_bilinear_operator((n, n), (n, n), matrix, offset,
                                         dtype=np.float64)
            out = (op @ img.ravel()).reshape(n, n)
            np.testing.assert_allclose(out, ref, rtol=0, atol=1e-12)
            if n_cases >= 50:
                break
        self.assertGreaterEqual(n_cases, 50)

    def test_transposed_columns(self):
        n_y, n_x = 9, 13
        matrix = np.array([[0.9, 0.1], [-0.1, 0.9]])
        offset = np.array([0.7, 1.1])
        op = build_bilinear_operator((n_y, n_x), (n_y, n_x), matrix, offset,
                                     dtype=np.float64)
        op_t = build_bilinear_operator((n_y, n_x), (n_y, n_x), matrix, offset,
                                       dtype=np.float64, input_transposed=True)
        img = np.random.default_rng(0).standard_normal((n_y, n_x))
        np.testing.assert_allclose(op @ img.ravel(), op_t @ img.T.ravel(),
                                   rtol=0, atol=1e-14)


class TestEndToEndEquivalence(unittest.TestCase):
    """Interaction and projection matrices computed with both methods."""

    @classmethod
    def setUpClass(cls):
        n = 96
        cls.n = n
        cls.nsubaps = 12
        yy, xx = np.mgrid[:n, :n]
        r = np.hypot(xx - n / 2 + 0.5, yy - n / 2 + 0.5)
        cls.pup_mask = ((r < n / 2) & (r > n / 8)).astype(np.float32)
        cls.dm_mask = (r < n / 2 + 4).astype(np.float32)
        rng = np.random.default_rng(0)
        x = (xx - n / 2) / (n / 2)
        y = (yy - n / 2) / (n / 2)
        modes = [rng.normal() * x ** (k % 4) * y ** (k // 4 % 4)
                 + rng.normal() * np.sin((k + 1) * x) * np.cos((k % 5) * y)
                 for k in range(30)]
        cls.dm_array = np.stack(modes, axis=2).astype(np.float32)

    def tearDown(self):
        set_interp_method('affine')

    def _im(self, method, **kw):
        set_interp_method(method)
        return cpuArray(synim_core.interaction_matrix(
            8.0, self.pup_mask, self.dm_array, self.dm_mask,
            wfs_nsubaps=self.nsubaps, wfs_fov_arcsec=4.0, **kw))

    def _assert_close(self, a, b, rtol=2e-5):
        self.assertEqual(a.shape, b.shape)
        scale = float(np.max(np.abs(a)))
        self.assertGreater(scale, 0)
        self.assertLessEqual(float(np.max(np.abs(a - b))) / scale, rtol)

    def test_interaction_matrix(self):
        cases = {
            'no_transform': dict(dm_height=0.0, dm_rotation=0.0, gs_pol_coo=(0, 0),
                                 gs_height=np.inf, wfs_rotation=0.0,
                                 wfs_translation=(0, 0), wfs_mag_global=1.0),
            'dm_only': dict(dm_height=5000.0, dm_rotation=2.0, gs_pol_coo=(30, 45),
                            gs_height=np.inf, wfs_rotation=0.0,
                            wfs_translation=(0, 0), wfs_mag_global=1.0),
            'wfs_rotation': dict(dm_height=0.0, dm_rotation=0.0, gs_pol_coo=(0, 0),
                                 gs_height=np.inf, wfs_rotation=4.0,
                                 wfs_translation=(0, 0), wfs_mag_global=1.0),
            'dm_and_wfs': dict(dm_height=8000.0, dm_rotation=1.0, gs_pol_coo=(20, 10),
                               gs_height=90e3, wfs_rotation=3.0,
                               wfs_translation=(0.5, 0.2), wfs_mag_global=0.98,
                               wfs_anamorphosis_45=1.02),
        }
        for name, kw in cases.items():
            with self.subTest(case=name, slope_method='derivatives'):
                self._assert_close(self._im('affine', **kw), self._im('sparse', **kw))
        # telsum: only cases that go through the combined workflow
        for name in ('wfs_rotation', 'dm_and_wfs'):
            with self.subTest(case=name, slope_method='telsum'):
                kw = dict(cases[name], slope_method='telsum')
                self._assert_close(self._im('affine', **kw), self._im('sparse', **kw))

    def test_interaction_matrices_multi_wfs(self):
        configs = [
            dict(name='a', nsubaps=self.nsubaps, fov_arcsec=4.0,
                 gs_pol_coo=(30, 0), gs_height=90e3),
            dict(name='b', nsubaps=self.nsubaps, fov_arcsec=4.0,
                 gs_pol_coo=(30, 120), gs_height=90e3, rotation=10.0,
                 translation=(0.3, -0.1)),
        ]
        results = {}
        for method in ('affine', 'sparse'):
            set_interp_method(method)
            im_dict, _ = synim_core.interaction_matrices_multi_wfs(
                8.0, self.pup_mask, self.dm_array, self.dm_mask, 8000.0, 1.0, configs)
            results[method] = {k: cpuArray(v) for k, v in im_dict.items()}
        for key in results['affine']:
            with self.subTest(wfs=key):
                self._assert_close(results['affine'][key], results['sparse'][key])

    def test_projection_matrix(self):
        idx = np.where(self.pup_mask > 0.5)
        base = self.dm_array[idx[0], idx[1], :10].astype(np.float64)
        base_inv = np.linalg.pinv(base).astype(np.float32)
        results = {}
        for method in ('affine', 'sparse'):
            set_interp_method(method)
            results[method] = cpuArray(synpm.projection_matrix(
                8.0, self.pup_mask, self.dm_array, self.dm_mask, base_inv,
                5000.0, 1.0, 0.0, (0, 0), (1, 1), (20, 10), np.inf))
        self._assert_close(results['affine'], results['sparse'])


if __name__ == '__main__':
    unittest.main()
