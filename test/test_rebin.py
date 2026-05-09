import unittest
import numpy as np
from synim.utils import rebin

class TestRebin(unittest.TestCase):

    def test_rebin_1d_to_2d_conversion(self):
        """Test that 1D arrays are converted to 2D before rebinning"""
        array_1d = np.array([1, 2, 3, 4])
        # After conversion to 2D: (4, 1), expanding to (8, 2)
        result = rebin(array_1d, (8, 2), method='average')
        self.assertEqual(result.shape, (8, 2))

    def test_rebin_2d_expansion(self):
        """Test expansion of 2D array"""
        array = np.array([[1, 2], [3, 4]])
        result = rebin(array, (4, 4), method='average')
        # np.tile replicates: [[1,2],[3,4]] -> [[1,2,1,2],[3,4,3,4],[1,2,1,2],[3,4,3,4]]
        expected = np.array([[1, 2, 1, 2],
                            [3, 4, 3, 4],
                            [1, 2, 1, 2],
                            [3, 4, 3, 4]])
        np.testing.assert_array_equal(result, expected)

    def test_rebin_2d_expansion_non_integer_factor(self):
        """Test that expansion with non-integer factors raises error"""
        array = np.array([[1, 2], [3, 4]])
        with self.assertRaises(ValueError):
            rebin(array, (5, 5), method='average')

    def test_rebin_3d_expansion(self):
        """Test expansion of 3D array preserves third dimension"""
        array = np.ones((2, 2, 3))
        result = rebin(array, (4, 4), method='average')
        self.assertEqual(result.shape, (4, 4, 3))
        np.testing.assert_array_equal(result[:, :, 0], np.ones((4, 4)))

    def test_rebin_2d_compression_sum(self):
        """Test compression with sum method"""
        array = np.array([[1, 2, 3, 4],
                         [5, 6, 7, 8],
                         [9, 10, 11, 12],
                         [13, 14, 15, 16]])
        result = rebin(array, (2, 2), method='sum')
        expected = np.array([[14, 22], [46, 54]])
        np.testing.assert_array_equal(result, expected)

    def test_rebin_2d_compression_average(self):
        """Test compression with average method"""
        array = np.array([[1, 2, 3, 4],
                         [5, 6, 7, 8],
                         [9, 10, 11, 12],
                         [13, 14, 15, 16]], dtype=float)
        result = rebin(array, (2, 2), method='average')
        expected = np.array([[3.5, 5.5], [11.5, 13.5]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_rebin_2d_compression_nanmean(self):
        """Test compression with nanmean method handling NaN values"""
        array = np.array([[1, 2, 3, 4],
                         [5, np.nan, 7, 8],
                         [9, 10, 11, 12],
                         [13, 14, 15, 16]], dtype=float)
        result = rebin(array, (2, 2), method='nanmean')
        # Top-left: mean of [1, 2, 5, nan] = (1+2+5)/3 = 8/3
        self.assertAlmostEqual(result[0, 0], 8/3, places=5)
        # Other quadrants should work normally
        self.assertAlmostEqual(result[0, 1], 5.5, places=5)

    def test_rebin_3d_compression_average(self):
        """Test compression of 3D array"""
        array = np.ones((4, 4, 2), dtype=float)
        array[:, :, 0] *= 2
        array[:, :, 1] *= 3
        result = rebin(array, (2, 2), method='average')
        self.assertEqual(result.shape, (2, 2, 2))
        np.testing.assert_array_almost_equal(result[:, :, 0], np.ones((2, 2)) * 2)
        np.testing.assert_array_almost_equal(result[:, :, 1], np.ones((2, 2)) * 3)

    def test_rebin_3d_compression_sum(self):
        """Test 3D compression with sum method"""
        array = np.ones((4, 4, 2), dtype=float)
        result = rebin(array, (2, 2), method='sum')
        # Each 2x2 block sums to 4
        np.testing.assert_array_almost_equal(result[:, :, 0], np.ones((2, 2)) * 4)
        np.testing.assert_array_almost_equal(result[:, :, 1], np.ones((2, 2)) * 4)

    def test_rebin_unsupported_method(self):
        """Test that unsupported method raises ValueError"""
        array = np.array([[1, 2], [3, 4]])
        with self.assertRaises(ValueError):
            rebin(array, (1, 1), method='median')

    def test_rebin_compression_non_divisible(self):
        """Test compression when dimensions aren't perfectly divisible.
        Checks the new anti-aliasing upscaling behavior that prevents spatial shifts."""
        array = np.ones((5, 5))
        result = rebin(array, (2, 2), method='average')

        # Shape should be correct
        self.assertEqual(result.shape, (2, 2))

        # With the new affine_transform upscaling (mode='constant', cval=0.0),
        # an array of ones will drop off at the edges.
        # The center of mass should remain perfectly centered.
        expected = np.array([[1.0, 0.6666667],
                             [0.6666667, 0.44444445]], dtype=np.float32)
        
        np.testing.assert_array_almost_equal(result, expected, decimal=5)

    def test_rebin_preserve_dtype_float(self):
        """Test that float dtype is preserved"""
        array = np.array([[1.5, 2.5], [3.5, 4.5]], dtype=np.float32)
        result = rebin(array, (1, 1), method='average')
        self.assertEqual(result.dtype, np.float32)

    def test_rebin_identity_operation(self):
        """Test rebinning to same size returns similar array"""
        array = np.array([[1, 2, 3, 4],
                         [5, 6, 7, 8],
                         [9, 10, 11, 12],
                         [13, 14, 15, 16]], dtype=float)
        result = rebin(array, (4, 4), method='average')
        np.testing.assert_array_almost_equal(result, array)

    def test_rebin_expansion_3d_different_channels(self):
        """Test 3D expansion maintains channel independence"""
        array = np.zeros((2, 2, 3))
        array[:, :, 0] = 1
        array[:, :, 1] = 2
        array[:, :, 2] = 3
        result = rebin(array, (4, 4), method='average')
        np.testing.assert_array_equal(result[:, :, 0], np.ones((4, 4)))
        np.testing.assert_array_equal(result[:, :, 1], np.ones((4, 4)) * 2)
        np.testing.assert_array_equal(result[:, :, 2], np.ones((4, 4)) * 3)

    def test_rebin_nanmean_all_nan_slice(self):
        """Test nanmean method with all-NaN slices"""
        array = np.array([[np.nan, np.nan, 1, 2],
                         [np.nan, np.nan, 3, 4],
                         [5, 6, 7, 8],
                         [9, 10, 11, 12]], dtype=float)
        with np.testing.suppress_warnings() as sup:
            sup.filter(RuntimeWarning, "Mean of empty slice")
            result = rebin(array, (2, 2), method='nanmean')
        # Top-left quadrant is all NaN, should result in NaN
        self.assertTrue(np.isnan(result[0, 0]))

    def test_rebin_1d_compression(self):
        """Test 1D array compression after conversion to 2D"""
        array_1d = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        # After conversion to 2D: (8, 1), compressing to (4, 1)
        result = rebin(array_1d, (4, 1), method='average')
        self.assertEqual(result.shape, (4, 1))
        expected = np.array([[1.5], [3.5], [5.5], [7.5]])
        np.testing.assert_array_almost_equal(result, expected)

    def test_rebin_expansion_factor_2(self):
        """Test expansion by factor of 2 using tile behavior"""
        array = np.array([[1, 2], [3, 4]])
        result = rebin(array, (4, 4), method='average')
        # np.tile((2,2)) creates: rows repeat, then cols repeat
        expected = np.array([[1, 2, 1, 2],
                            [3, 4, 3, 4],
                            [1, 2, 1, 2],
                            [3, 4, 3, 4]])
        np.testing.assert_array_equal(result, expected)

    def test_rebin_expansion_factor_3(self):
        """Test expansion by factor of 3"""
        array = np.array([[1, 2]])
        result = rebin(array, (3, 6), method='average')
        self.assertEqual(result.shape, (3, 6))
        # Each element repeated 3 times vertically, 3 times horizontally
        expected = np.array([[1, 2, 1, 2, 1, 2],
                            [1, 2, 1, 2, 1, 2],
                            [1, 2, 1, 2, 1, 2]])
        np.testing.assert_array_equal(result, expected)

    def test_rebin_compression_factor_2_sum(self):
        """Test 2x compression with sum"""
        array = np.arange(16).reshape(4, 4)
        result = rebin(array, (2, 2), method='sum')
        # Each 2x2 block: [[0,1],[4,5]]=10, [[2,3],[6,7]]=18, etc.
        expected = np.array([[10, 18], [42, 50]])
        np.testing.assert_array_equal(result, expected)

    def test_rebin_compression_factor_4(self):
        """Test 4x compression"""
        array = np.ones((8, 8))
        result = rebin(array, (2, 2), method='average')
        self.assertEqual(result.shape, (2, 2))
        np.testing.assert_array_equal(result, np.ones((2, 2)))

    def test_rebin_3d_compression_nanmean(self):
        """Test 3D array compression with nanmean"""
        array = np.ones((4, 4, 2), dtype=float)
        array[0, 0, 0] = np.nan
        result = rebin(array, (2, 2), method='nanmean')
        self.assertEqual(result.shape, (2, 2, 2))
        # First channel has one NaN in top-left 2x2 block
        self.assertAlmostEqual(result[0, 0, 0], 1.0, places=5)  # mean of 3 ones
        np.testing.assert_array_almost_equal(result[:, :, 1], np.ones((2, 2)))

    def test_rebin_mixed_compression_factors(self):
        """Test different compression factors for each dimension"""
        array = np.arange(24).reshape(6, 4)
        result = rebin(array, (3, 2), method='average')
        self.assertEqual(result.shape, (3, 2))
        # Each output pixel averages 2x2 input block

    def test_rebin_large_compression(self):
        """Test large compression factor"""
        array = np.ones((16, 16))
        result = rebin(array, (2, 2), method='sum')
        # Each output pixel sums 8x8=64 input pixels
        np.testing.assert_array_equal(result, np.ones((2, 2)) * 64)

    def test_rebin_asymmetric_expansion(self):
        """Test expansion with different factors per dimension"""
        array = np.array([[1, 2]])
        result = rebin(array, (2, 6), method='average')
        expected = np.array([[1, 2, 1, 2, 1, 2],
                            [1, 2, 1, 2, 1, 2]])
        np.testing.assert_array_equal(result, expected)

    def test_rebin_single_element(self):
        """Test rebinning single element array"""
        array = np.array([[5]])
        result_expand = rebin(array, (3, 3), method='average')
        np.testing.assert_array_equal(result_expand, np.ones((3, 3)) * 5)

    def test_rebin_3d_asymmetric_compression(self):
        """Test 3D with different compression per dimension"""
        array = np.ones((8, 4, 3))
        array[:, :, 1] *= 2
        array[:, :, 2] *= 3
        result = rebin(array, (4, 2), method='average')
        self.assertEqual(result.shape, (4, 2, 3))
        np.testing.assert_array_equal(result[:, :, 0], np.ones((4, 2)))
        np.testing.assert_array_equal(result[:, :, 1], np.ones((4, 2)) * 2)
        np.testing.assert_array_equal(result[:, :, 2], np.ones((4, 2)) * 3)

    def test_rebin_preserve_dtype_int(self):
        """Test that integer dtype is preserved in expansion"""
        array = np.array([[1, 2], [3, 4]], dtype=np.int32)
        result = rebin(array, (4, 4), method='average')
        self.assertEqual(result.dtype, np.int32)

    def test_rebin_zero_array(self):
        """Test rebinning array of zeros"""
        array = np.zeros((4, 4))
        result_compress = rebin(array, (2, 2), method='sum')
        np.testing.assert_array_equal(result_compress, np.zeros((2, 2)))
        result_expand = rebin(array, (8, 8), method='average')
        np.testing.assert_array_equal(result_expand, np.zeros((8, 8)))

    def test_rebin_nanmean_partial_nan(self):
        """Test nanmean with partial NaN in blocks"""
        array = np.array([[1, np.nan, 3, 4],
                        [np.nan, 2, 5, 6],
                        [7, 8, np.nan, 10],
                        [9, 10, 11, np.nan]], dtype=float)
        result = rebin(array, (2, 2), method='nanmean')
        # Top-left: mean of [1, nan, nan, 2] = (1+2)/2 = 1.5
        self.assertAlmostEqual(result[0, 0], 1.5, places=5)
