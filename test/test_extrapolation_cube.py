import unittest
import numpy as np
from synim.utils import calculate_extrapolation_indices_coeffs
from synim.utils import apply_extrapolation as _apply_extrapolation
from test import on_backend

# numpy arrays in and out, with the SynIM backend (CPU or GPU)
apply_extrapolation = on_backend(_apply_extrapolation)

def apply_extrapolation_legacy(data, edge_pixels, reference_indices, coefficients):
    """
    Applies linear extrapolation to edge pixels using precalculated indices and coefficients.

    Parameters:
        data (ndarray): Input array to extrapolate.
        edge_pixels (ndarray): Linear indices of edge pixels to extrapolate.
        reference_indices (ndarray): Indices of reference pixels.
        coefficients (ndarray): Coefficients for linear extrapolation.
        debug (bool): If True, displays debug information.
        problem_indices (list): Indices of problematic pixels to analyze.

    Returns:
        ndarray: Array with extrapolated pixels.
    """
    # Create a copy of the input array
    result = data.copy()
    flat_result = result.ravel()
    flat_data = data.ravel()

    # Create a mask for valid reference indices (>= 0)
    valid_ref_mask = reference_indices >= 0
    # Replace invalid indices with 0 to avoid indexing errors
    safe_ref_indices = np.where(valid_ref_mask, reference_indices, 0)
    # Get data values for all reference indices at once
    ref_data = flat_data[safe_ref_indices]  # Shape: (n_valid_edges, 8)
    # Zero out contributions from invalid references
    masked_coeffs = np.where(valid_ref_mask, coefficients, 0.0)
    # Compute all contributions at once and sum across reference positions
    contributions = masked_coeffs * ref_data  # Element-wise multiplication
    extrap_values = np.sum(contributions, axis=1)  # Sum across reference positions
    # Assign extrapolated values to edge pixels
    flat_result[edge_pixels] = extrap_values

    return result


class TestExtrapolationCube(unittest.TestCase):
    def setUp(self):
        # Create a simple mask with a square in the center
        self.n = 20
        self.mask = np.zeros((self.n, self.n), dtype=np.uint8)
        self.mask[5:15, 5:15] = 1

        # Create a 3D data cube with known values
        self.n_modes = 4
        self.data = np.zeros((self.n, self.n, self.n_modes), dtype=float)
        for i in range(self.n_modes):
            self.data[5:15, 5:15, i] = (i + 1) * np.arange(100).reshape(10, 10)

        # Calculate extrapolation indices and coefficients
        self.edge_pixels, self.reference_indices, self.coefficients = calculate_extrapolation_indices_coeffs(
            self.mask, debug=False
        )

    def test_vectorized_vs_loop(self):
        # --- Vectorized extrapolation ---
        data_vec = self.data.copy()
        data_vec_extr = apply_extrapolation(
            data_vec, self.edge_pixels, self.reference_indices, self.coefficients
        )

        # --- Loop (old) extrapolation ---
        data_loop = self.data.copy()
        for i in range(self.n_modes):
            data_loop[:, :, i] = apply_extrapolation_legacy(
                data_loop[:, :, i], self.edge_pixels, self.reference_indices, self.coefficients
            )

        difference = np.abs(data_vec_extr - data_loop)
        max_diff = np.max(difference)
        print(f"Max difference between vectorized and loop extrapolation: {max_diff:.2e}")

        # --- Compare ---
        np.testing.assert_allclose(data_vec_extr, data_loop, rtol=1e-10, atol=1e-12)

    def test_2d_consistency(self):
        # Test that 2D and 3D (single slice) give the same result
        data_2d = self.data[:, :, 0].copy()
        data_3d = self.data[:, :, [0]].copy()  # shape (n, n, 1)

        extr_2d = apply_extrapolation(
            data_2d, self.edge_pixels, self.reference_indices, self.coefficients
        )
        extr_3d = apply_extrapolation(
            data_3d, self.edge_pixels, self.reference_indices, self.coefficients
        )

        difference = np.abs(extr_3d[..., 0] - extr_2d)
        max_diff = np.max(difference)
        print(f"Max difference between 2D and 3D extrapolation: {max_diff:.2e}")

        # extr_3d is (n, n, 1), extr_2d is (n, n)
        np.testing.assert_allclose(extr_3d[..., 0], extr_2d, rtol=1e-10, atol=1e-12)
