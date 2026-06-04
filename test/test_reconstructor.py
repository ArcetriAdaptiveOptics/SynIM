"""Test reconstructor computation methods (MMSE vs pseudo-inverse).

Unit tests for compute_pseudoinverse_reconstructor and compute_mmse_reconstructor functions.
These tests validate the core functionality of both reconstructor computation methods.
"""

import unittest
import numpy as np
import specula

specula.init(device_idx=-1, precision=1)

from synim.params_utils import (
    compute_mmse_reconstructor,
    compute_pseudoinverse_reconstructor
)


class TestPseudoInverseReconstructor(unittest.TestCase):
    """Test the pseudo-inverse reconstructor function."""
    
    def test_pseudoinverse_basic(self):
        """Test that pseudo-inverse satisfies basic properties."""
        # Create a simple interaction matrix
        n_slopes = 100
        n_modes = 20
        
        # Random interaction matrix
        np.random.seed(42)
        A = np.random.randn(n_slopes, n_modes).astype(np.float32)
        
        # Compute pseudo-inverse reconstructor
        R = compute_pseudoinverse_reconstructor(A, verbose=False)
        
        # Check shape
        self.assertEqual(R.shape, (n_modes, n_slopes))
        
        # Check pseudo-inverse property: R @ A should be close to identity
        product = R @ A
        identity = np.eye(n_modes)
        
        # Should be close to identity (within numerical precision)
        np.testing.assert_allclose(product, identity, atol=1e-5)
        
    def test_pseudoinverse_overdetermined(self):
        """Test pseudo-inverse with overdetermined system (more equations than unknowns)."""
        # More slopes than modes (typical case)
        n_slopes = 200
        n_modes = 50
        
        np.random.seed(43)
        A = np.random.randn(n_slopes, n_modes).astype(np.float32)
        
        R = compute_pseudoinverse_reconstructor(A, verbose=False)
        
        # Reconstruct should satisfy A @ R @ A ≈ A
        reconstructed_A = A @ R @ A
        np.testing.assert_allclose(reconstructed_A, A, atol=1e-4)
        
    def test_pseudoinverse_reconstruction(self):
        """Test that pseudo-inverse can reconstruct known modes."""
        n_slopes = 100
        n_modes = 20
        
        np.random.seed(44)
        A = np.random.randn(n_slopes, n_modes).astype(np.float32)
        
        # Create a known mode vector
        true_modes = np.random.randn(n_modes).astype(np.float32)
        
        # Generate slopes from modes
        slopes = A @ true_modes
        
        # Compute reconstructor
        R = compute_pseudoinverse_reconstructor(A, verbose=False)
        
        # Reconstruct modes
        reconstructed_modes = R @ slopes
        
        # Should recover the original modes
        np.testing.assert_allclose(reconstructed_modes, true_modes, atol=1e-4)


class TestMMSEReconstructor(unittest.TestCase):
    """Test the MMSE reconstructor function."""
    
    def test_mmse_basic(self):
        """Test that MMSE reconstructor can be computed."""
        n_slopes = 100
        n_modes = 20
        
        np.random.seed(45)
        A = np.random.randn(n_slopes, n_modes).astype(np.float32)
        
        # Simple diagonal atmospheric covariance
        C_atm = np.eye(n_modes, dtype=np.float32)
        
        # Scalar noise variance
        noise_var = 1.0
        
        # Compute MMSE reconstructor
        R = compute_mmse_reconstructor(
            A, C_atm, 
            noise_variance=None,
            C_noise=noise_var * np.eye(n_slopes, dtype=np.float32),
            cinverse=False,
            verbose=False
        )
        
        # Check shape
        self.assertEqual(R.shape, (n_modes, n_slopes))
        
        # Check it's finite
        self.assertTrue(np.all(np.isfinite(R)))


if __name__ == '__main__':
    unittest.main()
