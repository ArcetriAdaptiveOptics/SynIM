import unittest
import numpy as np

from synim.registration.geometry import (
    AffineTransform, build_affine, decompose_affine, rotation_matrix,
)


class TestGeometryRoundTrip(unittest.TestCase):
    """Round-trip build_affine -> decompose_affine over a grid of parameters."""

    def _check_roundtrip(self, shift, rotation, magnification, anamorphosis_45,
                          anamorphosis_90=1.0):
        t = build_affine(shift, rotation, magnification, anamorphosis_45, anamorphosis_90)
        rec_shift, rec_rot, rec_mag, rec_anam45, rec_anam90 = decompose_affine(t)

        np.testing.assert_allclose(rec_shift, shift, atol=1e-9)
        self.assertAlmostEqual(
            ((rec_rot - rotation + 180) % 360) - 180, 0.0, places=6,
            msg=f"rotation mismatch: got {rec_rot}, expected {rotation}"
        )
        self.assertAlmostEqual(rec_mag, magnification, places=6)
        self.assertAlmostEqual(rec_anam45, anamorphosis_45, places=6)
        self.assertAlmostEqual(rec_anam90, anamorphosis_90, places=6)

    def test_identity(self):
        self._check_roundtrip((0.0, 0.0), 0.0, 1.0, 1.0)

    def test_shift_only(self):
        self._check_roundtrip((3.5, -2.1), 0.0, 1.0, 1.0)

    def test_rotation_only(self):
        for rot in [-179, -90, -45, -1, 1, 10, 45, 90, 135, 179]:
            self._check_roundtrip((0.0, 0.0), rot, 1.0, 1.0)

    def test_magnification_only(self):
        for mag in [0.5, 0.8, 1.0, 1.01, 1.5, 2.0]:
            self._check_roundtrip((0.0, 0.0), 0.0, mag, 1.0)

    def test_anamorphosis_45_only(self):
        for k in [0.5, 0.8, 1.0, 1.2, 1.5]:
            self._check_roundtrip((0.0, 0.0), 0.0, 1.0, k)

    def test_anamorphosis_90_only(self):
        for k in [0.5, 0.8, 1.0, 1.2, 1.5]:
            self._check_roundtrip((0.0, 0.0), 0.0, 1.0, 1.0, k)

    def test_anamorphosis_45_and_90_combined(self):
        for k45, k90 in [(0.8, 1.3), (1.2, 0.7), (0.9, 0.9), (1.4, 1.4)]:
            self._check_roundtrip((0.0, 0.0), 0.0, 1.0, k45, k90)
            self._check_roundtrip((0.0, 0.0), 37.0, 1.1, k45, k90)

    def test_combined_grid(self):
        rng = np.random.default_rng(0)
        for _ in range(200):
            shift = tuple(rng.uniform(-5, 5, size=2))
            rotation = rng.uniform(-179, 179)
            magnification = rng.uniform(0.5, 1.5)
            anam45 = rng.uniform(0.7, 1.3)
            anam90 = rng.uniform(0.7, 1.3)
            self._check_roundtrip(shift, rotation, magnification, anam45, anam90)


class TestGeometryComposition(unittest.TestCase):
    """Sanity checks on AffineTransform composition semantics."""

    def test_compose_is_function_composition(self):
        a = build_affine(shift=(1.0, 0.0), rotation=30.0, magnification=1.2)
        b = build_affine(shift=(0.0, 2.0), rotation=-10.0, magnification=0.9)

        p = np.array([1.3, -0.7])
        composed = a.compose(b)

        np.testing.assert_allclose(composed(p), a(b(p)), atol=1e-10)

    def test_inverse(self):
        a = build_affine(shift=(2.0, -1.0), rotation=17.0, magnification=1.3,
                          anamorphosis_45=1.1)
        identity = a.inverse().compose(a)
        np.testing.assert_allclose(identity.M, np.eye(2), atol=1e-10)
        np.testing.assert_allclose(identity.t, np.zeros(2), atol=1e-10)

    def test_pure_rotation_is_orthogonal(self):
        R = rotation_matrix(37.0)
        np.testing.assert_allclose(R @ R.T, np.eye(2), atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(R), 1.0, places=12)


if __name__ == "__main__":
    unittest.main()
