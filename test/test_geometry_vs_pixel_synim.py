import unittest
import numpy as np

from synim.utils import rotshiftzoom_array
from synim.registration.geometry import build_affine, build_dm_affine


def _centroid(img, threshold_frac=0.1):
    """Subpixel centroid of the thresholded blob, in (row, col) = (y, x)."""
    img = np.asarray(img, dtype=float)
    img = np.where(img > threshold_frac * img.max(), img, 0.0)
    ys, xs = np.indices(img.shape)
    total = img.sum()
    return np.array([(ys * img).sum() / total, (xs * img).sum() / total])


class TestGeometryVsPixelSynim(unittest.TestCase):
    """
    Cross-check the parameter-space `AffineTransform` model of
    `synim.registration.geometry` against the pixel-based
    `synim.utils.rotshiftzoom_array`, for the shift / rotation / isotropic
    magnification degrees of freedom (the paper's core mis-registration
    parametrization).

    `rotshiftzoom_array` works in (row, col) = (y, x) pixel coordinates
    with the array center as origin (see `test_rotshiftzoom.py`);
    `build_affine`/`AffineTransform` are convention-agnostic 2-vectors, so
    we simply feed them (y, x) pairs throughout this test to match.
    """

    def setUp(self):
        self.size = 200
        self.center = self.size / 2.0
        x = np.arange(self.size)
        y = np.arange(self.size)
        self._xx, self._yy = np.meshgrid(x, y)
        # Off-center blob: asymmetric offset so rotation is also observable.
        self.blob_center = np.array([self.center + 15.0, self.center + 8.0])  # (y, x)
        self.blob = np.exp(-(((self._yy - self.blob_center[0]) ** 2) / (2 * 6 ** 2)
                              + ((self._xx - self.blob_center[1]) ** 2) / (2 * 6 ** 2)))

    def _predicted_position(self, dm_translation=(0.0, 0.0), dm_rotation=0.0,
                             dm_magnification=1.0, wfs_translation=(0.0, 0.0),
                             wfs_rotation=0.0, wfs_magnification=1.0):
        dm = build_dm_affine(shift=dm_translation, rotation=dm_rotation,
                              magnification=dm_magnification)
        wfs = build_affine(shift=wfs_translation, rotation=wfs_rotation,
                            magnification=wfs_magnification)
        combined = wfs.compose(dm)  # DM step first, then WFS step
        p0 = self.blob_center - self.center
        return combined(p0) + self.center

    def _run(self, predicted_kwargs, **rszz_kwargs):
        transformed = rotshiftzoom_array(
            self.blob, output_size=(self.size, self.size), **rszz_kwargs
        )
        actual = _centroid(transformed)
        predicted = self._predicted_position(**predicted_kwargs)
        np.testing.assert_allclose(actual, predicted, atol=0.5,
                                    err_msg=f"actual={actual} predicted={predicted}")

    def test_dm_translation(self):
        t = (12.0, -7.0)
        self._run(dict(dm_translation=t), dm_translation=t)

    def test_wfs_translation(self):
        t = (-10.0, 6.0)
        self._run(dict(wfs_translation=t), wfs_translation=t)

    def test_dm_translation_scaled_by_dm_magnification(self):
        t = (20.0, -10.0)
        mag = 0.5
        self._run(dict(dm_translation=t, dm_magnification=mag),
                   dm_translation=t, dm_magnification=(mag, mag))

    def test_wfs_translation_independent_of_wfs_magnification(self):
        t = (20.0, 0.0)
        mag = 2.0
        self._run(dict(wfs_translation=t, wfs_magnification=mag),
                   wfs_translation=t, wfs_magnification=(mag, mag))

    def test_dm_translation_invariant_under_dm_rotation(self):
        t = (30.0, 0.0)
        rot = 90.0
        self._run(dict(dm_translation=t, dm_rotation=rot),
                   dm_translation=t, dm_rotation=rot)

    def test_wfs_translation_invariant_under_wfs_rotation(self):
        t = (20.0, 0.0)
        rot = 45.0
        self._run(dict(wfs_translation=t, wfs_rotation=rot),
                   wfs_translation=t, wfs_rotation=rot)

    def test_dm_translation_rotates_with_wfs(self):
        t = (30.0, 0.0)
        rot = 90.0
        self._run(dict(dm_translation=t, wfs_rotation=rot),
                   dm_translation=t, wfs_rotation=rot)

    def test_off_axis_orbit_arbitrary_rotation_and_magnification(self):
        dm_t = (20.0, 10.0)
        wfs_rot = 15.0
        wfs_mag = 1.05
        self._run(dict(dm_translation=dm_t, wfs_rotation=wfs_rot, wfs_magnification=wfs_mag),
                   dm_translation=dm_t, wfs_rotation=wfs_rot,
                   wfs_magnification=(wfs_mag, wfs_mag))

    def test_full_combination(self):
        dm_t = (8.0, -4.0)
        dm_rot = 12.0
        dm_mag = 0.92
        wfs_t = (-6.0, 3.0)
        wfs_rot = -20.0
        wfs_mag = 1.15
        self._run(
            dict(dm_translation=dm_t, dm_rotation=dm_rot, dm_magnification=dm_mag,
                 wfs_translation=wfs_t, wfs_rotation=wfs_rot, wfs_magnification=wfs_mag),
            dm_translation=dm_t, dm_rotation=dm_rot, dm_magnification=(dm_mag, dm_mag),
            wfs_translation=wfs_t, wfs_rotation=wfs_rot, wfs_magnification=(wfs_mag, wfs_mag),
        )


    def test_anamorphosis_convention(self):
        """
        `wfs_anamorphosis_45=k` in rotshiftzoom_array physically stretches
        the (1, -1) diagonal by 1/k (not k) and leaves the (1, 1) diagonal
        untouched; `build_affine(anamorphosis_45=...)` uses the physical
        stretch factor directly, so the two are reciprocal.
        """
        d = 20.0
        blob = (np.exp(-(((self._yy - (self.center + d)) ** 2
                           + (self._xx - (self.center + d)) ** 2) / (2 * 4 ** 2)))
                + np.exp(-(((self._yy - (self.center + d)) ** 2
                            + (self._xx - (self.center - d)) ** 2) / (2 * 4 ** 2))))

        k_code = 2.0
        out = rotshiftzoom_array(blob, wfs_anamorphosis_45=k_code,
                                  output_size=(self.size, self.size))

        from scipy.ndimage import label, center_of_mass
        mask = out > 0.3 * out.max()
        lab, n = label(mask)
        self.assertEqual(n, 2)
        coms = np.array(center_of_mass(out, lab, range(1, n + 1))) - self.center

        # (1, 1) diagonal point: untouched (offset (d, d))
        on_plus_diag = coms[np.argmin(np.abs(coms[:, 0] - coms[:, 1]))]
        np.testing.assert_allclose(on_plus_diag, (d, d), atol=0.5)

        # (1, -1) diagonal point: stretched by 1 / k_code (offset (d, -d))
        on_minus_diag = coms[np.argmax(np.abs(coms[:, 0] - coms[:, 1]))]
        np.testing.assert_allclose(on_minus_diag, (d / k_code, -d / k_code), atol=0.5)


if __name__ == "__main__":
    unittest.main()
