import unittest
import numpy as np

from synim.utils import shiftzoom_from_source_dm_params, polar_to_xy
from synim.registration.geometry import build_affine, decompose_affine
from synim.registration.model import GuideStar, DM, WFS, System, gs_parallax_transform


class TestGsParallaxVsExistingFormula(unittest.TestCase):
    """`gs_parallax_transform` must reduce exactly to
    `shiftzoom_from_source_dm_params` for a non-mis-registered guide star."""

    def test_matches_for_random_geometries(self):
        rng = np.random.default_rng(1)
        pixel_pitch = 0.2  # m / sub-aperture

        for _ in range(50):
            r_asec = rng.uniform(0, 60)
            theta_deg = rng.uniform(0, 360)
            dm_height = rng.uniform(0, 15000)
            # occasionally test the NGS (infinite height) case
            gs_height = np.inf if rng.random() < 0.2 else rng.uniform(dm_height + 1000, 100000)

            gs = GuideStar(name="gs", position=tuple(polar_to_xy(r_asec, np.deg2rad(theta_deg))),
                            height=gs_height)

            expected_shift, expected_zoom = shiftzoom_from_source_dm_params(
                (r_asec, theta_deg), gs_height, dm_height, pixel_pitch
            )
            transform = gs_parallax_transform(gs, dm_height, pixel_pitch)

            # `shiftzoom_from_source_dm_params` returns the raw shift meant
            # to be fed as `dm_translation` to `rotshiftzoom_array`, which
            # internally scales it by `dm_magnification` (see
            # `build_dm_affine`'s docstring/convention) - so the physical
            # shift is expected_shift * expected_zoom, matching transform.t.
            np.testing.assert_allclose(transform.t, np.array(expected_shift) * np.array(expected_zoom),
                                        atol=1e-9)
            np.testing.assert_allclose(np.diag(transform.M), expected_zoom, atol=1e-9)
            # pure scale, no rotation
            np.testing.assert_allclose(transform.M, np.diag(expected_zoom), atol=1e-9)

    def test_zero_altitude_is_always_identity(self):
        """No parallax/cone effect at the pupil conjugation, whatever the GS."""
        rng = np.random.default_rng(2)
        for _ in range(20):
            gs = GuideStar(name="gs", position=tuple(rng.uniform(-30, 30, size=2)),
                            height=rng.uniform(9000, 20000),
                            position_shift=tuple(rng.uniform(-1, 1, size=2)),
                            height_shift=rng.uniform(-500, 500))
            transform = gs_parallax_transform(gs, dm_height=0.0, pixel_pitch=0.2)
            np.testing.assert_allclose(transform.M, np.eye(2), atol=1e-12)
            np.testing.assert_allclose(transform.t, (0.0, 0.0), atol=1e-12)


class TestDmIsPupilCase(unittest.TestCase):
    """Sec. 2 of the paper: when a DM is the pupil, local == global for
    the WFS-DM registration (the WFS-pupil and DM-pupil groups merge)."""

    def test_local_equals_wfs_inverse_compose_dm_when_no_parallax(self):
        gs = GuideStar(name="ngs", position=(0.0, 0.0), height=np.inf)
        wfs = WFS(name="wfs1", guide_star=gs, shift=(1.3, -0.4), rotation=2.0,
                  magnification=1.01, anamorphosis_45=1.02)
        dm = DM(name="dm_pupil", height=0.0, shift=(0.5, 0.2), rotation=-1.5,
                magnification=0.99)

        system = System(pixel_pitch=0.2)
        system.add_wfs(wfs)
        system.add_dm(dm)

        expected = wfs.transform().inverse().compose(dm.transform())
        actual = system.local_transform("wfs1", "dm_pupil")

        np.testing.assert_allclose(actual.M, expected.M, atol=1e-10)
        np.testing.assert_allclose(actual.t, expected.t, atol=1e-10)

    def test_all_nominal_gives_identity_local_transform(self):
        gs = GuideStar(name="ngs", position=(0.0, 0.0), height=np.inf)
        wfs = WFS(name="wfs1", guide_star=gs)
        dm = DM(name="dm1", height=4500.0)

        system = System(pixel_pitch=0.2)
        system.add_wfs(wfs)
        system.add_dm(dm)

        shift, rotation, magnification, anam = system.local_params("wfs1", "dm1")
        np.testing.assert_allclose(shift, (0.0, 0.0), atol=1e-10)
        self.assertAlmostEqual(rotation, 0.0, places=8)
        self.assertAlmostEqual(magnification, 1.0, places=8)
        self.assertAlmostEqual(anam, 1.0, places=8)


class TestGsPositionErrorIsObservableLocally(unittest.TestCase):
    """A GS position error at a post-focal DM shows up as a local shift
    (entangled with an actual WFS-DM shift - Sec. 2/5 of the paper)."""

    def test_gs_shift_produces_expected_local_shift(self):
        pixel_pitch = 0.2
        dm_height = 6000.0
        gs_height = 90000.0
        pos_shift_asec = (2.0, -1.0)

        gs = GuideStar(name="lgs", position=(0.0, 0.0), height=gs_height,
                        position_shift=pos_shift_asec)
        wfs = WFS(name="wfs1", guide_star=gs)   # no own mis-registration
        dm = DM(name="dm2", height=dm_height)   # no own mis-registration

        system = System(pixel_pitch=pixel_pitch)
        system.add_wfs(wfs)
        system.add_dm(dm)

        shift, rotation, magnification, anam = system.local_params("wfs1", "dm2")

        mag_factor = gs_height / (gs_height - dm_height)
        expected_shift = -(np.array(pos_shift_asec) * (np.pi / 180 / 3600)
                            * dm_height / pixel_pitch) * mag_factor

        np.testing.assert_allclose(shift, expected_shift, atol=1e-8)
        self.assertAlmostEqual(rotation, 0.0, places=8)
        self.assertAlmostEqual(magnification, mag_factor, places=8)
        self.assertAlmostEqual(anam, 1.0, places=8)


class TestMultiWfsDmSystem(unittest.TestCase):
    """Smoke test: build a small MAVIS-like (3 WFS x 3 DM) system and make
    sure every pair produces a well-formed local transform."""

    def test_builds_and_decomposes_all_pairs(self):
        system = System(pixel_pitch=0.2)
        for i, (r, theta) in enumerate([(0.0, 0.0), (17.5, 0.0), (17.5, 120.0)]):
            gs = GuideStar(name=f"gs{i}", position=tuple(polar_to_xy(r, np.deg2rad(theta))),
                            height=90000.0 if r > 0 else np.inf)
            system.add_wfs(WFS(name=f"wfs{i}", guide_star=gs,
                                shift=(0.1 * i, -0.05 * i), rotation=0.5 * i))
        for j, height in enumerate([0.0, 6000.0, 13500.0]):
            system.add_dm(DM(name=f"dm{j}", height=height,
                              shift=(0.05 * j, 0.0), rotation=-0.2 * j))

        for wfs_name, dm_name in system.pairs():
            shift, rotation, magnification, anam = system.local_params(wfs_name, dm_name)
            self.assertTrue(np.all(np.isfinite(shift)))
            self.assertTrue(np.isfinite(rotation))
            self.assertGreater(magnification, 0.0)
            self.assertGreater(anam, 0.0)


if __name__ == "__main__":
    unittest.main()
