import os
import unittest
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless: just check the schematic renders, no display

from synim.utils import polar_to_xy
from synim.registration.model import GuideStar, DM, WFS, System
from synim.registration.reconstruction import ParameterSpec, jacobian
from synim.registration.analysis import svd_of_jacobian
from synim.registration.viz import (
    plot_altitude_schematic, plot_mis_registration_table,
    plot_dm_footprint, plot_dm_footprints, plot_system_overview,
    plot_mode_bar, plot_mode_gs_quiver,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _mavis_like_system():
    system = System(pixel_pitch=0.2)
    gs_specs = [(0.0, 0.0, np.inf), (17.5, 0.0, 90000.0), (17.5, 120.0, 90000.0)]
    for i, (r, theta, height) in enumerate(gs_specs):
        gs = GuideStar(name=f"gs{i}", position=tuple(polar_to_xy(r, np.deg2rad(theta))),
                        height=height)
        system.add_wfs(WFS(name=f"wfs{i}", guide_star=gs, shift=(0.1 * i, -0.05 * i),
                            rotation=0.5 * i))
    for j, height in enumerate([0.0, 6000.0, 13500.0]):
        system.add_dm(DM(name=f"dm{j}", height=height, shift=(0.05 * j, 0.0),
                          rotation=-0.2 * j))
    return system


class TestAltitudeSchematic(unittest.TestCase):
    def test_renders_mavis_like_system_without_error(self):
        system = _mavis_like_system()
        ax = plot_altitude_schematic(system)
        self.assertIsNotNone(ax)

        os.makedirs(OUT_DIR, exist_ok=True)
        ax.figure.savefig(os.path.join(OUT_DIR, "_test_altitude_schematic_mavis.png"),
                           bbox_inches="tight")

    def test_renders_single_wfs_single_dm_system(self):
        gs = GuideStar(name="ngs", position=(0.0, 0.0), height=np.inf)
        system = System(pixel_pitch=0.2)
        system.add_wfs(WFS(name="wfs0", guide_star=gs))
        system.add_dm(DM(name="dm0", height=0.0))
        ax = plot_altitude_schematic(system)
        self.assertIsNotNone(ax)

    def test_slice_axis_parameter(self):
        system = _mavis_like_system()
        ax = plot_altitude_schematic(system, slice_axis=1)
        self.assertIsNotNone(ax)

    def test_can_reuse_provided_axes(self):
        import matplotlib.pyplot as plt
        system = _mavis_like_system()
        fig, ax = plt.subplots()
        returned_ax = plot_altitude_schematic(system, ax=ax)
        self.assertIs(returned_ax, ax)


class TestMisRegistrationTable(unittest.TestCase):
    def test_renders_mavis_like_system_without_error(self):
        system = _mavis_like_system()
        ax = plot_mis_registration_table(system)
        self.assertIsNotNone(ax)

    def test_renders_system_with_no_dms_or_wfss(self):
        system = System(pixel_pitch=0.2)
        ax = plot_mis_registration_table(system)
        self.assertIsNotNone(ax)


class TestDmFootprint(unittest.TestCase):
    def test_renders_single_dm_without_error(self):
        system = _mavis_like_system()
        ax = plot_dm_footprint(system, "dm2", pupil_diameter=39.0)
        self.assertIsNotNone(ax)

    def test_renders_all_dms_without_error(self):
        system = _mavis_like_system()
        fig, axes = plot_dm_footprints(system, pupil_diameter=39.0)
        self.assertEqual(len(axes), 3)

    def test_ngs_footprint_equals_pupil_diameter(self):
        # An NGS's cone never converges: the footprint on any DM has the
        # same diameter as the pupil itself.
        gs = GuideStar(name="ngs", position=(0.0, 0.0), height=np.inf)
        system = System(pixel_pitch=0.2)
        system.add_wfs(WFS(name="wfs0", guide_star=gs))
        system.add_dm(DM(name="dm0", height=6000.0))
        ax = plot_dm_footprint(system, "dm0", pupil_diameter=39.0)
        circles = [p for p in ax.patches if p.get_edgecolor() is not None]
        # last-added circle is the WFS footprint (pupil reference circle
        # is added first)
        self.assertAlmostEqual(circles[-1].get_radius() * 2, 39.0)

    def test_technical_fov_draws_extra_circle(self):
        system = _mavis_like_system()
        ax_no_fov = plot_dm_footprint(system, "dm1", pupil_diameter=39.0)
        ax_with_fov = plot_dm_footprint(system, "dm1", pupil_diameter=39.0,
                                          technical_fov_arcsec=120.0)
        self.assertEqual(len(ax_with_fov.patches), len(ax_no_fov.patches) + 1)


class TestSystemOverview(unittest.TestCase):
    def test_renders_and_saves_mavis_like_system_table_only(self):
        system = _mavis_like_system()
        fig, axes = plot_system_overview(system)
        self.assertIn("schematic", axes)
        self.assertIn("table", axes)
        self.assertNotIn("footprints", axes)

    def test_renders_and_saves_mavis_like_system_with_footprints(self):
        system = _mavis_like_system()
        fig, axes = plot_system_overview(system, pupil_diameter=39.0)
        self.assertIn("footprints", axes)
        self.assertEqual(len(axes["footprints"]), 3)

        os.makedirs(OUT_DIR, exist_ok=True)
        fig.savefig(os.path.join(OUT_DIR, "_test_system_overview_mavis.png"),
                    bbox_inches="tight", dpi=150)


class TestModeVisualization(unittest.TestCase):
    def _system_and_mode(self):
        system = _mavis_like_system()
        pairs = system.pairs()
        specs = ([ParameterSpec("wfs", w, "shift", axis) for w in ["wfs0", "wfs1", "wfs2"]
                  for axis in (0, 1)]
                 + [ParameterSpec("gs", f"gs{i}", "position_shift", axis)
                    for i in range(3) for axis in (0, 1)])
        Lambda = jacobian(system, specs, pairs, dof=("shift_x", "shift_y"))
        _, _, Vt = svd_of_jacobian(Lambda)
        return system, specs, Vt[-1]  # worst-determined mode

    def test_plot_mode_bar_renders(self):
        _, specs, coeffs = self._system_and_mode()
        ax = plot_mode_bar(specs, coeffs)
        self.assertIsNotNone(ax)

    def test_plot_mode_bar_top_n_limits_bars(self):
        _, specs, coeffs = self._system_and_mode()
        ax = plot_mode_bar(specs, coeffs, top_n=3)
        self.assertEqual(len(ax.get_yticks()), 3)

    def test_plot_mode_gs_quiver_renders(self):
        system, specs, coeffs = self._system_and_mode()
        ax = plot_mode_gs_quiver(system, specs, coeffs)
        self.assertIsNotNone(ax)

    def test_plot_mode_gs_quiver_raises_without_gs_specs(self):
        system, _, _ = self._system_and_mode()
        specs = [ParameterSpec("wfs", "wfs0", "shift", 0)]
        with self.assertRaises(ValueError):
            plot_mode_gs_quiver(system, specs, [1.0])


if __name__ == "__main__":
    unittest.main()
