"""
Global <-> local mis-registration geometry for WFAO systems (MCAO/GLAO).

This implements the geometry of Sec. 2 of Agapito, Plantet & Heritier,
"SPRINT for WFAO systems", Proc. SPIE 13097, 130975P (2024): the relation
between the LOCAL mis-registration of each WFS-DM pair (shift, rotation,
magnification, anamorphosis of the sub-aperture grid vs the actuator grid)
and the GLOBAL mis-registration geometry of the system (each WFS and each
DM vs the pupil, and each guide star's actual vs nominal position/height).

Explicitly OUT OF SCOPE here (this is deliberately separate from SPRINT,
which already exists in SPECULA):
  - on-sky acquisition of interaction matrices,
  - probe signal / modulation-demodulation,
  - estimation of the LOCAL mis-registration parameters themselves.

This module takes local mis-registration parameters as GIVEN (measured
elsewhere, e.g. by SPRINT, with their own errors) and provides only:
  - the forward geometric map from global parameters to the local
    parameters a local estimator would measure (`System.local_params`);
  - the building blocks (`System.local_transform`) needed by
    `reconstruction.py` to invert that map from a set of (possibly noisy)
    local measurements back to the global geometry.
"""
from dataclasses import dataclass, field
import numpy as np

from .geometry import build_affine, build_dm_affine, decompose_affine

__all__ = ["GuideStar", "DM", "WFS", "System", "gs_parallax_transform"]

# NOTE on yml loading (System.from_params_manager, below): by design this
# reuses synim.params_utils's existing extraction helpers read-only and
# never modifies synim/params_manager.py or synim/params_utils.py. Only
# the WFS mis-registration fields (rotation, shift, magnification,
# anamorphosis_45) have an established key in that yml schema; DM
# shift/magnification and GS position/height error do not (there is no
# nominal-vs-actual split for them there), so they are left at their
# identity default and are meant to be set programmatically afterwards
# (e.g. via `reconstruction.apply_alpha`/`ParameterSpec`).

ARCSEC2RAD = np.pi / 180 / 3600


@dataclass
class GuideStar:
    """Nominal position/height of a guide star, plus its mis-registration
    (actual vs nominal), both to be treated as global unknowns."""
    name: str
    position: tuple = (0.0, 0.0)          # nominal (x, y) on sky, arcsec
    height: float = np.inf                # nominal height, m (np.inf for NGS)
    position_shift: tuple = (0.0, 0.0)    # actual - nominal position, arcsec
    height_shift: float = 0.0             # actual - nominal height, m (LGS only)

    @property
    def actual_position(self):
        return np.array(self.position, dtype=float) + np.array(self.position_shift, dtype=float)

    @property
    def actual_height(self):
        if np.isinf(self.height):
            return self.height
        return self.height + self.height_shift


@dataclass
class DM:
    """A DM's global mis-registration vs the pupil (shift, rotation,
    magnification, anamorphosis of the actuator grid - mirroring `WFS`).
    `height == 0` models a DM that IS the pupil (e.g. an ASM/DSM): combined
    with `gs_parallax_transform` being the identity at zero altitude, this
    makes local == global for this DM, exactly as noted in Sec. 2 of the
    paper."""
    name: str
    height: float = 0.0
    shift: tuple = (0.0, 0.0)
    rotation: float = 0.0
    magnification: float = 1.0
    anamorphosis_45: float = 1.0

    def transform(self):
        """Own mis-registration transform vs the pupil: nominal actuator
        grid -> true physical actuator position."""
        return build_dm_affine(self.shift, self.rotation, self.magnification,
                                self.anamorphosis_45)


@dataclass
class WFS:
    """A WFS's global mis-registration vs the pupil (shift, rotation,
    magnification, anamorphosis of the sub-aperture grid), plus the guide
    star it looks at."""
    name: str
    guide_star: GuideStar
    shift: tuple = (0.0, 0.0)
    rotation: float = 0.0
    magnification: float = 1.0
    anamorphosis_45: float = 1.0

    def transform(self):
        """Own mis-registration transform vs the pupil: nominal
        sub-aperture grid -> true physical sub-aperture position."""
        return build_affine(self.shift, self.rotation, self.magnification,
                             self.anamorphosis_45)


def gs_parallax_transform(guide_star, dm_height, pixel_pitch):
    """
    Forward transform mapping a DM's true physical actuator position (at
    `dm_height`, in the pupil/instrument frame) to the position of its
    footprint in the metapupil as seen by `guide_star` - i.e. the
    cone-effect shift + magnification also computed by
    `synim.utils.shiftzoom_from_source_dm_params`, generalized to a
    (possibly mis-registered) Cartesian guide star position/height so
    that `position_shift`/`height_shift` are directly usable as global
    unknowns.

    Reduces exactly to `shiftzoom_from_source_dm_params` when
    `position_shift == (0, 0)` and `height_shift == 0` - see
    `test_model_geometry.py::test_gs_parallax_matches_shiftzoom_from_source_dm_params`.
    At `dm_height == 0` this is always the identity (no parallax at the
    pupil conjugation, regardless of the guide star), which is what makes
    a `DM` with `height == 0` behave as "the DM is the pupil" in
    `System.local_transform`.
    """
    height = guide_star.actual_height
    if np.isinf(height):
        mag_factor = 1.0
    else:
        mag_factor = height / (height - dm_height)

    pos_m = guide_star.actual_position * ARCSEC2RAD * dm_height
    # Sign flip: same convention as shiftzoom_from_source_dm_params.
    shift_pix = -pos_m / pixel_pitch

    return build_dm_affine(shift=shift_pix, rotation=0.0, magnification=mag_factor)


def _read_shift_mag_anam(params):
    """
    Read (shift, magnification, anamorphosis_45) from a WFS or DM yml
    section, using the same key convention for both (see
    `System.from_params_manager`): `xShiftPhInPixel`/`yShiftPhInPixel` (or
    `translation`), `magnification`, `anamorph45`. All default to
    identity (0 shift, 1 magnification/anamorphosis) when absent, which is
    always the case today for DM sections (no established key there yet).
    """
    x_shift = params.get("xShiftPhInPixel", 0.0)
    y_shift = params.get("yShiftPhInPixel", 0.0)
    shift = tuple(params.get("translation", [x_shift, y_shift]))
    magnification = params.get("magnification", 1.0)
    anamorphosis_45 = params.get("anamorph45", 1.0)
    return shift, magnification, anamorphosis_45


@dataclass
class System:
    """
    A WFAO system: a set of WFSs (each with its own guide star) and a set
    of DMs, plus the pixel pitch used to express the GS-parallax shift in
    the same units as the WFS/DM shift mis-registration parameters.
    """
    wfss: dict = field(default_factory=dict)   # name -> WFS
    dms: dict = field(default_factory=dict)    # name -> DM
    pixel_pitch: float = 1.0

    def add_wfs(self, wfs):
        self.wfss[wfs.name] = wfs
        return wfs

    def add_dm(self, dm):
        self.dms[dm.name] = dm
        return dm

    def local_transform(self, wfs_name, dm_name):
        """
        Forward map for one WFS-DM pair: nominal DM actuator grid ->
        nominal WFS sub-aperture grid. Composes, DM-plane first:
          1. the DM's own mis-registration vs the pupil,
          2. the GS-parallax projection at this DM's height,
          3. the inverse of the WFS's own mis-registration vs the pupil
             (since we want the result expressed in the WFS's own nominal
             grid, not in true physical position).
        """
        wfs = self.wfss[wfs_name]
        dm = self.dms[dm_name]

        dm_transform = dm.transform()
        parallax = gs_parallax_transform(wfs.guide_star, dm.height, self.pixel_pitch)
        wfs_transform = wfs.transform()

        return wfs_transform.inverse().compose(parallax).compose(dm_transform)

    def local_params(self, wfs_name, dm_name):
        """Local (shift, rotation, magnification, anamorphosis_45) for one
        WFS-DM pair, i.e. what a local (SPRINT-like) estimator measures."""
        return decompose_affine(self.local_transform(wfs_name, dm_name))

    def pairs(self):
        """All (wfs_name, dm_name) pairs, in a fixed order."""
        return [(w, d) for w in self.wfss for d in self.dms]

    @classmethod
    def from_params_manager(cls, config, wfs_names=None, dm_names=None, pixel_pitch=None):
        """
        Build a System with NOMINAL geometry (WFS mis-registration, DM
        height/rotation, GS position/height) read from a SynIM
        `ParamsManager`-style configuration, reusing
        `synim.params_utils.extract_wfs_list`/`extract_dm_list`/
        `extract_source_coordinates`/`extract_source_height` read-only (no
        modification to `params_manager.py`/`params_utils.py`).

        DM shift/magnification and GS position/height error have no
        established key in that yml schema (see module docstring) and are
        left at their identity default (0 shift, 1 magnification) - set
        them programmatically afterwards, e.g. via
        `reconstruction.apply_alpha`/`ParameterSpec`, or by editing the
        returned `System`'s `WFS`/`DM`/`GuideStar` fields directly.

        Parameters
        ----------
        config : dict or object with a `.params` dict attribute
            A parsed yml configuration (e.g. `synim.params_utils.
            parse_params_file(path)`), or an already-built `ParamsManager`.
        wfs_names, dm_names : list of str, optional
            Restrict to these WFS/DM config keys (default: all found).
        pixel_pitch : float, optional
            Overrides `config['main']['pixel_pitch']` (default: 1.0 if
            neither is available).
        """
        from ..params_utils import (
            extract_wfs_list, extract_dm_list,
            extract_source_coordinates, extract_source_height,
        )
        from ..utils import polar_to_xy

        raw_config = getattr(config, "params", config)

        if pixel_pitch is None:
            pixel_pitch = raw_config.get("main", {}).get("pixel_pitch", 1.0)

        system = cls(pixel_pitch=pixel_pitch)

        for entry in extract_wfs_list(raw_config):
            if wfs_names is not None and entry["name"] not in wfs_names:
                continue
            wfs_key = entry["name"]
            wfs_params = entry["config"]

            rotation = wfs_params.get("rotation", wfs_params.get("rotAnglePhInDeg", 0.0))
            shift, magnification, anamorphosis_45 = _read_shift_mag_anam(wfs_params)

            gs_r, gs_theta_deg = extract_source_coordinates(raw_config, wfs_key)
            gs_position = tuple(polar_to_xy(gs_r, np.deg2rad(gs_theta_deg)))
            gs_height = extract_source_height(raw_config, wfs_key)

            guide_star = GuideStar(name=f"gs_{wfs_key}", position=gs_position, height=gs_height)
            system.add_wfs(WFS(name=wfs_key, guide_star=guide_star, shift=shift,
                                rotation=rotation, magnification=magnification,
                                anamorphosis_45=anamorphosis_45))

        for entry in extract_dm_list(raw_config):
            if dm_names is not None and entry["name"] not in dm_names:
                continue
            dm_params = entry["config"]
            dm_shift, dm_magnification, dm_anamorphosis_45 = _read_shift_mag_anam(dm_params)
            system.add_dm(DM(name=entry["name"],
                              height=dm_params.get("height", 0.0),
                              shift=dm_shift,
                              rotation=dm_params.get("rotation", 0.0),
                              magnification=dm_magnification,
                              anamorphosis_45=dm_anamorphosis_45))

        return system
