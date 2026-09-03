"""
Analytic 2D affine geometry for WFAO mis-registration modeling.

Convention
----------
Every optical element (a WFS relative to the pupil, a DM relative to the
pupil, the local WFS-DM pairing, ...) is described by a forward affine map

    y = M @ x + t

from a nominal (undistorted) 2D coordinate `x` (e.g. a sub-aperture or
actuator index expressed as a physical position) to the true, mis-registered
position `y`, both expressed in the same physical unit (e.g. sub-aperture
pitch).

`build_affine` builds `M` from the four named mis-registration parameters
used throughout this project and in `synim.utils.rotshiftzoom_array`:

    shift          (x, y) translation, in the same units as position
    rotation       degrees, standard image/matrix convention (+Y is "down",
                   a positive angle rotates +X towards +Y)
    magnification  isotropic scale factor (1.0 = no magnification); this is
                   the "one parameter for simplicity" of the SPIE paper
    anamorphosis_45  diagonal (45 deg) shear/stretch factor (1.0 = none),
                   same definition as `wfs_anamorphosis_45` in
                   `synim.utils.rotshiftzoom_array`

The composition order (magnification, then 45 deg anamorphosis, then
rotation, with the shift added last and therefore unaffected by the linear
part) was reverse-engineered from, and is validated against, the physical
behaviour of `rotshiftzoom_array` as pinned down by
`test/test_rotshiftzoom.py` (DM/WFS shift independent of/scaled by
magnification and rotation) - see
`test/test_geometry_vs_pixel_synim.py` for the numerical cross-check.

`AffineTransform` objects compose like functions (`a.compose(b)` means
"apply b, then a") so that the global mis-registration parameters of a WFS
and of a DM can be chained into the effective local WFS-DM transform that
SPRINT actually measures on sky.
"""
import numpy as np

__all__ = [
    "AffineTransform",
    "rotation_matrix",
    "anamorphosis_matrix",
    "build_affine",
    "build_dm_affine",
    "decompose_affine",
]


def rotation_matrix(angle_deg):
    """2x2 rotation matrix, standard matrix convention (+Y down, CCW positive)."""
    theta = np.deg2rad(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def anamorphosis_matrix(k):
    """
    Diagonal (45 deg) shear/stretch matrix.

    Eigenvalue 1 along the (1, 1) direction (the +45 deg diagonal) and
    eigenvalue k along the (1, -1) direction (the -45 deg diagonal); this
    matches `wfs_anamorphosis_45` in `synim.utils.rotshiftzoom_array`.
    k == 1.0 is the identity (no anamorphosis).
    """
    if k == 1.0:
        return np.eye(2)
    return np.array([[(1.0 + k) / 2.0, (1.0 - k) / 2.0],
                      [(1.0 - k) / 2.0, (1.0 + k) / 2.0]])


class AffineTransform:
    """
    A forward 2D affine map ``y = M @ x + t``.

    Composition follows function-composition semantics: ``a.compose(b)``
    returns the transform "first b, then a", i.e. ``a.compose(b)(x) ==
    a(b(x))``.
    """
    __slots__ = ("M", "t")

    def __init__(self, M=None, t=None):
        self.M = np.eye(2) if M is None else np.array(M, dtype=float)
        self.t = np.zeros(2) if t is None else np.array(t, dtype=float)

    @classmethod
    def identity(cls):
        return cls()

    def __call__(self, p):
        p = np.asarray(p, dtype=float)
        return self.M @ p + self.t

    def compose(self, other):
        """Return ``self`` after ``other``: apply `other` first, then `self`."""
        return AffineTransform(self.M @ other.M, self.M @ other.t + self.t)

    def inverse(self):
        Minv = np.linalg.inv(self.M)
        return AffineTransform(Minv, -Minv @ self.t)

    def __repr__(self):
        return f"AffineTransform(M={self.M.tolist()}, t={self.t.tolist()})"

    def __eq__(self, other):
        if not isinstance(other, AffineTransform):
            return NotImplemented
        return np.allclose(self.M, other.M) and np.allclose(self.t, other.t)


def build_affine(shift=(0.0, 0.0), rotation=0.0, magnification=1.0,
                  anamorphosis_45=1.0):
    """
    Build the forward `AffineTransform` for a WFS-side (or generic/local)
    mis-registration:

        y = R(rotation) @ Anam(anamorphosis_45) @ (magnification * x) + shift

    `shift` is added last (i.e. NOT affected by rotation, magnification or
    anamorphosis): it represents a pure detector shift, expressed directly
    in the common/output frame. This is also the right convention to
    describe a generic *local* WFS-DM transform (e.g. what SPRINT measures,
    or what `decompose_affine` recovers from a composed transform), since
    at that point "shift" has no remaining meaning of "pre- or
    post-magnification".

    For a DM element's mis-registration relative to the pupil, use
    `build_dm_affine` instead: there `shift` represents a GS footprint
    offset that gets scaled by the cone effect (`magnification`) together
    with everything else - see `test_geometry_vs_pixel_synim.py`.

    Parameters
    ----------
    shift : (x, y)
    rotation : float
        Degrees.
    magnification : float
        Isotropic scale factor (1.0 = none).
    anamorphosis_45 : float
        Diagonal (45 deg) shear/stretch factor (1.0 = none). Note: this is
        the *physical* stretch factor along the (1, -1) diagonal, which is
        the reciprocal of `wfs_anamorphosis_45` as passed to
        `synim.utils.rotshiftzoom_array` (verified empirically - see
        `test_geometry_vs_pixel_synim.py::test_anamorphosis_convention`).
    """
    S = np.eye(2) * float(magnification)
    A = anamorphosis_matrix(float(anamorphosis_45))
    R = rotation_matrix(float(rotation))
    M = R @ A @ S
    t = np.array(shift, dtype=float)
    return AffineTransform(M, t)


def build_dm_affine(shift=(0.0, 0.0), rotation=0.0, magnification=1.0,
                     anamorphosis_45=1.0):
    """
    Build the forward `AffineTransform` for a DM-side mis-registration
    relative to the pupil (or the GS-parallax/cone-effect map from
    `shiftzoom_from_source_dm_params`, which does not use `anamorphosis_45`):

        y = R(rotation) @ Anam(anamorphosis_45) @ (magnification * x) + magnification * shift

    Unlike `build_affine`, `shift` here IS scaled by `magnification` (it
    represents e.g. a GS footprint offset expressed at the DM plane, which
    is naturally subject to the same cone-effect scaling as everything
    else) but, like in `build_affine`, is still unaffected by `rotation`
    and `anamorphosis_45` (a pure shape distortion of the DM's own
    actuator grid does not move a shift already expressed as a physical
    offset). `rotshiftzoom_array` has no DM-side anamorphosis term (only
    `wfs_anamorphosis_45` on the WFS side); this is a registration-package
    extension for a mis-registration category the SPIE paper allows for
    ("X and Y magnification or higher order distortion can also be
    considered", Sec. 2) but the pixel-based pipeline does not yet
    implement, so there is no `rotshiftzoom_array` cross-check for it.

    Parameters
    ----------
    shift : (x, y)
    rotation : float
        Degrees.
    magnification : float
        Isotropic scale factor (1.0 = none).
    anamorphosis_45 : float
        Diagonal (45 deg) shear/stretch factor (1.0 = none).
    """
    S = np.eye(2) * float(magnification)
    A = anamorphosis_matrix(float(anamorphosis_45))
    R = rotation_matrix(float(rotation))
    M = R @ A @ S
    mag = float(magnification)
    t = mag * np.array(shift, dtype=float)
    return AffineTransform(M, t)


def decompose_affine(transform):
    """
    Invert `build_affine`: recover (shift, rotation, magnification,
    anamorphosis_45) from an `AffineTransform`.

    This is an exact, closed-form inverse under the assumption used
    throughout this project that magnification is isotropic (the "one
    parameter for simplicity" of the SPIE paper). Writing
    ``M = R(rotation) @ Anam(anamorphosis_45) @ (magnification * I)`` out
    with ``a = (1+k)/2``, ``b = (1-k)/2`` (``k`` = anamorphosis_45) gives:

        M00 + M11 = magnification * cos(rotation) * (1 + k)
        M10 - M01 = magnification * sin(rotation) * (1 + k)
        M01 + M10 = magnification * cos(rotation) * (1 - k)
        M11 - M00 = magnification * sin(rotation) * (1 - k)

    so the norms of the two pairs above are ``magnification * (1 + k)`` and
    ``magnification * |1 - k|`` (assuming k > -1, i.e. magnification * (1+k)
    > 0, which always holds for a physical anamorphosis); the sign lost in
    the second norm is recovered from whether ``(M01+M10, M11-M00)`` points
    along or against ``(M00+M11, M10-M01)`` (i.e. whether k is below or
    above 1). ``k == 1`` (no anamorphosis) is the well behaved limit where
    the second pair is identically zero.

    Returns
    -------
    shift : (x, y) ndarray
    rotation : float, degrees
    magnification : float
    anamorphosis_45 : float
    """
    M = transform.M
    P, S_ = M[0, 0] + M[1, 1], M[1, 0] - M[0, 1]  # mag * (1 + k) * (cos, sin)
    R_, Q = M[0, 1] + M[1, 0], M[1, 1] - M[0, 0]  # mag * (1 - k) * (cos, sin)

    plus = np.hypot(P, S_)                         # mag * (1 + k)
    minus_sign = np.sign(P * R_ + S_ * Q) or 1.0
    minus = minus_sign * np.hypot(R_, Q)           # mag * (1 - k), signed

    magnification = float((plus + minus) / 2.0)
    if magnification == 0.0:
        raise ValueError("Singular affine transform: cannot recover rotation/anamorphosis.")
    anamorphosis_45 = float((plus - minus) / (2.0 * magnification))
    rotation = float(np.degrees(np.arctan2(S_, P)))

    return np.array(transform.t, dtype=float), rotation, magnification, anamorphosis_45
