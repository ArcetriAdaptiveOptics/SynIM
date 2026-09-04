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

`build_affine` builds `M` from the five named mis-registration parameters
used throughout this project and in `synim.utils.rotshiftzoom_array`/
`synim.synim`:

    shift          (x, y) translation, in the same units as position
    rotation       degrees, standard image/matrix convention (+Y is "down",
                   a positive angle rotates +X towards +Y)
    magnification  isotropic scale factor (1.0 = no magnification); alone,
                   this is the "one parameter for simplicity" of the SPIE
                   paper
    anamorphosis_45  diagonal (45 deg) shear/stretch factor (1.0 = none),
                   same definition as `wfs_anamorphosis_45` in
                   `synim.utils.rotshiftzoom_array`
    anamorphosis_90  ratio between the Y and X magnification (1.0 = none,
                   i.e. isotropic), same definition as `wfs_anamorphosis_90`
                   in `synim.synim` (there folded directly into an
                   anisotropic `wfs_magnification = (mag, mag *
                   anamorphosis_90)` rather than kept as a separate factor)

Together, magnification, anamorphosis_45 and anamorphosis_90 span every
possible pure shape distortion (no rotation) of a 2x2 linear map - three
numbers for three degrees of freedom - so no further "shape" parameter is
needed.

The composition order (magnification and anamorphosis_90 together as an
anisotropic scale, then 45 deg anamorphosis, then rotation, with the shift
added last and therefore unaffected by the linear part) was
reverse-engineered from, and is validated against, the physical behaviour
of `rotshiftzoom_array`/`synim.synim` as pinned down by
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
                  anamorphosis_45=1.0, anamorphosis_90=1.0):
    """
    Build the forward `AffineTransform` for a WFS-side (or generic/local)
    mis-registration:

        S = magnification * diag(1, anamorphosis_90)
        y = R(rotation) @ Anam(anamorphosis_45) @ (S @ x) + shift

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
    anamorphosis_90 : float
        Ratio between the Y and X magnification (1.0 = none). This
        directly matches `wfs_anamorphosis_90` as used in `synim.synim`
        (``wfs_magnification = (mag, mag * anamorphosis_90)``) - no
        reciprocal correction needed, unlike `anamorphosis_45`.
    """
    S = np.diag([1.0, float(anamorphosis_90)]) * float(magnification)
    A = anamorphosis_matrix(float(anamorphosis_45))
    R = rotation_matrix(float(rotation))
    M = R @ A @ S
    t = np.array(shift, dtype=float)
    return AffineTransform(M, t)


def build_dm_affine(shift=(0.0, 0.0), rotation=0.0, magnification=1.0,
                     anamorphosis_45=1.0, anamorphosis_90=1.0):
    """
    Build the forward `AffineTransform` for a DM-side mis-registration
    relative to the pupil (or the GS-parallax/cone-effect map from
    `shiftzoom_from_source_dm_params`, which uses neither anamorphosis
    term):

        S = magnification * diag(1, anamorphosis_90)
        y = R(rotation) @ Anam(anamorphosis_45) @ (S @ x) + magnification * shift

    Unlike `build_affine`, `shift` here IS scaled by `magnification` (it
    represents e.g. a GS footprint offset expressed at the DM plane, which
    is naturally subject to the same cone-effect scaling as everything
    else) but, like in `build_affine`, is still unaffected by `rotation`
    and by the anamorphosis terms (a pure shape distortion of the DM's own
    actuator grid does not move a shift already expressed as a physical
    offset). `rotshiftzoom_array` has no DM-side anamorphosis term (only
    the WFS side has one); this is a registration-package extension for a
    mis-registration category the SPIE paper allows for ("X and Y
    magnification or higher order distortion can also be considered",
    Sec. 2) but the pixel-based pipeline does not implement on the DM
    side, so there is no `rotshiftzoom_array` cross-check for it.

    Parameters
    ----------
    shift : (x, y)
    rotation : float
        Degrees.
    magnification : float
        Isotropic scale factor (1.0 = none).
    anamorphosis_45 : float
        Diagonal (45 deg) shear/stretch factor (1.0 = none).
    anamorphosis_90 : float
        Ratio between the Y and X magnification (1.0 = none).
    """
    S = np.diag([1.0, float(anamorphosis_90)]) * float(magnification)
    A = anamorphosis_matrix(float(anamorphosis_45))
    R = rotation_matrix(float(rotation))
    M = R @ A @ S
    mag = float(magnification)
    t = mag * np.array(shift, dtype=float)
    return AffineTransform(M, t)


def decompose_affine(transform):
    """
    Invert `build_affine`: recover (shift, rotation, magnification,
    anamorphosis_45, anamorphosis_90) from an `AffineTransform`.

    This is an exact, closed-form inverse, built from the Gram matrix
    ``G = M.T @ M``, which is invariant under the left rotation `R` (i.e.
    it only sees the shape part ``N = Anam(k45) @ diag(mag, mag*k90)``,
    which the rotation does not affect: ``G = N.T @ N``). Using
    ``Anam(k45) @ Anam(k45) == Anam(k45**2)`` (Anam is a fixed-eigenbasis
    matrix, so squaring it squares its eigenvalues) gives, with
    ``p = (1 + k45**2) / 2`` and ``q = (1 - k45**2) / 2``:

        G00 = mag**2 * p
        G11 = (mag*k90)**2 * p
        G01 = mag**2 * k90 * q

    ``r = G01 / sqrt(G00*G11)`` therefore equals ``q / p``, giving
    ``k45**2 = (1-r)/(1+r)`` (the positive root: an anamorphosis factor is
    a ratio of two positive lengths). ``mag`` and ``k90`` follow from
    ``G00``/``G11`` and ``p``, and finally `rotation` from how `M`'s first
    column compares to the now fully known shape matrix `N`'s first
    column (see the source for the exact, short derivation).

    Returns
    -------
    shift : (x, y) ndarray
    rotation : float, degrees
    magnification : float
    anamorphosis_45 : float
    anamorphosis_90 : float
    """
    M = transform.M
    G = M.T @ M
    if G[0, 0] <= 0.0 or G[1, 1] <= 0.0:
        raise ValueError("Singular affine transform: cannot recover rotation/magnification/anamorphosis.")

    # |r| <= 1 always holds: G is a Gram matrix, hence positive
    # semi-definite, i.e. G01**2 <= G00*G11 (Cauchy-Schwarz).
    r = G[0, 1] / np.sqrt(G[0, 0] * G[1, 1])
    anamorphosis_45 = float(np.sqrt((1.0 - r) / (1.0 + r)))

    p = (1.0 + anamorphosis_45 ** 2) / 2.0
    d1 = float(np.sqrt(G[0, 0] / p))   # = magnification
    d2 = float(np.sqrt(G[1, 1] / p))   # = magnification * anamorphosis_90
    magnification = d1
    anamorphosis_90 = d2 / d1

    # N = Anam(anamorphosis_45) @ diag(d1, d2); its first column is
    # (a*d1, b*d1) with a=(1+k45)/2, b=(1-k45)/2. M's first column is
    # R(rotation) applied to that (known) vector, so rotation follows by
    # solving R(rotation) @ (a*d1, b*d1) = M[:, 0] for the angle (a
    # 2-equation, 1-unknown system, consistent by construction).
    a = (1.0 + anamorphosis_45) / 2.0
    b = (1.0 - anamorphosis_45) / 2.0
    tx, ty = M[0, 0], M[1, 0]
    cos_r = (a * tx + b * ty) / (p * d1)
    sin_r = (a * ty - b * tx) / (p * d1)
    rotation = float(np.degrees(np.arctan2(sin_r, cos_r)))

    return (np.array(transform.t, dtype=float), rotation, magnification,
            anamorphosis_45, anamorphosis_90)
