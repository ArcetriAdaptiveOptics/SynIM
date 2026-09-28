"""
SynIM interaction matrix vs a SPECULA Shack-Hartmann (SH processing object,
push-pull, centre of gravity over the subaperture pixels) with a
non-integer number of pixels per subaperture (90 and 100 pixels, 12
subapertures) and, as a control, an integer one (96 pixels).

The SH slopes include optical effects (diffraction, finite field of view)
that the geometric model of SynIM does not have: the comparison uses a
gain per mode fitted on the fully illuminated subapertures, and the
subapertures with less than half of the flux are not compared.

Relative rms difference measured with the modes of this test (max over
the modes; fully illuminated / partially illuminated with flux >= 0.5):
    96 px (integer, control)   derivatives 0.020 / 0.112, telsum 0.020 / 0.093
    90 and 100 px              derivatives 0.019 / 0.111, telsum 0.019 / 0.099
    before (upscaled rebin)    derivatives NaN (90 px) and 0.021 / 0.250 (100 px)
The residual in the partially illuminated subapertures is the same as with
an integer ratio: it comes from the optical effects, not from the rebinning.
"""
import contextlib
import io
import unittest

import numpy as np

import synim
import synim.synim as synim_core

try:
    import specula
    if specula.xp is None:
        specula.init(device_idx=synim.default_target_device_idx,
                     precision=synim.global_precision)
    from specula.data_objects.electric_field import ElectricField
    from specula.processing_objects.sh import SH
    HAVE_SPECULA = True
except ImportError:  # pragma: no cover
    HAVE_SPECULA = False

N_SUBAPS = 12
SUBAP_NPX = 8
AMPLITUDE_NM = 20.0


def _pupil(n, obstruction=0.25):
    c = (n - 1) / 2
    yy, xx = np.mgrid[:n, :n]
    r = np.hypot(xx - c, yy - c) / (n / 2)
    return ((r <= 1) & (r >= obstruction)).astype(np.float32)


def _modes(n):
    c = (n - 1) / 2
    yy, xx = np.mgrid[:n, :n]
    x = (xx - c) / (n / 2)
    y = (yy - c) / (n / 2)
    return np.stack([x, y, x * y, x ** 2 - y ** 2, x ** 3 - 2 * x * y ** 2],
                    axis=2).astype(np.float32)


def _sh_slopes(phase_nm, pup):
    """x, y slopes (centre of gravity / half the subaperture pixels) and flux."""
    n = pup.shape[0]
    ef = ElectricField(n, n, pixel_pitch=1.0 / n, S0=1, target_device_idx=-1)
    ef.A = pup
    ef.phaseInNm = phase_nm.astype(np.float32)
    ef.generation_time = 1
    sh = SH(wavelengthInNm=600, subap_wanted_fov=4.0, sensor_pxscale=0.5,
            subap_on_diameter=N_SUBAPS, subap_npx=SUBAP_NPX, target_device_idx=-1)
    sh.inputs['in_ef'].set(ef)
    with contextlib.redirect_stdout(io.StringIO()):
        sh.setup()
        sh.check_ready(1)
        sh.trigger()
        sh.post_trigger()
    image = np.asarray(specula.cpuArray(sh.outputs['out_i'].i), dtype=np.float64)
    blocks = image.reshape(N_SUBAPS, SUBAP_NPX, N_SUBAPS, SUBAP_NPX)
    flux = blocks.sum(axis=(1, 3))
    pix = np.arange(SUBAP_NPX) - (SUBAP_NPX - 1) / 2
    safe = np.where(flux > 0, flux, 1)
    cx = (blocks * pix[None, None, None, :]).sum(axis=(1, 3)) / safe
    cy = (blocks * pix[None, :, None, None]).sum(axis=(1, 3)) / safe
    return cx / (SUBAP_NPX / 2), cy / (SUBAP_NPX / 2), flux / flux.max()


@unittest.skipUnless(HAVE_SPECULA, 'SPECULA not available')
class TestShFractionalSpecula(unittest.TestCase):

    specula_cache = {}

    def _specula_im(self, n):
        if n not in self.specula_cache:
            pup, modes = _pupil(n), _modes(n)
            sx, sy = [], []
            for k in range(modes.shape[2]):
                px, py, flux = _sh_slopes(AMPLITUDE_NM / 2 * modes[:, :, k], pup)
                mx, my, _ = _sh_slopes(-AMPLITUDE_NM / 2 * modes[:, :, k], pup)
                sx.append((px - mx) / AMPLITUDE_NM)
                sy.append((py - my) / AMPLITUDE_NM)
            self.specula_cache[n] = (np.stack(sx, 2), np.stack(sy, 2), flux)
        return self.specula_cache[n]

    def _compare(self, n, slope_method):
        sx, sy, flux = self._specula_im(n)
        pup, modes = _pupil(n), _modes(n)
        im = synim.cpuArray(synim_core.interaction_matrix(
            1.0, pup, modes, np.ones((n, n), dtype=np.float32), 0.0, 0.0, N_SUBAPS, 4.0,
            (0.0, 0.0), np.inf, 0.0, (0.0, 0.0), 1.0, specula_convention=False,
            slope_method=slope_method))
        self.assertFalse(np.isnan(im).any())
        n2 = N_SUBAPS ** 2
        ix = im[:n2].reshape(N_SUBAPS, N_SUBAPS, -1)
        iy = im[n2:].reshape(N_SUBAPS, N_SUBAPS, -1)
        full = flux > 0.999
        partial = (flux >= 0.5) & ~full
        err_full, err_partial = [], []
        for k in range(modes.shape[2]):
            a = np.concatenate([sx[..., k][full], sy[..., k][full]])
            b = np.concatenate([ix[..., k][full], iy[..., k][full]])
            gain = np.dot(a, b) / np.dot(b, b)
            # same sign and similar optical gain for all the modes
            self.assertTrue(0.75 < gain < 0.95, f'mode {k}: gain {gain}')
            norm = np.sqrt(np.mean(a ** 2))

            def rms(mask):
                d = np.concatenate([(sx[..., k] - gain * ix[..., k])[mask],
                                    (sy[..., k] - gain * iy[..., k])[mask]])
                return np.sqrt(np.mean(d ** 2)) / norm
            err_full.append(rms(full))
            err_partial.append(rms(partial))
        return np.array(err_full), np.array(err_partial)

    def test_non_integer_ratio(self):
        for n in (90, 100, 96):
            for slope_method in ('derivatives', 'telsum'):
                with self.subTest(n=n, slope_method=slope_method):
                    err_full, err_partial = self._compare(n, slope_method)
                    self.assertLess(err_full.max(), 0.03, err_full)
                    self.assertLess(err_partial.max(), 0.13, err_partial)


if __name__ == '__main__':
    unittest.main()
