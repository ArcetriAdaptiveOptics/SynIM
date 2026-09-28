"""
MORFEO LGS interaction matrices: SynIM vs the SPECULA Shack-Hartmann.

For each LGS WFS and a few modes of dm1, the SPECULA slopes are computed
with the SH and ShSlopec processing objects (push-pull, noiseless
intensity, plain centre of gravity) and compared with the SynIM
interaction matrix computed with the same ParamsManager WFS parameters.
The modes are used as a phase in the pupil plane: SynIM is called without
DM height, DM rotation and guide star geometry (the SH receives the same
phase), so the comparison covers the WFS part (WFS rotation, shift and
magnification, subapertures and slopes), not the DM footprint.
MORFEO has 480 pixels and 68 subapertures (7.06 pixels per subaperture).

The SH slopes include optical effects that the geometric SynIM model does
not have: the comparison uses a gain per mode fitted on the fully
illuminated subapertures, and reports the relative rms difference on
    - fully illuminated subapertures (flux > 0.999)
    - partially illuminated subapertures with flux >= 0.5
    - partially illuminated subapertures with flux < 0.5
Run it on two SynIM branches to compare them, e.g.
    SYNIM_TEST_DEVICE_IDX=0 python -m unittest test.test_morfeo_intmat
The results are also saved to morfeo_intmat.txt (current directory).
"""
import inspect
import os
import unittest

import numpy as np

import synim
import synim.synim as synim_core
from synim.utils import rotshiftzoom_array

YAML_FILE = '/home/guido/pythonLib/SPECULA_scripts/morfeo/params_morfeo_calib.yml'
ROOT_DIR = '/raid1/guido/PASSATA/MAORYC'

DM_INDEX = 1
MODES = [0, 1, 2, 3, 5, 10, 30, 100]
AMPLITUDE_NM = 20.0          # push-pull amplitude of a mode with unit rms on the pupil
SLOPE_METHODS = ('derivatives', 'telsum')

AVAILABLE = os.path.exists(YAML_FILE) and os.path.exists(ROOT_DIR)
if AVAILABLE:
    import specula
    if specula.xp is None:
        specula.init(device_idx=synim.default_target_device_idx,
                     precision=synim.global_precision)
    from specula.data_objects.electric_field import ElectricField
    from specula.data_objects.subap_data import SubapData
    from specula.data_objects.pixels import Pixels
    from specula.processing_objects.sh import SH
    from specula.processing_objects.sh_slopec import ShSlopec
    from synim.params_manager import ParamsManager


def _scalar_kwargs(cls, config):
    """Entries of config that are scalar arguments of cls.__init__."""
    names = inspect.signature(cls.__init__).parameters
    return {k: v for k, v in config.items()
            if k in names and isinstance(v, (int, float, bool, str))
            and k not in ('target_device_idx', 'precision', 'data_dir')}


@unittest.skipUnless(AVAILABLE, f'MORFEO configuration not found at {YAML_FILE}')
class TestMorfeoIntmatSpecula(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pm = ParamsManager(YAML_FILE, root_dir=ROOT_DIR, verbose=False)
        cls.pixel_pitch = cls.pm.params['main']['pixel_pitch']
        cls.n_lgs = cls.pm._count_wfs('lgs')
        cls.results = {}

    def _specula(self, wfs_idx, pup, modes):
        """Push-pull slopes (n_slopes, n_modes) and flux per valid subaperture."""
        sh_cfg = self.pm.params[f'sh_lgs{wfs_idx}']
        slopec_cfg = self.pm.params[f'slopec_lgs{wfs_idx}']
        if slopec_cfg.get('interleave', False):
            self.skipTest('interleaved slopes not supported by this test')
        subapdata = SubapData.restore(self.pm.cm.filename('subapdata',
                                                          slopec_cfg['subapdata_object']))
        n = pup.shape[0]
        ef = ElectricField(n, n, pixel_pitch=self.pixel_pitch, S0=1)
        ef.A = pup
        sh = SH(**_scalar_kwargs(SH, sh_cfg))
        sh.inputs['in_ef'].set(ef)
        slopec = ShSlopec(subapdata)
        pixels = Pixels(*sh.outputs['out_i'].i.shape)
        slopec.inputs['in_pixels'].set(pixels)
        state = {'t': 0, 'setup': False}

        def measure(phase_nm):
            state['t'] += 1
            ef.phaseInNm = phase_nm
            ef.generation_time = state['t']
            if not state['setup']:
                sh.setup()
                state['setup'] = True
            sh.check_ready(state['t'])
            sh.trigger()
            sh.post_trigger()
            pixels.pixels = sh.outputs['out_i'].i
            pixels.generation_time = state['t']
            slopec.check_ready(state['t'])
            slopec.trigger()
            slopec.post_trigger()
            return (specula.cpuArray(slopec.outputs['out_slopes'].slopes).astype(np.float64),
                    specula.cpuArray(slopec.outputs['out_flux_per_subaperture'].value)
                    .astype(np.float64))

        columns = []
        for k in range(modes.shape[2]):
            push, flux = measure(AMPLITUDE_NM / 2 * modes[:, :, k])
            pull, _ = measure(-AMPLITUDE_NM / 2 * modes[:, :, k])
            columns.append((push - pull) / AMPLITUDE_NM)
        return np.stack(columns, axis=1), flux / flux.max()

    def _compare(self, wfs_idx):
        params = self.pm.prepare_interaction_matrix_params('lgs', wfs_idx, DM_INDEX, None)
        pup = synim.cpuArray(params['pup_mask']).astype(np.float32)
        n = pup.shape[0]
        modes = synim.xp.asarray(params['dm_array'][:, :, MODES], dtype=np.float32)
        dm_mask = synim.xp.asarray(params['dm_mask'], dtype=np.float32)
        if modes.shape[0] != n:
            # DM larger than the pupil (meta-pupil): central part, as SynIM
            # does for a DM at the ground and an on-axis source
            modes = rotshiftzoom_array(modes, output_size=(n, n))
            dm_mask = rotshiftzoom_array(dm_mask, output_size=(n, n))
            dm_mask[dm_mask < 0.5] = 0
        modes = synim.cpuArray(modes)
        dm_mask = synim.cpuArray(dm_mask)
        rms = np.sqrt(np.mean(modes[pup > 0] ** 2, axis=0))
        modes = modes / rms

        specula_im, flux = self._specula(wfs_idx, pup, modes)
        rows = {}
        for slope_method in SLOPE_METHODS:
            im = synim.cpuArray(synim_core.interaction_matrix(
                pup_diam_m=params['pup_diam_m'], pup_mask=pup, dm_array=modes,
                dm_mask=dm_mask,
                # phase in the pupil plane, as for the SH (see the module docstring)
                dm_height=0.0, dm_rotation=0.0, gs_pol_coo=(0.0, 0.0), gs_height=np.inf,
                wfs_nsubaps=params['wfs_nsubaps'], wfs_rotation=params['wfs_rotation'],
                wfs_translation=params['wfs_translation'],
                wfs_mag_global=params['wfs_magnification'],
                wfs_fov_arcsec=params['wfs_fov_arcsec'],
                idx_valid_sa=params['idx_valid_sa'], slope_method=slope_method))
            self.assertEqual(im.shape, specula_im.shape)
            n_nan = int(np.isnan(im).sum())
            im = np.nan_to_num(im)
            flux2 = np.tile(flux, 2)
            full = flux2 > 0.999
            high = (flux2 >= 0.5) & ~full
            low = flux2 < 0.5
            errors = []
            for k in range(len(MODES)):
                a, b = specula_im[:, k], im[:, k]
                gain = np.dot(a[full], b[full]) / np.dot(b[full], b[full])
                norm = np.sqrt(np.mean(a[full] ** 2))
                d = a - gain * b
                errors.append([gain] + [np.sqrt(np.mean(d[m] ** 2)) / norm if m.any() else np.nan
                                        for m in (full, high, low)])
            rows[slope_method] = (n_nan, np.array(errors), int(full.sum() // 2),
                                  int(high.sum() // 2), int(low.sum() // 2))
        return params['wfs_rotation'], rows

    def test_morfeo_lgs_intmat(self):
        lines = []
        for wfs_idx in range(1, self.n_lgs + 1):
            rotation, rows = self._compare(wfs_idx)
            for slope_method, (n_nan, errors, n_full, n_high, n_low) in rows.items():
                lines.append(f'LGS {wfs_idx} (rotation {rotation}), {slope_method}: NaN {n_nan},'
                             f' subapertures full {n_full}, >=0.5 {n_high}, <0.5 {n_low}')
                for mode, (gain, e_full, e_high, e_low) in zip(MODES, errors):
                    lines.append(f'    mode {mode:4d}: gain {gain:+.4f}  rel rms diff'
                                 f' full {e_full:.3f}  >=0.5 {e_high:.3f}  <0.5 {e_low:.3f}')
                self.results[(wfs_idx, slope_method)] = (n_nan, errors)
            print('\n'.join(lines[-2 * (len(MODES) + 1):]))

        with open('morfeo_intmat.txt', 'w') as f:
            f.write(f'SynIM {os.path.dirname(synim.__file__)}\n' + '\n'.join(lines) + '\n')

        for (wfs_idx, slope_method), (n_nan, errors) in self.results.items():
            with self.subTest(lgs=wfs_idx, slope_method=slope_method):
                self.assertEqual(n_nan, 0)
                # same sign and similar optical gain for all the modes
                gains = errors[:, 0]
                self.assertTrue((gains > 0).all(), gains)
                self.assertLess(gains.max() / gains.min(), 1.3, gains)
                # fully illuminated subapertures: geometric model close to the SH
                self.assertLess(errors[:, 1].max(), 0.1, errors[:, 1])


if __name__ == '__main__':
    unittest.main()
