"""
Default number of modes of a DM / layer component without 'nmodes':
- count_component_modes returns the number of modes in the files: columns of the
  M2C if there is one, otherwise rows of the influence function
- the count read from the FITS header equals the count obtained by loading the
  objects (SPECULA-saved file and raw PASSATA-like FITS)
- it is also the number of modes returned by load_influence_functions
- components without influence function files give None
- ParamsManager writes the total into params['nmodes'] when missing, leaves an
  explicit value unchanged, and the IM filename is the same as with the
  explicit total
"""
import copy
import os
import tempfile
import unittest

import numpy as np
import yaml
from astropy.io import fits

from specula.data_objects.ifunc import IFunc
from specula.data_objects.m2c import M2C
from specula.calib_manager import CalibManager

from synim.params_manager import ParamsManager
from synim.params_utils import (
    count_component_modes,
    generate_im_filename,
    load_influence_functions,
)

HERE = os.path.dirname(os.path.abspath(__file__))
YAML_FILE = os.path.join(HERE, 'params_scao_sh_test.yml')

N_PIX = 16          # side of the small pupil mask
N_MODES = 7         # modes of the modal basis
N_ACT = 12          # actuators of the zonal basis
N_M2C_MODES = 5     # modes of the M2C (< N_ACT)


def _mask():
    yy, xx = np.mgrid[:N_PIX, :N_PIX]
    return (np.hypot(xx - N_PIX / 2 + 0.5, yy - N_PIX / 2 + 0.5) < 6).astype(np.uint8)


class _TempCalibMixin:
    """Temporary calibration root with the ifunc/ and m2c/ subdirectories."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = self._tmp.name
        os.makedirs(os.path.join(self.root, 'ifunc'))
        os.makedirs(os.path.join(self.root, 'm2c'))
        self.cm = CalibManager(self.root)
        self.mask = _mask()
        self.n_valid = int(self.mask.sum())
        self.rng = np.random.default_rng(0)

    def _save_ifunc(self, tag, n_rows):
        arr = self.rng.standard_normal((n_rows, self.n_valid)).astype(np.float32)
        IFunc(ifunc=arr, mask=self.mask).save(self.cm.filename('ifunc', tag + '.fits'))
        return arr

    def _save_m2c(self, tag, n_act, n_modes):
        arr = self.rng.standard_normal((n_act, n_modes)).astype(np.float32)
        M2C(m2c=arr).save(self.cm.filename('m2c', tag + '.fits'))
        return arr


class TestCountComponentModes(_TempCalibMixin, unittest.TestCase):

    def test_specula_modal_ifunc(self):
        self._save_ifunc('modal', N_MODES)
        n = count_component_modes(self.cm, {'ifunc_tag': 'modal'})
        self.assertEqual(n, N_MODES)

    def test_passata_like_raw_fits(self):
        # Empty primary HDU, HDU 1 = [n_valid, n_modes] (unnamed), HDU 2 = mask
        arr = self.rng.standard_normal((self.n_valid, N_MODES)).astype(np.float32)
        path = self.cm.filename('ifunc', 'passata.fits')
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(arr), fits.ImageHDU(self.mask)]
                     ).writeto(path)

        restored = IFunc.restore(path)
        self.assertEqual(restored.influence_function.shape[0], N_MODES)
        np.testing.assert_allclose(restored.influence_function, arr.T)
        self.assertEqual(count_component_modes(self.cm, {'ifunc_tag': 'passata'}), N_MODES)

    def test_zonal_ifunc_with_m2c(self):
        self._save_ifunc('zonal', N_ACT)
        self._save_m2c('m2c_a', N_ACT, N_M2C_MODES)
        params = {'ifunc_tag': 'zonal', 'm2c_tag': 'm2c_a'}
        n = count_component_modes(self.cm, params)
        self.assertEqual(n, N_M2C_MODES)

        # Same number as the modes returned by load_influence_functions
        # (3D array [ny, nx, n_modes])
        dm_array, _ = load_influence_functions(self.cm, params, pixel_pupil=N_PIX)
        self.assertEqual(dm_array.shape[-1], n)
        dm_2d, _ = load_influence_functions(self.cm, params, pixel_pupil=N_PIX, as_2d=True)
        self.assertEqual(dm_2d.shape[0], n)

    def test_zonal_ifunc_without_m2c_counts_actuators(self):
        self._save_ifunc('zonal', N_ACT)
        self.assertEqual(count_component_modes(self.cm, {'ifunc_tag': 'zonal'}), N_ACT)

    def test_modal_ifunc_matches_load_influence_functions(self):
        self._save_ifunc('modal', N_MODES)
        params = {'ifunc_tag': 'modal'}
        n = count_component_modes(self.cm, params)
        dm_array, _ = load_influence_functions(self.cm, params, pixel_pupil=N_PIX)
        self.assertEqual(dm_array.shape[-1], n)

    def test_header_count_equals_loaded_count(self):
        # A non-square case where rows != columns on disk, to catch an axis swap
        self._save_ifunc('zonal', N_ACT)
        self._save_m2c('m2c_a', N_ACT, N_M2C_MODES)
        self._save_ifunc('modal', N_MODES)
        loaded_ifunc = IFunc.restore(self.cm.filename('ifunc', 'modal.fits')
                                     ).influence_function.shape[0]
        loaded_m2c = M2C.restore(self.cm.filename('m2c', 'm2c_a.fits')).m2c.shape[1]
        self.assertEqual(count_component_modes(self.cm, {'ifunc_tag': 'modal'}), loaded_ifunc)
        self.assertEqual(
            count_component_modes(self.cm, {'ifunc_tag': 'zonal', 'm2c_tag': 'm2c_a'}),
            loaded_m2c)

    def test_no_ifunc_returns_none(self):
        self.assertIsNone(count_component_modes(self.cm, {'type_str': 'zernike', 'nmodes': 5}))
        self.assertIsNone(count_component_modes(self.cm, {}))


class TestParamsManagerDefaultNmodes(_TempCalibMixin, unittest.TestCase):

    def setUp(self):
        super().setUp()
        self._save_ifunc('zonal', N_ACT)
        self._save_ifunc('modal', N_MODES)
        self._save_m2c('m2c_a', N_ACT, N_M2C_MODES)

    def _config(self, dm1=None, layer1=None):
        with open(YAML_FILE) as f:
            config = yaml.safe_load(f)
        del config['dm']
        config['dm1'] = dict(dm1 or {})
        config['dm1'].setdefault('height', 0.0)
        if layer1 is not None:
            config['layer1'] = dict(layer1)
        return config

    def _pm(self, config):
        return ParamsManager(copy.deepcopy(config), root_dir=self.root, verbose=False)

    def test_default_is_total_of_m2c(self):
        pm = self._pm(self._config(dm1={'ifunc_tag': 'zonal', 'm2c_tag': 'm2c_a'}))
        self.assertEqual(pm.params['dm1']['nmodes'], N_M2C_MODES)

    def test_default_is_total_of_modal_ifunc_for_dm_and_layer(self):
        config = self._config(dm1={'ifunc_tag': 'modal'},
                              layer1={'ifunc_tag': 'modal', 'height': 5000.0})
        pm = self._pm(config)
        self.assertEqual(pm.params['dm1']['nmodes'], N_MODES)
        self.assertEqual(pm.params['layer1']['nmodes'], N_MODES)

    def test_explicit_nmodes_unchanged(self):
        pm = self._pm(self._config(dm1={'ifunc_tag': 'zonal', 'm2c_tag': 'm2c_a',
                                        'nmodes': 3}))
        self.assertEqual(pm.params['dm1']['nmodes'], 3)

    def test_component_without_ifunc_unchanged(self):
        config = self._config(dm1={'type_str': 'zernike'})
        pm = self._pm(config)
        self.assertNotIn('nmodes', pm.params['dm1'])

    def test_im_filename_same_as_explicit_total(self):
        kwargs = dict(wfs_type='ngs', wfs_index=1, dm_index=1)
        implicit = self._config(dm1={'ifunc_tag': 'zonal', 'm2c_tag': 'm2c_a'})
        explicit = self._config(dm1={'ifunc_tag': 'zonal', 'm2c_tag': 'm2c_a',
                                     'nmodes': N_M2C_MODES})

        name_implicit = self._pm(implicit).generate_im_filename(**kwargs)
        name_explicit = self._pm(explicit).generate_im_filename(**kwargs)
        self.assertEqual(name_implicit, name_explicit)
        # the filename carries the number of modes (sanity: not a trivial equality)
        self.assertIn(f'{N_M2C_MODES}modes', name_implicit)
        # and without the default it would differ
        self.assertNotEqual(generate_im_filename(implicit, **kwargs), name_implicit)


if __name__ == '__main__':
    unittest.main()
