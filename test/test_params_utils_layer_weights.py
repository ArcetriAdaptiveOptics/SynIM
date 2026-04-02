import unittest

import numpy as np

from synim.params_utils import compute_layer_weights_from_turbulence


class TestLayerWeightsFromTurbulence(unittest.TestCase):
    def test_uses_airmass_scaled_heights_like_idl(self):
        # With zenith=60deg, airmass=2, turbulence heights are doubled before
        # nearest-layer assignment (IDL-compatible behavior).
        params = {
            'main': {'zenithAngleInDeg': 60.0},
            'atmo': {
                'heights': [800.0, 1600.0],
                'cn2': [0.5, 0.5],
            },
            'layer1': {'height': 1000.0},
            'layer2': {'height': 3000.0},
        }

        weights = compute_layer_weights_from_turbulence(
            params=params,
            component_indices=[1, 2],
            component_type='layer',
            verbose=False,
        )

        np.testing.assert_allclose(weights, np.array([0.5, 0.5]), rtol=0, atol=1e-12)

    def test_accepts_legacy_cn2_key_and_normalizes(self):
        params = {
            'main': {},
            'atmo': {
                'heights': [0.0, 1000.0],
                'Cn2': [2.0, 1.0],
            },
            'layer1': {'height': 0.0},
            'layer2': {'height': 1500.0},
        }

        weights = compute_layer_weights_from_turbulence(
            params=params,
            component_indices=[1, 2],
            component_type='layer',
            verbose=False,
        )

        np.testing.assert_allclose(weights.sum(), 1.0, rtol=0, atol=1e-12)
        np.testing.assert_allclose(weights, np.array([2.0 / 3.0, 1.0 / 3.0]), rtol=0, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
