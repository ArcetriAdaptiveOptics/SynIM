import unittest
import numpy as np
from synim.synim import compute_subaperture_illumination

class TestSAIllumination(unittest.TestCase):
    def test_sa_illumination_consistency(self):
        '''
        Test that the sub-aperture illumination computed by
        compute_subaperture_illumination matches an independent calculation
        for a simple case.
        '''
        pupil_size = 64
        pupil = np.zeros((pupil_size, pupil_size), dtype=np.float32)
        pupil[16:48, 20:44] = 1.0

        wfs_nsubaps = 8
        illum = compute_subaperture_illumination(
            pupil, wfs_nsubaps, wfs_rotation=0.0, wfs_translation=(0.0, 0.0),
            wfs_magnification=(1.0, 1.0), idx_valid_sa=None, verbose=False, specula_convention=False
        )

        # independently compute expected illumination
        rebinned = pupil.reshape(
            wfs_nsubaps, pupil_size//wfs_nsubaps, wfs_nsubaps, pupil_size//wfs_nsubaps
        )
        rebinned = rebinned.sum(axis=(1,3))
        rebinned /= rebinned.max()
        expected = rebinned.flatten()

        np.testing.assert_allclose(illum, expected, atol=1e-6)

    def test_sa_illumination_empty_pupil(self):
        '''
        Test that compute_subaperture_illumination raises a ValueError
        when the pupil mask is empty.
        '''
        pupil_size = 64
        pupil = np.zeros((pupil_size, pupil_size), dtype=np.float32)
        wfs_nsubaps = 8
        with self.assertRaises(ValueError):
            compute_subaperture_illumination(
                pupil,
                wfs_nsubaps,
                wfs_rotation=0.0,
                wfs_translation=(0.0, 0.0),
                wfs_magnification=(1.0, 1.0),
                idx_valid_sa=None,
                verbose=False,
                specula_convention=False
            )

    def test_sa_illumination_full_pupil(self):
        '''
        Test that compute_subaperture_illumination returns uniform illumination
        for a full pupil mask.
        '''
        pupil_size = 64
        pupil = np.ones((pupil_size, pupil_size), dtype=np.float32)
        wfs_nsubaps = 8
        illum = compute_subaperture_illumination(
            pupil, wfs_nsubaps, wfs_rotation=0.0, wfs_translation=(0.0, 0.0),
            wfs_magnification=(1.0, 1.0), idx_valid_sa=None, verbose=False, specula_convention=False
        )
        self.assertTrue(np.allclose(illum, illum[0]),
                        "Illumination not uniform for full pupil")
        self.assertTrue(np.all(illum > 0),
                        "Illumination zero for full pupil")

    def test_sa_illumination_single_pixel(self):
        '''
        Test that compute_subaperture_illumination correctly identifies
        a single illuminated sub-aperture when only one pixel in the pupil is set.
        '''
        pupil_size = 64
        pupil = np.zeros((pupil_size, pupil_size), dtype=np.float32)
        pupil[32, 32] = 1.0
        wfs_nsubaps = 8
        illum = compute_subaperture_illumination(
            pupil, wfs_nsubaps, wfs_rotation=0.0, wfs_translation=(0.0, 0.0),
            wfs_magnification=(1.0, 1.0), idx_valid_sa=None, verbose=False, specula_convention=False
        )
        self.assertTrue(np.count_nonzero(illum) == 1,
                        "Only one sub-aperture should be illuminated")
