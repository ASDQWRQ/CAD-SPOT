"""Geometry normalization must not change scores or physical-coordinate exports."""
import unittest
import numpy as np
import pandas as pd
from end_to_end import pitch_scaled_coordinates


class PitchTests(unittest.TestCase):
    def test_section_specific_pitch_and_original_preserved(self):
        d = pd.DataFrame(dict(sample_id=['a']*4+['b']*4,
            x=[0., 10., 0., 10., 0., 30., 0., 30.],
            y=[0., 0., 10., 10., 0., 0., 30., 30.],
            target_score=np.arange(8, dtype=float)))
        original = d.copy()
        transformed, factors = pitch_scaled_coordinates(d)
        pd.testing.assert_frame_equal(d, original)
        np.testing.assert_allclose(factors.full_resolution_pixels_per_pitch, [10, 30])
        np.testing.assert_allclose(transformed.loc[:3, ['x', 'y']],
                                   transformed.loc[4:, ['x', 'y']].to_numpy())
        np.testing.assert_array_equal(transformed.target_score, d.target_score)

    def test_integer_pixel_columns_are_promoted_only_in_transformed_copy(self):
        d = pd.DataFrame(dict(sample_id=['a']*4, x=[0, 11, 0, 11], y=[0, 0, 11, 11]))
        t, _ = pitch_scaled_coordinates(d)
        self.assertEqual(d.x.dtype.kind, 'i')
        self.assertEqual(t.x.dtype.kind, 'f')
        np.testing.assert_allclose(t.x, [0, 1, 0, 1])

    def test_duplicate_only_geometry_fails(self):
        d = pd.DataFrame(dict(sample_id=['a']*3, x=0., y=0.))
        with self.assertRaises(ValueError):
            pitch_scaled_coordinates(d)


if __name__ == '__main__':
    unittest.main()
