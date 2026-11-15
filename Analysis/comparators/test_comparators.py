"""Focused tests of joins, label parsing, graph and block randomization; not model surrogates."""
import unittest
import common
import numpy as np
import pandas as pd
from enrichment import block_groups, permute_labels, section_stats


class ComparatorTests(unittest.TestCase):
    def test_boolean_not_string_truthiness(self):
        np.testing.assert_array_equal(common.bool_labels(pd.Series(['False', 'True', '0', '1'])), [False, True, False, True])
        with self.assertRaises(ValueError):
            common.bool_labels(pd.Series(['unknown']))

    def test_first_ring_graph_does_not_bridge_gaps(self):
        coords = np.array([[0, 0], [1, 0], [0, 1], [100, 100], [101, 100], [100, 101]])
        g = common.graph(coords)
        self.assertEqual(g[:3, 3:].nnz, 0)
        self.assertEqual((g != g.T).nnz, 0)
        self.assertEqual(g.diagonal().sum(), 0)

    def test_block_permutation_preserves_counts_and_block_patterns(self):
        coords = np.array([(x, y) for y in range(24) for x in range(24)])
        groups, meta = block_groups(coords, 4)
        y = np.arange(len(coords)) % 7 == 0
        permuted = permute_labels(y, groups, np.random.default_rng(42))
        self.assertEqual(y.sum(), permuted.sum())
        self.assertGreater(meta['movable_fraction'], .9)
        for group in groups:
            before = sorted(tuple(row) for row in y[group])
            after = sorted(tuple(row) for row in permuted[group])
            self.assertEqual(before, after)

    def test_no_exchangeability_means_no_p_value(self):
        f = pd.DataFrame({'pxl_col_in_fullres': range(10), 'pxl_row_in_fullres': 0,
                          'cadspot_candidate': [True]*3+[False]*7,
                          'factor_01': np.linspace(.1, .9, 10), 'factor_02': np.linspace(.9, .1, 10)})
        result = section_stats(f, ['factor_01', 'factor_02'], 8, 9, 1)
        self.assertTrue(result.p_block_domain.isna().all())


if __name__ == '__main__':
    unittest.main()
