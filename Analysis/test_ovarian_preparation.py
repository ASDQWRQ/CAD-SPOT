"""Focused tests for Matrix Market and keyed Visium adapters."""
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from scipy import io, sparse
from prepare_ovarian import read_market, read_positions, harmonize_scores


class MarketTests(unittest.TestCase):
    def test_cohort_membership_is_coverage_only_and_exports_exclusion(self):
        from test_end_to_end import PreparationTests
        from end_to_end import prepare_section
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            for name,missing in [('one',()),('two',('ALDOA',))]:
                x,genes,pos=PreparationTests().fixture(missing=missing)
                prepare_section(x,genes,pos,name,folder/name)
            frames,coverage=harmonize_scores(folder,['one','two'])
            status=coverage.query("sample_id=='one' and panel=='target_score' and gene=='ALDOA'").status.iloc[0]
            self.assertEqual(status,'excluded_not_common_across_sections')
            one=pd.read_csv(folder/'one/score_parameters.tsv',sep='\t')
            two=pd.read_csv(folder/'two/score_parameters.tsv',sep='\t')
            for panel in one.panel.unique():
                left=one[(one.panel==panel)&one.status.eq('used')]
                right=two[(two.panel==panel)&two.status.eq('used')]
                self.assertEqual(left.gene.tolist(),right.gene.tolist())
                self.assertAlmostEqual(left.weight.sum(),1.)

    def test_duplicate_symbols_and_non_gene_features(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)
            pd.DataFrame([['id1','A','Gene Expression'],['id2','A','Gene Expression'],
                          ['id3','B','Gene Expression'],['id4','Ab','Antibody Capture']]).to_csv(path/'features.tsv',sep='\t',index=False,header=False)
            pd.DataFrame(['b2','b1']).to_csv(path/'barcodes.tsv',sep='\t',index=False,header=False)
            io.mmwrite(path/'matrix.mtx',sparse.csr_matrix([[1,2],[3,4],[5,6],[7,8]]))
            x,genes,ids=read_market(path/'matrix.mtx',path/'features.tsv',path/'barcodes.tsv')
            self.assertEqual(genes.tolist(),['A','B'])
            self.assertEqual(ids.tolist(),['b2','b1'])
            np.testing.assert_array_equal(x.toarray(),[[4,5],[6,6]])

    def test_keyed_positions_preserve_expression_order(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'positions.csv'
            pd.DataFrame([['b1',1,0,0,10,20],['b2',1,0,1,30,40]]).to_csv(path,index=False,header=False)
            d=read_positions(path,['b2','b1'])
            self.assertEqual(d.barcode.tolist(),['b2','b1'])
            self.assertEqual(d.x.tolist(),[40,20])
            with self.assertRaises(ValueError):read_positions(path,['missing'])

    def test_dimension_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)
            pd.DataFrame([['id1','A']]).to_csv(path/'features.tsv',sep='\t',index=False,header=False)
            pd.DataFrame(['b1']).to_csv(path/'barcodes.tsv',sep='\t',index=False,header=False)
            io.mmwrite(path/'matrix.mtx',sparse.csr_matrix([[1,2]]))
            with self.assertRaises(ValueError):read_market(path/'matrix.mtx',path/'features.tsv',path/'barcodes.tsv')


if __name__ == '__main__':unittest.main()
