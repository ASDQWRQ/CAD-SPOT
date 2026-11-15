"""Focused tests for the newly public count-to-score stage."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from end_to_end import PANELS, prepare_section


class PreparationTests(unittest.TestCase):
    def fixture(self, missing=()):
        genes=pd.Index(sorted((set(g for panel in PANELS.values() for g in panel) |
                               {f'BACKGROUND{i:03}' for i in range(150)})-set(missing)))
        rng=np.random.default_rng(72)
        counts=sparse.csr_matrix(rng.poisson(4,size=(12,len(genes)))+1,dtype=float)
        pos=pd.DataFrame(dict(barcode=[f'b{i}' for i in range(12)],in_tissue=1,
                              x=np.arange(12,dtype=float),y=np.zeros(12)))
        return counts,genes,pos

    def test_score_parameters_reconstruct_every_panel(self):
        counts,genes,pos=self.fixture()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            values,spots,_,_=prepare_section(counts,genes,pos,'s1',out)
            parameters=pd.read_csv(out/'score_parameters.tsv',sep='\t')
            for panel,sub in parameters.groupby('panel'):
                idx=genes.get_indexer(sub.gene)
                z=(values[:,idx].toarray()-sub.mean_logcp10k.to_numpy())/sub.sd_logcp10k.to_numpy()
                np.testing.assert_allclose(z@sub.weight.to_numpy(),spots[panel],atol=1e-12)

    def test_missing_genes_are_omitted_not_zero_imputed(self):
        counts,genes,pos=self.fixture(missing=['ALDOA','GAPDH'])
        with tempfile.TemporaryDirectory() as tmp:
            _,_,_,coverage=prepare_section(counts,genes,pos,'s1',Path(tmp))
            coverage=pd.DataFrame(coverage)
            self.assertEqual(coverage.query("panel=='target_score' and status=='used'").shape[0],22)
            self.assertEqual(coverage.query("panel=='reference_score' and status=='used'").shape[0],20)
            self.assertTrue(coverage[coverage.gene.isin(['ALDOA','GAPDH'])].status.eq('missing').all())

    def test_locked_qc_and_axis_alignment(self):
        counts,genes,pos=self.fixture()
        pos.loc[0,'in_tissue']=0
        pos.loc[1,'x']=np.nan
        counts=counts.tolil();counts[2,:]=0;counts=counts.tocsr();counts.eliminate_zeros()
        with tempfile.TemporaryDirectory() as tmp:
            values,spots,_,_=prepare_section(counts,genes,pos,'s1',Path(tmp))
            self.assertEqual(spots.barcode.tolist(),[f'b{i}' for i in range(3,12)])
            self.assertEqual(values.shape[0],9)

    def test_detection_fraction_uses_integer_count_denominator(self):
        counts,genes,pos=self.fixture()
        with tempfile.TemporaryDirectory() as tmp:
            _,_,meta,_=prepare_section(counts,genes,pos,'s1',Path(tmp))
            self.assertTrue(meta.detection_rate.between(0,1).all())
            np.testing.assert_array_equal(meta.detection_rate,np.ones(len(genes)))

    def test_proxy_panels_are_disjoint_from_scores(self):
        used=set(PANELS['target_score'])|set(PANELS['reference_score'])
        for key,panel in PANELS.items():
            if key not in ['target_score','reference_score']:
                self.assertFalse(used&set(panel))


if __name__=='__main__':unittest.main()
