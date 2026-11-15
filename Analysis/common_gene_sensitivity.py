"""Post hoc coverage check; does not replace the prespecified external run."""
from pathlib import Path
import sys

import pandas as pd
from scipy import sparse

from end_to_end import ROOT, TECH, COMP
import cadspot as core
import cadspot_null as null


def main():
    folder=ROOT/'External_EndToEnd/NYU_PDAC1'
    spots=pd.read_csv(folder/'score_table.tsv.gz',sep='\t')
    genes=pd.read_csv(folder/'genes.tsv',sep='\t').gene.tolist()
    values=sparse.load_npz(folder/'logcp10k.npz').T.tocsr()
    expression=null.ExpressionData(values,genes,spots)
    coverage=pd.read_csv(ROOT/'Primary/panel_coverage.tsv',sep='\t')
    for panel in ['target_score','reference_score']:
        selected=coverage[(coverage.panel==panel)&(coverage.status=='used')]
        lists=[tuple(g.gene) for _,g in selected.groupby('sample_id')]
        assert len(set(lists))==1
        spots[panel]=null.score_gene_set(expression,list(lists[0]),scoring='mean_gene_z_by_section')
    assoc=core.associations(spots,'target_score','reference_score',TECH,COMP)
    assoc['analysis']='post_hoc_primary_available_gene_coverage_external_sensitivity'
    assoc.to_csv(folder/'Results/common_gene_sensitivity_associations.tsv',sep='\t',index=False)
    _,a,b=core.fit_pair(spots,'target_score','reference_score',TECH+COMP)
    ranked=core.spatial_ranking(spots,a,b)
    pd.DataFrame([dict(n_spots=len(spots),target_genes=19,reference_genes=20,
        candidates=int(ranked.cadspot_candidate.sum()),
        interpretation='Coverage sensitivity only; primary external analysis remains full available panels')]).to_csv(
        folder/'Results/common_gene_sensitivity_counts.tsv',sep='\t',index=False)
    print(assoc.query("scope=='pooled'").to_string(index=False),flush=True)


if __name__=='__main__':main()
