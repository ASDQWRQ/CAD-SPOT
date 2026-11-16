"""Post hoc common-gene sensitivity, separate from declared cohort analyses."""
import pandas as pd
from scipy import sparse
from end_to_end import ROOT, TECH, COMP, tsv, json_out
import cadspot as core
import cadspot_null as null


def main():
    folders={'ovarian':ROOT/'Primary','pdac':ROOT/'PDAC_Transfer'}
    common={}
    for panel in ['target_score','reference_score']:
        sets=[]
        for folder in folders.values():
            c=pd.read_csv(folder/'panel_coverage.tsv',sep='\t')
            lists=[set(s.loc[s.status.eq('used'),'gene']) for _,s in c[c.panel.eq(panel)].groupby('sample_id')]
            sets.extend(lists)
        common[panel]=sorted(set.intersection(*sets))
        if len(common[panel])<3:raise ValueError('Insufficient common nonconstant panel coverage')
    rows=[]
    for cohort,folder in folders.items():
        spots=pd.read_csv(folder/'score_table.tsv.gz',sep='\t')
        genes=pd.read_csv(folder/'genes.tsv',sep='\t').gene.tolist()
        expression=null.ExpressionData(sparse.load_npz(folder/'logcp10k_all.npz').T.tocsr(),genes,spots)
        for panel,members in common.items():
            spots[panel]=null.score_gene_set(expression,members,scoring='mean_gene_z_by_section')
        associations=core.associations(spots,'target_score','reference_score',TECH,COMP)
        associations.insert(0,'cohort',cohort)
        associations['analysis']='post_hoc_common_target_reference_coverage'
        tsv(associations,folder/'Results/common_gene_sensitivity_associations.tsv')
        _,a,b=core.fit_pair(spots,'target_score','reference_score',TECH+COMP)
        ranked=core.spatial_ranking(spots,a,b)
        rows.append(dict(cohort=cohort,n_spots=len(spots),target_genes=len(common['target_score']),
            reference_genes=len(common['reference_score']),candidates=int(ranked.cadspot_candidate.sum()),
            composition_rho=float(associations.query("scope=='pooled' and layer=='composition_adjusted' and method=='score_residual_spearman'").spearman_rho.iloc[0])))
    output=ROOT/'Analysis/Cross_Cancer';output.mkdir(exist_ok=True)
    tsv(pd.DataFrame(rows),output/'common_gene_summary.tsv')
    json_out(common,output/'common_gene_membership.json')
    print(pd.DataFrame(rows).to_string(index=False),flush=True)


if __name__=='__main__':main()
