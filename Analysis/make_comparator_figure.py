"""Descriptive model overlap, not an accuracy benchmark or factor alignment."""
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
from make_figures import ROOT, OUT, panel_grid, save, plt


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--summary',type=Path,required=True);args=ap.parse_args()
    summary=args.summary
    settings=summary/'analysis_settings.json'
    if not settings.exists():raise FileNotFoundError(settings)
    enrichment=pd.read_csv(summary/'candidate_factor_enrichment.tsv',sep='\t')
    d=pd.read_csv(summary/'GSM8443449_PDAC-p1_common_eligible.tsv.gz',sep='\t')
    print('Common map columns',list(d.columns),flush=True)
    run_paths={
        'STAMP':ROOT/'Analysis/comparators/results/stamp/k10_seed20261115_retry_numpy2',
        'POPARI':ROOT/'Analysis/comparators/results/popari/k10_seed20261115_budget50',
        'SpaHDmap':ROOT/'Analysis/comparators/results/spahdmap/k10_seed20261115_production',
    }
    candidates=pd.read_csv(ROOT/'Source_Data/fresh_coordinates.tsv.gz',sep='\t',
                           usecols=['sample_id','barcode','x','y','cadspot_candidate'])
    sample='GSM8443449_PDAC-p1'
    base=candidates[candidates.sample_id==sample].merge(d[['sample_id','barcode']],on=['sample_id','barcode'],validate='one_to_one')
    fig,a=panel_grid()
    methods=['STAMP','POPARI','SpaHDmap']
    for ax,method in zip(a[:3],methods):
        factors=pd.read_csv(run_paths[method]/sample/'spot_factors.tsv.gz',sep='\t')
        plotted=base.merge(factors[['sample_id','barcode','dominant_domain']],on=['sample_id','barcode'],validate='one_to_one')
        codes=plotted.dominant_domain.str.extract(r'(\d+)$')[0].astype(int)-1
        spatial_colors=ax.scatter(plotted.x,plotted.y,c=codes,cmap='tab10',vmin=-.5,vmax=9.5,s=3,linewidths=0,rasterized=True)
        hot=plotted[plotted.cadspot_candidate.astype(str).str.lower().isin(['true','1'])]
        ax.scatter(hot.x,hot.y,facecolors='none',edgecolors='black',s=12,linewidths=.4,rasterized=True)
        ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values():spine.set_visible(False)
        ax.set_xlabel(method)
    colorbar=fig.colorbar(spatial_colors,ax=list(a[:3]),fraction=.022,pad=.03,ticks=range(10),label='Within-method factor')
    colorbar.ax.set_yticklabels(range(1,11))
    primary=enrichment[enrichment.block_width_pitches==8].copy()
    primary['plot_effect']=np.log2((primary.candidate_overlap+.5)/(primary.expected_overlap+.5))
    primary.loc[primary.n_domain==0,'plot_effect']=np.nan
    bound=max(1,float(np.ceil(primary.plot_effect.abs().max())))
    all_samples=sorted(primary.sample_id.unique())
    for ax,method in zip(a[3:],methods):
        sub=primary[primary.method==method]
        samples=all_samples
        factors=[f'factor_{i:02d}' for i in range(1,11)]
        matrix=sub.pivot(index='sample_id',columns='factor',values='plot_effect').reindex(index=samples,columns=factors)
        q=sub.pivot(index='sample_id',columns='factor',values='q_block_domain').reindex(index=samples,columns=factors)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#E8E8E8')
        im=ax.imshow(matrix.to_numpy(float),vmin=-bound,vmax=bound,cmap=cmap,aspect='auto')
        for row in range(len(samples)):
            for col in range(10):
                if q.iloc[row,col]<.05:ax.plot(col,row,'o',color='black',ms=2)
        ax.set_yticks(range(len(samples)),[s.rsplit('_',1)[-1].replace('PDAC-','') for s in samples])
        ax.set_xticks([0,2,4,6,9],[1,3,5,7,10],fontsize=8);ax.set_xlabel('Factor / topic')
    fig.colorbar(im,ax=list(a[3:]),fraction=.022,pad=.03,label='log2 overlap enrichment')
    save(fig,'Figure5_comparators')
    print(OUT/'Figure5_comparators.pdf',flush=True)


if __name__=='__main__':main()
