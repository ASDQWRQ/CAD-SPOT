"""Contextual overlap from preserved fits and current candidate labels."""
import numpy as np
import pandas as pd

from make_figures import ROOT, grid, save, plt


def main():
    summary=ROOT/'Comparators/current_summary';sample='GSM8443449_PDAC-p1'
    d=pd.read_csv(summary/f'{sample}_common_eligible.tsv.gz',sep='\t')
    candidates=pd.read_csv(ROOT/'PDAC_Transfer/Results/composition_coordinates.tsv.gz',sep='\t',
        usecols=['sample_id','barcode','x','y','cadspot_candidate'])
    base=candidates[candidates.sample_id.eq(sample)].merge(d[['sample_id','barcode']],
        on=['sample_id','barcode'],validate='one_to_one')
    runs={'STAMP':'stamp/k10_seed20261115_retry_numpy2',
          'POPARI':'popari/k10_seed20261115_budget50',
          'SpaHDmap':'spahdmap/k10_seed20261115_production'}
    fig,axes=grid()
    for ax,(method,run) in zip(axes[:3],runs.items()):
        factors=pd.read_csv(ROOT/f'Analysis/comparators/results/{run}/{sample}/spot_factors.tsv.gz',sep='\t')
        plotted=base.merge(factors[['sample_id','barcode','dominant_domain']],
            on=['sample_id','barcode'],validate='one_to_one')
        codes=plotted.dominant_domain.str.extract(r'(\d+)$')[0].astype(int)-1
        artist=ax.scatter(plotted.x,plotted.y,c=codes,cmap='tab10',vmin=-.5,vmax=9.5,
            s=3,linewidths=0,rasterized=True)
        hot=plotted[plotted.cadspot_candidate]
        ax.scatter(hot.x,hot.y,facecolors='none',edgecolors='black',s=12,linewidths=.4,rasterized=True)
        ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values():spine.set_visible(False)
        ax.set_xlabel(method)
    cb=fig.colorbar(artist,ax=list(axes[:3]),fraction=.022,pad=.03,ticks=range(10),label='Factor / topic')
    cb.ax.set_yticklabels(range(1,11))
    enrichment=pd.read_csv(summary/'candidate_factor_enrichment.tsv',sep='\t')
    primary=enrichment[enrichment.block_width_pitches.eq(8)].copy()
    primary['plot_effect']=np.log2((primary.candidate_overlap+.5)/(primary.expected_overlap+.5))
    primary.loc[primary.n_domain.eq(0),'plot_effect']=np.nan
    bound=max(1,float(np.ceil(primary.plot_effect.abs().max())))
    samples=sorted(primary.sample_id.unique());factors=[f'factor_{i:02d}' for i in range(1,11)]
    for ax,method in zip(axes[3:],runs):
        sub=primary[primary.method.eq(method)]
        values=sub.pivot(index='sample_id',columns='factor',values='plot_effect').reindex(index=samples,columns=factors)
        q=sub.pivot(index='sample_id',columns='factor',values='q_block_domain').reindex(index=samples,columns=factors)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#E8E8E8')
        im=ax.imshow(values.to_numpy(float),vmin=-bound,vmax=bound,cmap=cmap,aspect='auto')
        for row in range(len(samples)):
            for col in range(10):
                if q.iloc[row,col]<.05:ax.plot(col,row,'o',color='black',ms=2)
        ax.set_yticks(range(len(samples)),[f'p{i}' for i in range(1,6)])
        ax.set_xticks([0,2,4,6,9],[1,3,5,7,10]);ax.set_xlabel('Factor / topic')
    fig.colorbar(im,ax=list(axes[3:]),fraction=.022,pad=.03,label='log2 overlap enrichment')
    save(fig,axes,'FigureS2_contextual_comparisons')


if __name__=='__main__':main()
