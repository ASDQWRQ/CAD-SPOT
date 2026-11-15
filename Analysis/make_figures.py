"""Publication figures from computed revision outputs, with concise panel labels."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT/'Source_Data'
OUT = ROOT/'Additional_Figures'
BLUE, RED, GREEN, GRAY = '#0072B2','#D55E00','#009E73','#B8BEC4'
PALETTE = [BLUE,RED,GREEN,'#CC79A7','#56B4E9']
plt.rcParams.update({'font.family':'Arial','font.size':9.5,'axes.labelsize':9.5,
    'legend.fontsize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.frameon':False,'axes.spines.top':False,'axes.spines.right':False,
    'axes.linewidth':.7,'xtick.major.size':3,'ytick.major.size':3,'pdf.fonttype':42,'ps.fonttype':42,
    'savefig.facecolor':'white'})


def read(name):
    path=SRC/name
    if not path.exists() and path.with_name(path.name+'.gz').exists():path=path.with_name(path.name+'.gz')
    return pd.read_csv(path,sep='\t')


def panel_grid():
    fig,axes=plt.subplots(2,3,figsize=(8.0,5.47),layout='constrained')
    fig.set_constrained_layout_pads(w_pad=.10,h_pad=.10,wspace=.07,hspace=.10)
    for letter,ax in zip('ABCDEF',axes.ravel()):
        ax.text(-.12,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=13,ha='left',va='bottom',gid='panel_label',in_layout=False)
    return fig,axes.ravel()


def save(fig,name):
    OUT.mkdir(exist_ok=True)
    fig.canvas.draw()
    axes=fig.axes[:6]
    for i,ax in enumerate(axes):
        column=i%3;row=i//3
        left=min(axes[column].get_position().x0,axes[column+3].get_position().x0)
        top=max(axes[j].get_position().y1 for j in range(row*3,row*3+3))
        label=next(t for t in ax.texts if t.get_gid()=='panel_label')
        label.set_transform(fig.transFigure)
        label.set_position((max(.012,left-.035),min(.962,top+.018)))
    fig.savefig(OUT/f'{name}.pdf')
    fig.savefig(OUT/f'{name}.png',dpi=300)
    plt.close(fig)


def rho_label(ax,text):
    ax.text(.04,.96,text,transform=ax.transAxes,ha='left',va='top',fontsize=9,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=2))


def workflow(d):
    fig,a=panel_grid()
    a[0].set_xlim(0,1);a[0].set_ylim(0,1);a[0].axis('off')
    nodes=[(.05,.78,'Scores + covariates'),(.05,.51,'Residual models'),(.05,.24,'Six-neighbour ranking')]
    for x,y,text in nodes:
        a[0].add_patch(Rectangle((x,y),.90,.16,facecolor='#F0F3F5',edgecolor=GRAY,lw=.8))
        a[0].text(x+.45,y+.08,text,ha='center',va='center',fontsize=9)
    for y in [.76,.49]:
        a[0].annotate('',xy=(.5,y-.08),xytext=(.5,y),arrowprops=dict(arrowstyle='->',color='#555555',lw=.8))
    a[0].text(.5,.10,'Block CI  |  Random gene sets',ha='center',va='center',fontsize=8)
    available=[19,20];specified=[23,21]
    pos=np.arange(2)
    a[1].bar(pos,specified,color='#E5E8EB',width=.55,label='Specified')
    a[1].bar(pos,available,color=[BLUE,GREEN],width=.55,label='Available')
    a[1].set_xticks(pos,['Target','Reference']);a[1].set_ylabel('Genes');a[1].set_ylim(0,28)
    a[1].legend(loc='upper center',bbox_to_anchor=(.5,1.06),ncol=2)
    counts=d.groupby('sample_id',sort=True).size()
    a[2].barh(np.arange(5),counts,color=PALETTE,height=.58)
    a[2].set_yticks(range(5),[f'p{i}' for i in range(1,6)]);a[2].set_xlabel('Spots');a[2].invert_yaxis()
    for i,(sample,g) in enumerate(d.groupby('sample_id',sort=True)):
        rng=np.random.default_rng(0)
        sub=g.sample(min(800,len(g)),random_state=0)
        a[3].scatter(i+rng.uniform(-.18,.18,len(sub)),sub.target_score,s=1,alpha=.15,color=PALETTE[i],rasterized=True)
        a[3].plot([i-.22,i+.22],[g.target_score.median()]*2,color='#333333',lw=1)
    a[3].set_xticks(range(5),[f'p{i}' for i in range(1,6)]);a[3].set_ylabel('Input target score')
    ranks=np.column_stack([rankdata(d.target_score)/len(d),rankdata(d.reference_score)/len(d)])
    a[4].hexbin(ranks[:,0],ranks[:,1],gridsize=32,cmap='Blues',mincnt=1,rasterized=True)
    a[4].set(xlabel='Target rank',ylabel='Reference rank',xlim=(0,1),ylim=(0,1))
    cols=['target_score','reference_score','Epithelial_malignant_proxy','ECM_Stroma','Macrophage','T_cell']
    c=d[cols].corr(method='spearman').to_numpy()
    im=a[5].imshow(c,cmap='RdBu_r',vmin=-1,vmax=1)
    labels=['Target','Ref.','Epi.','ECM','Macro.','T cell']
    a[5].set_xticks(range(6),labels,rotation=55,ha='right');a[5].set_yticks(range(6),labels)
    fig.colorbar(im,ax=a[5],fraction=.046,pad=.03,label='Spearman rho')
    save(fig,'Figure1_workflow')


def decomposition(d):
    fig,a=panel_grid()
    assoc=read('fresh_associations.tsv')
    main=assoc[(assoc.scope=='pooled')&(assoc.method=='score_residual_spearman')]
    layer_order=['raw','technical_adjusted','composition_adjusted']
    vals=[main.loc[main.layer.eq(l),'spearman_rho'].iloc[0] for l in layer_order]
    a[0].plot(range(3),vals,'o-',color=BLUE,lw=1.4)
    a[0].set_xticks(range(3),['Raw','Technical','Composition'],rotation=25,ha='right')
    a[0].set(ylabel='Pooled Spearman rho',ylim=(0,.65))
    methods=['score_residual_spearman','partial_spearman','spline_score_residual_spearman']
    vals=[assoc[(assoc.scope=='pooled')&(assoc.layer=='composition_adjusted')&(assoc.method==m)].spearman_rho.iloc[0] for m in methods]
    a[1].scatter(vals,range(3),c=[BLUE,GREEN,RED],s=35)
    a[1].set_yticks(range(3),['Score residual','Partial Spearman','Spline residual']);a[1].set_xlim(.17,.22)
    a[1].set_xlabel('Pooled Spearman rho');a[1].invert_yaxis()
    for i,w in enumerate([1500,3000,4500]):
        ci=read(f'block_bootstrap_{w}_ci.tsv')
        row=ci[(ci.scope=='pooled')&(ci.layer=='composition_adjusted')].iloc[0]
        a[2].errorbar(row.spearman_rho,i,xerr=[[row.spearman_rho-row.ci_low],[row.ci_high-row.spearman_rho]],fmt='o',color=PALETTE[i],capsize=3)
    a[2].set_yticks(range(3),['1,500 px','3,000 px','4,500 px']);a[2].set_xlabel('Pooled rho (95% block CI)');a[2].set_xlim(.13,.25);a[2].invert_yaxis()
    ci=read('block_bootstrap_3000_ci.tsv')
    ci=ci[(ci.scope=='section')&(ci.layer=='composition_adjusted')].sort_values('section')
    for i,(_,row) in enumerate(ci.iterrows()):
        a[3].errorbar(row.spearman_rho,i,xerr=[[row.spearman_rho-row.ci_low],[row.ci_high-row.spearman_rho]],fmt='o',color=PALETTE[i],capsize=3)
    a[3].set_yticks(range(5),[f'p{i}' for i in range(1,6)]);a[3].set_xlabel('Section rho (95% block CI)');a[3].invert_yaxis();a[3].axvline(0,color=GRAY,lw=.8)
    a[4].hexbin(d.target_resid_z,d.reference_resid_z,gridsize=40,cmap='Blues',mincnt=1,rasterized=True)
    a[4].axhline(0,color=GRAY,lw=.6);a[4].axvline(0,color=GRAY,lw=.6)
    a[4].set(xlabel='Target residual Z',ylabel='Reference residual Z',xlim=(-3,5),ylim=(-3,5))
    for i,(s,g) in enumerate(d.groupby('sample_id',sort=True)):
        a[5].plot(np.sort(g.priority_robust_z),np.arange(1,len(g)+1)/len(g),color=PALETTE[i],lw=1.2,label=f'p{i+1}')
    a[5].axvline(3,color='#333333',ls='--',lw=.8);a[5].set(xlim=(-2,6),xlabel='Priority robust Z',ylabel='Cumulative fraction')
    a[5].legend(loc='lower right',ncol=2)
    save(fig,'Figure2_decomposition')


def candidate_maps(d):
    fig,a=panel_grid()
    manifest=json.loads((ROOT/'Data/GSE274103/manifest.json').read_text())
    for ax,m in zip(a[:5],manifest):
        g=d[d.sample_id.eq(m['sample_id'])]
        folder=ROOT/m['spatial_dir']
        sf=json.loads((folder/'scalefactors_json.json').read_text())
        image=np.asarray(Image.open(folder/'tissue_hires_image.png'))
        scale=sf['tissue_hires_scalef']
        ax.imshow(image,alpha=.80)
        selected=g[g.cadspot_candidate]
        ax.scatter(selected.x*scale,selected.y*scale,s=9,facecolors='none',edgecolors=RED,lw=.7,rasterized=True)
        ax.set_xlim(0,image.shape[1]);ax.set_ylim(image.shape[0],0);ax.set_aspect('equal');ax.axis('off')
    counts=d.groupby('sample_id',sort=True).cadspot_candidate.agg(['sum','size'])
    fractions=100*counts['sum']/counts['size']
    a[5].scatter(range(5),fractions,c=PALETTE,s=40)
    a[5].vlines(range(5),0,fractions,color=PALETTE,lw=1)
    a[5].set_xticks(range(5),[f'p{i}' for i in range(1,6)]);a[5].set_ylabel('Candidate spots (%)');a[5].set_ylim(0,6)
    a[5].legend(handles=[Line2D([0],[0],marker='o',mfc='none',mec=RED,ls='none',label='Candidate')],loc='upper right')
    save(fig,'Figure3_candidates')


def ranking(d):
    fig,a=panel_grid()
    comparison=read('fresh_ranking_baseline.tsv')
    eligible=d[(d.target_resid_z>0)&(d.reference_resid_z>0)]
    base=eligible.target_resid_z*eligible.reference_resid_z
    a[0].hexbin(rankdata(base)/len(base),rankdata(eligible.priority_score)/len(base),gridsize=30,cmap='Blues',mincnt=1,rasterized=True)
    a[0].plot([0,1],[0,1],color=GRAY,lw=.8,ls='--')
    a[0].set(xlabel='Product rank',ylabel='CAD-SPOT rank',xlim=(0,1),ylim=(0,1))
    subset=comparison[(comparison.scope=='section')&(comparison.comparator=='signed_residual_product')&(comparison.restriction=='positive_quadrant')&(comparison.selection=='top_k')].sort_values('section')
    a[1].scatter(subset.score_rank_spearman,range(5),c=PALETTE,s=35)
    a[1].set_yticks(range(5),[f'p{i}' for i in range(1,6)]);a[1].set(xlabel='Rank correlation',xlim=(0,1));a[1].invert_yaxis()
    for i,(s,sub) in enumerate(comparison[(comparison.scope=='section')&(comparison.comparator=='signed_residual_product')&(comparison.restriction=='positive_quadrant')&(comparison.selection=='top_fraction')].groupby('section',sort=True)):
        sub=sub.sort_values('top_fraction')
        a[2].plot(100*sub.top_fraction,100*sub.overlap_fraction,'o-',color=PALETTE[i],label=f'p{i+1}',lw=1.1,ms=3)
    a[2].set(xlabel='Top fraction (%)',ylabel='Overlap (%)',ylim=(0,100));a[2].legend(loc='lower right',ncol=2)
    signed=d.target_resid_z*d.reference_resid_z
    ix=signed.nlargest(int(np.ceil(.05*len(d)))).index
    bins=[int(((d.loc[ix,'target_resid_z']>0)&(d.loc[ix,'reference_resid_z']>0)).sum()),int(((d.loc[ix,'target_resid_z']<0)&(d.loc[ix,'reference_resid_z']<0)).sum())]
    a[3].bar(['Both positive','Both negative'],bins,color=[GREEN,GRAY],width=.55)
    a[3].set_ylabel('Top 5% signed-product spots');a[3].tick_params(axis='x',rotation=15)
    topn=max(1,int(np.ceil(.05*len(eligible))))
    full=set(eligible.priority_score.nlargest(topn).index);simple=set(base.nlargest(topn).index)
    groups=[sorted(full-simple),sorted(full&simple),sorted(simple-full)]
    values=[eligible.loc[ix,'neighbor_reference_resid_z'].to_numpy() for ix in groups]
    box=a[4].boxplot(values,tick_labels=['CAD only','Shared','Product only'],showfliers=False,patch_artist=True,widths=.5)
    for patch,color in zip(box['boxes'],[BLUE,GREEN,GRAY]):patch.set_facecolor(color);patch.set_alpha(.7)
    a[4].set_ylabel('Neighbour reference residual Z');a[4].tick_params(axis='x',rotation=20)
    components=d[['target_resid_z','reference_resid_z','neighbor_reference_resid_z','local_support','priority_score']]
    matrix=components.corr(method='spearman')
    im=a[5].imshow(matrix,cmap='RdBu_r',vmin=-1,vmax=1)
    labs=['Target','Ref.','Neighbour','Product','Priority']
    a[5].set_xticks(range(5),labs,rotation=60,ha='right');a[5].set_yticks(range(5),labs)
    fig.colorbar(im,ax=a[5],fraction=.046,pad=.03,label='Spearman rho')
    save(fig,'Figure4_ranking')


def raw_audit():
    fig,a=panel_grid()
    null=read('fresh_matched_null_500_rhos.tsv')
    summary=read('fresh_matched_null_500_summary.tsv').iloc[0]
    a[0].hist(null.rho,bins=24,color=GRAY,edgecolor='white')
    a[0].axvline(summary.observed_rho,color=RED,lw=1.5,label='Observed')
    a[0].axvline(summary.null_95th_percentile,color=BLUE,lw=1,ls='--',label='Null 95th')
    a[0].set(xlabel='Raw-matrix adjusted rho',ylabel='Random gene sets')
    a[0].legend(loc='upper center',bbox_to_anchor=(.5,1.10),ncol=2)
    matched=read('fresh_matched_sets_500.tsv')
    for ax,col,color in zip(a[1:4],['delta_log1p_mean','delta_detection','delta_log_length'],[BLUE,GREEN,RED]):
        ax.hist(matched[col],bins=25,color=color,alpha=.8,edgecolor='white')
        ax.axvline(0,color='#333333',lw=.7);ax.set_ylabel('Matched gene slots')
    a[1].set_xlabel('Difference: log1p mean');a[2].set_xlabel('Detection-fraction difference');a[3].set_xlabel('Exon-union length difference\n(log units)')
    scores=read('raw_matrix_pair_scores.tsv')
    for ax,label,color in [(a[4],'target',BLUE),(a[5],'reference',GREEN)]:
        ax.hexbin(scores[f'{label}_score'],scores[f'fresh_{label}_score'],gridsize=30,cmap='Blues' if label=='target' else 'Greens',mincnt=1,rasterized=True)
        ax.set_xlabel(f'Archived input {label}');ax.set_ylabel(f'Reconstructed {label}')
    save(fig,'FigureS1_raw_matrix_audit')


def external():
    path=SRC/'external_reanalysis_summary.tsv'
    if not path.exists():return
    d=read(path.name)
    if len(d)!=10:return
    fig,a=panel_grid()
    labels=[v.replace('GSM','').split('_',1)[1].replace('_Vis','').replace('NYU_','') for v in d.sample_id]
    for i,row in d.iterrows():
        a[0].plot([row.null_95th_percentile,row.observed_rho],[i,i],color=GRAY,lw=1)
        a[0].scatter(row.null_95th_percentile,i,c=BLUE,marker='|',s=60)
        a[0].scatter(row.observed_rho,i,c=RED if row.observed_rho>row.null_95th_percentile else GRAY,s=22)
    a[0].set_yticks(range(10),labels,fontsize=8);a[0].set_xlabel('Technical-adjusted rho');a[0].invert_yaxis()
    coverage=d[['n_target_available','n_reference_available']].to_numpy()
    im=a[1].imshow(coverage,aspect='auto',cmap='Blues',vmin=15,vmax=23)
    a[1].set_xticks([0,1],['Target','Reference']);a[1].set_yticks(range(10),labels,fontsize=8)
    for i in range(10):
        for j in range(2):a[1].text(j,i,str(coverage[i,j]),ha='center',va='center',fontsize=8.5,color='white' if coverage[i,j]>20 else 'black')
    a[2].scatter(d.n_spots,d.observed_rho,color=BLUE,s=28)
    a[2].set(xlabel='Spots',ylabel='Technical-adjusted rho')
    a[3].scatter(d.n_universe,d.null_95th_percentile,color=GREEN,s=28)
    a[3].set(xlabel='Matching-pool genes',ylabel='Matched-null 95th percentile')
    for i,w in enumerate([1500,3000,4500]):
        b=read(f'block_bootstrap_{w}_draws.tsv')
        sub=b[(b.scope=='pooled')&(b.layer=='composition_adjusted')]
        a[4].hist(sub.spearman_rho,bins=30,histtype='step',color=PALETTE[i],label=f'{w:,} px',lw=1)
    a[4].set(xlabel='Bootstrap pooled rho',ylabel='Resamples');a[4].legend(loc='upper left')
    old=read('archived_candidate_comparison.tsv').iloc[0]
    a[5].bar(['Shared','Old only','New only'],[old.overlap,old.archived_only,old.fresh_only],color=[GREEN,GRAY,BLUE],width=.55)
    a[5].tick_params(axis='x',rotation=25);a[5].set_ylabel('Coordinates')
    save(fig,'FigureS2_external_and_provenance')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--include-comparators',action='store_true');ap.parse_args()
    d=read('fresh_coordinates.tsv.gz')
    workflow(d);decomposition(d);candidate_maps(d);ranking(d);raw_audit();external()
