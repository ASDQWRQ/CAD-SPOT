"""Compact six-panel figures from the end-to-end follow-up outputs."""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'Additional_Figures'
BLUE, ORANGE, GREEN = '#0072B2', '#D55E00', '#009E73'
COLORS = [BLUE, ORANGE, GREEN, '#CC79A7', '#555555']
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,
    'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.frameon':False,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
    'pdf.fonttype':42,'savefig.facecolor':'white'})


def read(folder, name):
    return pd.read_csv(folder/name, sep='\t', keep_default_na=False, na_values=[''])


def grid():
    fig, axes = plt.subplots(2,3,figsize=(8.4,5.8),layout='constrained')
    fig.set_constrained_layout_pads(w_pad=.09,h_pad=.12,wspace=.10,hspace=.14)
    fig.get_layout_engine().set(rect=(.01,.01,.98,.94))
    return fig, axes.ravel()


def save(fig, axes, name):
    OUT.mkdir(exist_ok=True)
    fig.canvas.draw()
    for i, ax in enumerate(axes):
        row, col = divmod(i,3)
        left = min(axes[col].get_position().x0,axes[col+3].get_position().x0)
        top = max(axes[j].get_position().y1 for j in range(row*3,row*3+3))
        fig.text(max(.005,left-.038),min(.975,top+.022),'ABCDEF'[i],fontsize=12,fontweight='bold')
    fig.savefig(OUT/f'{name}.pdf')
    fig.savefig(OUT/f'{name}.png',dpi=300)
    plt.close(fig)


def association_rows(folder):
    d = read(folder/'Results','associations.tsv')
    return d[(d.scope=='pooled') & (d.method=='score_residual_spearman')]


def map_on_image(ax, spots, image, scale, values=None, candidates=False):
    im = Image.open(image).convert('RGB')
    ax.imshow(im)
    x, y = spots.x.to_numpy()*scale, spots.y.to_numpy()*scale
    if candidates:
        chosen = spots.cadspot_candidate.astype(bool).to_numpy()
        ax.scatter(x[chosen],y[chosen],s=10,facecolors='none',edgecolors=ORANGE,linewidths=.6)
    else:
        ax.scatter(x,y,c=values,s=5,cmap='viridis',linewidths=0,alpha=.8,rasterized=True)
    ax.set_xlim(max(0,x.min()-50),min(im.width,x.max()+50))
    ax.set_ylim(min(im.height,y.max()+50),max(0,y.min()-50))
    ax.set_aspect('equal'); ax.axis('off')


def primary_figures():
    folder = ROOT/'Primary'
    d = read(folder,'score_table.tsv.gz')
    r = read(folder/'Results','composition_coordinates.tsv.gz')
    groups = list(d.groupby('sample_id',sort=True))
    fig,a = grid()
    a[0].set(xlim=(0,1),ylim=(0,1)); a[0].axis('off')
    for y,text in [(.82,'Counts + coordinates'),(.56,'Scores + proxies'),(.30,'Residuals + priorities')]:
        a[0].add_patch(Rectangle((.03,y-.08),.94,.16,facecolor='#EEF2F4',edgecolor='#9FA8AE',lw=.7))
        a[0].text(.5,y,text,ha='center',va='center',fontsize=9)
    for y in [.70,.44]:
        a[0].annotate('',xy=(.5,y-.08),xytext=(.5,y),arrowprops=dict(arrowstyle='->',color='#555555'))
    a[0].text(.5,.10,'Block CI | Matched sets',ha='center',fontsize=8)
    a[1].bar([0,1],[23,21],color='#DFE4E6',width=.55,label='Intended')
    a[1].bar([0,1],[19,20],color=[BLUE,GREEN],width=.55,label='Available')
    a[1].set(xticks=[0,1],xticklabels=['Target','Reference'],ylabel='Genes',ylim=(0,28))
    a[1].legend(loc='upper center',bbox_to_anchor=(.5,1.14),ncol=2)
    a[2].barh(range(5),[len(g) for _,g in groups],color=COLORS,height=.55)
    a[2].set(yticks=range(5),yticklabels=[f'p{i}' for i in range(1,6)],xlabel='Spots');a[2].invert_yaxis()
    a[3].hexbin(d.target_score,d.reference_score,gridsize=34,mincnt=1,cmap='Blues',rasterized=True)
    a[3].set(xlabel='Target score',ylabel='Reference score')
    panels=['Epithelial_proxy','ECM_proxy','Contractile_proxy','Myeloid_proxy','Tcell_proxy']
    for i,col in enumerate(panels):
        vals=d[col].to_numpy()
        a[4].boxplot([vals],positions=[i],widths=.5,showfliers=False,patch_artist=True,
                     boxprops=dict(facecolor=COLORS[i],alpha=.6),medianprops=dict(color='black'))
    a[4].set(xticks=range(5),xticklabels=['Epi','ECM','Contract.','Myeloid','T cell'],ylabel='Proxy score')
    a[4].tick_params(axis='x',rotation=35)
    corr=d[['target_score','reference_score']+panels].corr(method='spearman')
    im=a[5].imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r')
    labels=['Target','Ref.','Epi','ECM','Contr.','Myeloid','T cell']
    a[5].set(xticks=range(7),yticks=range(7),xticklabels=labels,yticklabels=labels)
    a[5].tick_params(axis='x',rotation=60)
    fig.colorbar(im,ax=a[5],shrink=.75,pad=.03)
    save(fig,a,'Figure1_end_to_end')

    fig,a=grid()
    points=association_rows(folder)
    a[0].plot(range(3),points.spearman_rho,'o-',color=BLUE,lw=1)
    a[0].set(xticks=range(3),xticklabels=['Raw','Technical','Composition'],ylabel='Spearman rho',ylim=(0,.6))
    a[0].tick_params(axis='x',rotation=25)
    sensitivity=read(folder/'Results','associations.tsv').query("scope=='pooled' and layer=='composition_adjusted'")
    a[1].scatter(range(3),sensitivity.spearman_rho,c=[BLUE,GREEN,ORANGE],s=32)
    a[1].set(xticks=range(3),xticklabels=['OLS','Partial rank','Spline'],ylabel='Adjusted rho',ylim=(.25,.40))
    a[1].tick_params(axis='x',rotation=25)
    for i,width in enumerate([1500,3000,4500]):
        ci=read(folder/'Results',f'bootstrap_{width}_ci.tsv').query("scope=='pooled' and layer=='composition_adjusted'").iloc[0]
        a[2].errorbar(i,ci.spearman_rho,yerr=[[ci.spearman_rho-ci.ci_low],[ci.ci_high-ci.spearman_rho]],fmt='o',color=BLUE,capsize=3)
    a[2].set(xticks=range(3),xticklabels=['1,500','3,000','4,500'],xlabel='Block width (pixels)',ylabel='Adjusted rho',ylim=(.25,.40))
    ci=read(folder/'Results','bootstrap_3000_ci.tsv').query("scope=='section' and layer=='composition_adjusted'").reset_index(drop=True)
    for i,row in ci.iterrows():
        a[3].errorbar(row.spearman_rho,i,xerr=[[row.spearman_rho-row.ci_low],[row.ci_high-row.spearman_rho]],fmt='o',color=COLORS[i],capsize=2)
    a[3].set(yticks=range(5),yticklabels=[f'p{i}' for i in range(1,6)],xlabel='Adjusted rho');a[3].invert_yaxis()
    null=read(folder/'Results','matched_rhos.tsv')
    summary=read(folder/'Results','matched_summary.tsv').iloc[0]
    a[4].hist(null.rho,bins=25,color='#B5CDD9',edgecolor='white')
    a[4].axvline(summary.observed_rho,color=ORANGE,lw=1.5,label='Observed')
    a[4].axvline(summary.null_95th_percentile,color='#555555',ls='--',lw=1,label='Null 95%')
    a[4].set(xlabel='Adjusted rho',ylabel='Matched sets')
    a[4].legend(loc='upper center',bbox_to_anchor=(.5,1.17),ncol=2,columnspacing=.7,handlelength=1.2)
    counts=read(folder/'Results','candidate_counts.tsv').pivot(index='sample_id',columns='layer',values='n_candidates')
    for i,(_,row) in enumerate(counts.iterrows()):
        a[5].plot([i,i],[row.technical,row.composition],color='#9FA8AE',lw=1)
    a[5].scatter(range(5),counts.technical,color=BLUE,label='Technical',s=25)
    a[5].scatter(range(5),counts.composition,color=ORANGE,label='Composition',s=25,marker='s')
    a[5].set(xticks=range(5),xticklabels=[f'p{i}' for i in range(1,6)],ylabel='Candidate spots',ylim=(70,195))
    a[5].legend(loc='upper center',bbox_to_anchor=(.5,1.18),ncol=1)
    save(fig,a,'Figure2_associations_null')

    fig,a=grid()
    manifest=json.loads((folder/'input_manifest.json').read_text())
    for ax,source in zip(a[:5],manifest):
        sample=r[r.sample_id.eq(source['sample_id'])]
        spatial=Path(source['spatial_dir'])
        scale=json.loads((spatial/'scalefactors_json.json').read_text())['tissue_hires_scalef']
        map_on_image(ax,sample,spatial/'tissue_hires_image.png',scale,candidates=True)
    technical=read(folder/'Results','technical_coordinates.tsv.gz')[['sample_id','barcode','cadspot_candidate']]
    both=r[['sample_id','barcode','cadspot_candidate']].merge(technical,on=['sample_id','barcode'],suffixes=('_comp','_tech'),validate='one_to_one')
    matrix=np.zeros((2,2),int)
    for i in [0,1]:
        for j in [0,1]:
            matrix[i,j]=int(((both.cadspot_candidate_comp==bool(i))&(both.cadspot_candidate_tech==bool(j))).sum())
    a[5].imshow(np.log1p(matrix),cmap='Blues')
    for i in [0,1]:
        for j in [0,1]:a[5].text(j,i,f'{matrix[i,j]:,}',ha='center',va='center',color='white' if matrix[i,j]>1000 else 'black')
    a[5].set(xticks=[0,1],yticks=[0,1],xticklabels=['No','Yes'],yticklabels=['No','Yes'],xlabel='Technical candidate',ylabel='Composition candidate')
    save(fig,a,'Figure3_primary_coordinates')

    fig,a=grid()
    sim=read(ROOT/'Ranking_Ablation','all_replicates.tsv')
    for i,geometry in enumerate(['compact','dispersed']):
        ex=read(ROOT/'Ranking_Ablation',f'example_{geometry}.tsv.gz').query("sample_id=='s0'")
        a[i].scatter(ex.x,ex.y,c=np.where(ex.truth,GREEN,'#DFE4E6'),s=3,linewidths=0)
        a[i].set(xlabel='Array column',ylabel='Array row',aspect='equal');a[i].invert_yaxis()
    methodlabels=['Product','Legacy','Support']
    for i,method in enumerate(['positive_product','legacy','positive_support']):
        for j,geometry in enumerate(['compact','dispersed']):
            vals=sim.query('method==@method and geometry==@geometry').average_precision
            a[2].boxplot([vals],positions=[j+(i-1)*.22],widths=.18,showfliers=False,patch_artist=True,
                boxprops=dict(facecolor=COLORS[i],alpha=.7),medianprops=dict(color='black'))
        a[2].plot([],[],color=COLORS[i],label=methodlabels[i])
    a[2].set(xticks=[0,1],xticklabels=['Compact','Dispersed'],ylabel='Average precision',ylim=(.45,1.02))
    a[2].legend(loc='upper center',bbox_to_anchor=(.5,1.20),ncol=3,columnspacing=.6,handlelength=.9)
    for i,method in enumerate(['positive_product','legacy','positive_support']):
        means=[];sds=[]
        for geometry in ['compact','dispersed']:
            vals=sim.query('method==@method and geometry==@geometry').top_5pct_precision
            means.append(vals.mean());sds.append(vals.std())
        a[3].errorbar(np.arange(2)+(i-1)*.12,means,yerr=sds,fmt='o',color=COLORS[i],capsize=2)
    a[3].set(xticks=[0,1],xticklabels=['Compact','Dispersed'],ylabel='Top-5% precision',ylim=(.6,1.04))
    for i,method in enumerate(['positive_product','legacy','positive_support']):
        vals=100*sim.query("method==@method and geometry=='null'").candidate_fraction
        a[4].boxplot([vals],positions=[i],widths=.5,showfliers=False,patch_artist=True,
            boxprops=dict(facecolor=COLORS[i],alpha=.7),medianprops=dict(color='black'))
    a[4].set(xticks=range(3),xticklabels=methodlabels,ylabel='Null candidate spots (%)',ylim=(0,5))
    positive=r[(r.target_resid_z>0)&(r.reference_resid_z>0)].copy()
    fractions=[.01,.05,.10]
    for col,color,label in [('legacy_priority_score',ORANGE,'Legacy'),('local_support',GREEN,'Local support')]:
        overlaps=[]
        for f in fractions:
            k=int(np.ceil(f*len(positive)))
            first=set(positive.nlargest(k,'priority_score').index)
            other=set(positive.nlargest(k,col).index)
            overlaps.append(len(first&other)/k)
        a[5].plot(np.array(fractions)*100,overlaps,'o-',color=color,label=label)
    a[5].set(xlabel='Top-ranked fraction (%)',ylabel='Overlap with product',ylim=(0,1));a[5].legend(loc='lower right')
    save(fig,a,'FigureS1_neighborhood_ablation')


def external_figure():
    folder=ROOT/'External_EndToEnd/NYU_PDAC1'
    d=read(folder/'Results','composition_coordinates.tsv.gz')
    spatial=ROOT/'External/spatial/GSM6177618_NYU_PDAC1_Vis'
    scale=json.loads((spatial/'scalefactors_json.json').read_text())['tissue_hires_scalef']
    fig,a=grid()
    for ax in a[:4]:
        map_on_image(ax,d,spatial/'tissue_hires_image.png',scale,candidates=True)
    # The first view is unmarked H&E; other views are scored or nominated spots.
    a[0].collections[-1].remove()
    a[1].collections[-1].remove()
    a[2].collections[-1].remove()
    for ax,col in [(a[1],'target_score'),(a[2],'reference_score')]:
        artist=ax.scatter(d.x*scale,d.y*scale,c=d[col],s=7,cmap='viridis',linewidths=0,rasterized=True)
        fig.colorbar(artist,ax=ax,shrink=.72,pad=.02)
    rows=association_rows(folder)
    a[4].plot(range(3),rows.spearman_rho,'o-',color=BLUE,lw=1)
    a[4].set(xticks=range(3),xticklabels=['Raw','Technical','Composition'],ylabel='Spearman rho',ylim=(0,.7))
    a[4].tick_params(axis='x',rotation=25)
    controls=read(folder/'Results','matched_rhos.tsv')
    summary=read(folder/'Results','matched_summary.tsv').iloc[0]
    a[5].hist(controls.rho,bins=25,color='#B5CDD9',edgecolor='white')
    a[5].axvline(summary.observed_rho,color=ORANGE,lw=1.5,label='Observed')
    a[5].axvline(summary.null_95th_percentile,color='#555555',ls='--',lw=1,label='Null 95%')
    a[5].set(xlabel='Adjusted rho',ylabel='Matched sets')
    a[5].legend(loc='upper center',bbox_to_anchor=(.5,1.17),ncol=2,columnspacing=.7,handlelength=1.2)
    save(fig,a,'Figure5_external_execution')


if __name__=='__main__':
    primary_figures()
    external_figure()
