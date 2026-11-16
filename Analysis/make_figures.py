"""Ovarian worked-example figures, with PDAC explicitly a transfer assessment."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'Additional_Figures'
BLUE,ORANGE,GREEN='#0072B2','#D55E00','#009E73'
COLORS=[BLUE,ORANGE,GREEN,'#CC79A7','#555555','#56B4E9','#E69F00','#882255']
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'legend.fontsize':8,
    'xtick.labelsize':8,'ytick.labelsize':8,'legend.frameon':False,'axes.spines.top':False,
    'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'savefig.facecolor':'white'})


def read(folder,name):return pd.read_csv(folder/name,sep='\t',keep_default_na=False,na_values=[''])


def grid(rows=2,cols=3):
    fig,axes=plt.subplots(rows,cols,figsize=(8.4,5.8),layout='constrained')
    fig.set_constrained_layout_pads(w_pad=.09,h_pad=.12,wspace=.10,hspace=.14)
    fig.get_layout_engine().set(rect=(.01,.01,.98,.94))
    return fig,np.asarray(axes).ravel()


def save(fig,axes,name):
    OUT.mkdir(exist_ok=True);fig.canvas.draw()
    cols=axes[0].get_subplotspec().get_gridspec().ncols
    for i,ax in enumerate(axes):
        row,col=divmod(i,cols)
        left=min(axes[j].get_position().x0 for j in range(col,len(axes),cols))
        top=max(axes[j].get_position().y1 for j in range(row*cols,(row+1)*cols))
        fig.text(max(.005,left-.038),min(.975,top+.022),'ABCDEFGHI'[i],fontsize=12,fontweight='bold')
    fig.savefig(OUT/f'{name}.pdf');fig.savefig(OUT/f'{name}.png',dpi=300);plt.close(fig)


def layers(folder):
    return read(folder/'Results','associations.tsv').query("scope=='pooled' and method=='score_residual_spearman'")


def image_view(ax,spots,source,candidates=True):
    spatial=Path(source['spatial_dir'])
    scale=json.loads((spatial/'scalefactors_json.json').read_text())['tissue_hires_scalef']
    im=Image.open(spatial/'tissue_hires_image.png').convert('RGB');ax.imshow(im)
    x=spots.x.to_numpy()*scale;y=spots.y.to_numpy()*scale
    if candidates:
        hot=spots.cadspot_candidate.eq(True).to_numpy()
        ax.scatter(x[hot],y[hot],s=9,facecolors='none',edgecolors=ORANGE,linewidths=.55)
    ax.set_xlim(max(0,x.min()-35),min(im.width,x.max()+35))
    ax.set_ylim(min(im.height,y.max()+35),max(0,y.min()-35))
    ax.set_aspect('equal');ax.axis('off')


def overview():
    folder=ROOT/'Primary';d=read(folder,'score_table.tsv.gz');cov=read(folder,'panel_coverage.tsv')
    fig,a=grid();a[0].set(xlim=(0,1),ylim=(0,1));a[0].axis('off')
    for y,text in [(.82,'Counts + coordinates'),(.56,'Scores + covariates'),(.30,'Audits + priorities')]:
        a[0].add_patch(Rectangle((.03,y-.08),.94,.16,facecolor='#EEF2F4',edgecolor='#9FA8AE',lw=.7))
        a[0].text(.5,y,text,ha='center',va='center',fontsize=9)
    for y in [.70,.44]:a[0].annotate('',xy=(.5,y-.08),xytext=(.5,y),arrowprops=dict(arrowstyle='->',color='#555'))
    a[0].text(.5,.10,'Block CI | Matched sets',ha='center',fontsize=8)
    for j,(cohort,f) in enumerate([('OV',folder),('PDAC',ROOT/'PDAC_Transfer')]):
        c=read(f,'panel_coverage.tsv')
        values=[c[c.panel.eq(p)&c.status.eq('used')].gene.nunique() for p in ['target_score','reference_score']]
        for i,(v,color) in enumerate(zip(values,[BLUE,GREEN])):
            a[1].bar(j+(i-.5)*.30,v,width=.28,color=color,label=['Target','Reference'][i] if j==0 else None)
    a[1].set(xticks=[0,1],xticklabels=['OV','PDAC'],ylabel='Used genes',ylim=(0,28))
    a[1].legend(loc='upper center',bbox_to_anchor=(.5,1.14),ncol=2)
    counts=d.groupby('sample_id').size();a[2].barh(range(len(counts)),counts.values,color=BLUE,height=.6)
    a[2].set(yticks=range(len(counts)),yticklabels=[f'OV{i+1}' for i in range(len(counts))],xlabel='Spots');a[2].invert_yaxis()
    a[3].hexbin(d.target_score,d.reference_score,gridsize=34,mincnt=1,cmap='Blues',rasterized=True)
    a[3].set(xlabel='Target score',ylabel='Reference score')
    panels=['Epithelial_proxy','ECM_proxy','Contractile_proxy','Myeloid_proxy','Tcell_proxy']
    for i,p in enumerate(panels):a[4].boxplot([d[p]],positions=[i],widths=.5,showfliers=False,patch_artist=True,
        boxprops=dict(facecolor=COLORS[i],alpha=.6),medianprops=dict(color='black'))
    a[4].set(xticks=range(5),xticklabels=['Epi','ECM','Contract.','Myeloid','T cell'],ylabel='Proxy score');a[4].tick_params(axis='x',rotation=35)
    corr=d[['target_score','reference_score']+panels].corr(method='spearman')
    im=a[5].imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r');labels=['Target','Ref.','Epi','ECM','Contr.','Myeloid','T cell']
    a[5].set(xticks=range(7),yticks=range(7),xticklabels=labels,yticklabels=labels);a[5].tick_params(axis='x',rotation=60)
    fig.colorbar(im,ax=a[5],shrink=.75,pad=.03);save(fig,a,'Figure1_workflow_ovarian_inputs')


def forest(ax,ci,prefix):
    for i,r in enumerate(ci.itertuples()):
        ax.plot([r.ci_low,r.ci_high],[i,i],color=COLORS[i%len(COLORS)],lw=1.4)
        ax.plot(r.spearman_rho,i,'o',color=COLORS[i%len(COLORS)],ms=4)
    ax.set(yticks=range(len(ci)),yticklabels=[f'{prefix}{i+1}' for i in range(len(ci))],xlabel='Adjusted rho');ax.invert_yaxis()


def null_hist(ax,folder):
    n=read(folder/'Results','matched_rhos.tsv');s=read(folder/'Results','matched_summary.tsv').iloc[0]
    ax.hist(n.rho,bins=25,color='#B5CDD9',edgecolor='white')
    ax.axvline(s.observed_rho,color=ORANGE,lw=1.5,label='Observed')
    ax.axvline(s.null_95th_percentile,color='#555',ls='--',lw=1,label='Null 95%')
    ax.set(xlabel='Adjusted rho',ylabel='Matched sets')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,1.17),ncol=2,columnspacing=.7,handlelength=1.2)


def ovarian_audits():
    folder=ROOT/'Primary';fig,a=grid();p=layers(folder)
    a[0].plot(range(3),p.spearman_rho,'o-',color=BLUE,lw=1)
    a[0].set(xticks=range(3),xticklabels=['Raw','Technical','Composition'],ylabel='Spearman rho');a[0].tick_params(axis='x',rotation=25)
    s=read(folder/'Results','associations.tsv').query("scope=='pooled' and layer=='composition_adjusted'")
    a[1].scatter(range(3),s.spearman_rho,c=[BLUE,GREEN,ORANGE],s=32)
    a[1].set(xticks=range(3),xticklabels=['OLS','Partial rank','Spline'],ylabel='Adjusted rho');a[1].tick_params(axis='x',rotation=25)
    for i,width in enumerate([6,12,18]):
        r=read(folder/'Results',f'bootstrap_{width}_ci.tsv').query("scope=='pooled' and layer=='composition_adjusted'").iloc[0]
        a[2].plot([i,i],[r.ci_low,r.ci_high],color=BLUE,lw=1.5);a[2].plot(i,r.spearman_rho,'o',color=BLUE)
    a[2].set(xticks=range(3),xticklabels=['6','12','18'],xlabel='Block width (pitches)',ylabel='Adjusted rho')
    ci=read(folder/'Results','bootstrap_12_ci.tsv').query("scope=='section' and layer=='composition_adjusted'")
    forest(a[3],ci,'OV');null_hist(a[4],folder)
    c=read(folder/'Results','candidate_counts.tsv').pivot(index='sample_id',columns='layer',values='n_candidates')
    for i,r in enumerate(c.itertuples()):a[5].plot([i,i],[r.technical,r.composition],color='#9FA8AE',lw=1)
    a[5].scatter(range(len(c)),c.technical,color=BLUE,label='Technical',s=20)
    a[5].scatter(range(len(c)),c.composition,color=ORANGE,label='Composition',s=20,marker='s')
    a[5].set(xticks=range(len(c)),xticklabels=[f'{i+1}' for i in range(len(c))],xlabel='OV section',ylabel='Candidate spots')
    a[5].legend(loc='upper center',bbox_to_anchor=(.5,1.18),ncol=1)
    save(fig,a,'Figure4_ovarian_association_audits')


def ovarian_coordinates():
    folder=ROOT/'Primary';r=read(folder/'Results','composition_coordinates.tsv.gz')
    manifest=json.loads((folder/'input_manifest.json').read_text());fig,a=grid(2,4)
    for ax,source in zip(a,manifest):image_view(ax,r[r.sample_id.eq(source['sample_id'])],source)
    save(fig,a,'Figure5_ovarian_candidate_coordinates')


def pdac_transfer():
    f=ROOT/'PDAC_Transfer';r=read(f/'Results','composition_coordinates.tsv.gz')
    source=json.loads((f/'input_manifest.json').read_text())[0];fig,a=grid()
    image_view(a[0],r[r.sample_id.eq(source['sample_id'])],source)
    for name,folder,color in [('OV',ROOT/'Primary',BLUE),('PDAC',f,ORANGE)]:
        p=layers(folder);a[1].plot(range(3),p.spearman_rho,'o-',lw=1,color=color,label=name)
    a[1].set(xticks=range(3),xticklabels=['Raw','Technical','Composition'],ylabel='Spearman rho');a[1].tick_params(axis='x',rotation=25)
    a[1].legend(loc='upper center',bbox_to_anchor=(.5,1.14),ncol=2)
    ci=read(f/'Results','bootstrap_12_ci.tsv').query("scope=='section' and layer=='composition_adjusted'");forest(a[2],ci,'P')
    null_hist(a[3],f)
    common=read(ROOT/'Analysis/Cross_Cancer','common_gene_summary.tsv').set_index('cohort')
    for i,(name,folder,color) in enumerate([('ovarian',ROOT/'Primary',BLUE),('pdac',f,ORANGE)]):
        full=float(layers(folder).query("layer=='composition_adjusted'").spearman_rho.iloc[0])
        a[4].plot([i-.12,i+.12],[full,common.loc[name,'composition_rho']],'o-',color=color,lw=1)
    a[4].set(xticks=[-.12,.12,.88,1.12],xticklabels=['Full','Shared','Full','Shared'],ylabel='Adjusted rho')
    a[4].tick_params(axis='x',rotation=35)
    a[4].set_xlabel('OV                         PDAC')
    for i,(label,folder) in enumerate([('OV',ROOT/'Primary'),('PDAC',f)]):
        d=read(folder,'score_table.tsv.gz');c=read(folder/'Results','candidate_counts.tsv').groupby('layer').n_candidates.sum()/len(d)*100
        a[5].plot([i-.12,i+.12],[c.technical,c.composition],'-',color='#999',lw=1)
        a[5].plot(i-.12,c.technical,'o',color=BLUE,label='Technical' if i==0 else None)
        a[5].plot(i+.12,c.composition,'s',color=ORANGE,label='Composition' if i==0 else None)
    a[5].set(xticks=[0,1],xticklabels=['OV','PDAC'],ylabel='Candidate fraction (%)')
    a[5].legend(loc='upper center',bbox_to_anchor=(.5,1.18),ncol=1)
    save(fig,a,'Figure6_pdac_cross_cancer_transfer')


def neighborhood():
    folder=ROOT/'Ranking_Ablation';d=read(folder,'all_replicates.tsv');fig,a=grid()
    for ax,geometry in zip(a[:2],['compact','dispersed']):
        e=read(folder,f'example_{geometry}.tsv.gz');e=e[e.sample_id.eq('s0')]
        ax.scatter(e.x,e.y,c=np.where(e.truth,GREEN,'#DFE4E6'),s=3,linewidths=0,rasterized=True)
        ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='Array column',ylabel='Array row')
    methods=['positive_product','legacy','positive_support']
    for j,method in enumerate(methods):
        for i,geometry in enumerate(['compact','dispersed']):
            sub=d[d.method.eq(method)&d.geometry.eq(geometry)]
            a[2].boxplot([sub.average_precision],positions=[i+(j-1)*.2],widths=.16,showfliers=False,patch_artist=True,
                boxprops=dict(facecolor=[BLUE,ORANGE,GREEN][j],alpha=.65),medianprops=dict(color='black'))
            a[3].errorbar(i+(j-1)*.13,sub.top_5pct_precision.mean(),yerr=sub.top_5pct_precision.std(),fmt='o',color=[BLUE,ORANGE,GREEN][j],capsize=2)
        a[2].plot([],[],color=[BLUE,ORANGE,GREEN][j],label=['Product','Legacy','Support'][j])
        vals=100*d[d.method.eq(method)&d.geometry.eq('null')].candidate_fraction
        a[4].boxplot([vals],positions=[j],widths=.5,showfliers=False,patch_artist=True,
            boxprops=dict(facecolor=[BLUE,ORANGE,GREEN][j],alpha=.65),medianprops=dict(color='black'))
    for ax,ylabel in [(a[2],'Average precision'),(a[3],'Top-5% precision')]:ax.set(xticks=[0,1],xticklabels=['Compact','Dispersed'],ylabel=ylabel)
    a[2].legend(loc='upper center',bbox_to_anchor=(.5,1.18),ncol=3,columnspacing=.6,handlelength=1)
    a[4].set(xticks=range(3),xticklabels=['Product','Legacy','Support'],ylabel='Null candidate spots (%)')
    ranked=read(ROOT/'Primary/Results','composition_coordinates.tsv.gz')
    positive=ranked[(ranked.target_resid_z>0)&(ranked.reference_resid_z>0)]
    for col,color,label in [('legacy_priority_score',ORANGE,'Legacy'),('local_support',GREEN,'Local support')]:
        values=[]
        for frac in [.01,.05,.10]:
            k=max(1,int(np.ceil(frac*len(positive))))
            left=set(positive.nlargest(k,'priority_score').index);right=set(positive.nlargest(k,col).index)
            values.append(len(left&right)/k)
        a[5].plot([1,5,10],values,'o-',color=color,label=label,lw=1)
    a[5].set(xlabel='Top-ranked fraction (%)',ylabel='Overlap with product',ylim=(0,1));a[5].legend(loc='lower right')
    save(fig,a,'FigureS1_neighborhood_ablation')


if __name__=='__main__':
    overview();ovarian_audits();ovarian_coordinates();pdac_transfer();neighborhood()
