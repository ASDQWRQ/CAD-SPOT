"""Balanced six-panel utility figure with explicit ovarian/PDAC roles."""
import json
from pathlib import Path
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
METHODS=['raw_product','technical_product','cadspot_product','adjusted_target','adjusted_minimum','spline_product']
LABELS=['Raw','Tech','OLS','Target','Min','Spline']
BLUE,ORANGE,GREEN='#0072B2','#D55E00','#009E73'


def read(name):return pd.read_csv(ROOT/'Results'/name,sep='\t')


def figures():
    plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'legend.fontsize':8,
        'xtick.labelsize':8,'ytick.labelsize':8,'legend.frameon':False,'axes.spines.top':False,
        'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(8.4,5.8),layout='constrained');a=axes.ravel()
    fig.set_constrained_layout_pads(w_pad=.09,h_pad=.12,wspace=.12,hspace=.14)
    fig.get_layout_engine().set(rect=(.01,.01,.98,.94))
    summary=read('accuracy_summary.tsv');cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#E8E8E8')
    pure=['no_nuisance','linear_compact','linear_dispersed','noisy_proxy','nonlinear_compact','architecture_coupled','no_program']
    real=['ovarian_compact','ovarian_dispersed','pdac_compact','pdac_dispersed']
    for ax,experiment,order,labels in [
        (a[0],'simulation',pure,['None','Linear C','Linear D','Noisy','Nonlinear','Coupled','Null']),
        (a[1],'injection',real,['OV C','OV D','PDAC C','PDAC D'])]:
        v=summary[summary.experiment.eq(experiment)].pivot(index='scenario',columns='method',values='average_precision_mean').reindex(index=order,columns=METHODS)
        im=ax.imshow(v.to_numpy(),vmin=0,vmax=1,cmap=cmap,aspect='auto')
        ax.set_yticks(range(len(order)),labels);ax.set_xticks(range(6),LABELS,rotation=55)
        for i in range(len(order)):
            for j in range(6):
                x=v.iloc[i,j]
                ax.text(j,i,'-' if pd.isna(x) else f'{x:.2f}',ha='center',va='center',fontsize=6.5,
                    color='white' if pd.notna(x) and x<.65 else 'black')
    fig.colorbar(im,ax=list(a[:2]),fraction=.025,pad=.025,label='Average precision')
    diff=read('paired_accuracy_differences.tsv').query("experiment=='injection' and comparison=='cadspot_product-minus-raw_product'").set_index('scenario').loc[real]
    for i,r in enumerate(diff.itertuples()):
        a[2].plot([r.ci_low,r.ci_high],[i,i],color=BLUE,lw=1.5)
        a[2].plot(r.mean_difference,i,'o',color=BLUE)
    a[2].axvline(0,color='#999',ls='--',lw=.7)
    a[2].set(yticks=range(4),yticklabels=['OV C','OV D','PDAC C','PDAC D'],xlabel='AP gain over raw');a[2].invert_yaxis()
    perf=read('performance_summary.tsv');counts=read('performance_replicates.tsv').groupby('dataset').spots.first()
    datasets=counts.sort_values().index.tolist();x=counts.loc[datasets].to_numpy()
    for build,color in [('v2.1.0',ORANGE),('v2.2.0',BLUE)]:
        p=perf[perf.build.eq(build)].set_index('dataset').loc[datasets]
        a[3].plot(x,p.kernel_seconds_median,'o-',color=color,label=build,lw=1)
    a[3].set(xlabel='Spots',ylabel='Ranking time (s)',xticks=x,xticklabels=[f'{n:,}' for n in x])
    a[3].tick_params(axis='x',rotation=30);a[3].legend(loc='upper left')
    for i,(build,color) in enumerate([('v2.1.0',ORANGE),('v2.2.0',BLUE),('reference','#777777')]):
        p=perf[perf.build.eq(build)].set_index('dataset').loc[datasets]
        a[4].plot(np.arange(3)+(i-1)*.13,p.sampled_peak_rss_mib_median,'o',color=color,
            label='Reference' if build=='reference' else build)
    a[4].set(xticks=range(3),xticklabels=[f'{n:,}' for n in x],xlabel='Spots',ylabel='Peak tree RSS (MiB)')
    a[4].tick_params(axis='x',rotation=30)
    a[4].legend(loc='upper center',bbox_to_anchor=(.5,1.23),ncol=2,columnspacing=.7,handlelength=.7)
    from standalone_reference import rank_reference
    from benchmark import core,TECH,COMP
    d=pd.read_csv(ROOT/'Inputs/ovarian.tsv.gz',sep='\t')
    reference=rank_reference(d,TECH+COMP)
    engine=core();_,u,v=engine.fit_pair(d,'target_score','reference_score',TECH+COMP)
    ranked=engine.spatial_ranking(d,u,v)
    a[5].scatter(reference['priority_score'],ranked.priority_score,s=3,color=GREEN,alpha=.3,rasterized=True)
    bound=float(ranked.priority_score.max());a[5].plot([0,bound],[0,bound],color='#444',lw=.7,ls='--')
    a[5].set(xlabel='Independent OLS product',ylabel='CAD-SPOT priority')
    fig.canvas.draw()
    for i,ax in enumerate(a):
        row,col=divmod(i,3);left=min(a[col].get_position().x0,a[col+3].get_position().x0)
        top=max(a[j].get_position().y1 for j in range(row*3,row*3+3))
        fig.text(max(.005,left-.038),min(.975,top+.022),'ABCDEF'[i],fontsize=12,fontweight='bold')
    for suffix in ['.pdf','.png']:
        dest=ROOT/f'Results/Figure2_task_utility{suffix}'
        fig.savefig(dest,**({'dpi':300} if suffix=='.png' else {}))
        shutil.copy2(dest,ROOT.parent/f'Additional_Figures/{dest.name}')
    plt.close(fig)


def report():
    accuracy=read('accuracy_summary.tsv');diff=read('paired_accuracy_differences.tsv')
    lines=['# 卵巢癌主要实例 / PDAC迁移：任务级对照','',
        '本轮固定公式、阈值、六种基线和种子，仅更新真实评分背景及运行规模。PDAC此前已分析，不能称未接触的外部验证。','',
        '| 背景/几何 | Raw AP | Technical AP | CAD-SPOT AP | 相对Raw差值及95%配对重复抽样区间 |',
        '|---|---:|---:|---:|---|']
    for case in ['ovarian_compact','ovarian_dispersed','pdac_compact','pdac_dispersed']:
        s=accuracy.query("experiment=='injection'");s=s[s.scenario.eq(case)].set_index('method')
        r=diff[(diff.scenario.eq(case))&diff.comparison.eq('cadspot_product-minus-raw_product')].iloc[0]
        lines.append('| '+case+' | '+' | '.join(f'{s.loc[m,"average_precision_mean"]:.3f}' for m in METHODS[:3])+
            f' | {r.mean_difference:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}] |')
    lines += ['', '每种几何32个配对种子；每个队列前两张切片按barcode每三行取一行。注入5%坐标，位置选择不查看评分，幅度为1.2/1.1个原评分SD。区间针对这些评分背景与模拟重复，不是患者总体推断或病理学验证。','',
        '所有七种纯模拟场景、六种基线、无程序场景候选负担和不利结果均保留在TSV/补充表中。独立OLS参考重现默认优先级，最小值/样条基线可能更好；不根据结果重调公式。','',
        '运行时间只测排序内核，不包括下载、评分、bootstrap或匹配集。旧版/新版完整输出必须一致；精简参考实现少了审计/邻域导出，不能据此宣称功能等价速度优势。']
    (ROOT/'Results/实用增益对照报告.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':report();figures()
