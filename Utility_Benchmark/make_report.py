"""Generate a compact, balanced report from every prespecified result."""
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font

ROOT=Path(__file__).resolve().parent
METHODS=['raw_product','technical_product','cadspot_product','adjusted_target','adjusted_minimum','spline_product']
LABELS=['Raw','Tech','OLS','Target','Min','Spline']
BLUE,ORANGE,GREEN='#0072B2','#D55E00','#009E73'


def read(name):return pd.read_csv(ROOT/'Results'/name,sep='\t')


def figures():
    plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,
        'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.frameon':False,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(8.4,5.8),layout='constrained');a=axes.ravel()
    fig.set_constrained_layout_pads(w_pad=.09,h_pad=.12,wspace=.12,hspace=.14)
    fig.get_layout_engine().set(rect=(.01,.01,.98,.94))
    summary=read('accuracy_summary.tsv')
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#E8E8E8')
    pure=['no_nuisance','linear_compact','linear_dispersed','noisy_proxy','nonlinear_compact','architecture_coupled','no_program']
    real=['primary_compact','primary_dispersed','external_compact','external_dispersed']
    for ax,experiment,order,labels in [(a[0],'simulation',pure,['None','Linear C','Linear D','Noisy','Nonlinear','Coupled','Null']),
        (a[1],'injection',real,['Primary C','Primary D','External C','External D'])]:
        values=summary[summary.experiment.eq(experiment)].pivot(index='scenario',columns='method',values='average_precision_mean').reindex(index=order,columns=METHODS)
        im=ax.imshow(values.to_numpy(),vmin=0,vmax=1,cmap=cmap,aspect='auto')
        ax.set_yticks(range(len(order)),labels);ax.set_xticks(range(6),LABELS,rotation=55)
        for i in range(len(order)):
            for j in range(6):
                value=values.iloc[i,j]
                ax.text(j,i,'-' if pd.isna(value) else f'{value:.2f}',ha='center',va='center',fontsize=6.5,color='white' if pd.notna(value) and value<.65 else 'black')
    fig.colorbar(im,ax=list(a[:2]),fraction=.025,pad=.025,label='Average precision')
    diff=read('paired_accuracy_differences.tsv')
    diff=diff[(diff.experiment=='injection')&(diff.comparison=='cadspot_product-minus-raw_product')].set_index('scenario').loc[real]
    for i,row in enumerate(diff.itertuples()):
        a[2].errorbar(row.mean_difference,i,xerr=[[row.mean_difference-row.ci_low],[row.ci_high-row.mean_difference]],fmt='o',color=BLUE,capsize=2)
    a[2].axvline(0,color='#999',ls='--',lw=.7)
    a[2].set(yticks=range(4),yticklabels=['Primary C','Primary D','External C','External D'],xlabel='AP gain over raw');a[2].invert_yaxis()
    performance=read('performance_summary.tsv');datasets=['external','primary_10000','primary']
    for build,color in [('v2.1.0',ORANGE),('v2.1.1',BLUE)]:
        p=performance[performance.build.eq(build)].set_index('dataset').loc[datasets]
        a[3].plot(p.spots if 'spots' in p else [1759,10000,23392],p.kernel_seconds_median,'o-',color=color,label=build)
    a[3].set(xlabel='Spots',ylabel='Ranking time (s)',xticks=[1759,10000,23392],xticklabels=['1,759','10,000','23,392'])
    a[3].tick_params(axis='x',rotation=25);a[3].legend(loc='upper left')
    for i,(build,color) in enumerate([('v2.1.0',ORANGE),('v2.1.1',BLUE),('reference','#777777')]):
        p=performance[performance.build.eq(build)].set_index('dataset').loc[datasets]
        a[4].plot(np.arange(3)+(i-1)*.12,p.sampled_peak_rss_mib_median,'o',color=color,label='Reference' if build=='reference' else build)
    a[4].set(xticks=range(3),xticklabels=['1,759','10,000','23,392'],xlabel='Spots',ylabel='Peak tree RSS (MiB)')
    a[4].legend(loc='upper center',bbox_to_anchor=(.5,1.23),ncol=2,columnspacing=.7,handlelength=.7)
    from standalone_reference import rank_reference
    d=pd.read_csv(ROOT/'Inputs/primary.tsv.gz',sep='\t')
    reference=rank_reference(d,['log_nCount','log_nFeature','Epithelial_proxy','ECM_proxy','Contractile_proxy','Myeloid_proxy','Tcell_proxy'])
    base=Path(json.loads((ROOT/'Inputs/sources.json').read_text())[0]['source']).parent/'Results/composition_coordinates.tsv.gz'
    if not base.exists():base=ROOT.parent/'Primary/Results/composition_coordinates.tsv.gz'
    ranked=pd.read_csv(base,sep='\t');a[5].scatter(reference['priority_score'],ranked.priority_score,s=3,color=GREEN,alpha=.3,rasterized=True)
    bound=float(ranked.priority_score.max());a[5].plot([0,bound],[0,bound],color='#444',lw=.7,ls='--')
    a[5].set(xlabel='Independent OLS product',ylabel='CAD-SPOT priority')
    fig.canvas.draw()
    for i,ax in enumerate(a):
        row,col=divmod(i,3);left=min(a[col].get_position().x0,a[col+3].get_position().x0)
        top=max(a[j].get_position().y1 for j in range(row*3,row*3+3))
        fig.text(max(.005,left-.038),min(.975,top+.022),'ABCDEF'[i],fontsize=12,fontweight='bold')
    fig.savefig(ROOT/'Results/Figure4_task_utility.pdf');fig.savefig(ROOT/'Results/Figure4_task_utility.png',dpi=300);plt.close(fig)
    if (ROOT.parent/'Additional_Figures').exists():
        for suffix in ['.pdf','.png']:
            shutil.copy2(ROOT/f'Results/Figure4_task_utility{suffix}',ROOT.parent/f'Additional_Figures/Figure4_task_utility{suffix}')


def report():
    data=read('accuracy_summary.tsv');performance=read('performance_summary.tsv');differences=read('paired_accuracy_differences.tsv')
    lines=['# CAD-SPOT: 实用增益对照结果','',
        '本报告保留预先规定的全部结果。结论是“有条件的任务收益与可复现工程收益”，不是证明CAD-SPOT算法普遍优于现有方法。','',
        '## 1. 已知真值的模拟','',
        '| 场景 | 原始乘积 AP | 仅技术调整 AP | CAD-SPOT线性调整 AP | 调整后最小值 AP | 样条敏感性 AP |',
        '|---|---:|---:|---:|---:|---:|']
    for case in ['no_nuisance','linear_compact','linear_dispersed','noisy_proxy','nonlinear_compact','architecture_coupled']:
        sub=data[(data.experiment=='simulation')&(data.scenario==case)].set_index('method')
        values=[sub.loc[m,'average_precision_mean'] for m in ['raw_product','technical_product','cadspot_product','adjusted_minimum','spline_product']]
        lines.append('| '+case+' | '+' | '.join(f'{v:.3f}' for v in values)+' |')
    lines+=['','每个场景40个配对种子，每次两张24×28网格。AP在所有spot上评价，不仅保留正象限。无程序场景AP未定义，候选负担另见完整表；不填成0，也不推导临床假阳性率。','',
        '## 2. 真实空间背景上的盲注入','',
        '| 背景/几何 | 原始 AP | 仅技术调整 AP | CAD-SPOT AP | 相对原始的AP差值（95%配对重复抽样区间） |',
        '|---|---:|---:|---:|---|']
    for case in ['primary_compact','primary_dispersed','external_compact','external_dispersed']:
        sub=data[(data.experiment=='injection')&(data.scenario==case)].set_index('method')
        diff=differences[(differences.experiment=='injection')&(differences.scenario==case)&(differences.comparison=='cadspot_product-minus-raw_product')].iloc[0]
        lines.append('| '+case+' | '+' | '.join(f'{sub.loc[m,"average_precision_mean"]:.3f}' for m in ['raw_product','technical_product','cadspot_product'])+f' | {diff.mean_difference:+.3f} [{diff.ci_low:+.3f}, {diff.ci_high:+.3f}] |')
    lines+=['','每种背景/几何32个配对种子；注入位置不看表达，预算每切片5%，幅度1.2/1.1个原背景评分SD。主背景使用前两张切片按barcode固定每三行取一行，另有跨研究技术背景。这里只检验已知注入程序的回收，不是确认组织中的真实生物状态。外部疾病注释冲突仍然存在。','',
        '## 3. 实现速度及等价性','',
        '| spot数 | 旧版排序中位秒 | 优化版排序中位秒 | 旧/新比值 | 优化版峰值进程树RSS MiB |',
        '|---:|---:|---:|---:|---:|']
    for name,n in [('external',1759),('primary_10000',10000),('primary',23392)]:
        sub=performance[performance.dataset.eq(name)].set_index('build');old=sub.loc['v2.1.0'];new=sub.loc['v2.1.1']
        lines.append(f'| {n:,} | {old.kernel_seconds_median:.3f} | {new.kernel_seconds_median:.3f} | {old.kernel_seconds_median/new.kernel_seconds_median:.2f}× | {new.sampled_peak_rss_mib_median:.1f} |')
    lines+=['','每组5次单线程独立进程，顺序随机且成对。计时只针对排序内核，未计入矩阵下载、bootstrap及500随机集生成；完整进程墙钟时间另列。内存按启动器及递归子进程合计采样，非独占堆内存。首轮错误地只测启动器RSS，已保留并明确排除，不混入最终数据。','',
        '旧/新全部输出字段在23,392和1,759个spot上逐字段完全一致，三种排序模式的并列距离、重复坐标和缺失坐标测试也一致。速度收益来自减少逐行pandas赋值，没有改变邻域或排序公式。','',
        '## 4. 必须同时报告的不利结果','',
        '- 独立标准OLS加正向乘积与CAD-SPOT优先级数值误差小于1e-10，候选完全相同；不存在相对数学等价实现的精度优势。',
        '- 去掉干扰后，原始排序略优于组分调整；真实程序与被调整的架构相关时，调整明显损害回收。',
        '- 非线性场景中，默认线性调整不能改善排序，样条敏感性更好。不能宣称线性残差化普遍可靠。',
        '- 调整后最小值基线在四种真实背景注入中均略优于默认乘积；这批实验后不追改默认公式以制造胜出。',
        '- 精简参考实现不产生六邻域及审计字段，运行远快于完整CAD-SPOT；不能宣称CAD-SPOT比简单脚本更快。',
        '- 未进行用户研究、组织学验证或与Giotto/Squidpy同任务的完整性能比较。','',
        '## 可支持的论文主张','',
        'CAD-SPOT在已知组成相关干扰的条件下，相对原始或仅技术调整的排序改善局部注入程序的回收，并以可追踪的评分、模型敏感性及输出协议实现这一标准分析；优化版在保持全部输出不变的情况下缩短自身排序时间。其贡献属于经过测试的整合工作流，而非独有或普遍最优的排序算法。','',
        '## 不支持的主张','',
        '“比所有空间方法更准确”“比简单OLS脚本更快”“650个坐标已获生物学验证”“positive product是新算法”均未获支持。']
    (ROOT/'Results/实用增益对照报告.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    wb=Workbook();index=wb.active;index.title='Contents';index.append(['Sheet','Description'])
    for name,file in [('S13_Simulation','simulation_replicates.tsv'),('S14_Injection','injection_replicates.tsv'),
        ('S15_Paired_Differences','paired_accuracy_differences.tsv'),('S16_Performance','performance_replicates.tsv')]:
        frame=read(file);ws=wb.create_sheet(name);ws.append(list(frame.columns))
        for row in frame.itertuples(index=False,name=None):ws.append([None if pd.isna(v) else v.item() if hasattr(v,'item') else v for v in row])
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
        for cell in ws[1]:cell.font=Font(bold=True)
        for i,col in enumerate(frame.columns,1):ws.column_dimensions[ws.cell(1,i).column_letter].width=min(34,max(14,len(col)+1))
        index.append([name,file])
    eq=json.loads((ROOT/'Results/equivalence.json').read_text());ws=wb.create_sheet('S17_Equivalence')
    ws.append(list(eq['real_tables'][0]))
    for row in eq['real_tables']:ws.append(list(row.values()))
    index.append(['S17_Equivalence','Exact old/new output and independent reference checks'])
    wb.save(ROOT/'Results/Utility_Tables_S13_S17.xlsx')


if __name__=='__main__':report();figures()
