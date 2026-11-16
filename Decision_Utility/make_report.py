"""Summarize every extension result and draw compact decision-utility figures."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Analysis'))
from make_figures import grid, save, BLUE, ORANGE, GREEN, COLORS

R=ROOT/'Decision_Utility/Results'
METHODS=['raw_product','technical_product','cadspot_product','adjusted_target','adjusted_minimum','spline_product']
CASES=['ovarian_compact','ovarian_dispersed','pdac_compact','pdac_dispersed']
LABELS=['OV C','OV D','PDAC C','PDAC D']


def read(path):return pd.read_csv(path,sep='\t',keep_default_na=False,na_values=[''])


def main():
    fixed=read(R/'fixed_budget_summary.tsv')
    cost=read(R/'inspection_cost_paired.tsv')
    perturb=read(R/'perturbation_replicates.tsv.gz')
    qc=read(R/'perturbation_qc.tsv')
    all_audits=[]; all_pairs=[]; panel_coverage=[]
    summary={}
    for cohort in ['ovarian','pdac']:
        summary[cohort]={}
        for program in ['stress','ifn']:
            d=read(R/f'{cohort}_{program}_audit.coordinates.tsv.gz')
            audit_summary=read(R/f'{cohort}_{program}_audit.summary.tsv')
            audit_summary.insert(0,'cohort',cohort);audit_summary.insert(1,'program',program);all_audits.append(audit_summary)
            top=d.conditional_votes_top25.gt(0);common=d.conditional_votes_top25.eq(3)
            summary[cohort][program]=dict(any_three=int(top.sum()),all_three=int(common.sum()),
                agreement_fraction=float(common.sum()/top.sum()),
                composition_selected=int(d.composition_top25.sum()),
                technical_composition_changed=int(d.technical_composition_changed_top25.sum()))
        pair=ROOT/f'Decision_Utility/Program_Pair/{cohort}'
        assoc=read(pair/'associations.tsv');assoc.insert(0,'cohort',cohort);all_pairs.append(assoc)
        cov=read(pair/'panel_coverage.tsv');cov.insert(0,'cohort',cohort);panel_coverage.append(cov)
        manifest=json.loads((pair/'input_manifest.json').read_text())
        summary[cohort]['ifn_used']={p:len(v) for p,v in manifest['used_panels'].items()}
        pooled=assoc[assoc.scope.eq('pooled') & assoc.method.eq('score_residual_spearman')]
        summary[cohort]['ifn_rhos']=dict(zip(pooled.layer,pooled.spearman_rho))
        summary[cohort]['perturbations']={}
        for name in ['count_80pct','count_50pct','panel_loss_10pct']:
            q=qc[qc.cohort.eq(cohort)&qc.perturbation.eq(name)]
            sub=perturb[perturb.cohort.eq(cohort)&perturb.perturbation.eq(name)&perturb.model.eq('composition')&perturb.budget.eq(25)]
            summary[cohort]['perturbations'][name]=dict(requested=len(q),valid=int(q.status.eq('ok').sum()),
                mean_top25_retention=float(sub.retention.mean()),
                mean_below_original_qc=float(q.loc[q.status.eq('ok'),'below_original_qc'].mean()))
    pd.concat(all_audits).to_csv(R/'application_audit_summary.tsv',sep='\t',index=False)
    pd.concat(all_pairs).to_csv(R/'second_pair_associations.tsv',sep='\t',index=False)
    pd.concat(panel_coverage).to_csv(R/'second_pair_coverage.tsv',sep='\t',index=False)
    summary['fixed_budget25']=fixed[fixed.experiment.eq('injection')&fixed.budget.eq(25)].to_dict('records')
    summary['cost80']=cost[cost.experiment.eq('injection')].groupby(['scenario','method']).inspections.mean().reset_index().to_dict('records')
    summary['failures']=qc[qc.status.ne('ok')].to_dict('records')
    (R/'decision_evidence.json').write_text(json.dumps(summary,indent=2)+'\n')

    fig,a=grid()
    heat=fixed[fixed.experiment.eq('injection')&fixed.budget.eq(25)].pivot(index='scenario',columns='method',values='precision_mean').loc[CASES,METHODS]
    precision_table=heat.copy()
    im=a[0].imshow(heat,vmin=0,vmax=1,cmap='viridis',aspect='auto')
    for i in range(4):
        for j in range(6):a[0].text(j,i,f'{heat.iloc[i,j]:.2f}',ha='center',va='center',fontsize=7,color='white')
    a[0].set(yticks=range(4),yticklabels=LABELS,xticks=range(6),xticklabels=['Raw','Tech','OLS','Target','Min','Spline'])
    a[0].tick_params(axis='x',rotation=55);a[0].set_ylabel('Precision at 25')
    c=cost[cost.experiment.eq('injection')].groupby(['scenario','method']).inspections.mean().unstack().loc[CASES]
    for i in range(4):
        a[1].plot([c.iloc[i].raw_product,c.iloc[i].cadspot_product],[i,i],color='#999',lw=1)
    a[1].plot(c.raw_product,range(4),'o',color=BLUE,label='Raw')
    a[1].plot(c.cadspot_product,range(4),'s',color=ORANGE,label='OLS')
    a[1].set(yticks=range(4),yticklabels=LABELS,xlabel='Inspections for 80% recall');a[1].invert_yaxis()
    a[1].legend(loc='upper center',bbox_to_anchor=(.5,1.15),ncol=2)
    observed=read(R/'ovarian_stress_audit.coordinates.tsv.gz')
    sub=observed[(observed.technical_top25|observed.composition_top25)]
    for group,color in [('all_three',GREEN),('some_models',ORANGE)]:
        s=sub[sub.model_agreement.eq(group)]
        a[2].scatter(s.technical_audit_rank,s.composition_audit_rank,s=7,alpha=.6,color=color,label=['All three','Some models'][group=='some_models'],rasterized=True)
    a[2].axhline(25,color='#888',lw=.7,ls='--');a[2].axvline(25,color='#888',lw=.7,ls='--')
    a[2].set(xscale='log',yscale='log',xlabel='Technical rank',ylabel='Composition rank')
    a[2].legend(loc='upper center',bbox_to_anchor=(.5,1.18),ncol=1)
    names=['count_80pct','count_50pct','panel_loss_10pct']
    for cohort,color,offset in [('ovarian',BLUE,-.08),('pdac',ORANGE,.08)]:
        for i,name in enumerate(names):
            sub=perturb[perturb.cohort.eq(cohort)&perturb.perturbation.eq(name)&perturb.model.eq('composition')&perturb.budget.eq(25)]
            means=sub.groupby('replicate').retention.mean()
            a[3].errorbar(i+offset,means.mean(),yerr=means.std(),fmt='o',color=color,capsize=2,label=('OV' if cohort=='ovarian' else 'PDAC') if i==0 else None)
    a[3].set(xticks=range(3),xticklabels=['80% counts','50% counts','10% loss'],ylabel='Top-25 retention',ylim=(0,1.03));a[3].tick_params(axis='x',rotation=25)
    a[3].legend(loc='lower left')
    pair=pd.concat(all_pairs)
    for cohort,color in [('ovarian',BLUE),('pdac',ORANGE)]:
        d=pair[pair.cohort.eq(cohort)&pair.scope.eq('pooled')&pair.method.eq('score_residual_spearman')]
        a[4].plot(range(3),d.set_index('layer').loc[['raw','technical_adjusted','composition_adjusted'],'spearman_rho'],'o-',color=color,label='OV' if cohort=='ovarian' else 'PDAC',lw=1)
    a[4].set(xticks=range(3),xticklabels=['Raw','Technical','Composition'],ylabel='IFN-pair rho');a[4].tick_params(axis='x',rotation=25);a[4].legend(loc='upper right')
    for cohort,color in [('ovarian',BLUE),('pdac',ORANGE)]:
        a[5].plot([0,1],[summary[cohort][p]['agreement_fraction'] for p in ['stress','ifn']],'o-',color=color,label='OV' if cohort=='ovarian' else 'PDAC',lw=1)
    a[5].set(xticks=[0,1],xticklabels=['Stress pair','IFN pair'],ylabel='All-three / any-three',ylim=(0,1));a[5].legend(loc='upper right')
    save(fig,a,'Figure3_decision_utility')

    fig,a=grid()
    diff=read(R/'fixed_budget_paired_differences.tsv')
    for i,scenario in enumerate(CASES):
        s=diff[diff.scenario.eq(scenario)&diff.metric.eq('precision')&diff.comparison.eq('cadspot_product-minus-raw_product')].sort_values('budget')
        a[0].plot(s.budget,s.mean_difference,'o-',lw=1,color=COLORS[i],label=LABELS[i])
        a[0].fill_between(s.budget,s.ci_low,s.ci_high,color=COLORS[i],alpha=.1)
    a[0].axhline(0,color='#777',lw=.7);a[0].set(xlabel='Budget per section',ylabel='Precision gain over raw')
    a[0].legend(fontsize=7,loc='upper center',bbox_to_anchor=(.5,1.18),ncol=2)
    for ax,cohort in zip(a[1:3],['ovarian','pdac']):
        s=perturb[perturb.cohort.eq(cohort)&perturb.budget.eq(25)&perturb.perturbation.isin(names)]
        heat=s.groupby(['perturbation','model']).retention.mean().unstack().loc[names,['raw','technical','composition','spline']]
        ax.imshow(heat,vmin=0,vmax=1,cmap='viridis',aspect='auto')
        for i in range(3):
            for j in range(4):ax.text(j,i,f'{heat.iloc[i,j]:.2f}',ha='center',va='center',fontsize=7,color='white')
        ax.set(yticks=range(3),yticklabels=['80% counts','50% counts','10% loss'],xticks=range(4),xticklabels=['Raw','Tech','OLS','Spline'],ylabel=cohort.upper());ax.tick_params(axis='x',rotation=45)
    for cohort,color,offset in [('ovarian',BLUE,-.08),('pdac',ORANGE,.08)]:
        s=perturb[perturb.cohort.eq(cohort)&perturb.budget.eq(25)&perturb.model.eq('composition')&perturb.perturbation.str.startswith('omit_')]
        means=s.groupby('perturbation').retention.mean()
        a[3].plot(np.arange(len(means))+offset,means.values,'o',color=color,label='OV' if cohort=='ovarian' else 'PDAC')
    a[3].set(xticks=range(5),xticklabels=['Contr.','ECM','Epi','Myeloid','T cell'],ylabel='Top-25 retention',ylim=(0,1.03));a[3].tick_params(axis='x',rotation=45);a[3].legend(loc='lower right')
    # Fractions of successful scheduled replicates; failures remain in the denominator.
    for cohort,color,offset in [('ovarian',BLUE,-.1),('pdac',ORANGE,.1)]:
        values=[summary[cohort]['perturbations'][n]['valid']/20 for n in names]
        a[4].bar(np.arange(3)+offset,values,width=.18,color=color,label='OV' if cohort=='ovarian' else 'PDAC')
    a[4].set(xticks=range(3),xticklabels=['80% counts','50% counts','10% loss'],ylabel='Computable replicates',ylim=(0,1.1));a[4].tick_params(axis='x',rotation=25);a[4].legend(loc='upper center',bbox_to_anchor=(.5,1.16),ncol=2)
    for i,cohort in enumerate(['ovarian','pdac']):
        for j,kind in enumerate(['target_score','reference_score']):
            a[5].bar(i+(j-.5)*.25,summary[cohort]['ifn_used'][kind],width=.23,color=[BLUE,GREEN][j],label=['IFN','Antigen'][j] if i==0 else None)
    a[5].set(xticks=[0,1],xticklabels=['OV','PDAC'],ylabel='Second-pair used genes',ylim=(0,16));a[5].legend(loc='upper center',bbox_to_anchor=(.5,1.16),ncol=2,fontsize=7,columnspacing=.7,handlelength=1)
    save(fig,a,'FigureS3_budget_and_robustness')
    lines=['# Decision-utility results','', 'All primary budgets and perturbations were specified before outcome inspection.', '',
        'Top-25 mean precision (known score injections, not pathology):', '',precision_table.to_string(), '',
        'Complete numerical records are in decision_evidence.json and the accompanying tables.', '',
        'Count/gene-loss failures are retained explicitly. Model agreement is not correctness or FDR.']
    (R/'Decision_Utility_Report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({k:{v:s[v] for v in ['stress','ifn','ifn_rhos','perturbations']} for k,s in summary.items() if k in ['ovarian','pdac']},indent=2),flush=True)


if __name__=='__main__':main()
