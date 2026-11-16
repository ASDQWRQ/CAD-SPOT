"""Editable, ordered follow-up workbook; full coordinate tables remain TSV."""
from pathlib import Path
import sys

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

ROOT=Path(__file__).resolve().parents[1]


def read(path):return pd.read_csv(path,sep='\t',keep_default_na=False,na_values=[''])


def combined(name):
    frames=[]
    for label,folder in [('ovarian',ROOT/'Primary'),('pdac_transfer',ROOT/'PDAC_Transfer')]:
        d=read(folder/name);d.insert(0,'dataset',label);frames.append(d)
    return pd.concat(frames,ignore_index=True)


def main():
    out=ROOT/'Additional_Tables/Additional_file_1_Supplementary_Tables.xlsx'
    sheets=[]
    sheets.append(('S1_Panels',combined('panel_coverage.tsv')))
    import json
    rows=[]
    for folder in [ROOT/'Primary',ROOT/'PDAC_Transfer']:
        manifest=json.loads((folder/'input_manifest.json').read_text())
        for source in manifest:rows.append(json.loads((folder/source['sample_id']/'qc.json').read_text()))
    sheets.append(('S2_Input_QC',pd.DataFrame(rows)))
    sheets.append(('S3_Associations',combined('Results/associations.tsv')))
    intervals=[]
    for dataset,folder in [('ovarian',ROOT/'Primary'),('pdac_transfer',ROOT/'PDAC_Transfer')]:
        for path in sorted((folder/'Results').glob('bootstrap_*_ci.tsv')):
            frame=read(path);frame.insert(0,'dataset',dataset);intervals.append(frame)
    counts=combined('Results/candidate_counts.tsv')
    comparisons=combined('Results/composition_baseline_comparison.tsv')
    counts['record_type']='candidate_counts';comparisons['record_type']='ranking_comparison'
    legacy=read(ROOT/'Analysis/Cross_Cancer/legacy_ranking_comparison.tsv')
    sheets.append(('S4_Ranking_and_Counts',pd.concat([counts,comparisons,legacy],ignore_index=True)))
    sheets.append(('S5_Block_Intervals',pd.concat(intervals,ignore_index=True)))
    sheets.append(('S6_Matched_Summary',combined('Results/matched_summary.tsv')))
    sheets.append(('S7_Matched_Assignments',combined('Results/matched_assignments.tsv.gz')))
    sheets.append(('S8_Ranking_Simulation',read(ROOT/'Ranking_Ablation/all_replicates.tsv')))
    sheets.append(('S9_Paired_Differences',read(ROOT/'Ranking_Ablation/paired_differences.tsv')))
    sheets.append(('S10_Common_Gene_Check',combined('Results/common_gene_sensitivity_associations.tsv')))
    sheets.append(('S11_Model_Scope',read(ROOT/'Comparators/summary/primary_model_scope_audit.tsv')))
    sheets.append(('S12_Model_Enrichment',read(ROOT/'Comparators/current_summary/candidate_factor_enrichment.tsv')))
    utility=ROOT/'Utility_Benchmark/Results'
    for name,filename in [('S13_Utility_Simulation','simulation_replicates.tsv'),('S14_Real_Injection','injection_replicates.tsv'),
        ('S15_Paired_AP_Differences','paired_accuracy_differences.tsv'),('S16_Performance','performance_replicates.tsv')]:
        sheets.append((name,read(utility/filename)))
    sheets.append(('S17_Output_Equivalence',pd.DataFrame(json.loads((utility/'equivalence.json').read_text())['real_tables'])))
    decision=ROOT/'Decision_Utility/Results'
    for name,filename in [('S18_Fixed_Budget','fixed_budget_replicates.tsv.gz'),('S19_Inspection_Cost','inspection_cost_replicates.tsv.gz'),
        ('S20_Model_Agreement','application_audit_summary.tsv'),('S21_Input_Perturbation','perturbation_replicates.tsv.gz'),
        ('S22_Perturbation_Status','perturbation_qc.tsv')]:
        sheets.append((name,read(decision/filename)))
    coverage=read(decision/'second_pair_coverage.tsv');coverage['record_type']='coverage'
    second=read(decision/'second_pair_associations.tsv');second['record_type']='association'
    sheets.append(('S23_Second_Pair',pd.concat([coverage,second],ignore_index=True)))
    book=Workbook();contents=book.active;contents.title='Contents'
    contents.append(['Sheet','Rows','Description'])
    notes=[
        'Declared/used genes; proxies are not validated cell fractions',
        'Ovarian primary example and previously examined PDAC transfer; not untouched validation',
        'Raw, technical and composition-proxy models, pooled and separate sections',
        'Candidate counts and ranking comparisons; outlier heuristic, no coordinate-level FDR',
        'Within-observed-section spatial block percentile intervals; not patient-population CIs',
        'Observed and random sets share scoring and residual model',
        'Strict matching, all assignments, seed and caliper differences',
        'All 50 paired replicates for compact/dispersed/no-program geometries',
        'All method-minus-product differences; no selected seeds',
        'Post hoc shared target/reference coverage across cancers; not replacement main scores',
        'PDAC transfer versus preserved PDAC model-input and joint valid coverage',
        'PDAC-only preserved model overlap; not ovarian training or accuracy evidence',
        'All 40 paired replicates per prespecified utility scenario, including null/negative cases',
        'All 32 paired injections per background/geometry; not biological ground truth',
        'AP differences across replicates, not patient-population inference',
        'Five randomized subprocess repeats; final corrected process-tree RSS round only',
        'Exact old/new complete output; independent OLS-product agreement',
        'All fixed-budget section outcomes; aggregate within replicate for comparisons',
        'Inspections for 80 percent known-injection recovery; undefined with no positives',
        'Full-cohort model membership counts for both predeclared score pairs',
        'Every computable count/panel/proxy sensitivity; failures explicitly in S22',
        'All scheduled perturbations, constant-gene failures and original-QC changes',
        'Illustrative IFN/antigen panels and all associations, not validated signatures']
    for (name,frame),note in zip(sheets,notes):
        contents.append([name,len(frame),note]);ws=book.create_sheet(name)
        ws.append(list(frame.columns))
        for values in frame.itertuples(index=False,name=None):
            ws.append([None if pd.isna(v) else v.item() if hasattr(v,'item') else v for v in values])
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
        for cell in ws[1]:cell.font=Font(name='Arial',bold=True,size=10);cell.alignment=Alignment(wrap_text=True,vertical='top')
        for j,col in enumerate(frame.columns,1):ws.column_dimensions[ws.cell(1,j).column_letter].width=min(34,max(12,len(col)+1))
    contents.freeze_panes='A2';contents.column_dimensions['A'].width=28;contents.column_dimensions['B'].width=12;contents.column_dimensions['C'].width=105
    for cell in contents[1]:cell.font=Font(name='Arial',bold=True,size=10)
    book.save(out)
    print(out,flush=True)


if __name__=='__main__':main()
