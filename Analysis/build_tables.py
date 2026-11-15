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
    for label,folder in [('primary',ROOT/'Primary'),('external',ROOT/'External_EndToEnd/NYU_PDAC1')]:
        d=read(folder/name);d.insert(0,'dataset',label);frames.append(d)
    return pd.concat(frames,ignore_index=True)


def main():
    out=ROOT/'Additional_Tables/Additional_file_1_Supplementary_Tables.xlsx'
    sheets=[]
    sheets.append(('S1_Panels',combined('panel_coverage.tsv')))
    import json
    rows=[]
    for p in sorted((ROOT/'Primary').glob('GSM*/qc.json'))+[ROOT/'External_EndToEnd/NYU_PDAC1/qc.json']:
        rows.append(json.loads(p.read_text()))
    sheets.append(('S2_Input_QC',pd.DataFrame(rows)))
    sheets.append(('S3_Associations',combined('Results/associations.tsv')))
    intervals=[]
    for dataset,folder in [('primary',ROOT/'Primary'),('external',ROOT/'External_EndToEnd/NYU_PDAC1')]:
        for path in sorted((folder/'Results').glob('bootstrap_*_ci.tsv')):
            frame=read(path);frame.insert(0,'dataset',dataset);intervals.append(frame)
    counts=combined('Results/candidate_counts.tsv')
    comparisons=combined('Results/composition_baseline_comparison.tsv')
    counts['record_type']='candidate_counts';comparisons['record_type']='ranking_comparison'
    sheets.append(('S4_Ranking_and_Counts',pd.concat([counts,comparisons],ignore_index=True)))
    sheets.append(('S5_Block_Intervals',pd.concat(intervals,ignore_index=True)))
    sheets.append(('S6_Matched_Summary',combined('Results/matched_summary.tsv')))
    sheets.append(('S7_Matched_Assignments',combined('Results/matched_assignments.tsv.gz')))
    sheets.append(('S8_Ranking_Simulation',read(ROOT/'Ranking_Ablation/all_replicates.tsv')))
    sheets.append(('S9_Paired_Differences',read(ROOT/'Ranking_Ablation/paired_differences.tsv')))
    sheets.append(('S10_Common_Gene_Check',read(ROOT/'External_EndToEnd/NYU_PDAC1/Results/common_gene_sensitivity_associations.tsv')))
    sheets.append(('S11_Model_Scope',read(ROOT/'Comparators/summary/primary_model_scope_audit.tsv')))
    sheets.append(('S12_Model_Enrichment',read(ROOT/'Comparators/current_summary/candidate_factor_enrichment.tsv')))
    book=Workbook();contents=book.active;contents.title='Contents'
    contents.append(['Sheet','Rows','Description'])
    notes=[
        'Declared/used genes; proxies are not validated cell fractions',
        'Counts/spatial QC; sample disease-identity conflict described in source audit',
        'Raw, technical and composition-proxy models, pooled and separate sections',
        'Candidate counts and ranking comparisons; outlier heuristic, no coordinate-level FDR',
        'Within-observed-section spatial block percentile intervals; not patient-population CIs',
        'Observed and random sets share scoring and residual model',
        'Strict matching, all assignments, seed and caliper differences',
        'All 50 paired replicates for compact/dispersed/no-program geometries',
        'All method-minus-product differences; no selected seeds',
        'Post hoc external 19/20-gene coverage sensitivity; not replacement primary',
        'New primary versus old model-input and valid joint coverage',
        'Updated labels on reused trained models; contextual overlap, not accuracy']
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
