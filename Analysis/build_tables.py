"""Build ordered, editable supplementary tables from the actual revision outputs."""
import json
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'Source_Data'


def read(name):
    return pd.read_csv(SRC / name, sep='\t')


def build():
    target = set(read('target_genes.tsv').gene)
    reference = set(read('reference_genes.tsv').gene)
    available_t = set(read('target_available_genes.tsv').gene)
    available_r = set(read('reference_available_genes.tsv').gene)
    coverage = pd.DataFrame([
        dict(gene=g, panel=panel, intended=True, available_in_pdac_matrix=g in available,
             matrix_audit_action='included' if g in available else 'excluded_not_zero_imputed')
        for panel, genes, available in [('target', target, available_t), ('reference', reference, available_r)]
        for g in sorted(genes)
    ])
    intervals = pd.concat([read(f'block_bootstrap_{w}_ci.tsv') for w in (1500, 3000, 4500)], ignore_index=True)
    coordinates = read('fresh_coordinates.tsv.gz')
    candidates = coordinates.groupby('sample_id').agg(n_spots=('barcode', 'size'), n_candidates=('cadspot_candidate', 'sum')).reset_index()
    candidates['candidate_percent'] = 100 * candidates.n_candidates / candidates.n_spots
    candidates['source'] = 'fresh_OLS_score_residuals_six_neighbours_robust_Z_ge_3'
    comparator_dir = ROOT / 'Analysis/comparators'
    comparator_rows = []
    for path in sorted((comparator_dir / 'results').rglob('run.json')):
        record = json.loads(path.read_text(encoding='utf-8'))
        args = record.get('arguments', {})
        comparator_rows.append(dict(method=args.get('method', path.parent.parent.name),
            run=path.parent.name, status=record.get('status', ''), commit=record.get('commit', ''),
            seed=args.get('seed'), K=args.get('k'), device=args.get('device'),
            details=record.get('stop_reason', record.get('error', ''))))
    enrichment = []
    for path in sorted(comparator_dir.glob('summary_*/candidate_factor_enrichment.tsv')):
        d = pd.read_csv(path, sep='\t')
        d.insert(0, 'source_summary', path.parent.name)
        enrichment.append(d)
    final_summary = comparator_dir/'summary_all'
    comparator_summary = (pd.read_csv(final_summary/'fit_summary.tsv',sep='\t')
                          if (final_summary/'fit_summary.tsv').exists() else pd.DataFrame(comparator_rows))
    comparator_enrichment = (pd.read_csv(final_summary/'candidate_factor_enrichment.tsv',sep='\t')
                            if (final_summary/'candidate_factor_enrichment.tsv').exists()
                            else pd.concat(enrichment, ignore_index=True) if enrichment else pd.DataFrame({'status': ['not_yet_available']}))
    tables = [
        ('S1_GeneCoverage', 'Intended and measured gene sets; raw-matrix audit only', coverage),
        ('S2_InputAudit', 'Raw barcode, image-coordinate and count reconciliation', read('raw_barcode_coordinate_audit.tsv')),
        ('S3_AssociationLayers', 'Frozen-score models and sensitivity estimators', read('fresh_associations.tsv')),
        ('S4_BlockIntervals', 'Section-stratified spatial block bootstrap, all widths', intervals),
        ('S5_RankingComparisons', 'Signed residual-product and score-product comparisons on declared common universes', read('fresh_ranking_baseline.tsv')),
        ('S6_CandidateCounts', 'Fresh candidates from consistent OLS residuals', candidates),
        ('S7_MatchedNull', 'Observed and 500-control statistics from the separate 19/20-gene matrix audit', read('fresh_matched_null_500_summary.tsv')),
        ('S8_MatchedAssignments', 'All target-slot assignments and caliper differences for 500 PDAC sets', read('fresh_matched_sets_500.tsv')),
        ('S9_ExternalAudit', 'Ten human samples; technical-only models and 500 controls each', read('external_reanalysis_summary.tsv')),
        ('S10_ComparatorSummary', 'Completed and diagnostic run status; final reported checkpoints identified in README', comparator_summary),
        ('S11_ComparatorEnrichment', 'Candidate enrichment in learned factors/domains; approximate block permutations', comparator_enrichment),
        ('S12_ProvenanceChanges', 'Prior archived candidate list versus fresh residual computation', read('archived_candidate_comparison.tsv')),
    ]
    dest = ROOT / 'Additional_Tables'
    dest.mkdir(exist_ok=True)
    output = dest / 'Additional_file_1_Supplementary_Tables.xlsx'
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        pd.DataFrame([{'Table': name.split('_', 1)[0], 'Worksheet': name, 'Description': desc, 'Rows': len(df)}
                      for name, desc, df in tables]).to_excel(writer, sheet_name='Contents', index=False)
        for name, _, frame in tables:
            frame.to_excel(writer, sheet_name=name, index=False)
        for ws in writer.book.worksheets:
            ws.freeze_panes = 'A2'
            ws.auto_filter.ref = ws.dimensions
            ws.sheet_view.showGridLines = False
            for cell in ws[1]:
                cell.font = Font(name='Arial', size=10, bold=True, color='FFFFFF')
                cell.fill = PatternFill('solid', fgColor='344A56')
                cell.alignment = Alignment(wrap_text=True, vertical='center')
            ws.row_dimensions[1].height = 40
            for column in ws.columns:
                preview = [str(cell.value or '') for cell in list(column)[:150]]
                width = min(55, max(12, max(map(len, preview), default=12) + 2))
                ws.column_dimensions[get_column_letter(column[0].column)].width = width
            for row in ws.iter_rows(min_row=2):
                for cell in row:
                    cell.font = Font(name='Arial', size=10)
                    if isinstance(cell.value, float):cell.number_format = '0.000000'
    print(output, flush=True)


if __name__ == '__main__':
    build()
