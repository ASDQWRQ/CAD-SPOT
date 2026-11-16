#!/usr/bin/env python
"""Model-dependent top-budget coordinate audit; no automatic model selection."""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

import cadspot

MODELS = ('raw', 'technical', 'composition', 'spline')
CONDITIONAL = ('technical', 'composition', 'spline')
DEFAULT_BUDGETS = (10, 25, 50)


def model_scores(table, target, reference, technical, composition, sample_col='sample_id'):
    """Return positive products and residual pairs using existing core estimators."""
    scores, residuals, rhos = {}, {}, {}
    for model in MODELS:
        if model == 'raw':
            a, b = table[target].to_numpy(float), table[reference].to_numpy(float)
            rho = cadspot.correlation(a, b)
        else:
            cov = list(technical) if model == 'technical' else list(technical) + list(composition)
            estimator = 'spline_score_residual_spearman' if model == 'spline' else 'score_residual_spearman'
            rho, a, b = cadspot.fit_pair(table, target, reference, cov, sample_col, estimator)
        pair = pd.DataFrame({sample_col: table[sample_col], 'a': a, 'b': b})
        za = cadspot.z_by_sample(pair, 'a', sample_col).to_numpy(float)
        zb = cadspot.z_by_sample(pair, 'b', sample_col).to_numpy(float)
        scores[model] = cadspot.positive_product(za, zb)
        residuals[model] = np.column_stack([za, zb])
        rhos[model] = rho
    return scores, residuals, rhos


def budget_order(table, score, idx, id_col='barcode'):
    """Unique ordinal order; labels/truth are never used for tie breaking."""
    idx = np.asarray(idx, int)
    return idx[np.lexsort((table[id_col].to_numpy(dtype=str)[idx], -np.asarray(score)[idx]))]


def audit_table(table, target, reference, technical, composition, budgets=DEFAULT_BUDGETS,
                primary_budget=25, sample_col='sample_id', id_col='barcode', x_col='x', y_col='y'):
    budgets = tuple(sorted(set(budgets)))
    if not budgets or any(k <= 0 for k in budgets) or primary_budget not in budgets:
        raise ValueError('Positive budgets must include the primary budget')
    d = cadspot.validate_table(table, target, reference, technical, composition,
        sample_col, id_col, x_col, y_col, spatial=True)
    score, residuals, rhos = model_scores(d, target, reference, technical, composition, sample_col)
    out = d[[sample_col, id_col, x_col, y_col]].copy()
    valid = d[[x_col, y_col]].notna().all(axis=1).to_numpy()
    out['audit_spatial_status'] = np.where(valid, 'ok', 'excluded_missing_coordinates')
    for _, idx in d.groupby(sample_col, sort=True).indices.items():
        idx = np.asarray(idx)
        if valid[idx].sum() < 7:
            out.loc[idx[valid[idx]], 'audit_spatial_status'] = 'excluded_fewer_than_seven_positions'
    valid = out.audit_spatial_status.eq('ok').to_numpy()
    summary = []
    for model in MODELS:
        out[model + '_priority'] = score[model]
        out[model + '_target_z'] = residuals[model][:, 0]
        out[model + '_reference_z'] = residuals[model][:, 1]
        out[model + '_audit_rank'] = np.nan
        for _, idx in d.groupby(sample_col, sort=True).indices.items():
            idx = np.asarray(idx)[valid[idx]]
            ordered = budget_order(d, score[model], idx, id_col)
            out.loc[ordered, model + '_audit_rank'] = np.arange(1, len(ordered) + 1)
        for k in budgets:
            out[f'{model}_top{k}'] = out[model + '_audit_rank'].le(k) & (score[model] > 0)
    for k in budgets:
        votes = out[[f'{m}_top{k}' for m in CONDITIONAL]].sum(axis=1)
        out[f'conditional_votes_top{k}'] = votes
        out[f'agreement_top{k}'] = np.select([~valid, votes.eq(3), votes.gt(0)],
            ['excluded', 'all_three', 'some_models'], default='not_selected')
        out[f'technical_composition_changed_top{k}'] = out[f'technical_top{k}'] != out[f'composition_top{k}']
        out[f'ols_spline_changed_top{k}'] = out[f'composition_top{k}'] != out[f'spline_top{k}']
        for section, idx in d.groupby(sample_col, sort=True).indices.items():
            sub = out.iloc[idx]
            for model in MODELS:
                summary.append(dict(sample_id=section, budget=k, record='model_selection', model=model,
                    n_spots=len(idx), n_spatial_valid=int(valid[idx].sum()),
                    selected=int(sub[f'{model}_top{k}'].sum()), pooled_rho=rhos[model]))
            for category in ['all_three', 'some_models', 'not_selected', 'excluded']:
                summary.append(dict(sample_id=section, budget=k, record='model_agreement', model=category,
                    n_spots=len(idx), n_spatial_valid=int(valid[idx].sum()),
                    selected=int(sub[f'agreement_top{k}'].eq(category).sum()), pooled_rho=np.nan))
    out['primary_budget'] = primary_budget
    out['model_agreement'] = out[f'agreement_top{primary_budget}']
    return out, pd.DataFrame(summary)


def report_html(summary, manifest):
    primary = summary[summary.budget.eq(manifest['primary_budget'])]
    selections = primary[primary.record.eq('model_selection')].pivot(index='sample_id', columns='model', values='selected')
    agreement = primary[primary.record.eq('model_agreement')].pivot(index='sample_id', columns='model', values='selected')
    display = selections[list(MODELS)].join(agreement[['all_three', 'some_models']])
    display.columns = ['Raw', 'Technical', 'Composition', 'Spline', 'All three conditional', 'Some conditional']
    rows = primary[primary.record.eq('model_selection')].drop_duplicates('model')[['model', 'pooled_rho']]
    title = 'CAD-SPOT coordinate decision audit'
    return ('<!doctype html><html lang="en"><meta charset="utf-8"><title>' + title + '</title>'
        '<style>body{font:14px Arial,sans-serif;color:#111;max-width:1050px;margin:28px auto;padding:0 18px}'
        'h1{font-size:22px}h2{font-size:17px;margin-top:22px}p{line-height:1.45}'
        'table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:7px;border:1px solid #ccc;text-align:right}'
        'th:first-child,td:first-child{text-align:left}th{background:#edf1f3}'
        '@media print{body{margin:8mm;font-size:11px}h1{font-size:18px}table{font-size:10px}th,td{padding:4px}}</style>'
        '<h1>' + title + '</h1><p>Input: ' + html.escape(manifest['input_name']) +
        ' | Version: ' + html.escape(cadspot.__version__) + ' | Primary budget: ' + str(manifest['primary_budget']) +
        ' per section | Coordinates: ' + html.escape(manifest['coordinate_units']) + '</p>'
        '<h2>Priority membership</h2>' + display.to_html(escape=True, border=0) +
        '<h2>Pooled association record</h2>' + rows.to_html(index=False, escape=True, border=0, float_format=lambda x:f'{x:.4f}') +
        '<p>Membership requires a positive score and deterministic within-section rank within budget; lists can be shorter than the budget. '
        'All three means agreement of technical, composition and spline lists. Raw is a separate baseline. '
        'Agreement is not biological validation, coordinate-level FDR or a model-selection rule. '
        'Models condition on supplied covariates and may remove genuine composition-linked biology.</p>'
        '<p>Exact coordinates, all budget memberships, rank changes, covariate names and hashes are retained in the accompanying TSV/JSON files.</p></html>')


def export_audit(table, prefix, target, reference, technical, composition, input_name='in_memory',
                 input_sha256=None, coordinate_units='declared_by_user', budgets=DEFAULT_BUDGETS, primary_budget=25):
    out, summary = audit_table(table, target, reference, technical, composition, budgets, primary_budget)
    prefix = Path(prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    manifest = dict(version=cadspot.__version__, input_name=input_name, input_sha256=input_sha256,
        target=target, reference=reference, technical=list(technical), composition=list(composition),
        coordinate_units=coordinate_units, budgets=list(budgets), primary_budget=primary_budget,
        conditional_models=list(CONDITIONAL), positive_priority_required=True,
        rank_ties='ordinal descending priority, ascending barcode within section',
        meaning='descriptive model dependence, not FDR, accuracy or automatic model selection',
        core_sha256=cadspot.file_sha256(Path(cadspot.__file__)),
        audit_source_sha256=cadspot.file_sha256(Path(__file__)))
    out.to_csv(str(prefix)+'.coordinates.tsv.gz', sep='\t', index=False)
    summary.to_csv(str(prefix)+'.summary.tsv', sep='\t', index=False)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    Path(str(prefix)+'.report.html').write_text(report_html(summary, manifest), encoding='utf-8')
    return out, summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True); p.add_argument('--output-prefix', type=Path, required=True)
    p.add_argument('--target', required=True); p.add_argument('--reference', required=True)
    p.add_argument('--technical', nargs='*', default=[]); p.add_argument('--composition', nargs='*', default=[])
    p.add_argument('--budgets', type=int, nargs='+', default=list(DEFAULT_BUDGETS))
    p.add_argument('--primary-budget', type=int, default=25); p.add_argument('--coordinate-units', required=True)
    args = p.parse_args()
    export_audit(cadspot.read_table(args.input), args.output_prefix, args.target, args.reference,
        args.technical, args.composition, args.input.name, cadspot.file_sha256(args.input),
        args.coordinate_units, args.budgets, args.primary_budget)
    print(str(args.output_prefix)+'.report.html')


if __name__ == '__main__':
    main()
