"""Validate completed actual method exports and build one auditable comparison handoff."""
from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import common
import numpy as np
import pandas as pd


def parse_spahdmap_objectives(path: Path) -> pd.DataFrame:
    phase, rows = None, []
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('Step 2:'):
            phase = 'image_pretraining'
        elif line.startswith('Step 3:'):
            phase = 'graph_autoencoder'
        elif line.startswith('Step 5:'):
            phase = 'final_fusion'
        match = re.match(r'\[Iter:\s*(\d+)\s*/\s*(\d+)\]', line)
        if not match or phase is None:
            continue
        row = dict(phase=phase, iteration=int(match[1]), nominal_max_iterations=int(match[2]))
        for label, number in re.findall(r'([A-Za-z][A-Za-z ]*):\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line):
            label = label.strip().lower().replace(' ', '_')
            if label not in ['iter', 'epoch']:
                row[label] = float(number)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stamp', type=Path, default=common.ROOT / 'results/stamp/k10_seed20261115_retry_numpy2')
    p.add_argument('--popari', type=Path, default=common.ROOT / 'results/popari/k10_seed20261115_budget50')
    p.add_argument('--spahdmap', type=Path, default=common.ROOT / 'results/spahdmap/k10_seed20261115_production')
    p.add_argument('--inputs', type=Path, default=common.ROOT / 'inputs')
    p.add_argument('--out', type=Path, default=common.ROOT / 'summary_all')
    p.add_argument('--wait', action='store_true')
    p.add_argument('--timeout-minutes', type=float, default=45)
    args = p.parse_args()
    runs = dict(STAMP=args.stamp, POPARI=args.popari, SpaHDmap=args.spahdmap)
    deadline = time.monotonic()+args.timeout_minutes*60
    while True:
        statuses = {m: json.loads((root / 'run.json').read_text()) for m, root in runs.items()}
        if all(r['status'] == 'completed' for r in statuses.values()):
            break
        failed = {m: r['status'] for m, r in statuses.items() if r['status'] not in ['running', 'completed']}
        if failed:
            raise RuntimeError(f'Incomplete method evidence: {failed}')
        if not args.wait or time.monotonic() >= deadline:
            raise RuntimeError('Methods still running; no final comparison was claimed')
        time.sleep(5)
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.inputs / 'manifest.json').read_text())
    eligible = pd.read_csv(args.inputs / 'eligible_spots.tsv.gz', sep='\t')
    eligible['cadspot_candidate'] = common.bool_labels(eligible.cadspot_candidate)
    if len(eligible) != 23013 or int(eligible.cadspot_candidate.sum()) != 755:
        raise ValueError('Primary freeze mismatch; expected23013 eligible spots and755 candidates')
    fits, hashes = [], {}
    for method, root in runs.items():
        record = statuses[method]
        if record['arguments'].get('smoke'):
            raise ValueError('Smoke evidence is forbidden')
        if record['input_manifest_sha256'] != common.sha256(args.inputs / 'manifest.json'):
            raise ValueError('Input manifest differs between fitting and comparison')
        files = sorted(root.rglob('spot_factors.tsv.gz'))
        frame = pd.concat([pd.read_csv(f, sep='\t') for f in files], ignore_index=True)
        expected = eligible if method != 'SpaHDmap' else eligible.loc[eligible.sample_id == manifest['representative']]
        output_keys = set(map(tuple, frame[['sample_id', 'barcode']].to_numpy()))
        expected_keys = set(map(tuple, expected[['sample_id', 'barcode']].to_numpy()))
        if not output_keys <= expected_keys or (method != 'SpaHDmap' and output_keys != expected_keys):
            raise ValueError(f'Output universe mismatch: {method}')
        if method == 'SpaHDmap':
            raw = pd.read_csv(root / manifest['representative'] / 'raw_factor_scores.tsv.gz', sep='\t')
            if set(map(tuple, raw[['sample_id', 'barcode']].to_numpy())) != expected_keys:
                raise ValueError('SpaHDmap must retain all native raw vectors, including invalid ones')
        if frame.duplicated(['sample_id', 'barcode']).any():
            raise ValueError('Duplicate output join keys')
        factors = [c for c in frame if c.startswith('factor_')]
        values = frame[factors].to_numpy()
        if len(factors) != 10 or not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values.sum(axis=1), 1, atol=1e-7):
            raise ValueError(f'Invalid final factor export: {method}')
        valid_expected = expected.merge(frame[['sample_id', 'barcode']], on=['sample_id', 'barcode'], validate='one_to_one')
        fits.append(dict(method=method, sections=frame.sample_id.nunique(), n_spots=len(frame),
                         n_candidates=int(valid_expected.cadspot_candidate.sum()), n_input_spots=len(expected),
                         n_input_candidates=int(expected.cadspot_candidate.sum()), excluded_spots=len(expected)-len(frame),
                         excluded_candidates=int(expected.cadspot_candidate.sum()-valid_expected.cadspot_candidate.sum()), factor_count=len(factors),
                         seed=record['arguments']['seed'], official_commit=record['commit'],
                         fixed_spatial_budget=record.get('fixed_spatial_iteration_budget'),
                         convergence_claim=False if method == 'POPARI' else 'official_stopping_rule_not_global_optimum',
                         run_directory=str(root), checkpoint_only=method == 'POPARI'))
        for f in [root / 'run.json', *files]:
            hashes[str(f.relative_to(common.ROOT))] = common.sha256(f)
    subprocess.run([sys.executable, str(common.ROOT / 'enrichment.py'), '--inputs', str(args.inputs),
                    '--outputs', *[str(root) for root in runs.values()], '--out', str(args.out)], check=True)
    trajectory = parse_spahdmap_objectives(args.spahdmap / 'console.log')
    trajectory.to_csv(args.spahdmap / 'training_objective.tsv', sep='\t', index=False)
    if set(trajectory.phase) != {'image_pretraining', 'graph_autoencoder', 'final_fusion'}:
        raise ValueError('Missing a full SpaHDmap training phase')
    for phase, target in [('graph_autoencoder', 5000), ('final_fusion', 2000)]:
        if trajectory.loc[trajectory.phase == phase, 'iteration'].max() < target:
            raise ValueError(f'Incomplete official phase: {phase}')
    pd.DataFrame(fits).to_csv(args.out / 'fit_summary.tsv', sep='\t', index=False)
    result = pd.read_csv(args.out / 'candidate_factor_enrichment.tsv', sep='\t')
    result['domain_q_lt_0_05'] = result.q_block_domain < .05
    result['loading_q_lt_0_05'] = result.q_block_loading < .05
    sensitivity = result.groupby(['method', 'sample_id', 'block_width_pitches'], as_index=False).agg(
        n_factors=('factor', 'size'), enriched_domains=('domain_q_lt_0_05', 'sum'),
        enriched_loadings=('loading_q_lt_0_05', 'sum'),
        n_valid_domain_tests=('q_block_domain','count'), n_valid_loading_tests=('q_block_loading','count'),
        candidate_count=('n_candidates', 'first'),
        eligible_count=('n_eligible', 'first'), common_scope=('common_scope', 'first'), fit_note=('fit_note', 'first'))
    for count, valid in [('enriched_domains','n_valid_domain_tests'),('enriched_loadings','n_valid_loading_tests')]:
        sensitivity[count]=sensitivity[count].astype(float)
        sensitivity.loc[sensitivity[valid]==0,count]=np.nan
    sensitivity.to_csv(args.out / 'block_sensitivity_counts.tsv', sep='\t', index=False)
    sensitivity.loc[sensitivity.block_width_pitches == 8].to_csv(args.out / 'primary_enrichment_counts.tsv', sep='\t', index=False)
    for f in [*common.ROOT.glob('*.py'), common.ROOT / 'README.md', args.inputs / 'selected_genes.tsv', args.inputs / 'eligible_spots.tsv.gz',
              *args.out.glob('*.tsv')]:
        hashes[str(f.relative_to(common.ROOT))] = common.sha256(f)
    core = ['numpy', 'pandas', 'scipy', 'scikit-learn', 'torch']
    (common.ROOT / 'requirements-core.txt').write_text('\n'.join(f'{name}=={metadata.version(name)}' for name in core)+'\n', encoding='utf-8')
    staged = sorted({(d.metadata['Name'], d.version) for d in metadata.distributions(path=[str(common.DEPS)])})
    (common.ROOT / 'requirements-staged.txt').write_text('\n'.join(f'{name}=={version}' for name, version in staged)+'\n', encoding='utf-8')
    common.write_json(args.out / 'comparison_audit.json', dict(status='completed', primary_candidates=755, n_eligible=23013,
                      representative=manifest['representative'], no_eligible_attrition=False, smoke_evidence=False,
                      representative_joint_valid_spots=int(result.loc[result.sample_id == manifest['representative'], 'n_eligible'].iloc[0]),
                      representative_joint_valid_candidates=int(result.loc[result.sample_id == manifest['representative'], 'n_candidates'].iloc[0]),
                      popari_convergence_claim=False, popari_budget=50, methods=fits,
                      spatial_test='9999 equal-occupancy spatial block permutations; widths8(primary)/4/12; method-wise BH',
                      interpretation='Contextual overlap/enrichment, not winner ranking or independent biological validation',
                      utc_finalized=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), sha256=hashes))
    print(f'Final validated comparison: {args.out}', flush=True)


if __name__ == '__main__':
    main()
