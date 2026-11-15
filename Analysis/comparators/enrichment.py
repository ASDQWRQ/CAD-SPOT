"""Within-section candidate/domain overlap with occupancy-matched spatial block permutations."""
from __future__ import annotations

import argparse
from collections import defaultdict
import itertools
import json
from pathlib import Path
import common
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control


def block_groups(coords: np.ndarray, width: float) -> tuple[list[np.ndarray], dict]:
    pitch = float(np.median(cKDTree(coords).query(coords, k=2)[0][:, 1]))
    if not np.isfinite(pitch) or pitch <= 0:
        raise ValueError('Invalid spot pitch')
    cell = np.floor((coords - coords.min(axis=0)) / (pitch * width)).astype(int)
    blocks = defaultdict(list)
    for i, key in enumerate(map(tuple, cell)):
        blocks[key].append(i)
    by_size = defaultdict(list)
    for indices in blocks.values():
        indices = np.asarray(indices)
        indices = indices[np.lexsort((coords[indices, 0], coords[indices, 1]))]
        by_size[len(indices)].append(indices)
    groups = [np.stack(group) for group in by_size.values() if len(group) >= 2]
    movable = sum(g.size for g in groups)
    meta = dict(n_blocks=len(blocks), n_exchangeable_blocks=sum(len(g) for g in groups),
                movable_fraction=movable/len(coords), block_width_pitches=width,
                spot_pitch=pitch, permutation_null='equal-occupancy grid-block exchange; within-block spatial ordering preserved')
    return groups, meta


def permute_labels(labels: np.ndarray, groups: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    shuffled = labels.copy()
    for g in groups:
        shuffled[g.ravel()] = labels[g[rng.permutation(len(g))].ravel()]
    return shuffled


def section_stats(frame: pd.DataFrame, factors: list[str], width: float, n_perm: int, seed: int) -> pd.DataFrame:
    coords = frame[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy()
    y = frame.cadspot_candidate.to_numpy(dtype=bool)
    weights = frame[factors].to_numpy(dtype=float)
    n, nc = len(y), int(y.sum())
    dominant = weights.argmax(axis=1)
    indicators = np.eye(len(factors))[dominant]
    groups, meta = block_groups(coords, width)
    matrix = np.column_stack([indicators, weights])
    observed = y.astype(float) @ matrix
    null_exceed = np.zeros(matrix.shape[1], dtype=int)
    null_total = np.zeros(matrix.shape[1], dtype=float)
    testable = 0 < nc < n and meta['n_exchangeable_blocks'] >= 8 and meta['movable_fraction'] >= 0.5
    if testable:
        rng = np.random.default_rng(seed)
        for _ in range(n_perm):
            yp = permute_labels(y, groups, rng)
            permuted_stat = yp.astype(float) @ matrix
            null_exceed += permuted_stat >= (observed-1e-12)
            null_total += permuted_stat
        pvals = (null_exceed+1)/(n_perm+1)
    else:
        pvals = np.full(matrix.shape[1], np.nan)
    rows = []
    for j, name in enumerate(factors):
        inside = dominant == j
        nd, overlap = int(inside.sum()), int((inside & y).sum())
        expected = nc*nd/n
        b, c, d = nd-overlap, nc-overlap, n-nd-nc+overlap
        rows.append(dict(factor=name, n_eligible=n, n_candidates=nc, n_domain=nd, candidate_overlap=overlap,
                 expected_overlap=expected, fold_enrichment=overlap/expected if expected else np.nan,
                 block_null_expected_overlap=null_total[j]/n_perm if testable else np.nan,
                 odds_ratio_haldane=(overlap+.5)*(d+.5)/((b+.5)*(c+.5)),
                 candidate_fraction_in_domain=overlap/nd if nd else np.nan,
                 candidate_recall=overlap/nc if nc else np.nan,
                 jaccard=overlap/(nd+nc-overlap) if nd+nc-overlap else np.nan,
                 candidate_mean_loading=float(weights[y, j].mean()) if nc else np.nan,
                 noncandidate_mean_loading=float(weights[~y, j].mean()) if nc<n else np.nan,
                 mean_loading_difference=float(weights[y, j].mean()-weights[~y, j].mean()) if 0<nc<n else np.nan,
                 p_block_domain=pvals[j], p_block_loading=pvals[len(factors)+j],
                 test_status='approximate_block_null' if testable else 'insufficient_exchangeable_blocks_or_labels',
                 n_permutations=n_perm if testable else 0, **meta))
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, default=common.ROOT / 'inputs')
    p.add_argument('--outputs', type=Path, nargs='+', required=True, help='Method run roots; supports common-universe intersection')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--n-perm', type=int, default=9999)
    p.add_argument('--seed', type=int, default=20261115)
    p.add_argument('--block-widths', type=float, nargs='+', default=[8, 4, 12])
    p.add_argument('--allow-partial', action='store_true', help='Explicitly analyze valid completed sections from an unfinished/failed run')
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    eligible = pd.read_csv(args.inputs / 'eligible_spots.tsv.gz', sep='\t')
    eligible['cadspot_candidate'] = common.bool_labels(eligible.cadspot_candidate)
    loaded = []
    provenance = {}
    attrition, rows, overlaps = [], [], []
    for root in args.outputs:
        status_path = root / 'run.json'
        if not status_path.is_file():
            raise ValueError(f'Missing official run provenance: {status_path}')
        status = json.loads(status_path.read_text())
        if status['status'] != 'completed' and not args.allow_partial:
            raise ValueError(f'Run not completed: {root}; use --allow-partial only with explicit reporting')
        files = sorted(root.rglob('spot_factors.tsv.gz'))
        if not files:
            raise ValueError(f'No completed exports in {root}')
        frame = pd.concat([pd.read_csv(f, sep='\t') for f in files], ignore_index=True)
        if common.bool_labels(frame.smoke_test).any():
            raise ValueError('Smoke-test exports are not final enrichment evidence')
        if frame.duplicated(['sample_id', 'barcode']).any() or frame.method.nunique() != 1:
            raise ValueError('Root must contain one method, one seed, no duplicate section/barcode')
        factors = [c for c in frame if c.startswith('factor_')]
        if not np.isfinite(frame[factors]).all().all():
            raise ValueError('Nonfinite model output')
        method = str(frame.method.iloc[0])
        budget = status.get('fixed_spatial_iteration_budget')
        provenance[method] = dict(fit_note=f'fixed{budget}_spatial_updates_NOT_claimed_converged' if budget else 'official_training_pipeline_and_logged_stopping_rule',
                                  fit_spatial_iteration_budget=budget, fit_run_path=str(root),
                                  fit_source_commit=status.get('commit'), fit_status=status['status'])
        loaded.append((method, frame, factors))
    all_sections = sorted(set.union(*(set(f.sample_id) for _, f, _ in loaded)))
    for sample in all_sections:
        available = [(m, f.loc[f.sample_id == sample].copy(), factors) for m, f, factors in loaded if sample in set(f.sample_id)]
        u = eligible.loc[eligible.sample_id == sample].copy()
        keys = set(u.barcode)
        for _, f, _ in available:
            keys &= set(f.barcode)
        u = u.loc[u.barcode.isin(keys)].sort_values('barcode').reset_index(drop=True)
        if not len(u):
            raise ValueError(f'Empty common eligible universe for {sample}')
        scope = '+'.join(m for m, _, _ in available)
        u.to_csv(args.out / f'{sample}_common_eligible.tsv.gz', sep='\t', index=False)
        aligned = []
        for method, f, factors in available:
            joined = u.merge(f, on=['sample_id', 'barcode'], validate='one_to_one')
            native = joined[factors].to_numpy()
            domain = native.argmax(axis=1)
            candidate = joined.cadspot_candidate.to_numpy(dtype=bool)
            freq_all = np.bincount(domain, minlength=len(factors))/len(domain)
            freq_cand = np.bincount(domain[candidate], minlength=len(factors))/max(int(candidate.sum()), 1)
            def entropy(freq: np.ndarray) -> float:
                positive = freq[freq > 0]
                return float(-np.sum(positive*np.log(positive))/np.log(len(factors)))
            attrition.append(dict(sample_id=sample, method=method, common_scope=scope, n_frozen_eligible=len(eligible.loc[eligible.sample_id == sample]),
                                  n_method_output=len(f), n_common=len(joined), n_candidates_common=int(joined.cadspot_candidate.sum()),
                                  n_candidates_frozen=int(eligible.loc[eligible.sample_id == sample, 'cadspot_candidate'].sum()),
                                  n_spots_excluded_from_common=len(eligible.loc[eligible.sample_id == sample])-len(joined),
                                  n_candidates_excluded_from_common=int(eligible.loc[eligible.sample_id == sample, 'cadspot_candidate'].sum()-joined.cadspot_candidate.sum()),
                                  n_domains_observed=int((freq_all>0).sum()), n_domains_with_candidates=int((freq_cand>0).sum()),
                                  dominant_domain_entropy_all=entropy(freq_all), dominant_domain_entropy_candidates=entropy(freq_cand),
                                  candidate_fraction=int(candidate.sum())/len(candidate), **provenance[method]))
            for width in args.block_widths:
                result = section_stats(joined, factors, width, args.n_perm, args.seed)
                result.insert(0, 'sample_id', sample)
                result.insert(0, 'method', method)
                result['common_scope'] = scope
                result['n_frozen_eligible'] = len(eligible.loc[eligible.sample_id == sample])
                result['n_frozen_candidates'] = int(eligible.loc[eligible.sample_id == sample, 'cadspot_candidate'].sum())
                result['excluded_spots_from_common'] = result['n_frozen_eligible']-len(joined)
                result['excluded_candidates_from_common'] = result['n_frozen_candidates']-int(joined.cadspot_candidate.sum())
                for key, value in provenance[method].items():
                    result[key] = value
                rows.append(result)
            aligned.append((method, joined[factors].to_numpy().argmax(axis=1)))
        for (m1, d1), (m2, d2) in itertools.combinations(aligned, 2):
            for i in np.unique(d1):
                for j in np.unique(d2):
                    left, right = d1 == i, d2 == j
                    both = int((left & right).sum())
                    overlaps.append(dict(sample_id=sample, method1=m1, domain1=int(i+1), method2=m2, domain2=int(j+1),
                                    n_common=len(u), intersection=both, jaccard=both/int((left | right).sum()),
                                    candidate_overlap=int((left & right & u.cadspot_candidate.to_numpy()).sum())))
    result = pd.concat(rows, ignore_index=True)
    for width in args.block_widths:
        for method in result.method.unique():
            ix = (result.block_width_pitches == width) & (result.method == method)
            values = result.loc[ix, ['p_block_domain', 'p_block_loading']].to_numpy()
            finite = np.isfinite(values)
            adjusted = np.full_like(values, np.nan)
            if finite.any():
                adjusted[finite] = false_discovery_control(values[finite])
            result.loc[ix, 'q_block_domain'] = adjusted[:, 0]
            result.loc[ix, 'q_block_loading'] = adjusted[:, 1]
    result.to_csv(args.out / 'candidate_factor_enrichment.tsv', sep='\t', index=False)
    pd.DataFrame(attrition).to_csv(args.out / 'within_section_summary.tsv', sep='\t', index=False)
    pd.DataFrame(overlaps).to_csv(args.out / 'pairwise_domain_overlap.tsv', sep='\t', index=False)
    common.write_json(args.out / 'analysis_settings.json', dict(arguments=vars(args),
                      primary_block_width=8, alternative_widths=[4, 12], p_alternative='one-sided enrichment',
                      multiple_testing='BH jointly over domain and loading tests per method and block-width across sections',
                      descriptive_expected_overlap='candidate fraction times domain size; not the conditional block-null expectation',
                      biological_scope='sections retained separately; independent STAMP factors not aligned',
                      eligible_sha256=common.sha256(args.inputs / 'eligible_spots.tsv.gz')))
    print(result[['method', 'sample_id', 'n_eligible', 'n_candidates', 'movable_fraction', 'test_status']].drop_duplicates().to_string(index=False))


if __name__ == '__main__':
    main()
