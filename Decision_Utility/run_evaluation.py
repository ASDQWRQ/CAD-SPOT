"""Fixed-budget recovery and perturbation records; settings in PROTOCOL.md."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'Software'))
sys.path.insert(0, str(ROOT/'Utility_Benchmark'))
import cadspot
import decision_audit as audit
import benchmark as old

OUT = ROOT/'Decision_Utility/Results'
BUDGETS = (10, 25, 50)
TECH, COMP = old.TECH, old.COMP


def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=[''])


def save(d, name):
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(d).to_csv(OUT/name, sep='\t', index=False)


def instances():
    for case in old.CASES:
        for rep in range(40):
            d, truth = old.generate(881001+rep, case)
            yield 'simulation', case, rep, d, truth, ['depth'], ['composition', 'stroma']
    for background in ['ovarian', 'pdac']:
        d = read(ROOT/f'Utility_Benchmark/Inputs/{background}.tsv.gz')
        names = sorted(d.sample_id.unique())[:2]
        d = pd.concat([sub.sort_values('barcode').iloc[::3] for name, sub in d.groupby('sample_id')
                       if name in names], ignore_index=True)
        d[['target_score', 'reference_score']] = old.standardize_sections(d, d[['target_score', 'reference_score']].to_numpy())
        for geometry in ['compact', 'dispersed']:
            for rep in range(32):
                rng = np.random.default_rng(771001+rep); truth = np.zeros(len(d), bool)
                for idx in d.groupby('sample_id', sort=True).indices.values():
                    k = int(np.ceil(.05*len(idx)))
                    if geometry == 'compact':
                        center = int(rng.choice(idx)); xy = d[['x', 'y']].to_numpy()
                        distances = ((xy[idx]-xy[center])**2).sum(axis=1)
                        chosen = idx[np.lexsort((d.barcode.to_numpy()[idx], distances))[:k]]
                    else:
                        chosen = rng.choice(idx, k, replace=False)
                    truth[chosen] = True
                table = d.copy(); table.target_score += 1.2*truth; table.reference_score += 1.1*truth
                yield 'injection', f'{background}_{geometry}', rep, table, truth, TECH, COMP


def known_signal():
    records, costs = [], []
    expected = pd.concat([read(ROOT/f'Utility_Benchmark/Results/{v}_replicates.tsv')
                          for v in ['simulation', 'injection']], ignore_index=True)
    expected = expected.set_index(['experiment', 'scenario', 'replicate', 'method'])
    for iteration, (experiment, scenario, rep, d, truth, tech, comp) in enumerate(instances()):
        scores = old.scores(d, tech, comp, cadspot)
        for method, values in scores.items():
            ap = old.metrics(d, truth, values)['average_precision']
            prior = expected.loc[(experiment, scenario, rep, method), 'average_precision']
            if not (pd.isna(ap) and pd.isna(prior)):
                assert np.isclose(ap, prior, atol=1e-12), (scenario, rep, method)
            for section, idx in d.groupby('sample_id', sort=True).indices.items():
                order = audit.budget_order(d, values, idx)
                ntrue = int(truth[idx].sum())
                for k in BUDGETS:
                    selected = order[:k]
                    recovered = int(truth[selected].sum())
                    records.append(dict(experiment=experiment, scenario=scenario, replicate=rep,
                        section=section, method=method, budget=k, inspected=len(selected), n_true=ntrue,
                        recovered=recovered, precision=recovered/len(selected),
                        recall=recovered/ntrue if ntrue else np.nan,
                        positive_score_selected=int((values[selected] > 0).sum())))
                required = int(np.ceil(.8*ntrue))
                needed = int(np.searchsorted(np.cumsum(truth[order]), required)+1) if ntrue else np.nan
                costs.append(dict(experiment=experiment, scenario=scenario, replicate=rep,
                    section=section, method=method, n_true=ntrue, target_recall=.8, inspections=needed))
        if iteration % 40 == 0:
            print('fixed budget', iteration, experiment, scenario, flush=True)
    save(records, 'fixed_budget_replicates.tsv.gz'); save(costs, 'inspection_cost_replicates.tsv.gz')
    d = pd.DataFrame(records)
    group = ['experiment', 'scenario', 'replicate', 'method', 'budget']
    paired = d.groupby(group)[['precision', 'recall']].mean().reset_index()
    summary = paired.groupby(['experiment', 'scenario', 'method', 'budget'])[['precision', 'recall']].agg(['mean', 'std'])
    summary.columns = ['_'.join(v) for v in summary.columns]
    save(summary.reset_index(), 'fixed_budget_summary.tsv')
    costs = pd.DataFrame(costs).groupby(['experiment', 'scenario', 'replicate', 'method']).inspections.mean().reset_index()
    save(costs, 'inspection_cost_paired.tsv')
    diffs = []; rng = np.random.default_rng(20261117)
    comparisons = [(paired, 'precision', ['experiment', 'scenario', 'budget']),
                   (costs, 'inspections', ['experiment', 'scenario'])]
    for table, metric, groups in comparisons:
        for keys, sub in table.groupby(groups, sort=True):
            keys = keys if isinstance(keys, tuple) else (keys,)
            wide = sub.pivot(index='replicate', columns='method', values=metric)
            for baseline in [v for v in old.METHODS if v != 'cadspot_product']:
                delta = (wide.cadspot_product-wide[baseline]).dropna().to_numpy()
                draws = rng.choice(delta, (10000, len(delta)), replace=True).mean(axis=1) if len(delta) else np.array([np.nan])
                diffs.append(dict(zip(groups, keys), metric=metric,
                    comparison='cadspot_product-minus-'+baseline, n_valid=len(delta),
                    mean_difference=delta.mean() if len(delta) else np.nan,
                    ci_low=np.quantile(draws, .025), ci_high=np.quantile(draws, .975)))
    save(diffs, 'fixed_budget_paired_differences.tsv')
    print('fixed-budget outputs', len(records), len(costs), flush=True)


def cohort_inputs(folder):
    d = read(folder/'score_table.tsv.gz')
    selected = sorted(d.sample_id.unique())[:2]
    parts = []
    for sample in selected:
        f = folder/sample
        spots = read(f/'score_table.tsv.gz')
        counts = sparse.load_npz(f/'counts.npz').astype(np.int64)
        genes = pd.Index(read(f/'genes.tsv').gene)
        cov = read(f/'panel_coverage.tsv')
        panels = {p:sub.loc[sub.status.eq('used'), 'gene'].tolist() for p, sub in cov.groupby('panel')}
        parts.append((spots, counts, genes, panels))
    return parts


def rescore(spots, counts, genes, panels):
    d = spots.copy()
    lib = np.asarray(counts.sum(axis=1)).ravel()
    detect = np.diff(counts.indptr)
    if (lib <= 0).any():
        raise ValueError('Zero library after perturbation; no silent row exclusion')
    x = counts.multiply((10000/lib)[:, None]).tocsr().astype(float)
    x.data = np.log1p(x.data)
    for panel, members in panels.items():
        values = x[:, genes.get_indexer(members)].toarray()
        sd = values.std(axis=0)
        if (sd <= 1e-12).any():
            failed = [g for g, value in zip(members, sd) if value <= 1e-12]
            raise ValueError('Fixed panel constant: '+str(spots.sample_id.iloc[0])+'/'+panel+'/'+','.join(failed))
        d[panel] = ((values-values.mean(axis=0))/sd).mean(axis=1)
    d['nCount'] = lib; d['nFeature'] = detect
    d['log_nCount'] = np.log1p(lib); d['log_nFeature'] = np.log1p(detect)
    return d, int(((lib < 100) | (detect < 100)).sum())


def measure_perturbation(baseline, changed, label, perturbation, rep, records, strata):
    assert list(zip(baseline.sample_id, baseline.barcode)) == list(zip(changed.sample_id, changed.barcode))
    for section, idx in baseline.groupby('sample_id', sort=True).indices.items():
        a, b = baseline.iloc[idx], changed.iloc[idx]
        for model in audit.MODELS:
            rho = cadspot.correlation(a[model+'_audit_rank'], b[model+'_audit_rank'], rank=False)
            for k in BUDGETS:
                # Evaluate fixed-budget rank retention even if positive audit membership is shorter.
                left = a[model+'_audit_rank'].le(k).to_numpy(); right = b[model+'_audit_rank'].le(k).to_numpy()
                shared, union = int((left & right).sum()), int((left | right).sum())
                records.append(dict(cohort=label, perturbation=perturbation, replicate=rep,
                    section=section, model=model, budget=k, baseline_selected=int(left.sum()),
                    retained=shared, retention=shared/left.sum(), jaccard=shared/union, rank_rho=rho))
        for k in BUDGETS:
            for group in ['all_three', 'some_models']:
                mask = a[f'composition_top{k}'] & a[f'agreement_top{k}'].eq(group)
                n = int(mask.sum())
                retained = int((mask.to_numpy() & b[f'composition_top{k}'].to_numpy()).sum())
                strata.append(dict(cohort=label, perturbation=perturbation, replicate=rep,
                    section=section, budget=k, agreement=group, baseline_n=n, retained=retained,
                    retention=retained/n if n else np.nan,
                    baseline_mean_composition_rank=a.loc[mask, 'composition_audit_rank'].mean()))


def perturbations():
    records, strata, qc, provenance = [], [], [], []
    for label, folder in [('ovarian', ROOT/'Primary'), ('pdac', ROOT/'PDAC_Transfer')]:
        parts = cohort_inputs(folder)
        original = pd.concat([p[0] for p in parts], ignore_index=True)
        for spots, counts, genes, panels in parts:
            rebuilt, _ = rescore(spots, counts, genes, panels)
            np.testing.assert_allclose(rebuilt[list(panels)], spots[list(panels)], atol=1e-10)
            provenance.append(dict(cohort=label, section=spots.sample_id.iloc[0],
                count_sha256=cadspot.file_sha256(folder/spots.sample_id.iloc[0]/'counts.npz')))
        baseline, _ = audit.audit_table(original, 'target_score', 'reference_score', TECH, COMP)
        for rep in range(20):
            for name, retention in [('count_80pct', .8), ('count_50pct', .5), ('panel_loss_10pct', None)]:
                rng = np.random.default_rng(2026111700+rep)
                frames = []; below = 0
                drops = {}
                if retention is None:
                    for panel, members in parts[0][3].items():
                        if len(members)-int(np.ceil(.1*len(members))) < 3:
                            raise ValueError('Prespecified gene loss leaves fewer than three genes')
                        drops[panel] = set(rng.choice(members, int(np.ceil(.1*len(members))), replace=False))
                status, failure = 'ok', ''
                try:
                    for spots, counts, genes, panels in parts:
                        current = counts.copy()
                        if retention is not None:
                            current.data = rng.binomial(current.data, retention)
                            current.eliminate_zeros()
                        used = {p:[g for g in v if g not in drops.get(p, set())] for p, v in panels.items()}
                        d, nbelow = rescore(spots, current, genes, used); below += nbelow; frames.append(d)
                    d = pd.concat(frames, ignore_index=True)
                    changed, _ = audit.audit_table(d, 'target_score', 'reference_score', TECH, COMP)
                    measure_perturbation(baseline, changed, label, name, rep, records, strata)
                except ValueError as error:
                    if not str(error).startswith(('Fixed panel constant:', 'Zero library after perturbation')):
                        raise
                    status, failure = 'not_computable', str(error)
                qc.append(dict(cohort=label, perturbation=name, replicate=rep, seed=2026111700+rep,
                    n_spots=len(original), below_original_qc=below if status == 'ok' else np.nan, qc_rows_removed=0,
                    status=status, failure_reason=failure,
                    dropped_genes=json.dumps({p:sorted(v) for p, v in drops.items()}, sort_keys=True)))
            if rep % 5 == 0:
                print('perturbation', label, rep, flush=True)
        for omitted in COMP:
            changed, _ = audit.audit_table(original, 'target_score', 'reference_score', TECH, [p for p in COMP if p != omitted])
            measure_perturbation(baseline, changed, label, 'omit_'+omitted, 0, records, strata)
        save(records, 'perturbation_replicates.tsv.gz'); save(strata, 'stratum_retention.tsv.gz'); save(qc, 'perturbation_qc.tsv')
    (OUT/'perturbation_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


def main():
    p = argparse.ArgumentParser(); p.add_argument('stage', choices=['budget', 'perturbations']); args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == 'budget': known_signal()
    else: perturbations()
    print('Completed', args.stage, flush=True)


if __name__ == '__main__':
    main()
