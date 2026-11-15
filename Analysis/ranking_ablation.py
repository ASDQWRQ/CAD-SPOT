"""Paired ranking ablation with known compact/dispersed coupling, no tuning."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from sklearn.metrics import average_precision_score
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'Software'))
import cadspot as core

SEED = 20261115
N_REPS = 50
METHODS = ['positive_product', 'legacy', 'positive_support']


def standardize(x):
    return (x - np.mean(x)) / (np.std(x) or 1)


def generate(seed, geometry):
    rng = np.random.default_rng(seed)
    nr, nc = 40, 45
    yy, xx = np.mgrid[:nr, :nc]
    frames, truths = [], []
    for section in range(2):
        architecture = standardize((xx / nc + gaussian_filter(rng.normal(size=(nr, nc)), 4)).ravel())
        depth = .7 * architecture + .7 * standardize(rng.normal(size=nr*nc))
        stroma = -.4 * architecture + standardize(gaussian_filter(rng.normal(size=(nr, nc)), 3).ravel())
        noise_a = .8 * standardize(rng.normal(size=nr*nc)) + .2 * standardize(gaussian_filter(rng.normal(size=(nr, nc)), 1).ravel())
        noise_b = .8 * standardize(rng.normal(size=nr*nc)) + .2 * standardize(gaussian_filter(rng.normal(size=(nr, nc)), 1).ravel())
        cx, cy = rng.uniform(14, 30), rng.uniform(12, 27)
        truth = ((xx-cx)**2 + (yy-cy)**2 <= 7.5**2).ravel()
        if geometry == 'dispersed':
            truth = rng.permutation(truth)
        elif geometry == 'null':
            truth[:] = False
        signal = truth.astype(float)
        frame = pd.DataFrame(dict(sample_id=f's{section}', barcode=[f'b{i:04d}' for i in range(nr*nc)],
            x=xx.ravel().astype(float), y=yy.ravel().astype(float),
            depth=depth, architecture=architecture, stroma=stroma,
            target_score=1.2*architecture + .3*depth + .2*stroma + 2.0*signal + noise_a,
            reference_score=architecture + .2*depth + .5*stroma + 1.8*signal + noise_b))
        frames.append(frame)
        truths.append(truth)
    return pd.concat(frames, ignore_index=True), np.concatenate(truths)


def compute_scores(df):
    _, ra, rb = core.fit_pair(df, 'target_score', 'reference_score', ['depth', 'architecture', 'stroma'])
    a, b, h = np.zeros(len(df)), np.zeros(len(df)), np.zeros(len(df))
    for _, group in df.groupby('sample_id', sort=True):
        idx = group.index.to_numpy()
        a[idx], b[idx] = standardize(ra[idx]), standardize(rb[idx])
        xy = group[['x', 'y']].to_numpy()
        nn = NearestNeighbors(n_neighbors=7).fit(xy)
        distance, neighbors = nn.kneighbors(xy)
        ids = group.barcode.to_numpy()
        local = np.zeros(len(group))
        # Use the same all-boundary-ties and deterministic ID tie policy as core.
        for i in range(len(group)):
            candidates = nn.radius_neighbors(xy[i:i+1], radius=np.nextafter(distance[i, -1], np.inf), return_distance=False)[0]
            candidates = candidates[candidates != i]
            d2 = ((xy[candidates]-xy[i])**2).sum(axis=1)
            chosen = candidates[np.lexsort((ids[candidates], d2))[:6]]
            local[i] = b[idx[chosen]].mean()
        h[idx] = standardize(local)
    positive = (a > 0) & (b > 0)
    values = dict(positive_product=np.maximum(a, 0)*np.maximum(b, 0),
                  legacy=(a+b+h+a*h)/4,
                  positive_support=np.maximum(a, 0)*np.maximum(b, 0)*np.maximum(h, 0))
    return values, positive, h


def evaluate(df, truth, scores, positive, h, method):
    score = scores[method]
    eligible = positive & ((h > 0) if method != 'positive_product' else True)
    # Identical positive-quadrant universe for the continuous-score comparison;
    # candidate-rule differences are separately measured below.
    comparison_score = np.where(positive, score, min(-1., float(score.min())-1))
    ap = average_precision_score(truth, comparison_score) if truth.any() else np.nan
    budget = int(np.ceil(.05*len(df)))
    order = np.lexsort((df.barcode.to_numpy(), df.sample_id.to_numpy(), -comparison_score))[:budget]
    selected = np.zeros(len(df), dtype=bool)
    scale_fallbacks = 0
    for _, group in df.groupby('sample_id', sort=True):
        idx = group.index.to_numpy()
        s = score[idx]
        scale = float(1.4826*np.median(np.abs(s-np.median(s))))
        if scale <= 0:
            scale = float(np.std(s) or 1)
            scale_fallbacks += 1
        selected[idx] = ((s-np.median(s))/scale >= 3) & eligible[idx]
    return dict(average_precision=ap,
        top_5pct_precision=float(truth[order].mean()) if truth.any() else np.nan,
        top_5pct_recall=float(truth[order].sum()/truth.sum()) if truth.any() else np.nan,
        candidate_count=int(selected.sum()), candidate_fraction=float(selected.mean()),
        candidate_precision=float(truth[selected].mean()) if selected.any() and truth.any() else np.nan,
        candidate_recall=float(truth[selected].sum()/truth.sum()) if truth.any() else np.nan,
        sd_fallback_sections=scale_fallbacks, n_spots=len(df), n_true=int(truth.sum()))


def main():
    out = ROOT/'Ranking_Ablation'
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for geometry in ['compact', 'dispersed', 'null']:
        for rep in range(N_REPS):
            seed = SEED + rep
            df, truth = generate(seed, geometry)
            scores, positive, h = compute_scores(df)
            for method in METHODS:
                rows.append(dict(geometry=geometry, replicate=rep, seed=seed, method=method,
                                 **evaluate(df, truth, scores, positive, h, method)))
            if rep == 0:
                example = df.copy()
                example['truth'] = truth
                for method in METHODS:
                    example[method] = scores[method]
                example.to_csv(out/f'example_{geometry}.tsv.gz', sep='\t', index=False)
            if rep % 10 == 0:
                print(geometry, rep, flush=True)
        pd.DataFrame(rows).to_csv(out/'all_replicates.tsv', sep='\t', index=False)
    result = pd.DataFrame(rows)
    metrics = ['average_precision', 'top_5pct_precision', 'top_5pct_recall',
               'candidate_fraction', 'candidate_precision', 'candidate_recall']
    summary = result.groupby(['geometry', 'method'])[metrics].agg(['mean', 'std', 'median'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary.reset_index().to_csv(out/'summary.tsv', sep='\t', index=False)
    paired = []
    for geometry, group in result.groupby('geometry'):
        for method in METHODS[1:]:
            for metric in metrics:
                wide = group.pivot(index='replicate', columns='method', values=metric)
                delta = (wide[method]-wide['positive_product']).dropna()
                paired.append(dict(geometry=geometry, comparison=f'{method}-positive_product', metric=metric,
                    n_valid=len(delta), mean_difference=float(delta.mean()) if len(delta) else np.nan,
                    median_difference=float(delta.median()) if len(delta) else np.nan,
                    min_difference=float(delta.min()) if len(delta) else np.nan,
                    max_difference=float(delta.max()) if len(delta) else np.nan,
                    n_positive=int(delta.gt(0).sum()), n_negative=int(delta.lt(0).sum())))
    pd.DataFrame(paired).to_csv(out/'paired_differences.tsv', sep='\t', index=False)
    (out/'design.json').write_text(json.dumps(dict(seed=SEED, replicates_per_geometry=N_REPS,
        geometries=['compact','dispersed','null'], n_sections=2, grid=[40,45],
        signal_amplitude_target=2.0, signal_amplitude_reference=1.8,
        score_model='Known linear architecture/depth/stroma plus injected shared positive program and noise',
        comparisons='Paired seeds; all methods tested on every generated dataset',
        interpretation='Controlled algorithm behavior only, not independent tissue validation or general FDR control'), indent=2)+'\n')
    print(summary.to_string(), flush=True)


if __name__ == '__main__':
    main()
