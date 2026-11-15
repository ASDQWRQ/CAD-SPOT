"""Postprocess already-completed official SpaHDmap training; never rerun a model."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import common
import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from finalize_comparison import parse_spahdmap_objectives


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, default=common.ROOT / 'results/spahdmap/k10_seed20261115_production')
    p.add_argument('--manifest', type=Path, default=common.ROOT / 'inputs/manifest.json')
    args = p.parse_args()
    record = json.loads((args.run / 'run.json').read_text())
    if record['arguments'].get('smoke'):
        raise ValueError('Cannot promote smoke training to final evidence')
    manifest = json.loads(args.manifest.read_text())
    s = next(s for s in manifest['sections'] if s['sample_id'] == manifest['representative'])
    native = args.run / 'official_outputs' / s['sample_id'] / 'SpaHDmap/SpaHDmap_spot_score.npy'
    checkpoint = args.run / 'official_outputs/models/trained_model.pth'
    if not native.is_file() or not checkpoint.is_file():
        raise FileNotFoundError('Final official fusion checkpoint and native scores are required')
    trajectory = parse_spahdmap_objectives(args.run / 'console.log')
    for phase, target in [('graph_autoencoder', 5000), ('final_fusion', 2000)]:
        if trajectory.loc[trajectory.phase == phase, 'iteration'].max() < target:
            raise ValueError('Full training phase not demonstrated in original objective log')
    trajectory.to_csv(args.run / 'training_objective.tsv', sep='\t', index=False)
    a = ad.read_h5ad(s['h5ad'])
    if common.sha256(Path(s['h5ad'])) != record['input_h5ad_sha256'][s['sample_id']]:
        raise ValueError('Training input changed')
    z = np.load(native)
    valid = common.export_scores(args.run / s['sample_id'], a, z, 'SpaHDmap', record['arguments']['seed'])
    labels = KMeans(n_clusters=record['arguments']['k'], n_init=20, random_state=record['arguments']['seed']).fit_predict(z[valid])
    pd.DataFrame({'sample_id': s['sample_id'], 'barcode': a.obs.barcode.to_numpy()[valid],
                  'kmeans_domain': labels+1}).to_csv(args.run / s['sample_id'] / 'downstream_domains.tsv', sep='\t', index=False)
    candidate = common.bool_labels(a.obs.cadspot_candidate)
    coverage = dict(n_input_spots=a.n_obs, n_valid_spots=int(valid.sum()), n_zero_vectors=int((~valid).sum()),
                    n_input_candidates=int(candidate.sum()), n_valid_candidates=int(candidate[valid].sum()),
                    n_excluded_candidates=int(candidate[~valid].sum()), native_score_sha256=common.sha256(native),
                    final_checkpoint_sha256=common.sha256(checkpoint), no_pseudocount=True,
                    scope='Only positive-sum final vectors normalized; invalid vectors retained raw and unassigned',
                    row_alignment='Original prepared h5ad order; official preprocessing skipped expression filtering and only swapped coordinates')
    common.write_json(args.run / 'coverage.json', coverage)
    common.write_json(args.run / 'initial_export_failure.json', record)
    record.update(status='completed', model_training_status='completed', export_status='completed_with_explicit_unassigned_zero_vectors',
                  initial_export_error=record.pop('error', None), export_coverage=coverage,
                  recovery='Export only from already-trained official final native scores; no training rerun',
                  utc_export_recovered=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    common.write_json(args.run / 'run.json', record)
    print(json.dumps(coverage, indent=2), flush=True)


if __name__ == '__main__':
    main()
