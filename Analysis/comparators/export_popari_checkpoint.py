"""Freeze an actual official POPARI checkpoint at an approved fixed spatial-update budget."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import time
import common
import anndata as ad
import numpy as np
import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--iteration', type=int, default=50)
    p.add_argument('--timeout-minutes', type=float, default=40)
    args = p.parse_args()
    if args.iteration % 25:
        p.error('The active training wrapper saves checkpoints every25 iterations')
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True, exist_ok=True)
    common.write_json(args.out / 'run.json', dict(status='waiting_for_checkpoint', arguments=vars(args)))
    deadline = time.monotonic()+args.timeout_minutes*60
    checkpoint = args.run / 'spatial_checkpoint.h5ad'
    history_path = args.run / 'training_objective.tsv'
    while time.monotonic() < deadline:
        parent = json.loads((args.run / 'run.json').read_text())
        if parent['arguments'].get('smoke'):
            raise ValueError('Cannot promote a smoke-test checkpoint to production evidence')
        if parent['status'] not in ['running', 'completed']:
            raise RuntimeError(f'Parent training cannot reach checkpoint: {parent["status"]}')
        if history_path.exists() and checkpoint.exists():
            try:
                history = pd.read_csv(history_path, sep='\t')
            except (pd.errors.EmptyDataError, pd.errors.ParserError):
                time.sleep(1)
                continue
            last = int(history.iteration.max())
            if last >= args.iteration+25:
                raise RuntimeError('Requested checkpoint window already passed; refuse mislabeled later model')
            # Iteration51 cannot be logged until the iteration50 checkpoint is closed.
            if args.iteration < last < args.iteration+25:
                frozen = args.out / 'popari_model.h5ad'
                shutil.copyfile(checkpoint, frozen)
                after = pd.read_csv(history_path, sep='\t')
                if int(after.iteration.max()) >= args.iteration+25:
                    raise RuntimeError('Checkpoint window crossed while copying; refuse ambiguous iteration label')
                a = ad.read_h5ad(frozen)
                if 'X' not in a.obsm or a.obsm['X'].shape[1] != parent['arguments']['k']:
                    raise ValueError('Official checkpoint does not contain expected learned weights')
                for sample in common.SAMPLES:
                    section = a[a.obs.sample_id == sample].copy()
                    if not section.n_obs:
                        raise ValueError(f'Missing section {sample}')
                    common.export_scores(args.out / sample, section, section.obsm['X'], 'POPARI', parent['arguments']['seed'])
                    pd.DataFrame(section.uns['M'][sample], index=section.var_names,
                                 columns=[f'factor_{i+1:02d}' for i in range(parent['arguments']['k'])]).to_csv(
                                     args.out / sample / 'gene_factor_weights.tsv', sep='\t')
                history = history.loc[history.iteration <= args.iteration]
                history.to_csv(args.out / 'training_objective.tsv', sep='\t', index=False)
                common.write_json(args.out / 'convergence.json', dict(iterations=args.iteration,
                                  reason='main_approved_fixed_iteration_budget', convergence_claim=False,
                                  last_objective=history.iloc[-1].to_dict()))
                record = dict(parent)
                record.update(status='completed', source_run=str(args.run), source_run_status_at_export=parent['status'],
                              fixed_spatial_iteration_budget=args.iteration, convergence_claim=False,
                              export_arguments=vars(args), checkpoint_sha256=common.sha256(frozen),
                              interpretation=f'Actual POPARI spatial model at{args.iteration} updates after5 nonspatial warmup; fixed-budget contextual comparison, not claimed converged',
                              utc_exported=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
                common.write_json(args.out / 'run.json', record)
                print(f'Frozen actual POPARI {args.iteration}-iteration checkpoint: {args.out}', flush=True)
                return
        time.sleep(2)
    raise TimeoutError('Approved POPARI checkpoint not yet available; no surrogate result exported')


if __name__ == '__main__':
    main()
