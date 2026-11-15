"""Execute actual official POPARI, STAMP, or full SpaHDmap, with provenance."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import common
import anndata as ad
import numpy as np
import pandas as pd


def load_section(section: dict) -> ad.AnnData:
    a = ad.read_h5ad(section['h5ad'])
    # Candidate labels cannot enter fitting, initialization or model selection.
    a.obs = a.obs[['sample_id', 'barcode', 'raw_count', 'array_row', 'array_col']].copy()
    return a


def run_popari(args: argparse.Namespace, sections: list[dict], out: Path) -> None:
    import torch
    import awkward as ak
    from popari.model import Popari
    dtype = torch.float64 if args.popari_dtype == 'float64' else torch.float32
    arrays = [load_section(s) for s in sections]
    for a in arrays:
        a.X = common.log_normalize(a.X, a.obs.raw_count.to_numpy())
        a.obsp['adjacency_matrix'] = a.obsp['spatial_connectivities'].copy()
        g = a.obsp['adjacency_matrix'].tocsr()
        a.obsm['adjacency_list'] = ak.Array([g.indices[g.indptr[i]:g.indptr[i+1]].tolist() for i in range(a.n_obs)])
    model = Popari(K=args.k, datasets=arrays, replicate_names=[s['sample_id'] for s in sections],
                   random_state=args.seed, initialization_method='svd', hierarchical_levels=1,
                   spatial_affinity_groups='disjoint', torch_context={'device': args.device, 'dtype': dtype},
                   initial_context={'device': 'cpu', 'dtype': dtype}, verbose=0)
    for iteration in range(2 if args.smoke else 5):
        model.estimate_parameters(update_spatial_affinities=False)
        model.estimate_weights(use_neighbors=False)
        print(f'POPARI initialization iteration {iteration+1}', flush=True)
    model.save_results(str(out / 'initialization_checkpoint.h5ad'))
    objectives = []
    stable = 0
    previous_embeddings = [a.obsm['X'].copy() for a in model.datasets]
    for iteration in range(2 if args.smoke else args.popari_iterations):
        model.estimate_parameters()
        model.estimate_weights()
        loss = np.asarray(model.nll(use_spatial=True)).item()
        if not np.isfinite(loss):
            raise FloatingPointError('Nonfinite POPARI spatial objective')
        relative_objective_change = abs(loss-objectives[-1]['negative_log_likelihood'])/max(abs(loss), 1) if objectives else None
        changes = [float(np.max(np.abs(a.obsm['X']-old))/max(np.max(np.abs(old)), 1e-12))
                   for a, old in zip(model.datasets, previous_embeddings)]
        loading_change = max(changes)
        previous_embeddings = [a.obsm['X'].copy() for a in model.datasets]
        if iteration+1 >= args.popari_min_iterations and relative_objective_change is not None and relative_objective_change < args.popari_objective_tol and loading_change < args.popari_loading_tol:
            stable += 1
        else:
            stable = 0
        objectives.append(dict(iteration=iteration+1, negative_log_likelihood=loss,
                               relative_objective_change=relative_objective_change, max_relative_loading_change=loading_change,
                               consecutive_stable_iterations=stable))
        pd.DataFrame(objectives).to_csv(out / 'training_objective.tsv', sep='\t', index=False)
        print(f'POPARI spatial iteration {iteration+1}: nll={loss:.6g}, loading_change={loading_change:.4g}, stable={stable}', flush=True)
        if (iteration+1) % 25 == 0:
            model.save_results(str(out / 'spatial_checkpoint.h5ad'))
        if not args.smoke and stable >= args.popari_patience:
            break
    common.write_json(out / 'convergence.json', dict(iterations=len(objectives),
                      reason='objective_and_loading_stability' if stable >= args.popari_patience else 'fixed_iteration_limit',
                      convergence_claim=stable >= args.popari_patience, last_objective=objectives[-1]))
    model.synchronize_datasets()
    model.save_results(str(out / 'popari_model.h5ad'))
    for a in model.datasets:
        common.export_scores(out / a.name, a, a.obsm['X'], 'POPARI', args.seed, smoke=args.smoke)


def run_stamp(args: argparse.Namespace, sections: list[dict], out: Path) -> None:
    import pyro
    from sctm.stamp import STAMP
    for s in sections:
        common.seed_all(args.seed)
        pyro.set_rng_seed(args.seed)
        a = load_section(s)
        model = STAMP(a, n_topics=args.k, n_layers=1, hidden_size=128, mode='sign', gene_likelihood='nb')
        model.train(max_epochs=2 if args.smoke else args.stamp_epochs,
                    min_epochs=0 if args.smoke else 100, early_stop=not args.smoke,
                    batch_size=args.batch_size, device=args.device)
        destination = out / s['sample_id']
        destination.mkdir(parents=True, exist_ok=True)
        history = np.asarray(model.training_history)
        pd.DataFrame({'epoch': np.arange(1, len(history)+1), 'mean_batch_scaled_negative_elbo': history}).to_csv(
            destination / 'training_objective.tsv', sep='\t', index=False)
        if not np.isfinite(history).all():
            raise FloatingPointError('Nonfinite STAMP objective; no final evidence exported')
        scores = model.get_cell_by_topic().loc[a.obs_names].to_numpy()
        common.export_scores(destination, a, scores, 'STAMP', args.seed, probability=True, smoke=args.smoke)
        model.save(str(destination / 'stamp_model.pt'))
        model.get_feature_by_topic().to_csv(destination / 'gene_topic_weights.tsv', sep='\t')
        print(f'STAMP complete: {s["sample_id"]}', flush=True)


def run_spahdmap(args: argparse.Namespace, sections: list[dict], out: Path) -> None:
    import SpaHDmap as hdmap
    from sklearn.cluster import KMeans
    if len(sections) != 1:
        raise ValueError('SpaHDmap scope is exactly one prespecified representative section')
    s = sections[0]
    a = load_section(s)
    sf = json.loads(Path(s['scalefactors']).read_text())
    scale = float(sf['tissue_hires_scalef'])
    radius = float(sf['spot_diameter_fullres']) * scale / 2
    # Official STData subtracts 1 internally; add 1 to zero-based PNG coordinates.
    a.obsm['spatial'] = a.obsm['spatial'] * scale + 1.0
    a.X = common.log_normalize(a.X, a.obs.raw_count.to_numpy())
    section = hdmap.prepare_stdata(section_name=s['sample_id'], adata=a,
               image_path=s['image_path'], radius=radius, scale_rate=1,
               swap_coord=True, select_hvgs=False, image_type='HE', create_mask=True)
    mapper = hdmap.Mapper(section=section, results_path=str(out / 'official_outputs'),
                         rank=args.k, ratio_pseudo_spots=args.pseudo_ratio, verbose=True)
    mapper.args.num_workers = 0
    mapper.args.batch_size = args.batch_size
    mapper.args.batch_size_train = args.batch_size
    if args.smoke:
        mapper.args.total_iter_pretrain = 2
        mapper.args.total_iter_gcn = 2
        mapper.args.total_iter_train = 2
        mapper.args.fix_iter_train = 1
        mapper.args.rec_iter = 1
    common.write_json(out / 'spahdmap_settings.json', dict(mapper_args=vars(mapper.args),
                      source_image=s['image_path'], fullres_to_hires_scale=scale, hires_radius=radius,
                      coordinate_conversion='zero-based fullres XY * scale + 1; official swap to row/col then subtract1',
                      source='registered Space Ranger hires PNG; not original full-resolution pathology TIFF'))
    section.show()
    import matplotlib.pyplot as plt
    plt.savefig(out / 'registered_overlay.png', dpi=120)
    plt.close('all')
    mapper.run_SpaHDmap(save_score=True, save_model=True, load_model=False, visualize=False)
    final = section.scores['SpaHDmap_spot']
    if final is None:
        raise ValueError('Official final fusion scores missing; NMF is not a substitute')
    valid = common.export_scores(out / s['sample_id'], section.adata, final, 'SpaHDmap', args.seed, smoke=args.smoke)
    # Downstream domains use fixed K on final embeddings. Do not label these as native topics.
    domains = KMeans(n_clusters=args.k, n_init=20, random_state=args.seed).fit_predict(final[valid])
    pd.DataFrame({'sample_id': s['sample_id'], 'barcode': section.adata.obs['barcode'].to_numpy()[valid],
                  'kmeans_domain': domains + 1}).to_csv(out / s['sample_id'] / 'downstream_domains.tsv', sep='\t', index=False)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('method', choices=['popari', 'stamp', 'spahdmap'])
    p.add_argument('--manifest', type=Path, default=common.ROOT / 'inputs' / 'manifest.json')
    p.add_argument('--out', type=Path)
    p.add_argument('--sample', help='Optional exact section ID; omission fits five sections except representative-only SpaHDmap')
    p.add_argument('--k', type=int, default=10)
    p.add_argument('--seed', type=int, default=20261115)
    p.add_argument('--device', choices=['cpu', 'cuda:0'], default='cuda:0')
    p.add_argument('--batch-size', type=int, default=128)
    p.add_argument('--popari-iterations', type=int, default=200)
    p.add_argument('--popari-dtype', choices=['float32', 'float64'], default='float64')
    p.add_argument('--popari-min-iterations', type=int, default=50)
    p.add_argument('--popari-patience', type=int, default=10)
    p.add_argument('--popari-objective-tol', type=float, default=1e-4)
    p.add_argument('--popari-loading-tol', type=float, default=1e-3)
    p.add_argument('--stamp-epochs', type=int, default=800)
    p.add_argument('--pseudo-ratio', type=int, default=5)
    p.add_argument('--smoke', action='store_true', help='Two-step execution test, never final manuscript evidence')
    args = p.parse_args()
    if args.k < 2:
        p.error('K must be at least 2')
    if args.method == 'spahdmap' and args.device == 'cpu':
        # Official Mapper auto-selects CUDA; hide devices before torch import.
        os.environ['CUDA_VISIBLE_DEVICES'] = ''
    import matplotlib
    matplotlib.use('Agg')
    common.seed_all(args.seed)
    manifest = json.loads(args.manifest.read_text())
    sections = manifest['sections']
    sample = args.sample or (manifest['representative'] if args.method == 'spahdmap' else None)
    if sample:
        sections = [s for s in sections if s['sample_id'] == sample]
    if not sections:
        p.error('No matching section')
    out = args.out or common.ROOT / ('smoke' if args.smoke else 'results') / args.method / f'k{args.k}_seed{args.seed}'
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Refusing overwrite of prior results: {out}')
    out.mkdir(parents=True, exist_ok=True)
    class Tee:
        def __init__(self, original: object, file: object):
            self.original, self.file = original, file
        def write(self, value: str) -> None:
            self.original.write(value)
            self.file.write(value)
            self.file.flush()
        def flush(self) -> None:
            self.original.flush()
            self.file.flush()
        def isatty(self) -> bool:
            return False
    logfile = (out / 'console.log').open('w', encoding='utf-8')
    sys.stdout = Tee(sys.stdout, logfile)
    sys.stderr = Tee(sys.stderr, logfile)
    repo = dict(popari='popari', stamp='scTM', spahdmap='SpaHDmap')[args.method]
    sys.path.insert(0, str(common.ROOT / 'vendor' / repo))
    log = common.RunLog(out / 'run.json', args, repo)
    if args.method == 'popari' and args.popari_iterations == 50 and not args.smoke:
        log.data['fixed_spatial_iteration_budget'] = 50
        log.data['convergence_claim'] = False
    log.data['input_manifest_sha256'] = common.sha256(args.manifest)
    log.data['input_h5ad_sha256'] = {s['sample_id']: common.sha256(Path(s['h5ad'])) for s in sections}
    common.write_json(log.path, log.data)
    common.execute_logged(log, lambda: dict(popari=run_popari, stamp=run_stamp, spahdmap=run_spahdmap)[args.method](args, sections, out))


if __name__ == '__main__':
    main()
