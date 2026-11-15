"""Read original 10X H5 and positions, join archived labels, freeze a common universe."""
from __future__ import annotations

import argparse
from pathlib import Path
import common
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def read_h5(path: Path) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    with h5py.File(path) as handle:
        m = handle['matrix']
        x = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
        barcodes = m['barcodes'][:].astype(str)
        features = m['features']
        genes = features['name'][:].astype(str)
        if 'feature_type' in features:
            keep = features['feature_type'][:].astype(str) == 'Gene Expression'
            x, genes = x[:, keep], genes[keep]
    # Sum duplicate symbols using an explicit sparse feature mapping.
    unique, inverse = np.unique(genes, return_inverse=True)
    mapping = sparse.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), inverse)), shape=(len(genes), len(unique)))
    return (x @ mapping).tocsr(), barcodes, unique


def positions(path: Path) -> pd.DataFrame:
    names = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    frame = pd.read_csv(path)
    if 'barcode' not in frame:
        frame = pd.read_csv(path, header=None, names=names)
    if frame['barcode'].duplicated().any():
        raise ValueError('Duplicate coordinate barcode')
    return frame.set_index('barcode')


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, default=common.ROOT.parents[1] / 'Data/GSE274103/raw')
    p.add_argument('--spatial-root', type=Path, default=common.ROOT.parents[1] / 'Data/GSE274103/prepared')
    p.add_argument('--candidates', type=Path, default=common.ROOT.parents[1] / 'Source_Data/fresh_coordinates.tsv.gz')
    p.add_argument('--out', type=Path, default=common.ROOT / 'inputs')
    p.add_argument('--n-genes', type=int, default=3000)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    labels = pd.read_csv(args.candidates, sep='\t')
    if labels.duplicated(['sample_id', 'barcode']).any():
        raise ValueError('Duplicate archived sample/barcode')
    datasets, audits, sections = [], [], []
    for sample in common.SAMPLES:
        count_path = args.raw / f'{sample}_filtered_feature_bc_matrix.h5'
        matches = list((args.spatial_root / sample).rglob('tissue_positions*.csv'))
        if len(matches) != 1:
            raise ValueError(f'Expected exactly one positions file for {sample}: {matches}')
        pos = positions(matches[0])
        x, barcodes, genes = read_h5(count_path)
        raw = pd.DataFrame({'barcode': barcodes, 'raw_index': np.arange(len(barcodes))})
        subset = labels.loc[labels.sample_id == sample, ['sample_id', 'barcode', 'cadspot_candidate']].copy()
        universe = subset.merge(raw, on='barcode', how='left', validate='one_to_one').merge(pos.reset_index(), on='barcode', how='left', validate='one_to_one')
        universe['present_raw'] = universe.raw_index.notna()
        universe['finite_coordinates'] = np.isfinite(universe[['pxl_col_in_fullres', 'pxl_row_in_fullres']]).all(axis=1)
        universe['raw_count'] = 0.0
        universe['detected_genes'] = 0
        present = universe.present_raw.to_numpy()
        indices = universe.loc[present, 'raw_index'].astype(int).to_numpy()
        universe.loc[present, 'raw_count'] = np.asarray(x[indices].sum(axis=1)).ravel()
        universe.loc[present, 'detected_genes'] = x[indices].getnnz(axis=1)
        universe['eligible'] = universe.present_raw & universe.finite_coordinates & universe.in_tissue.eq(1) & (universe.raw_count > 0) & (universe.detected_genes >= 3)
        universe['cadspot_candidate'] = common.bool_labels(universe.cadspot_candidate)
        audits.append(universe)
        obs = universe.loc[universe.eligible].copy().reset_index(drop=True)
        a = ad.AnnData(x[obs.raw_index.astype(int).to_numpy()].astype(np.float32), obs=obs, var=pd.DataFrame(index=genes))
        a.obs_names = obs.barcode.to_numpy()
        a.obsm['spatial'] = obs[['pxl_col_in_fullres', 'pxl_row_in_fullres']].to_numpy(dtype=float)
        datasets.append(a)
        directory = matches[0].parent
        sections.append(dict(sample_id=sample, h5ad=str((args.out / f'{sample}.h5ad').resolve()),
                             spatial_dir=str(directory.resolve()), image_path=str((directory / 'tissue_hires_image.png').resolve()),
                             scalefactors=str((directory / 'scalefactors_json.json').resolve()),
                             n_eligible=a.n_obs, raw_h5=str(count_path), raw_sha256=common.sha256(count_path)))
    common_genes = sorted(set.intersection(*(set(a.var_names) for a in datasets)))
    normalized = sparse.vstack([common.log_normalize(a[:, common_genes].X, a.obs.raw_count.to_numpy()) for a in datasets]).tocsr()
    mean = np.asarray(normalized.mean(axis=0)).ravel()
    var = np.asarray(normalized.power(2).mean(axis=0)).ravel() - mean**2
    expressed = np.asarray((normalized > 0).sum(axis=0)).ravel() >= 3
    ranks = np.lexsort((np.asarray(common_genes), -var))
    selected = [common_genes[i] for i in ranks if expressed[i]][:args.n_genes]
    if len(selected) != args.n_genes:
        raise ValueError('Insufficient expressed common genes')
    pd.DataFrame({'gene': selected}).to_csv(args.out / 'selected_genes.tsv', sep='\t', index=False)
    for a, section in zip(datasets, sections):
        a = a[:, selected].copy()
        a.obsp['spatial_connectivities'] = common.graph(a.obsm['spatial'])
        a.uns['preparation'] = 'raw counts; common top log-normalized-variance genes; candidate-blind; first-ring radius graph'
        a.write_h5ad(section['h5ad'], compression='gzip')
    audit = pd.concat(audits, ignore_index=True)
    audit.to_csv(args.out / 'universe_audit.tsv.gz', sep='\t', index=False)
    audit.loc[audit.eligible].to_csv(args.out / 'eligible_spots.tsv.gz', sep='\t', index=False)
    summary = audit.groupby('sample_id').agg(n_evaluated=('barcode', 'size'), n_eligible=('eligible', 'sum'), n_candidate_evaluated=('cadspot_candidate', 'sum'))
    summary.to_csv(args.out / 'universe_summary.tsv', sep='\t')
    representative = next((s['sample_id'] for s in sections if Path(s['image_path']).is_file() and Path(s['scalefactors']).is_file()), None)
    common.write_json(args.out / 'manifest.json', dict(sections=sections, representative=representative,
                      factor_count=10, seed=20261115, n_genes=len(selected), gene_rule='pooled normalized log-expression variance; labels unused',
                      candidate_table=str(args.candidates), candidate_sha256=common.sha256(args.candidates),
                      eligible_rule='archived evaluated intersection raw in-tissue finite coordinates nonzero counts >=3 genes'))
    print(summary.to_string())
    print(f'Representative (fixed sample order and image availability): {representative}')


if __name__ == '__main__':
    main()
