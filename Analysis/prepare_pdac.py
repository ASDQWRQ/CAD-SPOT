"""Recover GEO counts, registered images and score provenance by barcode."""
import argparse
import gzip
import hashlib
import json
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

TARGET = 'BAX BAK1 BID PMAIP1 BNIP3 BNIP3L MCL1 BCL2L1 HK2 LDHA PDK1 SLC2A1 CA9 ENO1 ALDOA PGK1 PGAM1 PKM GPI TPI1 VEGFA HIF1A NDRG1'.split()
REFERENCE = 'ADM AK4 ALDOC ANGPTL4 ANKZF1 BHLHE40 CA12 CCNG2 DDIT4 ENO2 FAM162A GAPDH GCK HK1 KDM3A LDHC MXI1 P4HA1 P4HA2 PFKL SLC2A3'.split()


def decode(a):
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in a])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    args = ap.parse_args()
    root = args.root
    raw = root / 'Data/GSE274103/raw'
    prepared = root / 'Data/GSE274103/prepared'
    prepared.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(root / 'Source_Data/archived_coordinate_table.tsv.gz', sep='\t')
    matrices, obs, audits, manifests = [], [], [], []
    gene_names = None
    for sample, sub in df.groupby('sample_id', sort=True):
        st = prepared / sample
        st.mkdir(exist_ok=True)
        archive = raw / f'{sample}_spatial.tar.gz'
        with tarfile.open(archive) as tf:
            tf.extractall(st, filter='data')
        spatial_dir = next(p.parent for p in st.rglob('scalefactors_json.json'))
        pos = pd.read_csv(spatial_dir / 'tissue_positions.csv')
        if 'barcode' not in pos:
            pos = pd.read_csv(spatial_dir / 'tissue_positions.csv', header=None)
        pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
        pos = pos.set_index('barcode').loc[sub.barcode]
        h5file = raw / f'{sample}_filtered_feature_bc_matrix.h5'
        with h5py.File(h5file) as f:
            m = f['matrix']
            x = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
            barcodes = decode(m['barcodes'][:])
            symbols = decode(m['features']['name'][:])
            stable = decode(m['features']['id'][:])
        if gene_names is None:
            gene_names = pd.Index(sorted(set(symbols)))
        transform = sparse.coo_matrix((np.ones(len(symbols)), (np.arange(len(symbols)), gene_names.get_indexer(symbols))), shape=(len(symbols), len(gene_names))).tocsr()
        order = pd.Index(barcodes).get_indexer(sub.barcode)
        if (order < 0).any():
            raise ValueError(f'Missing archived barcode in {sample}')
        x = (x[order] @ transform).tocsr().astype(np.float32)
        sparse.save_npz(st / 'counts.npz', x)
        fullcol = pos.pxl_col_in_fullres.to_numpy(float)
        fullrow = pos.pxl_row_in_fullres.to_numpy(float)
        direct = np.max(np.abs(sub.x.to_numpy() - fullcol)) == 0 and np.max(np.abs(sub.y.to_numpy() - fullrow)) == 0
        swapped = np.max(np.abs(sub.x.to_numpy() - fullrow)) == 0 and np.max(np.abs(sub.y.to_numpy() - fullcol)) == 0
        sub = sub.copy()
        sub['pxl_col_in_fullres'] = fullcol
        sub['pxl_row_in_fullres'] = fullrow
        sub['array_row'] = pos.array_row.to_numpy()
        sub['array_col'] = pos.array_col.to_numpy()
        library = np.asarray(x.sum(axis=1)).ravel()
        detected = np.diff(x.indptr)
        audits.append(dict(sample_id=sample, n_spots=len(sub), all_barcodes_present=True,
                           archived_xy_matches_fullres_colrow=bool(direct), archived_xy_matches_fullres_rowcol=bool(swapped),
                           ncount_max_abs_difference=float(np.max(np.abs(library-sub.nCount_Spatial))),
                           nfeature_max_abs_difference=int(np.max(np.abs(detected-sub.nFeature_Spatial)))))
        sub.to_csv(st / 'spots.tsv', sep='\t', index=False)
        matrices.append(x)
        obs.append(sub)
        manifests.append(dict(sample_id=sample, counts=str((st/'counts.npz').relative_to(root)),
                              spots=str((st/'spots.tsv').relative_to(root)), spatial_dir=str(spatial_dir.relative_to(root)),
                              coordinates='Full-resolution image pixels; x=column, y=row',
                              scaling='Multiply full-resolution pixels by tissue_hires_scalef for the hires image'))
    x = sparse.vstack(matrices, format='csr')
    spots = pd.concat(obs, ignore_index=True)
    sparse.save_npz(prepared / 'counts_all.npz', x)
    pd.Series(gene_names, name='gene').to_csv(prepared / 'genes.tsv', sep='\t', index=False)
    spots.to_csv(prepared / 'spots_all.tsv.gz', sep='\t', index=False)
    lib = np.asarray(x.sum(axis=1)).ravel()
    logx = x.multiply((1e4 / lib)[:, None]).tocsr()
    logx.data = np.log1p(logx.data)
    sparse.save_npz(prepared / 'logcp10k_all.npz', logx)
    mean = np.asarray(logx.mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(0, np.asarray(logx.power(2).mean(axis=0)).ravel()-mean**2))
    detection = np.asarray((x>0).sum(axis=0)).ravel()/x.shape[0]
    metadata = pd.DataFrame(dict(gene=gene_names, mean_expression=mean, detection_rate=detection, sd_log_expression=sd))
    metadata.to_csv(prepared / 'gene_expression_metadata.tsv', sep='\t', index=False)
    scoring = []
    for label, genes, archived in [('target', TARGET, 'target_score'), ('reference', REFERENCE, 'reference_score')]:
        idx = gene_names.get_indexer(genes)
        missing = [g for g, i in zip(genes, idx) if i < 0]
        idx = idx[idx >= 0]
        standardized = (logx[:, idx].toarray()-mean[idx])/np.where(sd[idx]>0, sd[idx], 1)
        score = standardized.mean(axis=1)
        spots[f'recomputed_{label}_mean_gene_z'] = score
        scoring.append(dict(score=label, n_specified=len(genes), n_available=len(idx), missing_genes=';'.join(missing), archived_vs_recomputed_spearman=float(spearmanr(score, spots[archived]).statistic),
                            archived_vs_recomputed_pearson=float(np.corrcoef(score, spots[archived])[0,1])))
    spots.to_csv(prepared / 'spots_all.tsv.gz', sep='\t', index=False)
    pd.DataFrame({'gene': TARGET}).to_csv(root / 'Source_Data/target_genes.tsv', sep='\t', index=False)
    pd.DataFrame({'gene': REFERENCE}).to_csv(root / 'Source_Data/reference_genes.tsv', sep='\t', index=False)
    pd.DataFrame(audits).to_csv(root / 'Source_Data/raw_barcode_coordinate_audit.tsv', sep='\t', index=False)
    pd.DataFrame(scoring).to_csv(root / 'Source_Data/score_reconstruction_audit.tsv', sep='\t', index=False)
    (root / 'Data/GSE274103/manifest.json').write_text(json.dumps(manifests, indent=2)+'\n')
    print(pd.DataFrame(audits).to_string(index=False), flush=True)
    print(pd.DataFrame(scoring).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
