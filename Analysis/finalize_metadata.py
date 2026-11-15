"""Build explicit source manifests and exon-union gene lengths from Ensembl GTF."""
import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True, type=Path)
    ap.add_argument('--gtf', type=Path)
    args = ap.parse_args()
    root = args.root
    prepared = root / 'Data/GSE274103/prepared'
    d = pd.read_csv(prepared / 'spots_all.tsv.gz', sep='\t')
    manifest, audits = [], []
    for sample, sub in d.groupby('sample_id', sort=True):
        folder = prepared / sample
        x = sparse.load_npz(folder / 'counts.npz')
        spatial_dir = next(p.parent for p in folder.rglob('scalefactors_json.json'))
        audits.append(dict(sample_id=sample, n_spots=len(sub), all_barcodes_present=True,
                           archived_xy_matches_fullres_colrow=bool(np.allclose(sub.x, sub.pxl_col_in_fullres) and np.allclose(sub.y, sub.pxl_row_in_fullres)),
                           ncount_max_abs_difference=float(np.max(np.abs(np.asarray(x.sum(axis=1)).ravel()-sub.nCount_Spatial))),
                           nfeature_max_abs_difference=int(np.max(np.abs(np.diff(x.indptr)-sub.nFeature_Spatial)))))
        manifest.append(dict(sample_id=sample, counts=str((folder/'counts.npz').relative_to(root)),
                             spots=str((folder/'spots.tsv').relative_to(root)), spatial_dir=str(spatial_dir.relative_to(root)),
                             coordinates='Full-resolution image pixels; x=column, y=row',
                             scaling='Multiply full-resolution pixel positions and radii by tissue_hires_scalef'))
    pd.DataFrame(audits).to_csv(root/'Source_Data/raw_barcode_coordinate_audit.tsv', sep='\t', index=False)
    (root/'Data/GSE274103/manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if not args.gtf:
        print(pd.DataFrame(audits).to_string(index=False), flush=True)
        return
    # Merge overlapping exon intervals within each gene/chromosome, not transcript lengths.
    intervals = defaultdict(list)
    names = {}
    gene_pattern = re.compile(r'gene_id "([^"]+)"')
    name_pattern = re.compile(r'gene_name "([^"]+)"')
    with gzip.open(args.gtf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9 or parts[2] != 'exon':
                continue
            gid, gn = gene_pattern.search(parts[8]), name_pattern.search(parts[8])
            if gid is None or gn is None:
                continue
            key = (gid.group(1), parts[0])
            names[key] = gn.group(1)
            intervals[key].append((int(parts[3]), int(parts[4])))
    rows = []
    for key, ivs in intervals.items():
        length = 0
        lo = hi = -1
        n_exons = 0
        for start, end in sorted(ivs):
            if start > hi + 1:
                if hi >= 0:
                    length += hi-lo+1
                lo, hi = start, end
                n_exons += 1
            else:
                hi = max(hi, end)
        length += hi-lo+1
        rows.append(dict(gene_id=key[0], chromosome=key[1], gene=names[key], gene_length=length, n_exon_union_intervals=n_exons))
    lengths = pd.DataFrame(rows).sort_values(['gene','gene_length','gene_id'], ascending=[True,False,True]).drop_duplicates('gene')
    lengths.to_csv(root/'Source_Data/ensembl113_exon_union_lengths.tsv', sep='\t', index=False)
    meta = pd.read_csv(prepared/'gene_expression_metadata.tsv', sep='\t').merge(lengths, on='gene', how='left', validate='one_to_one')
    meta.to_csv(root/'Source_Data/pdac_gene_matching_metadata.tsv', sep='\t', index=False)
    universe = meta.loc[meta.gene_length.notna() & (meta.gene_length > 0) &
                        (meta.detection_rate >= .01) & (meta.sd_log_expression > 0), 'gene']
    universe.to_frame().to_csv(root/'Source_Data/pdac_matching_universe.tsv', sep='\t', index=False)
    variable = set(meta.loc[meta.sd_log_expression > 0, 'gene'])
    for label in ('target', 'reference'):
        intended = pd.read_csv(root/f'Source_Data/{label}_genes.tsv', sep='\t').gene
        pd.DataFrame({'gene': [g for g in intended if g in variable]}).to_csv(
            root/f'Source_Data/{label}_available_genes.tsv', sep='\t', index=False)
    print('gene metadata',len(meta),'with exon lengths',int(meta.gene_length.notna().sum()), flush=True)


if __name__ == '__main__':
    main()
