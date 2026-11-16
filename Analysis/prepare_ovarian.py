"""Rebuild the ovarian worked example from public Matrix Market/Visium files."""
import argparse
import json
from pathlib import Path
import tarfile
import zipfile

import numpy as np
import pandas as pd
from scipy import io, sparse

from end_to_end import ROOT, PANELS, TECH, COMP, prepare_section, moments, tsv, json_out, core


def read_market(matrix, features, barcodes):
    feature = pd.read_csv(features, sep='\t', header=None, dtype=str)
    barcode = pd.read_csv(barcodes, sep='\t', header=None, dtype=str)[0].to_numpy()
    if feature.shape[1] < 2:
        raise ValueError('Matrix Market features require gene IDs and symbols')
    x = io.mmread(str(matrix)).tocsr()
    if x.shape != (len(feature), len(barcode)):
        raise ValueError('Matrix/feature/barcode dimension mismatch')
    used = feature.iloc[:, 2].eq('Gene Expression').to_numpy() if feature.shape[1] >= 3 else np.ones(len(feature), bool)
    names = feature.loc[used].iloc[:, 1].to_numpy()
    if pd.isna(names).any() or pd.Index(barcode).has_duplicates:
        raise ValueError('Missing gene names or duplicated expression barcodes')
    genes = pd.Index(sorted(set(names)))
    transform = sparse.coo_matrix((np.ones(len(names)), (np.arange(len(names)), genes.get_indexer(names))),
        shape=(len(names), len(genes))).tocsr()
    x = (x[used].T @ transform).tocsr().astype(float)
    if not np.isfinite(x.data).all() or (x.data < 0).any() or not np.allclose(x.data, np.round(x.data)):
        raise ValueError('Expected finite nonnegative integer UMI counts')
    x.eliminate_zeros()
    return x, genes, barcode


def extract_spatial(path, destination):
    destination.mkdir(parents=True, exist_ok=True)
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            if archive.testzip() is not None:
                raise ValueError('Corrupt spatial ZIP')
            base = destination.resolve()
            for name in archive.namelist():
                if not (destination/name).resolve().is_relative_to(base):
                    raise ValueError('Unsafe spatial archive member')
            archive.extractall(destination)
    else:
        with tarfile.open(path) as archive:
            archive.extractall(destination, filter='data')
    files = list(destination.rglob('tissue_positions*'))
    if len(files) != 1:
        raise ValueError(f'Expected one tissue positions file in {path}, found {files}')
    return files[0]


def read_positions(path, barcodes):
    columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'y', 'x']
    d = pd.read_csv(path, header=None)
    if str(d.iloc[0, 0]) == 'barcode':
        d = d.iloc[1:].reset_index(drop=True)
    if d.shape[1] != 6:
        raise ValueError('Expected six Visium coordinate columns')
    d.columns = columns
    if d.barcode.duplicated().any():
        raise ValueError('Duplicate spatial barcodes')
    for col in columns[1:]:
        d[col] = pd.to_numeric(d[col], errors='raise')
    d = d.set_index('barcode').reindex(barcodes).reset_index()
    if d.in_tissue.isna().any():
        raise ValueError('Expression barcodes missing from spatial file')
    return d


def harmonize_scores(destination, samples):
    """Fix cohort-wide gene membership using coverage only, before association."""
    all_coverage = pd.concat([pd.read_csv(destination/s/'panel_coverage.tsv', sep='\t')
                             for s in samples], ignore_index=True)
    memberships = {}
    for panel, intended in PANELS.items():
        used = [set(sub.loc[sub.status.eq('used'), 'gene']) for _, sub in
                all_coverage[all_coverage.panel.eq(panel)].groupby('sample_id')]
        common = set.intersection(*used)
        memberships[panel] = [g for g in intended if g in common]
        if len(common) < 3:
            raise ValueError(f'Fewer than three shared variable genes in {panel}')
    frames, coverage_frames = [], []
    for sample in samples:
        folder = destination/sample
        spots = pd.read_csv(folder/'score_table.tsv.gz', sep='\t')
        parameters = pd.read_csv(folder/'score_parameters.tsv', sep='\t')
        genes = pd.Index(pd.read_csv(folder/'genes.tsv', sep='\t').gene)
        values = sparse.load_npz(folder/'logcp10k.npz')
        for panel, common in memberships.items():
            panel_rows = parameters.panel.eq(panel)
            excluded = panel_rows & parameters.status.eq('used') & ~parameters.gene.isin(common)
            parameters.loc[excluded, 'status'] = 'excluded_not_common_across_sections'
            parameters.loc[panel_rows, 'weight'] = 0.
            used = panel_rows & parameters.status.eq('used')
            parameters.loc[used, 'weight'] = 1/len(common)
            sub = parameters.loc[used]
            spots[panel] = ((values[:, genes.get_indexer(sub.gene)].toarray() - sub.mean_logcp10k.to_numpy()) /
                            sub.sd_logcp10k.to_numpy()) @ sub.weight.to_numpy()
        coverage = parameters[['sample_id', 'panel', 'gene', 'status']]
        tsv(spots, folder/'score_table.tsv.gz')
        tsv(parameters, folder/'score_parameters.tsv')
        tsv(coverage, folder/'panel_coverage.tsv')
        frames.append(spots)
        coverage_frames.append(coverage)
    json_out(dict(rule='intersection_of_available_nonconstant_genes_before_any_association',
                  used_panels=memberships), destination/'cohort_panel_membership.json')
    return frames, pd.concat(coverage_frames, ignore_index=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, default=ROOT/'Data/GSE211956')
    args = p.parse_args()
    raw = args.source/'raw'
    destination = ROOT/'Primary'
    destination.mkdir(exist_ok=True)
    frames, matrices, axes, coverages, sds, sources = [], [], [], [], [], []
    for gsm in range(6506110, 6506118):
        matrix = next(raw.glob(f'GSM{gsm}_*matrix.mtx.gz'))
        sample = matrix.name.removesuffix('_matrix.mtx.gz')
        features = raw/f'{sample}_features.tsv.gz'
        barcodes = raw/f'{sample}_barcodes.tsv.gz'
        candidates = list(raw.glob(f'GSM{gsm}_*spatial*'))
        if len(candidates) != 1:
            raise ValueError(f'Expected one public spatial archive for {sample}: {candidates}')
        archive = candidates[0]
        pospath = extract_spatial(archive, args.source/'prepared'/sample)
        x, genes, ids = read_market(matrix, features, barcodes)
        positions = read_positions(pospath, ids)
        print('Preparing ovarian', sample, x.shape, flush=True)
        logx, spots, meta, coverage = prepare_section(x, genes, positions, sample, destination/sample)
        if axes and not axes[0].equals(genes):
            raise ValueError('Section gene axes differ; align explicitly before pooling')
        axes.append(genes)
        matrices.append(logx)
        frames.append(spots)
        sds.append(meta.sd_log_expression.to_numpy())
        coverages.extend(coverage)
        sources.append(dict(sample_id=sample, cohort='ovarian_primary', accession='GSE211956',
            container_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211956/suppl/GSE211956_RAW.tar',
            counts=str(matrix), counts_sha256=core.file_sha256(matrix), counts_archive_member=matrix.name,
            features_sha256=core.file_sha256(features), barcodes_sha256=core.file_sha256(barcodes),
            spatial=str(archive), spatial_sha256=core.file_sha256(archive), spatial_archive_member=archive.name,
            tissue_positions=str(pospath), spatial_dir=str(pospath.parent),
            coordinate_units='full_resolution_pixels', orientation='x=column, y=row'))
    frames, harmonized_coverage = harmonize_scores(destination, [s['sample_id'] for s in sources])
    logx = sparse.vstack(matrices, format='csr')
    genes = axes[0]
    spots = pd.concat(frames, ignore_index=True)
    mean, sd = moments(logx)
    metadata = pd.DataFrame(dict(gene=genes, mean_expression=mean, sd_log_expression=sd,
        min_section_sd=np.min(np.stack(sds), axis=0),
        detection_rate=np.asarray((logx>0).sum(axis=0)).ravel()/len(spots)))
    sparse.save_npz(destination/'logcp10k_all.npz', logx)
    tsv(pd.DataFrame({'gene': genes}), destination/'genes.tsv')
    tsv(spots, destination/'score_table.tsv.gz')
    tsv(metadata, destination/'gene_metadata.tsv.gz')
    tsv(harmonized_coverage, destination/'panel_coverage.tsv')
    json_out(sources, destination/'input_manifest.json')
    json_out(dict(dataset='GSE211956', role='primary_ovarian_worked_example', panels=PANELS,
        technical=TECH, composition=COMP, scoring='mean_gene_z_by_section', normalization='log1p_CP10K',
        qc='in_tissue=1,total_UMI>=100,detected_unique_symbols>=100,finite_coordinates',
        membership='cohort_intersection_of_available_nonconstant_genes_before_association',
        frozen_score_table_used=False, prior_prediction_or_deconvolution_used=False), destination/'pipeline_config.json')
    print('Ovarian reconstruction complete:', len(spots), 'spots', flush=True)


if __name__ == '__main__':
    main()
