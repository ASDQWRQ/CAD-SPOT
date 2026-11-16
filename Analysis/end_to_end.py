"""Count-matrix-to-coordinate follow-up; no archived scores are inputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tarfile
import time
import urllib.request

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'Software'))
import cadspot as core
import cadspot_null as null

TARGET = 'BAX BAK1 BID PMAIP1 BNIP3 BNIP3L MCL1 BCL2L1 HK2 LDHA PDK1 SLC2A1 CA9 ENO1 ALDOA PGK1 PGAM1 PKM GPI TPI1 VEGFA HIF1A NDRG1'.split()
REFERENCE = 'ADM AK4 ALDOC ANGPTL4 ANKZF1 BHLHE40 CA12 CCNG2 DDIT4 ENO2 FAM162A GAPDH GCK HK1 KDM3A LDHC MXI1 P4HA1 P4HA2 PFKL SLC2A3'.split()
PANELS = dict(target_score=TARGET, reference_score=REFERENCE,
    Epithelial_proxy='EPCAM KRT7 KRT8 KRT18 KRT19'.split(),
    ECM_proxy='COL1A1 COL1A2 COL3A1 DCN LUM'.split(),
    Contractile_proxy='ACTA2 TAGLN FAP PDGFRB'.split(),
    Myeloid_proxy='LST1 TYROBP FCER1G AIF1 CD68'.split(),
    Tcell_proxy='CD3D CD3E TRAC CD247 LCK'.split())
TECH = ['log_nCount', 'log_nFeature']
COMP = list(PANELS)[2:]
SEED = 20261115


def tsv(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False)


def json_out(value, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, default=str) + '\n', encoding='utf-8')


def decode(values):
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def read_h5(path):
    with h5py.File(path) as f:
        m = f['matrix']
        x = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
        names = decode(m['features']['name'][:])
        barcodes = decode(m['barcodes'][:])
    genes = pd.Index(sorted(set(names)))
    transform = sparse.coo_matrix((np.ones(len(names)), (np.arange(len(names)), genes.get_indexer(names))), shape=(len(names), len(genes))).tocsr()
    x = (x @ transform).tocsr().astype(np.float64)
    x.eliminate_zeros()
    return x, genes, barcodes


def moments(x):
    mean = np.asarray(x.mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(0, np.asarray(x.power(2).mean(axis=0)).ravel() - mean**2))
    return mean, sd


def prepare_section(x, genes, positions, sample_id, output):
    """positions must have the same row ordering as x, verified by caller."""
    output.mkdir(parents=True, exist_ok=True)
    if len(positions) != x.shape[0]:
        raise ValueError('Counts/positions row mismatch')
    lib = np.asarray(x.sum(axis=1)).ravel()
    detect = np.diff(x.indptr)
    keep = (lib >= 100) & (detect >= 100) & positions.in_tissue.eq(1).to_numpy()
    keep &= np.isfinite(positions[['x', 'y']].to_numpy(dtype=float)).all(axis=1)
    qc = dict(sample_id=sample_id, n_input=x.shape[0], n_kept=int(keep.sum()),
              n_umi_below_100=int((lib < 100).sum()), n_detected_below_100=int((detect < 100).sum()),
              n_out_of_tissue=int(positions.in_tissue.ne(1).sum()),
              n_missing_coordinates=int(positions[['x', 'y']].isna().any(axis=1).sum()))
    x = x[keep].tocsr()
    lib, detect = lib[keep], detect[keep]
    spots = positions.loc[keep].reset_index(drop=True).copy()
    spots.insert(0, 'sample_id', sample_id)
    spots['nCount'] = lib
    spots['nFeature'] = detect
    spots['log_nCount'] = np.log1p(lib)
    spots['log_nFeature'] = np.log1p(detect)
    logx = x.multiply((1e4 / lib)[:, None]).tocsr()
    logx.data = np.log1p(logx.data)
    mean, sd = moments(logx)
    detection = np.asarray((x > 0).sum(axis=0)).ravel() / x.shape[0]
    metadata = pd.DataFrame(dict(gene=genes, mean_expression=mean, sd_log_expression=sd, detection_rate=detection))
    coverage, parameters = [], []
    for panel, members in PANELS.items():
        if panel in COMP and set(members) & (set(TARGET) | set(REFERENCE)):
            raise ValueError('Proxy panel overlaps target/reference')
        idx = genes.get_indexer(members)
        available = [g for g, i in zip(members, idx) if i >= 0 and sd[i] > 1e-12]
        if len(available) < 3:
            raise ValueError(f'Insufficient panel coverage {sample_id}/{panel}: {available}')
        used = genes.get_indexer(available)
        spots[panel] = ((logx[:, used].toarray() - mean[used]) / sd[used]).mean(axis=1)
        for g, i in zip(members, idx):
            status = 'used' if g in available else ('missing' if i < 0 else 'constant')
            coverage.append(dict(sample_id=sample_id, panel=panel, gene=g, status=status))
            parameters.append(dict(sample_id=sample_id, panel=panel, gene=g, status=status,
                                   mean_logcp10k=float(mean[i]) if i >= 0 else np.nan,
                                   sd_logcp10k=float(sd[i]) if i >= 0 else np.nan,
                                   weight=1 / len(available) if status == 'used' else 0))
    sparse.save_npz(output / 'counts.npz', x)
    sparse.save_npz(output / 'logcp10k.npz', logx)
    tsv(pd.DataFrame({'gene': genes}), output / 'genes.tsv')
    tsv(spots, output / 'score_table.tsv.gz')
    tsv(metadata, output / 'gene_metadata.tsv.gz')
    tsv(pd.DataFrame(coverage), output / 'panel_coverage.tsv')
    tsv(pd.DataFrame(parameters), output / 'score_parameters.tsv')
    json_out(qc, output / 'qc.json')
    return logx, spots, metadata, coverage


def prepare_primary(source_root, destination=None):
    destination = Path(destination) if destination is not None else ROOT / 'Primary'
    matrices, frames, coverages, sources, sds = [], [], [], [], []
    genes = None
    for gsm in range(8443449, 8443454):
        path = next((source_root / 'Data/GSE274103/raw').glob(f'GSM{gsm}_*matrix.h5'))
        sample_id = path.name.replace('_filtered_feature_bc_matrix.h5', '')
        print('Preparing', sample_id, flush=True)
        x, names, barcodes = read_h5(path)
        if genes is not None and not genes.equals(names):
            raise ValueError('Primary sections have different gene axes; explicit alignment required')
        genes = names
        cached = list((source_root / 'Data/GSE274103/prepared' / sample_id).rglob('tissue_positions.csv'))
        if cached:
            spatial_dir = cached[0].parent
        else:
            extracted = ROOT / 'Data/GSE274103/prepared' / sample_id
            extracted.mkdir(parents=True, exist_ok=True)
            with tarfile.open(source_root / f'Data/GSE274103/raw/{sample_id}_spatial.tar.gz') as archive:
                archive.extractall(extracted, filter='data')
            spatial_dir = next(extracted.rglob('tissue_positions.csv')).parent
        pos = pd.read_csv(spatial_dir / 'tissue_positions.csv')
        if 'barcode' not in pos:
            pos = pd.read_csv(spatial_dir / 'tissue_positions.csv', header=None)
        pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'y', 'x']
        if pos.barcode.duplicated().any():
            raise ValueError('Duplicate position barcode')
        pos = pos.set_index('barcode').reindex(barcodes).reset_index()
        if pos.in_tissue.isna().any():
            raise ValueError('Missing raw barcode in public tissue positions')
        folder = destination / sample_id
        logx, spots, metadata, coverage = prepare_section(x, names, pos, sample_id, folder)
        matrices.append(logx)
        frames.append(spots)
        coverages.extend(coverage)
        sds.append(metadata.sd_log_expression.to_numpy())
        archive = source_root / f'Data/GSE274103/raw/{sample_id}_spatial.tar.gz'
        sources.append(dict(sample_id=sample_id, counts=str(path), counts_sha256=core.file_sha256(path),
            container_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274103/suppl/GSE274103_RAW.tar',
            counts_archive_member=path.name,
            spatial=str(archive), spatial_sha256=core.file_sha256(archive),
            spatial_archive_member=archive.name,
            tissue_positions=str(spatial_dir / 'tissue_positions.csv'), spatial_dir=str(spatial_dir),
            coordinate_units='full_resolution_pixels', orientation='x=column, y=row'))
    logx = sparse.vstack(matrices, format='csr')
    spots = pd.concat(frames, ignore_index=True)
    mean, sd = moments(logx)
    counts_detection = np.asarray((logx > 0).sum(axis=0)).ravel() / logx.shape[0]
    metadata = pd.DataFrame(dict(gene=genes, mean_expression=mean, sd_log_expression=sd,
        min_section_sd=np.min(np.stack(sds), axis=0), detection_rate=counts_detection))
    sparse.save_npz(destination / 'logcp10k_all.npz', logx)
    tsv(pd.DataFrame({'gene': genes}), destination / 'genes.tsv')
    tsv(spots, destination / 'score_table.tsv.gz')
    tsv(metadata, destination / 'gene_metadata.tsv.gz')
    tsv(pd.DataFrame(coverages), destination / 'panel_coverage.tsv')
    json_out(sources, destination / 'input_manifest.json')
    json_out(dict(panels=PANELS, technical=TECH, composition=COMP, scoring='mean_gene_z_by_section',
        normalization='log1p_CP10K', proxy_interpretation='transcriptional_proxies_not_cell_fractions',
        qc='in_tissue=1,total_UMI>=100,detected_unique_symbols>=100,finite_coordinates',
        frozen_score_table_used=False), destination / 'pipeline_config.json')
    print('Primary prepared:', len(spots), 'spots', flush=True)


def pitch_scaled_coordinates(spots):
    """Transform geometry only; exported capture positions remain image pixels."""
    from sklearn.neighbors import NearestNeighbors
    transformed = spots.copy().astype({'x': float, 'y': float})
    records = []
    for section, group in spots.groupby('sample_id', sort=True):
        xy = group[['x', 'y']].to_numpy(float)
        if len(group) < 2 or not np.isfinite(xy).all():
            raise ValueError('Pitch scaling requires two finite positions per section')
        distances = NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)[0][:, 1]
        pitch = float(np.median(distances))
        if not pitch > 0:
            raise ValueError('Nonpositive median nearest-other-spot pitch')
        transformed.loc[group.index, ['x', 'y']] = xy / pitch
        records.append(dict(sample_id=section, full_resolution_pixels_per_pitch=pitch,
            bootstrap_coordinate_units='section_specific_spot_pitch', x_origin=0, y_origin=0))
    return transformed, pd.DataFrame(records)


def run_dataset(folder, reps=1000, widths=None, n_sets=500, annotation_root=None,
                seed=SEED, bootstrap_in_pitches=False):
    folder = Path(folder)
    results = folder / 'Results'
    results.mkdir(exist_ok=True)
    spots = pd.read_csv(folder / 'score_table.tsv.gz', sep='\t')
    units = 'full_resolution_pixels'
    tick = time.perf_counter()
    associations = core.associations(spots, 'target_score', 'reference_score', TECH, COMP)
    tsv(associations, results / 'associations.tsv')
    ranking_frames = []
    for layer, covs in [('technical', TECH), ('composition', TECH + COMP)]:
        _, a, b = core.fit_pair(spots, 'target_score', 'reference_score', covs)
        ranked = core.spatial_ranking(spots, a, b, ranking_method='positive_product',
            target_scores=spots.target_score.to_numpy(), reference_scores=spots.reference_score.to_numpy())
        ranked['coordinate_units'] = units
        tsv(ranked, results / f'{layer}_coordinates.tsv.gz')
        tsv(core.ranking_comparison(ranked), results / f'{layer}_baseline_comparison.tsv')
        rows = ranked.groupby('sample_id').agg(n_spots=('barcode', 'size'),
                    n_candidates=('cadspot_candidate', 'sum')).reset_index()
        rows['layer'] = layer
        ranking_frames.append(rows)
    tsv(pd.concat(ranking_frames), results / 'candidate_counts.tsv')
    bootstrap_spots = spots
    if bootstrap_in_pitches:
        bootstrap_spots, transform = pitch_scaled_coordinates(spots)
        tsv(transform, results / 'bootstrap_coordinate_transforms.tsv')
        if widths is None:
            widths = [6, 12, 18]
        pitch = None
    elif widths is None:
        # Pixel scales differ between studies; widths are declared multiples of
        # median nearest-neighbor pitch, not transferred absolute pixel values.
        from sklearn.neighbors import NearestNeighbors
        pitches = []
        for _, sub in spots.groupby('sample_id'):
            distances = NearestNeighbors(n_neighbors=2).fit(sub[['x', 'y']]).kneighbors()[0]
            pitches.append(float(np.median(distances[:, 0])))
        pitch = float(np.median(pitches))
        widths = [6*pitch, 12*pitch, 18*pitch]
    else:
        pitch = None
    if reps:
        for width in widths:
            label = f'{width:.3f}'.rstrip('0').rstrip('.')
            print('Bootstrap', folder.name, width, flush=True)
            ci, draws, audit, blocks = core.block_bootstrap(bootstrap_spots, 'target_score', 'reference_score',
                TECH, COMP, block_size=width, n_boot=reps, seed=seed)
            tsv(ci, results / f'bootstrap_{label}_ci.tsv')
            tsv(draws, results / f'bootstrap_{label}_draws.tsv.gz')
            tsv(audit, results / f'bootstrap_{label}_audit.tsv.gz')
            tsv(blocks[['sample_id', 'barcode', 'block_x', 'block_y']], results / f'blocks_{label}.tsv.gz')
    if n_sets:
        genes = pd.read_csv(folder / 'genes.tsv', sep='\t').gene.tolist()
        logx = sparse.load_npz(folder / ('logcp10k_all.npz' if (folder / 'logcp10k_all.npz').exists() else 'logcp10k.npz'))
        expression = null.ExpressionData(logx.T.tocsr(), genes, spots)
        coverage = pd.read_csv(folder / 'panel_coverage.tsv', sep='\t')
        available = {}
        for name, intended in [('target_score', TARGET), ('reference_score', REFERENCE)]:
            # Score membership must be identical across sections for a pooled
            # matched-set audit. Fail explicitly if coverage differs.
            used = [tuple(sub.loc[sub.status.eq('used'), 'gene']) for _, sub in coverage[coverage.panel.eq(name)].groupby('sample_id')]
            if len(set(used)) != 1:
                raise ValueError('Panel membership differs across sections; audit separately')
            available[name] = list(used[0])
        meta = pd.read_csv(folder / 'gene_metadata.tsv.gz', sep='\t')
        annotation_root = annotation_root or ROOT
        length_path = ROOT / 'Provenance/ensembl113_exon_union_lengths.tsv'
        if not length_path.exists():
            length_path = annotation_root / 'Source_Data/ensembl113_exon_union_lengths.tsv'
        lengths = pd.read_csv(length_path, sep='\t')
        meta = meta.merge(lengths[['gene', 'gene_length']], on='gene', how='left', validate='one_to_one')
        variable = meta.min_section_sd.gt(1e-12) if 'min_section_sd' in meta else meta.sd_log_expression.gt(1e-12)
        universe = meta.loc[variable & meta.gene_length.gt(0) & meta.detection_rate.ge(.01), 'gene'].tolist()
        exclusions = sorted({g for p in COMP for g in PANELS[p]})
        matched = null.generate_matched_sets(expression, meta, available['target_score'], available['reference_score'],
            universe, exclusions=exclusions, n_sets=n_sets, seed=seed,
            mean_tolerance=.25, detection_tolerance=.05, length_tolerance=.25)
        summary, rhos, scores = null.calibrate_gene_sets(expression, spots,
            available['target_score'], available['reference_score'], matched, TECH, COMP,
            scoring='mean_gene_z_by_section')
        main = associations.query("scope == 'pooled' and layer == 'composition_adjusted' and method == 'score_residual_spearman'").spearman_rho.iloc[0]
        if not np.isclose(summary.observed_rho.iloc[0], main, atol=1e-10):
            raise ValueError('Null observed statistic differs from primary: scoring mismatch')
        for frame, name in [(matched, 'matched_assignments.tsv.gz'), (summary, 'matched_summary.tsv'),
                            (rhos, 'matched_rhos.tsv'), (meta, 'matching_metadata.tsv.gz')]:
            tsv(frame, results / name)
        tsv(pd.DataFrame({'gene': universe}), results / 'matching_universe.tsv')
    core.write_manifest(results / 'run_manifest.json', dict(input_sha256=core.file_sha256(folder / 'score_table.tsv.gz'),
        primary_estimator='score_residual_spearman', scoring='mean_gene_z_by_section', ranking='positive_product',
        technical=TECH, composition=COMP, bootstrap_reps=reps, bootstrap_widths=widths,
        median_nearest_neighbor_pitch=pitch, coordinate_units=units, n_sets=n_sets, seed=seed,
        bootstrap_coordinate_units='section_specific_spot_pitch' if bootstrap_in_pitches else units,
        frozen_score_table_used=False, runtime_seconds=time.perf_counter()-tick))
    print(associations.query("scope == 'pooled'").to_string(index=False), flush=True)
    print(pd.concat(ranking_frames).to_string(index=False), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['download-primary', 'prepare-primary', 'download-pdac', 'prepare-pdac', 'prepare-external', 'run'], required=True)
    p.add_argument('--source-root', type=Path, default=ROOT)
    p.add_argument('--folder', type=Path, default=ROOT / 'Primary')
    p.add_argument('--reps', type=int, default=1000)
    p.add_argument('--widths', nargs='+', type=float)
    p.add_argument('--n-sets', type=int, default=500)
    p.add_argument('--seed', type=int, default=SEED)
    p.add_argument('--bootstrap-in-pitches', action='store_true')
    p.add_argument('--prepared-external', type=Path)
    args = p.parse_args()
    if args.stage in ('download-primary', 'download-pdac'):
        destination = args.source_root / 'Data'
        destination.mkdir(parents=True, exist_ok=True)
        path = destination / 'GSE274103_RAW.tar'
        if not path.exists():
            from download_ovarian import download_ranges
            download_ranges('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274103/suppl/GSE274103_RAW.tar', path)
        raw = destination / 'GSE274103/raw'
        raw.mkdir(parents=True, exist_ok=True)
        with tarfile.open(path) as archive:
            archive.extractall(raw, filter='data')
    elif args.stage in ('prepare-primary', 'prepare-pdac'):
        prepare_primary(args.source_root, args.folder)
    elif args.stage == 'prepare-external':
        source = args.prepared_external
        if source is None:
            raise ValueError('--prepared-external required')
        x = sparse.load_npz(source / 'counts.npz').tocsr().astype(np.float64)
        genes = pd.Index(pd.read_csv(source / 'genes.tsv', sep='\t').gene)
        positions = pd.read_csv(source / 'positions.tsv', sep='\t').drop(columns=['sample_id'], errors='ignore')
        barcodes = pd.read_csv(source / 'barcodes.tsv', sep='\t').barcode
        if list(positions.barcode) != list(barcodes):
            raise ValueError('Prepared external barcode order mismatch')
        prepare_section(x, genes, positions, source.name, args.folder)
    else:
        run_dataset(args.folder, args.reps, args.widths, args.n_sets, args.source_root,
                    args.seed, args.bootstrap_in_pitches)


if __name__ == '__main__':
    main()
