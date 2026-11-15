"""Recompute the documented review analyses; paths are relative to this package."""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'Software'))
import cadspot
import cadspot_null as cn

TECHNICAL = ['log_nCount', 'log_nFeature']
COMPOSITION = ['Epithelial_malignant_proxy','ECM_Stroma','CAF','Macrophage','T_cell']
SEED = 20261115
EXTERNAL_IDS = [6177599,6177601,6177603,6177607,6177609,6177612,6177614,6177617,6177618,6177623]


def run_cli(script, arguments):
    command = [sys.executable, str(ROOT/'Software'/script)] + [str(v) for v in arguments]
    print('Running', script, arguments[0:2], flush=True)
    subprocess.run(command, check=True)


def core(stage, reps):
    source = ROOT/'Source_Data'
    base = ['--input',source/'archived_coordinate_table.tsv.gz','--target','target_score','--reference','reference_score',
            '--technical',*TECHNICAL,'--composition',*COMPOSITION,'--coordinate-units','full_resolution_pixels',
            '--seed',SEED]
    if stage == 'core':
        run_cli('cadspot.py', base + ['--ranking-method','legacy','--top-k',100,
            '--coordinate-output',source/'fresh_coordinates.tsv.gz',
            '--ranking-comparison-output',source/'fresh_ranking_baseline.tsv',
            '--archived-comparison-output',source/'archived_candidate_comparison.tsv',
            '--output',source/'fresh_associations.tsv','--manifest-output',source/'fresh_run_manifest.json'])
    else:
        for width in (1500,3000,4500):
            run_cli('cadspot.py',base + ['--block-size',width,'--bootstrap-reps',reps,
                '--bootstrap-output',source/f'block_bootstrap_{width}_ci.tsv',
                '--bootstrap-draws-output',source/f'block_bootstrap_{width}_draws.tsv',
                '--bootstrap-audit-output',source/f'block_bootstrap_{width}_audit.tsv',
                '--block-map-output',source/f'block_map_{width}.tsv',
                '--output',source/f'bootstrap_{width}_point_estimates.tsv'])


def null():
    source = ROOT/'Source_Data'
    data = ROOT/'Data/GSE274103/prepared'
    base = ['--expression',data/'logcp10k_all.npz','--expression-spots',data/'spots_all.tsv.gz',
            '--genes',data/'genes.tsv','--matrix-orientation','spots-by-genes','--expression-scale','log1p_CP10K',
            '--spot-key-cols','sample_id','barcode','--target-genes',source/'target_available_genes.tsv',
            '--reference-genes',source/'reference_available_genes.tsv',
            '--intended-target-genes',source/'target_genes.tsv','--intended-reference-genes',source/'reference_genes.tsv']
    run_cli('cadspot_null.py',['generate',*base,'--metadata',source/'pdac_gene_matching_metadata.tsv',
        '--universe',source/'pdac_matching_universe.tsv','--mean-tolerance',.25,'--detection-tolerance',.05,
        '--length-tolerance',.25,'--n-sets',500,'--seed',SEED,'--output',source/'fresh_matched_sets_500.tsv'])
    run_cli('cadspot_null.py',['calibrate',*base,'--spot-table',source/'archived_coordinate_table.tsv.gz',
        '--matched-sets',source/'fresh_matched_sets_500.tsv','--technical',*TECHNICAL,'--composition',*COMPOSITION,
        '--scoring','mean_gene_z','--null-output',source/'fresh_matched_null_500_rhos.tsv',
        '--score-output',source/'raw_matrix_pair_scores.tsv','--output',source/'fresh_matched_null_500_summary.tsv'])


def external(external_raw, n_sets):
    import h5py
    from prepare_pdac import TARGET, REFERENCE, decode
    source = ROOT/'Source_Data'
    folder = source/'external_reanalysis'
    folder.mkdir(exist_ok=True)
    lengths = pd.read_csv(source/'ensembl113_exon_union_lengths.tsv',sep='\t')
    summaries = []
    for gsm in EXTERNAL_IDS:
        h5path = next(external_raw.glob(f'GSM{gsm}_*.h5'))
        sample = h5path.name.split('_processed_')[0]
        print('External',sample,flush=True)
        dest = folder/sample
        dest.mkdir(exist_ok=True)
        with h5py.File(h5path) as f:
            m = f['matrix']
            x = sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr()
            names, barcodes = decode(m['features']['name'][:]), decode(m['barcodes'][:])
        lib = np.asarray(x.sum(axis=1)).ravel()
        detected = np.diff(x.indptr)
        keep = (lib>0)&(detected>=3)
        x, lib, detected, barcodes = x[keep],lib[keep],detected[keep],barcodes[keep]
        genes = pd.Index(sorted(set(names)))
        transform = sparse.coo_matrix((np.ones(len(names)),(np.arange(len(names)),genes.get_indexer(names))),shape=(len(names),len(genes))).tocsr()
        x = (x@transform).tocsr().astype(float)
        logx = x.multiply((1e4/lib)[:,None]).tocsr()
        logx.data = np.log1p(logx.data)
        mean = np.asarray(logx.mean(axis=0)).ravel()
        sd = np.sqrt(np.maximum(0,np.asarray(logx.power(2).mean(axis=0)).ravel()-mean**2))
        detection = np.asarray((x>0).sum(axis=0)).ravel()/x.shape[0]
        meta = pd.DataFrame(dict(gene=genes,mean_expression=mean,detection_rate=detection,sd_log_expression=sd)).merge(lengths,on='gene',how='left',validate='one_to_one')
        valid = set(meta.loc[meta.sd_log_expression>0,'gene'])
        target, ref = [g for g in TARGET if g in valid], [g for g in REFERENCE if g in valid]
        universe = meta.loc[meta.gene_length.notna() & (meta.gene_length>0) & (meta.detection_rate>=.01)&(meta.sd_log_expression>0),'gene'].tolist()
        spots = pd.DataFrame(dict(sample_id=sample,barcode=barcodes,log_nCount=np.log1p(lib),log_nFeature=np.log1p(detected)))
        expression = cn.ExpressionData(logx.T.tocsr(),genes.tolist(),spots)
        for label, values in [('target',target),('reference',ref),('universe',universe)]:
            pd.DataFrame({'gene':values}).to_csv(dest/f'{label}_genes.tsv',sep='\t',index=False)
        meta.to_csv(dest/'matching_metadata.tsv.gz',sep='\t',index=False)
        spots.to_csv(dest/'technical_covariates.tsv',sep='\t',index=False)
        # The same prespecified strict calipers are used for every sample.
        matched = cn.generate_matched_sets(expression,meta,target,ref,universe,n_sets=n_sets,seed=SEED,
                                           mean_tolerance=.25,detection_tolerance=.05,length_tolerance=.25)
        summary, controls, scores = cn.calibrate_gene_sets(expression,spots,target,ref,matched,TECHNICAL,
                                                        composition=[],scoring='mean_gene_z')
        for frame, name in [(matched,'matched_sets.tsv'),(controls,'null_rhos.tsv'),(scores,'scores.tsv'),(summary,'summary.tsv')]:
            frame.to_csv(dest/name,sep='\t',index=False)
        row = summary.iloc[0].to_dict()
        row.update(sample_id=sample,n_spots=len(spots),n_universe=len(universe),
                   n_target_available=len(target),n_reference_available=len(ref),
                   technical_covariates='log1p_total_UMI;log1p_detected_features',composition_covariates='not_supplied',
                   analysis_layer='technical_adjusted',coordinates='not_supplied_in_selected_expression_files',
                   spatial_inference=False,seed=SEED,n_random_sets=n_sets,
                   missing_target=';'.join(g for g in TARGET if g not in valid),
                   missing_reference=';'.join(g for g in REFERENCE if g not in valid))
        summaries.append(row)
        pd.DataFrame(summaries).to_csv(source/'external_reanalysis_summary.tsv',sep='\t',index=False)
        (dest/'provenance.json').write_text(json.dumps(dict(source_h5=str(h5path),source_sha256=cadspot.file_sha256(h5path),
            normalization='log1p_CP10K',scoring='mean_gene_z',mean_tolerance=.25,detection_tolerance=.05,length_tolerance=.25,
            matching_universe='available,length_annotated,detected_in_at_least_1_percent,positive_variance',
            length_source='Ensembl113_GRCh38_exon_union',excluded_genes='target_and_reference',seed=SEED),indent=2)+'\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stages',nargs='+',choices=['core','bootstrap','null','external'],required=True)
    ap.add_argument('--reps',type=int,default=1000)
    ap.add_argument('--n-external-sets',type=int,default=500)
    ap.add_argument('--external-raw',type=Path,default=ROOT/'Data/GSE203612/raw')
    args=ap.parse_args()
    for stage in args.stages:
        if stage in ('core','bootstrap'):
            core(stage,args.reps)
        elif stage=='null':
            null()
        else:
            external(args.external_raw,args.n_external_sets)
