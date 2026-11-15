"""Prespecified utility experiments with all outcomes retained."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from sklearn.metrics import average_precision_score

from standalone_reference import fit_residuals,rank_reference,standardize_sections

ROOT=Path(__file__).resolve().parent
DEFAULT_DATA=ROOT.parent if (ROOT.parent/'Primary').exists() else Path(r'D:\CADSPOT_BMC_INTEGRATED_20261115')
METHODS=['raw_product','technical_product','cadspot_product','adjusted_target','adjusted_minimum','spline_product']
CASES=['no_nuisance','linear_compact','linear_dispersed','noisy_proxy','nonlinear_compact','architecture_coupled','no_program']
TECH=['log_nCount','log_nFeature']
COMP=['Epithelial_proxy','ECM_proxy','Contractile_proxy','Myeloid_proxy','Tcell_proxy']


def core(name='Software'):
    path=ROOT/name/'cadspot.py';spec=importlib.util.spec_from_file_location('core_'+name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def write_json(data,path):path.write_text(json.dumps(data,indent=2,default=str)+'\n',encoding='utf-8')


def z(x):
    x=np.asarray(x,float);return (x-x.mean())/(x.std() or 1.)


def generate(seed,case):
    rng=np.random.default_rng(seed);ny,nx=24,28;yy,xx=np.mgrid[:ny,:nx];frames=[];truths=[]
    for section in range(2):
        arch=z((xx/nx+.35*gaussian_filter(rng.normal(size=(ny,nx)),2)).ravel())
        depth=.6*arch+.8*z(rng.normal(size=ny*nx))
        stroma=-.4*arch+z(gaussian_filter(rng.normal(size=(ny,nx)),2).ravel())
        noisy=1. if case=='noisy_proxy' else .25
        measured=arch+noisy*rng.normal(size=ny*nx)
        cx,cy=rng.uniform(7,21),rng.uniform(6,18)
        truth=(((xx-cx)**2+(yy-cy)**2)<=4.5**2).ravel()
        if case=='linear_dispersed':truth=rng.permutation(truth)
        elif case=='architecture_coupled':truth=arch>=np.quantile(arch,.90)
        elif case=='no_program':truth[:]=False
        nuisance=1.5*arch
        if case=='nonlinear_compact':nuisance=2.*z(np.maximum(arch-.25,0)**2)
        technical=.6*depth+.25*stroma
        if case=='no_nuisance':nuisance=np.zeros(len(arch));technical=np.zeros(len(arch))
        noise_a=.8*z(rng.normal(size=ny*nx))+.2*z(gaussian_filter(rng.normal(size=(ny,nx)),1).ravel())
        noise_b=.8*z(rng.normal(size=ny*nx))+.2*z(gaussian_filter(rng.normal(size=(ny,nx)),1).ravel())
        frames.append(pd.DataFrame(dict(sample_id=f's{section}',barcode=[f'b{i:04d}' for i in range(ny*nx)],
            x=xx.ravel(),y=yy.ravel(),depth=depth,composition=measured,stroma=stroma,
            target_score=nuisance+technical+1.2*truth+noise_a,
            reference_score=.9*nuisance+.8*technical+1.1*truth+noise_b)))
        truths.append(truth)
    return pd.concat(frames,ignore_index=True),np.concatenate(truths)


def scores(table,technical,composition,engine):
    raw=standardize_sections(table,table[['target_score','reference_score']].to_numpy())
    tr=standardize_sections(table,fit_residuals(table,technical))
    independently_fitted=fit_residuals(table,technical+composition)
    _,ca,cb=engine.fit_pair(table,'target_score','reference_score',technical+composition)
    fitted=np.column_stack([ca,cb])
    np.testing.assert_allclose(fitted,independently_fitted,atol=1e-10,rtol=1e-10)
    cr=standardize_sections(table,fitted)
    _,sa,sb=engine.fit_pair(table,'target_score','reference_score',technical+composition,method='spline_score_residual_spearman')
    sp=standardize_sections(table,np.column_stack([sa,sb]))
    product=lambda v:np.maximum(v[:,0],0)*np.maximum(v[:,1],0)
    return dict(raw_product=product(raw),technical_product=product(tr),cadspot_product=product(cr),
        adjusted_target=cr[:,0],adjusted_minimum=np.minimum(cr[:,0],cr[:,1]),spline_product=product(sp))


def metrics(table,truth,score):
    k=int(np.ceil(.05*len(table)))
    order=np.lexsort((table.barcode.astype(str),table.sample_id.astype(str),-score))[:k]
    selected=np.zeros(len(table),bool)
    for idx in table.groupby('sample_id',sort=True).indices.values():
        p=score[idx];center=np.median(p);scale=1.4826*np.median(np.abs(p-center))
        if scale<=0:scale=p.std() or 1.
        selected[idx]=((p-center)/scale>=3)&(p>0)
    return dict(average_precision=float(average_precision_score(truth,score)) if truth.any() else np.nan,
        top5_precision=float(truth[order].mean()) if truth.any() else np.nan,
        top5_recall=float(truth[order].sum()/truth.sum()) if truth.any() else np.nan,
        candidate_fraction=float(selected.mean()),n_spots=len(table),n_true=int(truth.sum()))


def prepare_inputs(data_root):
    out=ROOT/'Inputs';out.mkdir(exist_ok=True);records=[]
    for label,source in [('primary',data_root/'Primary/score_table.tsv.gz'),
                         ('external',data_root/'External_EndToEnd/NYU_PDAC1/score_table.tsv.gz')]:
        dest=out/f'{label}.tsv.gz';shutil.copy2(source,dest)
        records.append(dict(dataset=label,source=str(source),sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))
    d=pd.read_csv(out/'primary.tsv.gz',sep='\t')
    d.iloc[np.sort(np.random.default_rng(661001).choice(len(d),10000,replace=False))].to_csv(out/'primary_10000.tsv.gz',sep='\t',index=False)
    write_json(records,out/'sources.json')


def simulations(engine):
    rows=[]
    for case in CASES:
        for rep in range(40):
            d,truth=generate(881001+rep,case)
            for method,values in scores(d,['depth'],['composition','stroma'],engine).items():
                rows.append(dict(experiment='simulation',scenario=case,replicate=rep,seed=881001+rep,method=method,**metrics(d,truth,values)))
            if rep%10==0:print('simulation',case,rep,flush=True)
        pd.DataFrame(rows).to_csv(ROOT/'Results/simulation_replicates.tsv',sep='\t',index=False)


def injections(engine):
    rows=[]
    for background in ['primary','external']:
        d=pd.read_csv(ROOT/f'Inputs/{background}.tsv.gz',sep='\t')
        if background=='primary':
            names=sorted(d.sample_id.unique())[:2]
            d=pd.concat([sub.sort_values('barcode').iloc[::3] for name,sub in d.groupby('sample_id') if name in names],ignore_index=True)
        d[['target_score','reference_score']]=standardize_sections(d,d[['target_score','reference_score']].to_numpy())
        for geometry in ['compact','dispersed']:
            for rep in range(32):
                rng=np.random.default_rng(771001+rep);truth=np.zeros(len(d),bool)
                for idx in d.groupby('sample_id',sort=True).indices.values():
                    k=int(np.ceil(.05*len(idx)))
                    if geometry=='compact':
                        center=int(rng.choice(idx));xy=d[['x','y']].to_numpy()
                        distances=((xy[idx]-xy[center])**2).sum(axis=1)
                        chosen=idx[np.lexsort((d.barcode.to_numpy()[idx],distances))[:k]]
                    else:chosen=rng.choice(idx,k,replace=False)
                    truth[chosen]=True
                table=d.copy();table.target_score+=1.2*truth;table.reference_score+=1.1*truth
                for method,values in scores(table,TECH,COMP,engine).items():
                    rows.append(dict(experiment='injection',scenario=f'{background}_{geometry}',replicate=rep,
                        seed=771001+rep,method=method,**metrics(table,truth,values)))
                if rep%10==0:print('injection',background,geometry,rep,flush=True)
            pd.DataFrame(rows).to_csv(ROOT/'Results/injection_replicates.tsv',sep='\t',index=False)


def summarize():
    frames=[pd.read_csv(ROOT/f'Results/{name}_replicates.tsv',sep='\t') for name in ['simulation','injection']]
    data=pd.concat(frames,ignore_index=True)
    summary=data.groupby(['experiment','scenario','method'])[['average_precision','top5_precision','top5_recall','candidate_fraction']].agg(['mean','std'])
    summary.columns=['_'.join(col) for col in summary.columns]
    summary.reset_index().to_csv(ROOT/'Results/accuracy_summary.tsv',sep='\t',index=False)
    rng=np.random.default_rng(551001);rows=[]
    for (experiment,scenario),sub in data.groupby(['experiment','scenario'],sort=True):
        wide=sub.pivot(index='replicate',columns='method',values='average_precision')
        for baseline in [m for m in METHODS if m!='cadspot_product']:
            delta=(wide.cadspot_product-wide[baseline]).dropna().to_numpy()
            draws=rng.choice(delta,size=(10000,len(delta)),replace=True).mean(axis=1) if len(delta) else np.array([np.nan])
            rows.append(dict(experiment=experiment,scenario=scenario,comparison='cadspot_product-minus-'+baseline,
                n_valid=len(delta),mean_difference=float(delta.mean()) if len(delta) else np.nan,
                ci_low=float(np.quantile(draws,.025)),ci_high=float(np.quantile(draws,.975)),
                wins=int((delta>0).sum()),losses=int((delta<0).sum())))
    pd.DataFrame(rows).to_csv(ROOT/'Results/paired_accuracy_differences.tsv',sep='\t',index=False)


def check_equivalence():
    old,new=core('Software_baseline'),core();records=[]
    for background in ['external','primary']:
        d=pd.read_csv(ROOT/f'Inputs/{background}.tsv.gz',sep='\t')
        _,a,b=old.fit_pair(d,'target_score','reference_score',TECH+COMP)
        before=old.spatial_ranking(d,a,b,target_scores=d.target_score,reference_scores=d.reference_score)
        after=new.spatial_ranking(d,a,b,target_scores=d.target_score,reference_scores=d.reference_score)
        pd.testing.assert_frame_equal(before,after,check_exact=True)
        reference=rank_reference(d,TECH+COMP)
        maxerr=max(float(np.max(np.abs(reference[col]-after[col]))) for col in ['target_resid','reference_resid','target_resid_z','reference_resid_z','priority_score'])
        assert maxerr<1e-10
        assert np.array_equal(reference['cadspot_candidate'],after.cadspot_candidate)
        records.append(dict(dataset=background,spots=len(d),full_old_new_exact=True,
            reference_max_abs_error=maxerr,reference_same_candidates=True,candidates=int(after.cadspot_candidate.sum())))
        print('equivalent',background,len(d),flush=True)
    # Ties, duplicates, shuffled identifiers and missing coordinates, all modes.
    d,truth=generate(881001,'linear_compact');d=d.iloc[:160].copy()
    d.loc[1,['x','y']]=d.loc[0,['x','y']].to_numpy();d.loc[2,'x']=np.nan
    _,a,b=old.fit_pair(d,'target_score','reference_score',['depth','composition','stroma'])
    for mode in old.RANKING_METHODS:
        pd.testing.assert_frame_equal(old.spatial_ranking(d,a,b,ranking_method=mode),
            new.spatial_ranking(d,a,b,ranking_method=mode),check_exact=True)
    write_json(dict(real_tables=records,edge_cases='exact_all_three_ranking_modes',status='passed'),ROOT/'Results/equivalence.json')


def worker(args):
    table=pd.read_csv(ROOT/f'Inputs/{args.dataset}.tsv.gz',sep='\t')
    if args.build=='reference':
        start=time.perf_counter();out=rank_reference(table,TECH+COMP);elapsed=time.perf_counter()-start
        count=int(out['cadspot_candidate'].sum())
    else:
        engine=core('Software_baseline' if args.build=='v2.1.0' else 'Software')
        _,a,b=engine.fit_pair(table,'target_score','reference_score',TECH+COMP)
        start=time.perf_counter();out=engine.spatial_ranking(table,a,b);elapsed=time.perf_counter()-start
        count=int(out.cadspot_candidate.sum())
    print(json.dumps(dict(kernel_seconds=elapsed,spots=len(table),candidates=count)),flush=True)


def performance():
    import psutil
    jobs=[(dataset,build,rep) for dataset in ['external','primary_10000','primary'] for build in ['v2.1.0','v2.1.1','reference'] for rep in range(5)]
    order=np.random.default_rng(991001).permutation(len(jobs));rows=[]
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    for number,i in enumerate(order):
        dataset,build,rep=jobs[i];start=time.perf_counter();peak=0
        proc=subprocess.Popen([sys.executable,'-B',str(Path(__file__).resolve()),'worker','--dataset',dataset,'--build',build],cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        ps=psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                members=[ps]+ps.children(recursive=True);rss=0
                for member in members:
                    try:rss+=member.memory_info().rss
                    except psutil.NoSuchProcess:pass
                peak=max(peak,rss)
            except psutil.NoSuchProcess:pass
            time.sleep(.02)
        stdout,stderr=proc.communicate();wall=time.perf_counter()-start
        if proc.returncode:raise RuntimeError(stderr)
        values=json.loads(stdout.strip().splitlines()[-1])
        rows.append(dict(dataset=dataset,build=build,replicate=rep,execution_order=number,
            process_wall_seconds=wall,sampled_peak_rss_mib=peak/1024**2,**values))
        pd.DataFrame(rows).to_csv(ROOT/'Results/performance_replicates.tsv',sep='\t',index=False)
        print('performance',number+1,'/',len(jobs),dataset,build,round(values['kernel_seconds'],3),flush=True)
    data=pd.DataFrame(rows);summary=data.groupby(['dataset','build'])[['kernel_seconds','process_wall_seconds','sampled_peak_rss_mib']].agg(['median','min','max'])
    summary.columns=['_'.join(v) for v in summary.columns]
    summary.reset_index().to_csv(ROOT/'Results/performance_summary.tsv',sep='\t',index=False)
    write_json(dict(python=sys.version,platform=platform.platform(),processor=platform.processor(),
        logical_cpus=psutil.cpu_count(),physical_cpus=psutil.cpu_count(logical=False),memory_gib=psutil.virtual_memory().total/1024**3,
        openblas_threads=1,omp_threads=1,sampling_interval_seconds=.02,
        memory_measure='Sampled aggregate RSS of the launcher plus recursive worker process tree; shared pages may be counted more than once',
        caveat='Single-machine randomized serial subprocess measurements; not controlled cross-platform or competitor-platform speed evidence'),ROOT/'Results/performance_environment.json')


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','accuracy','equivalence','performance','worker','summarize'])
    p.add_argument('--data-root',type=Path,default=DEFAULT_DATA);p.add_argument('--dataset');p.add_argument('--build');args=p.parse_args()
    (ROOT/'Results').mkdir(exist_ok=True)
    if args.stage=='prepare':prepare_inputs(args.data_root)
    elif args.stage=='accuracy':simulations(core());injections(core());summarize()
    elif args.stage=='equivalence':check_equivalence()
    elif args.stage=='performance':performance()
    elif args.stage=='worker':worker(args)
    else:summarize()


if __name__=='__main__':main()
