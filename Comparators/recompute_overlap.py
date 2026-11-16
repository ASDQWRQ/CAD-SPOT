"""Recompute new-candidate overlap with preserved candidate-blind model fits."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
COMP=ROOT/'Analysis/comparators'
RUNS=[COMP/'results/stamp/k10_seed20261115_retry_numpy2',
      COMP/'results/popari/k10_seed20261115_budget50',
      COMP/'results/spahdmap/k10_seed20261115_production']


def main():
    keys=['sample_id','barcode']
    primary=pd.read_csv(ROOT/'PDAC_Transfer/Results/composition_coordinates.tsv.gz',sep='\t')
    original=pd.read_csv(ROOT/'Comparators/inputs/eligible_spots.tsv.gz',sep='\t')
    universe=original[keys+['pxl_col_in_fullres','pxl_row_in_fullres']].copy()
    current=universe.merge(primary[keys+['cadspot_candidate','x','y']],on=keys,validate='one_to_one')
    if len(current)!=len(universe):raise ValueError('Current primary no longer covers model inputs')
    if not (current.x.eq(current.pxl_col_in_fullres)&current.y.eq(current.pxl_row_in_fullres)).all():
        raise ValueError('Physical coordinate mismatch')
    inputs=ROOT/'Comparators/current_inputs';inputs.mkdir(exist_ok=True)
    current.drop(columns=['x','y']).to_csv(inputs/'eligible_spots.tsv.gz',sep='\t',index=False)
    universe.to_csv(inputs/'model_input_keys.tsv.gz',sep='\t',index=False)
    out=ROOT/'Comparators/current_summary'
    command=[sys.executable,str(COMP/'enrichment.py'),'--inputs',str(inputs),
             '--outputs',*[str(p) for p in RUNS],'--out',str(out),
             '--n-perm','9999','--seed','20261115','--block-widths','8','4','12']
    subprocess.run(command,check=True)
    payload=dict(comparison_cohort='PDAC_GSE274103_transfer',primary_spots=len(primary),primary_candidates=int(primary.cadspot_candidate.sum()),
        model_input_spots=len(current),model_input_candidates=int(current.cadspot_candidate.sum()),
        comparison_type='reused_candidate_blind_model_outputs_new_labels_not_retrained',
        command=command,spatial_permutations=9999,seed=20261115,
        methods=['STAMP','POPARI','SpaHDmap'],popari='fixed_budget_50_not_converged',
        spahdmap='representative_section_only_zero_vectors_unassigned')
    (ROOT/'Comparators/current_analysis.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(payload,indent=2),flush=True)


if __name__=='__main__':main()
