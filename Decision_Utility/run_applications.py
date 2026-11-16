"""Full-cohort audits and a predeclared, illustrative non-metabolic pair."""
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'Software'))
import cadspot
from decision_audit import export_audit

TECH = ['log_nCount', 'log_nFeature']
COMP = ['Epithelial_proxy', 'ECM_proxy', 'Contractile_proxy', 'Myeloid_proxy', 'Tcell_proxy']
PANELS = {
    'target_score': 'IFIT1 IFIT2 IFIT3 ISG15 MX1 MX2 OAS1 OAS2 OAS3 RSAD2 IFI6 IFI44 IFI44L'.split(),
    'reference_score': 'B2M HLA-A HLA-B HLA-C TAP1 TAP2 TAPBP PSMB8 PSMB9 PSMB10 ERAP1 NLRC5'.split()}


def main():
    if set(PANELS['target_score']) & set(PANELS['reference_score']):
        raise ValueError('Second target/reference panels overlap')
    for label, folder in [('ovarian', ROOT/'Primary'), ('pdac', ROOT/'PDAC_Transfer')]:
        original = pd.read_csv(folder/'score_table.tsv.gz', sep='\t')
        prior_coverage = pd.read_csv(folder/'panel_coverage.tsv', sep='\t')
        proxy_genes = set(prior_coverage.loc[prior_coverage.panel.isin(COMP), 'gene'])
        if proxy_genes & set(sum(PANELS.values(), [])):
            raise ValueError('New programs overlap proxy genes')
        prefix = ROOT/f'Decision_Utility/Results/{label}_stress_audit'
        export_audit(original, prefix, 'target_score', 'reference_score', TECH, COMP,
            input_name=folder.name+'/score_table.tsv.gz',
            input_sha256=cadspot.file_sha256(folder/'score_table.tsv.gz'), coordinate_units='full_resolution_pixels')
        cached = []
        for sample in sorted(original.sample_id.unique()):
            f = folder/sample; x = sparse.load_npz(f/'logcp10k.npz')
            genes = pd.Index(pd.read_csv(f/'genes.tsv', sep='\t').gene)
            values = {}; available = {}
            for panel, members in PANELS.items():
                present = [g for g in members if g in genes]
                v = x[:, genes.get_indexer(present)].toarray()
                values[panel] = (present, v)
                available[panel] = set(g for g, sd in zip(present, v.std(axis=0)) if sd > 1e-12)
            cached.append((sample, values, available))
        used = {p:[g for g in members if all(g in c[2][p] for c in cached)] for p, members in PANELS.items()}
        if any(len(v) < 3 for v in used.values()):
            raise ValueError('Insufficient common nonconstant second-pair genes')
        output = ROOT/f'Decision_Utility/Program_Pair/{label}'; output.mkdir(parents=True, exist_ok=True)
        frames, parameters = [], []
        for sample, values, _ in cached:
            d = original[original.sample_id.eq(sample)].copy()
            for panel, intended in PANELS.items():
                present, v = values[panel]
                idx = [present.index(g) for g in used[panel]]
                mean, sd = v.mean(axis=0), v.std(axis=0)
                d[panel] = ((v[:, idx]-mean[idx])/sd[idx]).mean(axis=1)
                for g in intended:
                    i = present.index(g) if g in present else None
                    status = 'used' if g in used[panel] else 'missing' if i is None else 'excluded_constant_or_not_common'
                    parameters.append(dict(sample_id=sample, panel=panel, gene=g, status=status,
                        mean_logcp10k=float(mean[i]) if i is not None else np.nan,
                        sd_logcp10k=float(sd[i]) if i is not None else np.nan,
                        weight=1/len(used[panel]) if status == 'used' else 0))
            frames.append(d)
        d = pd.concat(frames, ignore_index=True)
        d.to_csv(output/'score_table.tsv.gz', sep='\t', index=False)
        pd.DataFrame(parameters).to_csv(output/'score_parameters.tsv', sep='\t', index=False)
        pd.DataFrame(parameters)[['sample_id','panel','gene','status']].to_csv(output/'panel_coverage.tsv', sep='\t', index=False)
        cadspot.associations(d, 'target_score', 'reference_score', TECH, COMP).to_csv(output/'associations.tsv', sep='\t', index=False)
        export_audit(d, ROOT/f'Decision_Utility/Results/{label}_ifn_audit', 'target_score', 'reference_score', TECH, COMP,
            input_name='Decision_Utility/Program_Pair/'+label+'/score_table.tsv.gz',
            input_sha256=cadspot.file_sha256(output/'score_table.tsv.gz'), coordinate_units='full_resolution_pixels')
        manifest = dict(cohort=label, n_sections=d.sample_id.nunique(), n_spots=len(d), intended_panels=PANELS,
            used_panels=used, interpretation='illustrative panels; not complete Reactome sets or validated activation signatures',
            biological_motivation=['https://reactome.org/content/detail/R-HSA-909733', 'https://reactome.org/content/detail/R-HSA-983170'],
            retrieved='2026-11-16', coverage='common nonconstant genes across cohort, without association-based selection',
            upstream_source=str(folder/'input_manifest.json'), upstream_sha256=cadspot.file_sha256(folder/'input_manifest.json'))
        (output/'input_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print(label, 'second-pair coverage', {p:len(v) for p, v in used.items()}, flush=True)


if __name__ == '__main__':
    main()
