# Coordinate decision audit and fixed-budget evaluation

The extension answers which positions enter a limited inspection budget and how
that choice depends on model/input assumptions. It does not automatically choose
the largest-rho model, claim consensus correctness, or supply coordinate FDR.

## One-command score-table audit

Run from the release root with `Software/requirements.txt`:

```bash
python Software/decision_audit.py --input Primary/score_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_proxy ECM_proxy Contractile_proxy Myeloid_proxy Tcell_proxy --coordinate-units full_resolution_pixels --output-prefix my_ovarian_audit
```

Outputs: `.coordinates.tsv.gz`, `.summary.tsv`, `.manifest.json`, `.report.html`.
HTML is self-contained and printable; submitted examples also have PDF previews.
Budgets default to 10/25/50, with 25 primary. Change `--budgets` and
`--primary-budget` explicitly. Ranks descend by product, with barcode tie breaking;
the unique `*_audit_rank` differs from core minimum-tie `priority_rank`.
Membership requires rank within budget and strictly positive priority. Original
outlier candidates are separate. Raw is a baseline, not part of agreement.
`all_three` means technical/composition/spline membership; `some_models` is the
remaining union; `not_selected` is outside it. Every budget/model flag is kept.
Missing positions and sections with fewer than seven positions are excluded
explicitly. Agreement does not identify biological truth or a winning model.

## Reproduce the extension

```bash
python Decision_Utility/run_evaluation.py budget
python Decision_Utility/run_applications.py
python Decision_Utility/run_evaluation.py perturbations
python Decision_Utility/make_report.py
```

`budget` uses supplied Utility_Benchmark tables and every prior seed, asserting
AP equality against earlier outcomes. All six baselines and three budgets are
kept. Sections are averaged within replicate for paired inference. Inspection
cost concerns known injected labels, not measured laboratory savings.

`run_applications` rebuilds the illustrative IFN/antigen pair from prepared public
logCP10K sections, then exports both pair/cohort audits. `perturbations` needs
public `counts.npz` prepared through Analysis/README.md; large matrices are not
in the small supplement. Count-thinning membership stays fixed: six ovarian
50%-count repetitions become noncomputable when CA9 loses variation. Seeds and
reasons remain in `perturbation_qc.tsv`, without silent repair or omission.

Supplied score inputs replay full-cohort reports without retrieving raw files.
For the second pair substitute
`Decision_Utility/Program_Pair/ovarian/score_table.tsv.gz` in the audit command.
Its exact weights are under Program_Pair. Protocol, outcomes, dropped genes,
failures and source records are retained. Agreement and input stability measure
different properties and cannot be treated as interchangeable.
