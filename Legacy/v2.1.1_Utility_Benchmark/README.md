# Task-aligned CAD-SPOT utility benchmark

This supplement tests specific claims instead of assuming an overall winner.
Read `PROTOCOL.md` for the prespecified scenarios, seeds, comparators and metrics.
All cases, including unfavorable results, are retained.

## Reproduction

Run from this directory in Python 3.10+ with NumPy, pandas, SciPy, scikit-learn,
psutil, Matplotlib and openpyxl. Use single-thread BLAS for the recorded timing
protocol. The two source score tables are public-derived inputs, not raw reads.

```bash
python -B benchmark.py prepare --data-root ..
python -B benchmark.py accuracy
python -B benchmark.py equivalence
python -B benchmark.py performance
python -B make_report.py
```

`--data-root` points to the directory containing `Primary/` and
`External_EndToEnd/NYU_PDAC1/`. Included `Inputs/` already permits all subsequent
stages without downloading count matrices. Preparation overwrites only these
derived local input copies; result stages overwrite their own tables.

- `Software_baseline/`: unchanged CAD-SPOT v2.1.0 for timing/equality, not an
  alternative primary analysis or historical frozen score table.
- `Software/`: optimized v2.1.1. It changes only buffered neighbor-output writes.
  Neighbor queries, tie handling, all formulas and candidate rules are unchanged.
- `standalone_reference.py`: independently implemented standard OLS + positive
  product; no CAD-SPOT imports. Equal mathematical results are expected, not an
  accuracy win for CAD-SPOT.
- `Results/*_replicates.tsv`: every prespecified outcome, not selected seeds.
- `Results/paired_accuracy_differences.tsv`: paired AP differences and percentile
  intervals across simulation/injection replicates, conditional on these designs.
- `Results/equivalence.json`: exact complete-output old/new equality and
  independent-reference agreement on both real score tables.
- `Results/performance_*`: randomized serial single-machine performance; the
  minimal reference has fewer outputs and is not an equivalent-function speed
  comparator. Reference timing includes its fit and rank, whereas old/new kernel
  timing isolates spatial ranking after the same residual fit.

Seven pure scenarios have 40 paired seeds each. Real backgrounds have 32 seeds
for each of four background/geometry combinations. Injections are score-level,
not gene-count or biological labels; locations are chosen without inspecting
scores. The separate-study disease-identity conflict remains unresolved.

## Scope of conclusions

Conditional local-program recovery can improve over raw/technical products when
broad composition drives co-enrichment. It need not improve if composition
captures the desired biological program, if nuisance is absent, or if the linear
model is misspecified. The adjusted-minimum and spline baselines are retained
even where they outperform the default. The default was not retuned after seeing
these results.

The optimization can support a speed claim only relative to v2.1.0 with identical
outputs. It is not a speed comparison with Giotto, Squidpy or latent spatial
models, and the minimal OLS reference is much faster. No user study, tissue ground
truth, coordinate-level FDR or universal superiority is established.

The first performance round incorrectly sampled the Python launcher alone for
RSS on Windows. It is retained under `Results/invalid_launcher_memory_round/`,
excluded from final summaries, and explicitly documented. Final RSS sums the
launcher and recursive worker tree and can double-count shared pages.
