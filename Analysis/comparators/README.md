# Preserved candidate-blind spatial-program fits

These optional analyses are contextual overlap tests, not an accuracy benchmark.
The model fits used a fixed 23,013-spot input scope. The current primary has
23,392 spots and 650 candidates; 23,013/646 lie in model-input scope. Recompute
the current new-label intersection/enrichment from the repository root:

```bash
python -B Comparators/recompute_overlap.py
python -B Analysis/make_comparator_figure.py
```

`Comparators/current_summary/` is active. Inherited fields containing `frozen`
refer to the relabelled fixed model input, not the full new primary.

## Actual fits

- STAMP: separate section fits, K=10, seed 20261115, 800 maximum epochs, official
  early stop after 100; completed epochs 140/127/177/150/154.
- POPARI: five-section joint fit, K=10, 5 nonspatial warm-up and 50 spatial
  updates, float64. Relative loading change 0.007573 exceeds the 0.001 criterion:
  this is a budget-limited checkpoint, NOT a converged fit.
- SpaHDmap: first section only, official image/graph/fusion sequence; image
  early stop at 2600, 5000 graph and 2000 fusion steps. Its 343 all-zero native
  vectors remain unassigned without pseudocounts. The current intersection is
  4461 spots/108 candidates; 13 new candidates fall among the zero vectors.

All methods used 3000 pooled top-log-expression-variance genes, without candidate
labels to choose genes or tune settings. Independent STAMP topics are not aligned
across sections. Exact source commits, compatibility patches, settings and
stopping records are retained in run.json and source_versions.json. Package
freezes document optional execution environments, not core installation needs.

## Training from public inputs

Small factor exports and model-input keys are included; checkpoints, external
source trees and runtime environments are not. Restore sources with
`setup_sources.py`, pinned commits and recorded patches, in separate compatible
environments. Stage primary data as in `Analysis/README.md`. From this directory:

```bash
python setup_sources.py
python prepare_inputs.py --raw ../../Data/GSE274103/raw --spatial-root ../../Data/GSE274103/prepared --candidates ../../Comparators/current_inputs/eligible_spots.tsv.gz --out inputs
python run_methods.py stamp --batch-size 512 --out results/stamp/reproduction
python run_methods.py popari --popari-dtype float64 --popari-iterations 50 --out results/popari/reproduction
python run_methods.py spahdmap --batch-size 4 --out results/spahdmap/reproduction
python enrichment.py --outputs results/stamp/reproduction results/popari/reproduction results/spahdmap/reproduction --out summary_reproduction
```

Use fresh output directories. These commands document the training protocol,
not models refitted on the expanded primary universe. Full setup needs optional
dependencies/assets, including SpaHDmap's color.csv. GPU/hardware differences
may prevent bitwise reproducibility.

The overlap null uses 9999 equal-occupancy block-label permutations, widths 8
median pitches (main), 4/12 (sensitivities). BH correction combines finite domain
and loading tests across sections per method/width. Fewer than 8 blocks or less
than half the spots movable gives undefined P/Q, not no association. Only the
first section has three-method scope. Overlap does not establish independent
biology or method superiority.
