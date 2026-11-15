# CAD-SPOT task-aligned utility benchmark

Protocol fixed before running the experiments in this directory, 15 Nov 2026.
The purpose is to test claims, not select favorable outcomes. All scenarios,
replicates and negative cases will be retained. This is not a clinical validation
or a claim that CAD-SPOT outperforms Giotto, Squidpy or latent spatial models.

## Hypotheses and comparisons

1. Conditional-coordinate recovery: does composition-adjusted ranking recover
   a known injected local program better than raw and technical-only products
   when broad architecture drives both scores? Primary metric: average precision
   over every spot. Secondary: fixed top-5% precision/recall and outlier burden.
2. Boundary cases: no nuisance, noisy composition proxies, nonlinear nuisance,
   no injected program and a biologically architecture-coupled program. These
   may favor baselines and must remain visible. Residualization cannot recover
   an unconditional biological program if that program is represented by a
   covariate being removed.
3. Positive-quadrant ranking is not a novel estimator: an independently written
   OLS + positive-product reference should reproduce the default. Equality is
   expected and will not be described as an accuracy gain.
4. Optimize implementation without changing the algorithm: replace repeated
   per-row pandas assignments in six-neighbor annotation with section-wise
   buffered assignment, retaining the same neighbor queries/tie policy. Compare
   v2.1.0 with the optimized build on identical ranking inputs and outputs.
5. Real expression-derived background: test recovery of known score-level
   injections, with locations selected without inspecting expression scores.
   This measures injected-program recovery, not biological ground truth.

## Locked design

- Pure simulations: 40 paired seeds starting at 881001 (not the prior ranking
  ablation seeds), two 24 x 28 grids, 7 scenarios: no nuisance, linear compact,
  linear dispersed, noisy proxy, nonlinear compact, architecture-coupled biology,
  no-program. Program amplitude target/reference 1.2/1.1, radius 4.5 grid units.
- Methods: raw product, technical-only product, composition OLS product
  (CAD-SPOT default), composition target-only, composition minimum, and additive
  spline-residual product (an analysis sensitivity, not the CLI default rank).
  Identical inputs/covariates and known labels within each replicate. No label
  tuning; all-spot AP, not AP restricted to a favorable eligible quadrant.
- Real backgrounds: every third barcode-ordered spot from the first two primary
  sections, and all retained separate-study sample spots. 32 paired seeds per
  compact/dispersed geometry, starting at 771001; exactly ceil(5%) injected per
  section, target/reference amplitude 1.2/1.1 section SDs. Source disease metadata
  are not used to define truth; the external identity conflict remains unresolved.
- Summaries: every replicate, mean/SD; paired AP difference intervals from 10,000
  resamples of independent simulation/injection replicates, not spots. No P-value
  shopping or inference that these simulation intervals generalize to tissues.
- Performance: five single-thread subprocess repeats per build and dataset,
  randomized order seed 991001. External 1,759 spots, a deterministic 10,000-spot
  primary subset, and the full 23,392-spot primary. Time the ranking kernel, not
  downloads/bootstrap/null calibration; separately record total process wall time
  and sampled peak process RSS. Include a lower-functionality product reference
  for context, not a functionally equivalent speed contest.
- Verify old/new full output equality on real tables and tie/duplicate/missing
  fixtures. No performance result will be called improvement unless outputs
  agree and the measured paired comparison supports it.

No wet-lab validation, user study or speed comparison with unrelated spatial
models is part of this protocol. Those claims remain unsupported.

## Measurement correction before final performance run

The first timing round measured the Python launcher RSS rather than its worker
child on this Windows installation. Those memory values are invalid; the entire
round is retained in `Results/invalid_launcher_memory_round/` and is not pooled
into final statistics. The final round repeats the same 45 prespecified jobs,
order and settings while summing launcher/recursive-child RSS. Shared pages may
be counted more than once, so this is sampled process-tree RSS, not unique heap
allocation. Ranking times in the first round are retained for audit only.
