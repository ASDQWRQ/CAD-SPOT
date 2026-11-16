# Task-level utility: ovarian worked example and PDAC transfer

`PROTOCOL.md` and `../Study_Design.md` specify the retrospective design. The
default formula, thresholds and six baselines are fixed. Ovarian outcomes were
new; PDAC and the pure simulation scenarios were previously explored.

```bash
python -B Utility_Benchmark/benchmark.py prepare --data-root .
python -B Utility_Benchmark/benchmark.py accuracy
python -B Utility_Benchmark/benchmark.py equivalence
python -B Utility_Benchmark/benchmark.py performance
python -B Utility_Benchmark/make_report.py
```

Run from the release root. Dependencies include the scientific Python stack,
psutil for timing/memory and Matplotlib for figures. Run performance with no
other numerical jobs active. Inputs are explicit ovarian/PDAC score tables;
generation never imports candidate labels or historical image-model targets.

`accuracy` reruns seven 40-seed score-level scenarios and four 32-seed real
background/geometry combinations against six fixed baselines. AP uses all spots.
No-program AP remains undefined. Real labels are injected score perturbations,
not independently validated biological tissue regions. Every replicate and all
unfavorable comparisons are retained. Paired AP intervals concern these
backgrounds, not patients.

`equivalence` checks complete old/new ranking output equality and independent
OLS-product agreement. `performance` measures five randomized serial subprocess
repeats for ovarian 10,000 spots, full ovarian and full PDAC. Kernel time excludes
downloads, preprocessing, bootstrap and matching. Aggregate process-tree RSS can
double-count shared pages. The reduced-output reference does less work and is
not a functionally equivalent speed comparator. These tests support conditional
task utility and reproducible implementation, not a new estimator, universal
superiority or clinical/biological validation.
