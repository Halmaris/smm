# Additional simulations: equal tuning budgets and the smoothness boundary

This optional follow-up design is implemented but excluded from the current
manuscript. No inputs or results are bundled. Training-size/channel extensions
are not implemented here.

## Scientific questions

1. At the same total tuning budget, does jointly tuning both penalties improve
   held-out performance relative to choosing among simpler submodels by CV?
2. Are conclusions on the smooth, low-rank scenario sensitive to the upper
   endpoint of the adjacency-penalty domain?

These are targeted sensitivity experiments, not an additional broad benchmark.

## Frozen inputs and search

Reuse original multiclass development samples, TEST samples, labels and five CV
folds, replications **1–10**. No new data are generated for these comparisons.
Both TPE seeds, **91001 and 91002**, are applied to each scenario–replication block.
They measure search variability and are not additional independent data samples.

| Scenario | Classes | Rank | Roughness | Development observations/class |
|---|---:|---:|---|---:|
| s09 | 3 | 2 | Smooth | 100 |
| s13 | 3 | 2 | Rough | 100 |
| s15 | 3 | 8 | Rough | 100 |

The original generator uses 8 rows, 32 time points and 1,000 TEST observations per
class. `paper_design`, `mc_plan`, `mc_seed` and `mc_inputs` retain its definitions.
The new experiments read the original input files rather than relying on
regeneration across R/platform versions.

The common protocol is Optuna 4.5.0 TPE: 20 startup trials, multivariate/grouped
sampling and no constant liar. Domains are `Cbar ∈ [1e-2, 1e3]`, positive
`tau ∈ [1e-4, 100]`, positive `gamma ∈ [1e-5, 1e4]`, all logarithmic, with explicit
zero branches for allowed penalties. Fixed lexicographically ordered anchors use
`Cbar = .1, 1, 10, 100` and each allowed penalty `0, .1, 1`. Anchor counts are
4/12/12/36 for M00/M10/M01/M11; anchors count toward the budget and are never pruned.

Selection uses **pooled out-of-fold balanced accuracy**, not the unweighted mean
of fold scores. A proposal is pruned only if, even with every unprocessed
prediction correct, its final BA is below the best complete proposal by more than
`1e-12`. Ties continue. Pruned intermediate scores are reported to TPE as in the
original implementation. Within a variant, BA ties within `1e-12` are broken by
`Cbar`, `tau`, `gamma`, then trial number.

## Stage 1: equal-budget comparison

For every scenario–replication–seed combination:

- Joint model: M11, **100 proposals**.
- Portfolio: M00/M10/M01, **34/33/33 proposals**, total **100**.
- Select the portfolio member with highest development CV BA; exact ties use
  M00, then M10, then M01.
- Freeze both choices, then fit/evaluate M11 and the selected member on TRAIN/TEST.

Total: **240 tunings, 12,000 proposals, at most 60,000 CV fits, 120 final refits**.
Anchors are included, not added to these budgets. The two alternatives have equal
proposal budgets, not necessarily equal wall-clock cost or equal numbers of
completed fold fits. Pruning savings can be calculated directly from `trials$folds`.

## Stage 2: adjacency-domain sensitivity

Use s09, the same 10 replications and two seeds. Repeat M01 (33 proposals) and
M11 (100 proposals), changing only the positive gamma endpoint to **1e5**.
The objective, data, folds, other domains and solver tolerances remain unchanged.
Compare expanded versus original domain within each method. The enlarged domain
changes the sampled positive gamma values; it is not a nested extension of the
original list of candidates.

Total: **40 additional tunings, 2,660 proposals, at most 13,300 CV fits**.
There are **40 expanded-domain final refits**. Original M11 assessments are reused.
Original M01 assessments are reused when M01 won the core portfolio. Otherwise
its already-frozen core selection needs one final refit: **0–20 additional refits**.
The combined plan therefore needs **160–180 final refits**, not unconditionally
160. No tuning is repeated merely to obtain this missing assessment.

## Execution

Supply an input root outside the code directory containing:

```text
INPUT_ROOT/blocks/s09_r001/inputs.rds
INPUT_ROOT/blocks/s09_r002/inputs.rds
...
```

Each file is a list with `train$X`, `train$y`, `train$folds`, `test$X`, `test$y`.
The original archive can contain additional fields. A small input-file checksum
binds boundary results to the same frozen inputs as their core result.

Print the plan without fitting:

```r
source('R/load.R')
smm <- load_smm('.', compile = FALSE)
smm$additional_plan('core')
smm$additional_plan('boundary')
```

The following are **manual execution examples**, not commands run during package
preparation. Set the one-thread environment from the README first.

```sh
Rscript --vanilla scripts/run_additional.R --execute core s09 1 91001 /path/to/INPUT_ROOT /path/to/OUTPUT_ROOT
Rscript --vanilla scripts/run_additional.R --execute boundary s09 1 91001 /path/to/INPUT_ROOT /path/to/OUTPUT_ROOT
```

One command executes one scenario–replication–seed block sequentially. The boundary
command requires the matching completed core result. Select other blocks from
`additional_plan`; there is intentionally no background launcher or worker pool.
Existing output files are not overwritten. An interrupted block has no automatic
checkpoint recovery; diagnose the interruption before explicitly rerunning it.

Each completed block writes one RDS object containing tuning traces, frozen
selections and held-out assessments, plus a comparison CSV. Fitted coefficient
objects are not retained. Outputs belong in an external output directory, not in
the source archive.

## Analysis

After completing the corresponding full stage:

```sh
Rscript --vanilla scripts/summarize_additional.R /path/to/OUTPUT_ROOT/core /path/to/SUMMARY/core
Rscript --vanilla scripts/summarize_additional.R /path/to/OUTPUT_ROOT/boundary /path/to/SUMMARY/boundary
```

The summary rejects missing/duplicated planned blocks. Compute each seed-specific
paired TEST difference, then average the two seed differences **within each
replication**. Each scenario thus contributes **10 independent replication means**.

Report BA and macro-F1 differences, ordinary paired 95% t intervals across those
replication means, and two-sided exploratory t tests. The core BA family has
three scenario contrasts with Holm correction; macro-F1 is a separate secondary
family. Boundary analysis has two method contrasts per metric. The confidence
intervals are pointwise, not simultaneous Holm-adjusted intervals. With only ten
replications, these are exploratory sensitivity estimates, not an equivalence or
noninferiority design. A zero estimated variance is reported with an unavailable
t-test p-value rather than fabricated significance.

The code stores all three simpler-model CV traces. TEST performance never selects
the portfolio, a seed, an upper penalty endpoint or a method. A better result in
one sensitivity setting does not retrospectively replace the original benchmark.
