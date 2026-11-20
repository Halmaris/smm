# Reproducing the experiments

Commands below run from `code/`, with the environment in [ENVIRONMENT.md](ENVIRONMENT.md).
Inputs and outputs stay in local directories; neither is distributed with the code.
Dataset sources and matrix layouts are in [DATA_SOURCES.md](DATA_SOURCES.md).

## Experiment map

| Experiment | Implementation and configuration | Saved-result analysis |
| --- | --- | --- |
| Temporal benchmark | `fit_method`, `fit_ovr`, `fit_direct`; `configs/quality.json` records the direct-method jobs | `quality` for each metric |
| Grayscale and RGB images | Spatial models; three separate adjacency grids for RGB; `configs/quality.json`, `rgb_direct.json`, `beans.json`, `retina.json` | `predictions`, then `quality` within a common method set |
| Binary simulations | `simulate_smm`; `configs/binary.json` records the direct-method jobs | `binary` |
| Multiclass ablations | `paper_design`, `mc_inputs`, `run_ablation`; `configs/multiclass.json` | `multiclass` |
| Training times | Selected TRAIN refits, one worker and one numerical thread; one warm-up followed by five measured repetitions for the 85 direct/RGB pairs | `timing` |

The README gives model-fitting and tuning examples. JSON files record historical
job seeds, candidate grids and input checksums. Shared grids are stored once under
`templates`, with job identifiers selecting the corresponding grid. These files
are experiment metadata; the public functions use their documented R arguments.
Use the saved inputs, including their folds, when rerunning a historical block.
For independent reproduction from a dataset release, follow the preprocessing
and split instructions in DATA_SOURCES.md. Checksums identify the serialized
historical inputs, including fold assignments; reserializing inputs can change a checksum.

```sh
python3 scripts/verify_inputs.py configs/quality.json /path/to/quality/data
python3 scripts/verify_inputs.py configs/rgb_direct.json /path/to/rgb/data
```

The original multiclass design has 16 scenarios, 50 replications and four variants,
with 100 TPE proposals per variant. M00/M10/M01/M11 allow neither penalty, nuclear,
adjacency, or both. Fixed anchors are included in each 100-proposal budget;
remaining proposals use TPE and optimistic-bound pruning. The seed rule and
independent streams are listed in `configs/multiclass.json`.

To execute one block explicitly, supply `INPUT_ROOT/blocks/s01_r001/inputs.rds`
containing `train$X`, `train$y`, `train$folds`, `test$X` and `test$y`:

```sh
Rscript --vanilla scripts/run_multiclass.R --execute s01 1 /path/to/INPUT_ROOT /path/to/results
```

The command uses the compact model implementation in this repository. Historical
runtime measurements refer to their recorded solver implementations and hardware.
The optional follow-up experiments are documented separately in [SIMULATIONS.md](SIMULATIONS.md).

## Temporal OvR timing initializer

The timed temporal OvR refits used the archived [initializer](../timing/library_initialization.R).
It checks the nuclear-norm optimality condition for a zero start class by class.
For linear SVM starts it uses inverse-frequency class weights, cost `Cbar / n`,
tolerance `1e-8`, no internal scaling or shrinking, and a 256 MB LIBSVM cache.
The compact OvR fitter uses a zero start. To obtain the original ADMM starting
state for a training partition, run from `code/`:

```r
source('timing/context.R')
source('timing/library_initialization.R')
context <- smm_fit_context(X_train, y_train)
start <- smm_library_initial(context, list(Cbar = C, tau = tau, gamma = gamma))
```

Supply training observations as a list of channel-by-time matrices, or as a
channel-by-time-by-observation array, with one class label per observation.
Use the selected values of `C`, `tau` and `gamma` for that refit. Create a fresh
context for each independent refit. The return value contains the primal and dual
starting variables for the historical solver. The two files use `digest`,
`Matrix` and `e1071` from the project environment.

## Analysing saved results

These commands read CSV files and fit no models. Use a new or empty output directory.
Outputs are numerical summaries, not publication tables or figures.

```sh
Rscript --vanilla scripts/analyse_saved.R predictions /path/to/predictions.csv /path/to/results/metrics
Rscript --vanilla scripts/analyse_saved.R quality /path/to/BA_scores.csv /path/to/results/BA
Rscript --vanilla scripts/analyse_saved.R quality /path/to/F1_scores.csv /path/to/results/F1
Rscript --vanilla scripts/analyse_saved.R binary /path/to/binary.csv /path/to/results/binary
Rscript --vanilla scripts/analyse_saved.R multiclass /path/to/replicate_results.csv /path/to/results/multiclass
Rscript --vanilla scripts/analyse_saved.R timing /path/to/measured_times.csv /path/to/results/timing
```

| Mode | Required CSV columns |
| --- | --- |
| `predictions` | `dataset,method,id,truth,prediction`; one TEST prediction per observation and method, matching observations across methods |
| `quality` | `dataset` followed by numeric method columns; one row per dataset and one metric per file |
| `binary` | `scenario,replication,method,balanced_accuracy,macro_f1`; replications 1–50 per scenario and method |
| `multiclass` | `scenario,replication,method,BA,macro_F1`; s01–s16, replications 1–50, M00/M10/M01/M11 |
| `timing` | `dataset,method,repeat,seconds`; measured repetitions only, with warm-ups in a separate file |

BA and macro-F1 inputs are proportions in [0,1]. Prediction labels must include
all assessed classes in `truth`. Quality analysis returns mean scores, mean ranks,
Friedman/Iman–Davenport tests and Bergmann–Hommel pairwise comparisons. Its small
rank-inference helper is compiled with Rcpp; it does not load the model solver.
Multiclass analysis uses paired replication differences, pointwise 95% t intervals
in percentage points and Holm adjustment across the 32 primary conditional BA
comparisons. Binary analysis returns per-scenario means and 95% t intervals.
Timing analysis returns medians and repetition counts without pooling different
protocols. A constant paired difference has no estimated t-test variance and is
reported with an unavailable p-value.
