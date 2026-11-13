# Structured support matrix machines

Research code for classifying matrix-valued observations with independently
selectable nuclear and temporal/spatial adjacency penalties. Includes
one-versus-rest and direct multiclass models, baseline classifiers, simulation
generators, TPE tuning with deterministic pruning, and paired analyses.

This repository contains **code and documentation only**. Datasets, fitted models,
results, manuscript sources and figures are not included. Data are obtained
separately under each provider's terms.

## Installation

The native solver currently targets **macOS 14 or later on Apple Silicon** and
uses Apple's Accelerate framework. The checked environment is CRAN R 4.5.1 and
CPython 3.12.14, with Xcode Command Line Tools and GNU gettext headers.
See [prerequisites](code/docs/ENVIRONMENT.md#prerequisites) for the required paths.

From the repository root:

```sh
cd code
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
Rscript --vanilla scripts/restore_environment.R

export R_LIBS_USER="$PWD/.renv/library"
export RETICULATE_PYTHON="$PWD/.venv/bin/python"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

.venv/bin/python -m pip check
Rscript --vanilla scripts/check_environment.R
```

`renv.lock` and `requirements.txt` pin the R and Python dependencies. The R
installer selects the checked compiler, SDK and gettext paths for its own
process, restores packages into a local library and leaves system settings
unchanged. Installation and the environment check fit no models. Repeat the
exports in each new terminal; see [environment details](code/docs/ENVIRONMENT.md).

## Local data

Supply observations as a list of equally sized numeric matrices, or a
`p × q × n` array, with class labels in `y`. Temporal matrices have channels in
rows and ordered measurements in columns. RGB inputs stack three image planes
vertically; their adjacency graphs have no edges between channels.

Keep the released TEST partition separate and supply fixed stratified fold
identifiers for tuning. Scaling is learned within each training partition.
[DATA_SOURCES.md](code/docs/DATA_SOURCES.md) lists dataset releases, preprocessing
and the expected input format. Store downloaded inputs under an ignored local
`data/` directory and generated outputs under `results/`.

## Run

From `code/`, after setting the environment above:

```r
source('R/load.R')
smm <- load_smm('.')

# X and y are the development observations and labels; X_test is held out.
context <- smm$smm_context(X, y, structure = 'temporal')
model <- smm$fit_direct(context, Cbar = 10, tau = .1, gamma = 1)
prediction <- smm$predict_model(model, X_test)
```

For spatial inputs use `structure = 'spatial'`; add `channels = 3` for RGB.
`fit_ovr` fits one-versus-rest models, while `fit_method` and `predict_method`
provide a common interface to the comparators. Loading the API compiles the
native backend; model fitting starts only when a fitting function is called.

To tune a direct model using development observations and fixed five-fold splits:

```r
train <- list(X = X, y = y, folds = smm$stratified_folds(y, seed = 91001L))
choice <- smm$tune_direct(train, method = 'M11', seed = 91001L, budget = 100L)
```

The [research API](code/README.md) documents these functions. `paper_design`,
`mc_plan` and `mc_inputs` define the original multiclass simulation; `run_ablation`
runs one supplied four-variant block. The optional equal-budget and penalty-domain
experiments have separate [designs and execution commands](code/docs/SIMULATIONS.md).
Full simulations are substantial computations and require explicit execution.

## Reproducing the paper

[REPRODUCING.md](code/docs/REPRODUCING.md) maps the experiments to code and commands.
Historical job seeds, candidate grids and input checksums are in `code/configs/`.
Saved predictions, quality scores, simulation results and timing measurements can
be analysed with `code/scripts/analyse_saved.R`, without fitting models.

## Methods and configuration

- Nuclear regularization promotes low-rank coefficient matrices. The adjacency
  penalty controls temporal or spatial roughness; either penalty can be zero.
- M00, M10, M01 and M11 allow neither penalty, only nuclear, only adjacency, or
  both penalties, respectively.
- Direct multiclass models use a coupled hinge loss, uniform observation costs
  and no intercept. One-versus-rest models use balanced binary losses and free
  intercepts.
- Tuning uses Optuna 4.5.0 TPE, fixed anchors included in the proposal budget,
  pooled five-fold balanced accuracy and optimistic-bound pruning. TEST never
  selects hyperparameters.
- The native ADMM backend uses FP32 updates with FP64 numerical checks. See
  [METHODS.md](code/docs/METHODS.md) for objectives and numerical criteria.

## Tests

From `code/`, with the same environment:

```sh
.venv/bin/python tests/test_tpe.py
Rscript --vanilla tests/test_math.R
Rscript --vanilla tests/test_tuning.R
```

These check algebra, artificial-score searches and mock workflows without fitting
classifiers. `tests/test_models.R` separately fits small synthetic models.
[VALIDATION.md](code/docs/VALIDATION.md) describes the checks and tested versions.

## Repository contents

| Path | Contents |
| --- | --- |
| `code/R/` | Inputs, penalties, classifiers, tuning, generators and paired analyses |
| `code/src/` | Native FP32 ADMM solver |
| `code/python/` | Optuna TPE interface |
| `code/scripts/` | Environment setup and explicit experiment/analysis commands |
| `code/tests/` | Algebraic, mock and small synthetic solver tests |
| `code/configs/` | Historical experiment settings, seeds and input checksums |
| `code/docs/` | Methods, data sources, simulations and environment instructions |

Historical launchers, checkpoint-recovery services and publication-table/figure
generators are excluded. The compact implementation preserves the model
definitions and search protocol; it is not an archive of every historical timing
implementation. Review `git status --short` before committing local changes.

## Author and license

Software author: **Tomasz Górecki** ([Halmaris](https://github.com/Halmaris)).

Original code and documentation are released under the [MIT License](LICENSE).
Dependencies retain their own licenses; distribution of compiled code linked with
GPL-licensed R/Rcpp must also comply with the applicable GPL terms.
Software citation metadata are in [CITATION.cff](CITATION.cff).
