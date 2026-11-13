# Research API

Run the examples from `code/`. Installation and environment variables are in the
[root README](../README.md); [ENVIRONMENT.md](docs/ENVIRONMENT.md) explains the
lockfiles and local libraries. The API uses one native FP32 ADMM backend with FP64
numerical checks; functions document their purpose and parameters in English.

## Classification

```r
source('R/load.R')
smm <- load_smm('.')

# X is a list of equally sized matrices; y contains the training labels.
context <- smm$smm_context(X, y, structure = 'temporal')
model <- smm$fit_direct(context, Cbar = 10, tau = .1, gamma = 1)
prediction <- smm$predict_model(model, X_test)
```

A `p × q × n` array is also accepted. Scaling is learned from the supplied training
partition. Never reuse a context across CV folds. RGB stacks three image planes
vertically and uses `structure = 'spatial', channels = 3`; adjacency does not cross
channel boundaries. See [input formats and sources](docs/DATA_SOURCES.md).

`fit_ovr` uses balanced binary losses and free intercepts. `fit_direct` uses the
coupled multiclass hinge, uniform observation costs and no intercept. `fit_method`
and `predict_method` expose all comparators; [METHODS.md](docs/METHODS.md) gives
the objectives. Loading the API compiles in R's temporary directory and fits
nothing until a fitting function is called.

## Tuning and simulations

```r
# train contains X, y and fixed stratified fold identifiers 1..5.
choice <- smm$tune_direct(train, method = 'M11', seed = 91001L, budget = 100L)
```

`tune_direct` uses TPE, fixed anchors, pooled five-fold balanced accuracy and
optimistic-bound pruning. M00/M10/M01/M11 allow neither, nuclear, adjacency or both
penalties. `select_grid` is a separate helper for supplied grids, not the paper's
TPE protocol.

`paper_design`, `mc_plan` and `mc_inputs` define the multiclass generator;
`run_ablation` runs one supplied four-variant block. `simulate_smm` implements
the binary generator. See function documentation for seeds and parameters.
The optional follow-up designs and explicit execution commands are in
[SIMULATIONS.md](docs/SIMULATIONS.md). Inspect the plan without fits:

```r
source('R/load.R')
smm <- load_smm('.', compile = FALSE)
smm$additional_plan('core')
smm$additional_plan('boundary')
```

## Layout

- `R/data.R`, `R/numerics.R`: inputs, scaling, adjacency, proximal operations.
- `R/models.R`, `R/baselines.R`, `src/`: classifiers and the native solver.
- `R/generators.R`: simulated inputs and numerical metrics.
- `R/tuning.R`, `python/`: TPE and pruning.
- `R/experiments.R`, `R/analysis.R`: experiment designs and paired estimates.
- `scripts/`: environment setup, experiment and saved-result analysis commands.
- `tests/`: algebra, mock workflows and small synthetic solver checks.

No publication tables, figures, background controllers, automatic retries or
checkpoint-recovery machinery are included. Tests and their limits are described
in [VALIDATION.md](docs/VALIDATION.md). This code preserves the model definitions
and search protocol, not every historical timing implementation.
