# Validation of the research refactoring

Validation performed on 10 November 2026, macOS / Apple Silicon:

| Component | Tested version |
|---|---|
| R | 4.5.1 |
| Rcpp | 1.1.2 |
| Matrix | 1.7-3 |
| reticulate | 1.43.0 |
| e1071 | 1.7-16 |
| Python | 3.12.14 |
| Optuna | 4.5.0 |

## Checks

- Temporal and disconnected RGB penalty inverses agree with independently
  assembled FP64 matrices; both primal and Woodbury paths are tested.
- A constant image in each channel has zero spatial roughness even when channel
  constants differ. No cross-channel graph edges are created.
- Nuclear soft thresholding, simplex projection, canonical intercepts and the
  optimistic pooled-BA bound pass explicit small algebraic examples.
- Small synthetic direct fits cover all four penalty variants and analytically
  certified zero solutions. Temporal and RGB fits satisfy the unchanged residual
  and gap criteria, including a small test at the expanded `gamma = 1e5` endpoint. All named method APIs fit/predict on a small synthetic case.
- Native FP32 inverse and factor actions agree with their FP64 counterparts to
  `1e-5` on the tested small RGB problem.
- Four small direct-model fits match the frozen multiclass R wrapper exactly
  (maximum coefficient difference zero, identical update counts and objectives).
  This comparison uses the shared numerical kernel and does not establish
  bitwise equality on every dataset or platform.
- Reduced-size synthetic generator outputs, class labels and fold assignments
  are identical to the frozen generator.
- 180 artificial-score TPE proposals across M00/M10/M01/M11 reproduce the frozen
  implementation's parameters, pruning decisions and selected choices.
- The R/Python interface completes all five folds with deterministic mock models.
- Mock experiment tests cover CV-only portfolio choice, boundary reuse, the
  otherwise missing original M01 refit, and averaging two seeds within each of
  ten replications before inference.
- Plans enumerate 240/40 tunings and 12,000/2,660 proposals for core/boundary stages.

`tests/test_math.R` and `tests/test_tpe.py` check algebra and workflow without
research fitting. `tests/test_tuning.R` uses mock classifiers.
`tests/test_models.R` performs only small synthetic validation fits. The optional
historical comparisons above were run against the local frozen sources, which
are deliberately not bundled as a second implementation.

No new research campaign, empirical-data refit or timing study was run. Full-scale
performance, long-running numerical edge cases and exact historical timing
reproduction are outside the validation performed here. Nonconvergence still
raises an error; tests do not justify relaxing solver criteria.

## Project layout

The code-only package is `code/` within the main project. Historical scientific
inputs/results are stored separately under `data/` and `results/`; frozen source
snapshots are in `results/provenance/`. Operational queues, checkpoint caches,
duplicate environments and manuscript previews are not part of the active code.

## Relocation check, 12 November 2026

The algebra tests, artificial-score TPE tests and R/Python mock tuning test passed
again during GitHub packaging. No synthetic solver fits or research refits were
run in this check. The native solver source was unchanged. Dataset-release links
and public-file selection were checked separately from numerical validation.

## Environment lock verification, 13 November 2026

The active R and Python sources were scanned for dependencies. `renv.lock`
records R 4.5.1 and 17 required packages, including renv and transitive dependencies.
`requirements.txt` records Optuna 4.5.0 and its 10 transitive dependencies for
CPython 3.12.14 on macOS arm64.

A clean temporary Python virtual environment was installed from PyPI; `pip check`
and imports passed. The initial R check copied matching installed packages into
a local library and verified the reticulate/Optuna bridge.

A subsequent R check downloaded and built all 17 locked packages from CRAN
sources in an empty temporary library, without the shared renv cache or copying
installed R packages. It identified conflicting host Xcode/C++ headers and a
missing gettext include path. `scripts/restore_environment.R` now selects the
checked Command Line Tools/SDK, gettext headers and R-bundled Fortran runtime
libraries automatically through a temporary Makevars file.

The updated script was tested with the documented `Rscript --vanilla` command,
without external compiler overrides, in a path containing spaces. Fresh
restoration, exact package versions, all 17 local namespace imports and API
loading with `compile = FALSE` passed. Repeating restoration retained the
synchronized library. Temporary Makevars removal was checked both after success
and an injected installation error. Global/user configuration files were
unchanged, and the disposable validation environment was removed.

No model fits, tuning runs or SMM backend compilation were performed in these
environment checks; compilation concerned R package dependencies only.

These lockfiles define the environment checked on this date; they do not claim
to reconstruct every transitive package version used in historical campaigns.

## Saved-result tools

`Rscript --vanilla tests/test_saved_analysis.R` checks classification metrics,
paired replication alignment and incomplete blocks without fitting models.
The added CSV analyser was also checked against the saved manuscript inputs:
multiclass effects, intervals and Holm p-values, and temporal BA rank comparisons.
The checksum utility was tested with matching, changed and missing files.
The new multiclass runner was syntax-checked; no new model fits were run.
