# Reproducing the research environment

The two dependency files live in `code/`: `renv.lock` for R and
`requirements.txt` for Python. They cover the active research API, simulation
scripts and tests. Data, results, archived code and private manuscript tools are
not installed by these commands.

## Prerequisites

- macOS 14 or later on Apple Silicon (arm64); the solver links to Accelerate.
- The CRAN macOS arm64 distribution of R 4.5.1 and CPython 3.12.14 (arm64).
  Python must provide a shared library for reticulate; use `python3.12`, not the
  macOS system Python.
- Xcode Command Line Tools under `/Library/Developer/CommandLineTools`
  (`xcode-select --install` if absent). The script selects their compiler and SDK
  for its own process; it does not change the system's selected developer tools.
- GNU gettext headers: `libintl.h` in `/opt/homebrew/opt/gettext/include` or
  `/opt/R/arm64/include`. The script detects these locations automatically,
  preferring the first. If neither exists, install gettext separately, for
  example with `brew install gettext` on a machine already using Homebrew.
- The Fortran runtime libraries bundled with CRAN R, `libgfortran.5.dylib` and
  `libquadmath.0.dylib`, in `R.home('lib')`. The script uses these matching
  libraries, not a separate system Fortran installation.
- Internet access to CRAN and PyPI.

The interpreters and native build prerequisites are not installed by the script.
It checks the native build prerequisites before downloading R packages and stops
with a message if they are missing. The configuration above is the checked
macOS / Apple Silicon setup.

## Restore

From the project root, use the following commands. If your Python 3.12.14 binary
has a different name, replace only `python3.12` with its path.

```sh
cd code
python3.12 --version
Rscript --vanilla -e 'cat(R.version.string, "\n")'
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

The R script bootstraps renv 1.2.3 into `.renv/library` and restores the locked
versions from CRAN source packages. Only base R and the project library are
visible during restoration; copying from global libraries and the shared renv
cache is disabled. An empty `.renv/` therefore gives a fresh installation,
including the locked recommended packages. A repeated run retains matching
packages already in this project's library. The Python environment is `.venv/`.

Repeat the four `export` commands above in each new terminal before using the
code. There is no automatic startup hook: the exported `R_LIBS_USER` also works
with the existing `Rscript --vanilla` commands. Setting `RETICULATE_PYTHON`
selects the Python environment explicitly. Both environment directories are
ignored by Git; Google Drive synchronization is independent of `.gitignore`.

## Automatic build configuration

No compiler exports or edits to `~/.R/Makevars` are needed before the restoration
command. The script creates a temporary `Makevars` inside `.renv/` and removes it
on return, including when restoration raises an error. It sets:

- the Command Line Tools `clang`/`clang++` and their macOS SDK;
- `-O2`, C11 and C++17 defaults, with matching C++11/14/17/20 compiler entries;
- the detected gettext include directory alongside `/opt/R/arm64/include`;
- the two Fortran runtime libraries from the current R installation.

For this R process and its children, it ignores user/site Makevars and clears
`CPATH`, `C_INCLUDE_PATH`, `CPLUS_INCLUDE_PATH` and `LIBRARY_PATH`. This avoids the
conflicting Xcode/C++ headers and missing `libintl.h` found in the fresh-install
check. Package compilation uses one job. These settings expire with the script;
no shell profile, global R configuration, system compiler selection or installed
global package is changed.

All renv libraries, downloads and metadata stay in `.renv/`, which is ignored by
Git. The temporary compiler settings apply to package restoration, not later
model fitting or compilation of the research backend. Verification is described
in [VALIDATION.md](VALIDATION.md).

## Dependency scope

Direct R dependencies are Rcpp 1.1.2, Matrix 1.7-3, reticulate 1.43.0, dplyr 1.1.4, tidyr 1.3.1 and e1071
1.7-16. The lockfile also includes their required dependencies and renv: 35
packages in total; base R packages are supplied by R 4.5.1.

Optuna 4.5.0 is the sole direct third-party Python dependency. Its 10 transitive
dependencies are pinned as well, including NumPy and SQLAlchemy. No unrelated
packages from the author's Python installation are included.

The lock includes the tidyverse packages used by the analysis scripts. It preserves
the documented direct research dependencies, but are not a recovered complete
snapshot of every historical campaign environment.

## Checks without fitting

`scripts/check_environment.R` verifies the local locked R versions, imports the
research API with `compile = FALSE`, and imports Optuna through reticulate. It
neither fits a model nor creates an optimization study. `pip check` verifies
Python dependency compatibility.

The broader `scripts/run_checks.sh` includes `tests/test_models.R`, which fits
small synthetic models. Do not run it when only an environment check is wanted.

To refresh the locks intentionally, first update and validate the packages, then
snapshot only the research dependencies:

```r
renv::snapshot(project = '.', library = '.renv/library',
               packages = c('Rcpp', 'Matrix', 'reticulate', 'e1071', 'dplyr', 'tidyr', 'renv'),
               lockfile = 'renv.lock', prompt = FALSE)
```

For Python, resolve dependencies in a clean virtual environment, run the import
check, and use `python -m pip freeze` to record its runtime packages. Keep the
bootstrap renv version in `scripts/restore_environment.R` aligned with the lock.

References: [renv snapshot](https://rstudio.github.io/renv/reference/snapshot.html),
[renv restore](https://rstudio.github.io/renv/reference/restore.html),
[reticulate Python configuration](https://rstudio.github.io/reticulate/articles/versions.html).
