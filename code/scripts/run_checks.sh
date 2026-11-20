#!/bin/sh
# Run validation from any working directory; Python defaults to the local venv.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export RETICULATE_PYTHON="${RETICULATE_PYTHON:-$root/.venv/bin/python}"
"$RETICULATE_PYTHON" "$root/tests/test_tpe.py"
Rscript --vanilla "$root/tests/test_math.R"
Rscript --vanilla "$root/tests/test_tuning.R"
Rscript --vanilla "$root/tests/test_models.R"
(cd "$root" && Rscript --vanilla tests/test_saved_analysis.R)
