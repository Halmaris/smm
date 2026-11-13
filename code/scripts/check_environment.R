# Check locked R versions and the R/Python bridge; do not compile or fit models.
# Run from code/ after exporting R_LIBS_USER and RETICULATE_PYTHON as documented.
library_path <- normalizePath('.renv/library', mustWork = TRUE)
lock <- renv::lockfile_read('renv.lock')
stopifnot(as.character(getRversion()) == lock$R$Version)
for (name in names(lock$Packages)) {
  description <- packageDescription(name, lib.loc = library_path)
  stopifnot(description$Version == lock$Packages[[name]]$Version)
  stopifnot(normalizePath(find.package(name)) ==
              normalizePath(file.path(library_path, name)))
}
for (name in c('Rcpp', 'Matrix', 'reticulate', 'e1071'))
  stopifnot(requireNamespace(name, quietly = TRUE))

source('R/load.R')
api <- load_smm('.', compile = FALSE)
stopifnot(is.function(api$fit_direct), is.function(api$tune_direct))
bridge <- reticulate::import_from_path('tpe', normalizePath('python'))
optuna <- reticulate::import('optuna')
stopifnot(optuna$`__version__` == '4.5.0')
sys <- reticulate::import('sys')
stopifnot(startsWith(sys$version, '3.12.14 '))
cat('PASS: locked R packages, research API and Optuna import; no model fits.\n')
