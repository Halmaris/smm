# Shared imports for the research scripts and tests.
suppressPackageStartupMessages({
  library(Rcpp, include.only = c('sourceCpp', 'cppFunction'))
  library(Matrix, include.only = c('Diagonal', 'sparseMatrix', 'Cholesky'))
  library(e1071, include.only = 'svm')
  library(reticulate, include.only = c('import', 'import_from_path'))
  library(tidyr)
  library(dplyr)
})
