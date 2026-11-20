#' Load the research functions into a separate environment.
#' @param root Repository directory, relative or absolute.
#' @param compile Whether to compile the single native FP32 backend.
#' @return An environment containing the research API; loading does not fit models.
load_smm <- function(root = '.', compile = TRUE) {
  root <- normalizePath(root, mustWork = TRUE)
  api <- new.env(parent = globalenv())
  api$.repo <- root
  for (file in c('data.R', 'numerics.R', 'models.R', 'baselines.R', 'generators.R',
                 'tuning.R', 'experiments.R', 'saved_analysis.R', 'analysis.R'))
    sys.source(file.path(root, 'R', file), envir = api)
  if (compile) api$compile_solver(file.path(root, 'src', 'admm_fp32.cpp'))
  api
}
