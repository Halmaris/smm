# Run from code/ with Rscript --vanilla scripts/restore_environment.R.
# This script installs packages; it never compiles the SMM backend or fits models.

#' Configure the checked macOS build tools for this R process and its children.
#' @param directory Project-local directory for the temporary Makevars file.
#' @return The Makevars path, to be removed by the caller after restoration.
configure_macos_build <- function(directory) {
  tools <- '/Library/Developer/CommandLineTools'
  cc <- file.path(tools, 'usr/bin/clang')
  cxx <- file.path(tools, 'usr/bin/clang++')
  sdk <- file.path(tools, 'SDKs/MacOSX.sdk')
  gettext <- c('/opt/homebrew/opt/gettext/include', '/opt/R/arm64/include')
  gettext <- gettext[file.exists(file.path(gettext, 'libintl.h'))]
  runtime <- file.path(R.home('lib'), c('libgfortran.5.dylib', 'libquadmath.0.dylib'))
  if (!all(file.exists(c(cc, cxx, sdk, runtime))) || !length(gettext))
    stop('Missing Command Line Tools, gettext headers or CRAN R runtime libraries; ',
         'see docs/ENVIRONMENT.md, Prerequisites.')

  makevars <- tempfile('Makevars-restore-', tmpdir = directory)
  include <- unique(c('/opt/R/arm64/include', gettext[1L]))
  writeLines(c(
    sprintf('CC = %s -arch arm64 -std=gnu11', shQuote(cc)),
    sprintf('CXX = %s -arch arm64 -std=gnu++17', shQuote(cxx)),
    paste0('CXX', c('11', '14', '17', '20'), ' = ', shQuote(cxx), ' -arch arm64'),
    paste('CPPFLAGS =', paste0('-I', shQuote(include), collapse = ' ')),
    'CFLAGS = -O2', 'CXXFLAGS = -O2',
    paste0('CXX', c('11', '14', '17', '20'), 'FLAGS = -O2'),
    paste('FLIBS =', paste(shQuote(runtime), collapse = ' '))), makevars)

  # Ignore conflicting shell include paths and user/site Makevars for this run only.
  Sys.unsetenv(c('CPATH', 'LIBRARY_PATH', 'CPLUS_INCLUDE_PATH', 'C_INCLUDE_PATH'))
  Sys.setenv(R_MAKEVARS_USER = makevars, R_MAKEVARS_SITE = '/dev/null', SDKROOT = sdk)
  message('Build tools: ', tools, '; gettext headers: ', gettext[1L])
  makevars
}

#' Restore locked packages into an isolated project library from CRAN sources.
#' @param project Directory containing renv.lock (normally code/).
#' @return Invisibly, the restored library path.
restore_environment <- function(project = '.') {
  project <- normalizePath(project, mustWork = TRUE)
  if (Sys.info()[['sysname']] != 'Darwin' || R.version$arch != 'aarch64')
    stop('This environment targets CRAN R on macOS / Apple Silicon.')
  stopifnot(file.exists(file.path(project, 'renv.lock')))
  state <- file.path(project, '.renv')
  library_path <- file.path(state, 'library')
  dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
  makevars <- configure_macos_build(state)
  on.exit(unlink(makevars), add = TRUE)

  Sys.setenv(RENV_PATHS_ROOT = state, RENV_PATHS_RENV = state,
             RENV_PATHS_LIBRARY = library_path,
             RENV_PATHS_SANDBOX = file.path(state, 'sandbox'),
             RENV_CONFIG_CACHE_ENABLED = 'FALSE',
             RENV_CONFIG_INSTALL_SHORTCUTS = 'FALSE',
             RENV_CONFIG_USER_LIBRARY = 'FALSE',
             RENV_CONFIG_SANDBOX_ENABLED = 'TRUE',
             RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE',
             RENV_CONFIG_PAK_ENABLED = 'FALSE',
             RENV_SANDBOX_LOCKING_ENABLED = 'FALSE',
             OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1',
             MKL_NUM_THREADS = '1', VECLIB_MAXIMUM_THREADS = '1')
  options(repos = c(CRAN = 'https://cloud.r-project.org'), pkgType = 'source',
          Ncpus = 1L, renv.config.install.jobs = 1L, renv.sandbox.priority = 'base')

  # Bootstrap locally; renv needs no non-base R dependencies.
  bootstrap_version <- '1.2.3'
  description <- file.path(library_path, 'renv', 'DESCRIPTION')
  if (!file.exists(description) || read.dcf(description)[1L, 'Version'] != bootstrap_version)
    install.packages(paste0('https://cloud.r-project.org/src/contrib/Archive/renv/renv_',
                            bootstrap_version, '.tar.gz'),
                     repos = NULL, type = 'source', lib = library_path)
  library('renv', lib.loc = library_path)
  lock <- renv::lockfile_read(file.path(project, 'renv.lock'))
  stopifnot(as.character(getRversion()) == lock$R$Version,
            lock$Packages$renv$Version == bootstrap_version)

  # Expose only base R and this project: also rebuild recommended packages locally.
  renv::load(project, quiet = TRUE)
  .libPaths(c(library_path, .libPaths()))
  renv::restore(project = project, library = library_path,
                lockfile = file.path(project, 'renv.lock'), prompt = FALSE)
  for (name in names(lock$Packages)) {
    stopifnot(packageDescription(name, lib.loc = library_path)$Version ==
                lock$Packages[[name]]$Version,
              normalizePath(find.package(name)) ==
                normalizePath(file.path(library_path, name)))
  }
  cat('R packages restored. Set R_LIBS_USER to code/.renv/library before running R.\n')
  invisible(library_path)
}

restore_environment()
