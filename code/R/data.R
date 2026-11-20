# Design matrices and training context ----------------------------------

#' Pack a list of p-by-q observations into an n-by-(p*q) design matrix.
#' @param X List of equally sized numeric matrices, or a p-by-q-by-n array.
#' @return A list containing the vectorized observations and matrix dimensions.
.smm_pack <- function(X) {
  if (is.array(X) && length(dim(X)) == 3L)
    X <- lapply(seq_len(dim(X)[3]), function(i) return(matrix(X[, , i], nrow = dim(X)[1], ncol = dim(X)[2])))
  if (!is.list(X) || !length(X) || !all(vapply(X, is.matrix, TRUE)))
    stop('X must be a nonempty list of matrices or a three-dimensional array.')
  dims <- dim(X[[1]])
  stopifnot(all(vapply(X, function(x) return(is.numeric(x) && identical(dim(x), dims)), TRUE)))
  Z <- t(vapply(X, as.vector, numeric(prod(dims))))
  if (!all(is.finite(Z))) stop('Predictors must be finite and complete.')
  return(list(X = Z, p = dims[1], q = dims[2]))
}

#' Prepare training-only scaling and the structural operator for one partition.
#' @param X Matrix-valued observations, in the format accepted by .smm_pack.
#' @param y Class labels, one per observation.
#' @param structure Temporal second differences or spatial first differences.
#' @param standardize Whether to center and scale each row over observations and columns.
#' @param channels Number of disconnected colour planes stacked along rows (1 or 3).
#' @return A training context shared only within this training partition.
smm_context <- function(X, y, structure = 'temporal', standardize = TRUE, channels = 1L) {
  packed <- .smm_pack(X)
  stopifnot(length(y) == nrow(packed$X), !anyNA(y), structure %in% c('temporal', 'spatial'))
  classes <- if (is.factor(y)) levels(droplevels(y)) else sort(unique(y))
  stopifnot(length(classes) >= 2L, channels %in% c(1L, 3L), packed$p %% channels == 0L,
    structure == 'spatial' || channels == 1L)
  context <- new.env(parent = emptyenv())
  context$p <- packed$p; context$q <- packed$q; context$y <- y
  context$labels <- list(classes = classes, factor = is.factor(y))
  context$preprocessing <- .smm_preprocess(packed$X, packed$p, packed$q, standardize)
  context$Z <- context$preprocessing$X
  context$center <- colMeans(context$Z)
  context$Zc <- sweep(context$Z, 2L, context$center, '-')
  context$structure <- structure; context$row_blocks <- as.integer(channels)
  context$cache <- new.env(parent = emptyenv())
  if (structure == 'spatial') {
    laplacian <- function(n) return(if (n < 2L) matrix(0, n, n) else crossprod(diff(diag(n))))
    context$row_laplacian <- kronecker(diag(channels), laplacian(packed$p / channels))
    context$column_laplacian <- laplacian(packed$q)
    context$roughness_fn <- function(W) {
      height <- nrow(W) / channels
      return(sum(vapply(seq_len(channels), function(j) {
        plane <- W[seq_len(height) + (j - 1L) * height, , drop = FALSE]
        return(sum(diff(plane)^2) + sum(diff(t(plane))^2))
      }, 0)))
    }
  } else {
    D <- if (packed$q < 3L) matrix(0, 0, packed$q) else diff(diag(packed$q), differences = 2L)
    context$roughness_fn <- function(W) return(sum((W %*% t(D))^2))
  }
  return(context)
}

#' Reuse a matrix factorization within one training partition.
#' @param context Partition-specific training context.
#' @param key String identifying an operator and its parameters.
#' @param build Function that constructs the requested value.
.smm_cached <- function(context, key, build) {
  if (exists(key, context$cache, inherits = FALSE)) return(get(key, context$cache))
  # Keep at most two large factorizations of each type, as in the study.
  prefix <- if (startsWith(key, 'solver:')) 'solver:' else if (startsWith(key, 'inverse:'))
    'inverse:' else if (startsWith(key, 'native_inverse32:')) 'native_inverse32:' else NULL
  if (!is.null(prefix)) {
    order <- context$cache_order
    if (is.null(order)) order <- character()
    related <- order[startsWith(order, prefix)]
    if (length(related) >= 2L) {
      rm(list = related[1], envir = context$cache)
      order <- order[order != related[1]]
    }
    context$cache_order <- c(order, key)
  }
  value <- build(); assign(key, value, context$cache); return(value)
}

# Cross-validation and pruning ------------------------------------------

#' Create a fixed stratified cross-validation partition.
#' @param y Class labels, with at least folds observations in every class.
#' @param seed Integer RNG seed for the partition.
#' @param folds Number of folds; five in the research protocol.
stratified_folds <- function(y, seed, folds = 5L) {
  stopifnot(all(table(y) >= folds))
  set.seed(seed)
  answer <- integer(length(y))
  for (label in unique(y)) {
    ids <- which(y == label)
    answer[ids] <- sample(rep(seq_len(folds), length.out = length(ids)))
  }
  return(answer)
}

#' Bound pooled balanced accuracy after an incomplete set of validation folds.
#' @param truth Full training labels.
#' @param prediction Out-of-fold predictions; NA denotes an unprocessed observation.
#' @return The optimistic BA achievable if all remaining predictions are correct.
ba_upper_bound <- function(truth, prediction) {
  return(mean(vapply(unique(as.character(truth)), function(label) {
    ids <- as.character(truth) == label
    return(1 - sum(prediction[ids] != label, na.rm = TRUE) / sum(ids))
  }, 0)))
}
