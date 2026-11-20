# Original training-partition preprocessing and helpers for the temporal OvR timing initializer.
.smm_scalar <- function(x, name, lower = 0, strict = FALSE, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      if (strict) x <= lower else x < lower) {
    stop(name, ' must be a finite scalar ', if (strict) '> ' else '>= ', lower,
         call. = FALSE)
  }
  if (integer && x != floor(x)) stop(name, ' must be an integer.', call. = FALSE)
  invisible(x)
}

.smm_pack <- function(X, p = NULL, q = NULL, vectorization = 'column',
                      single = FALSE) {
  if (!is.null(p)) .smm_scalar(p, 'p', 0, TRUE, TRUE)
  if (!is.null(q)) .smm_scalar(q, 'q', 0, TRUE, TRUE)
  if (is.list(X) && !is.data.frame(X)) {
    if (!length(X) || !all(vapply(X, is.matrix, logical(1)))) {
      stop('X must be a nonempty list of numeric matrices.', call. = FALSE)
    }
    dims <- dim(X[[1L]])
    if (!all(vapply(X, function(z) is.numeric(z) && identical(dim(z), dims),
                    logical(1)))) {
      stop('All matrices in X must be numeric and have the same dimensions.',
           call. = FALSE)
    }
    inferred <- dims
    Z <- t(matrix(vapply(X, as.vector, numeric(prod(dims))), nrow = prod(dims)))
  } else if (is.array(X) && length(dim(X)) == 3L && is.numeric(X)) {
    inferred <- dim(X)[1:2]
    Z <- t(matrix(X, nrow = prod(inferred)))
  } else if (is.matrix(X) && is.numeric(X)) {
    if (is.null(p) || is.null(q)) {
      stop('Supply p and q for a flattened matrix X.', call. = FALSE)
    }
    inferred <- c(p, q)
    if (single && identical(as.integer(dim(X)), as.integer(c(p, q)))) {
      Z <- matrix(as.vector(X), nrow = 1L)
    } else {
      if (ncol(X) != p * q) stop('ncol(X) must equal p * q.', call. = FALSE)
      Z <- X
      if (vectorization == 'row') {
        order <- as.vector(matrix(seq_len(p * q), p, q, byrow = TRUE))
        Z <- Z[, order, drop = FALSE]
      }
    }
  } else {
    stop('X must be a list of matrices, a p-by-q-by-n array, or an n-by-(p*q) matrix.',
         call. = FALSE)
  }
  if ((!is.null(p) && p != inferred[1L]) ||
      (!is.null(q) && q != inferred[2L])) {
    stop('The supplied p or q disagrees with the input matrices.', call. = FALSE)
  }
  if (!nrow(Z) || !ncol(Z) || !all(is.finite(Z))) {
    stop('X must be nonempty and contain finite values only.', call. = FALSE)
  }
  storage.mode(Z) <- 'double'
  list(X = Z, p = inferred[1L], q = inferred[2L])
}

.smm_labels <- function(y, n, binary = TRUE) {
  if (!is.atomic(y) || !is.null(dim(y)) || length(y) != n || anyNA(y) ||
      (is.numeric(y) && any(!is.finite(y)))) {
    stop('y must be a finite label vector with one label per observation.',
         call. = FALSE)
  }
  classes <- if (is.factor(y)) levels(droplevels(y)) else sort(unique(y))
  if (length(classes) < 2L || (binary && length(classes) != 2L)) {
    stop(if (binary) 'Exactly two observed classes are required.' else
           'At least two observed classes are required.', call. = FALSE)
  }
  list(classes = classes, factor = is.factor(y),
       signed = ifelse(match(y, classes) == 1L, -1, 1))
}

.smm_preprocess <- function(X, p, q, standardize) {
  center <- rep(0, p)
  scale <- rep(1, p)
  if (standardize) {
    values <- matrix(t(X), nrow = p)
    center <- rowMeans(values)
    centered <- sweep(values, 1L, center, '-')
    scale <- sqrt(rowSums(centered^2) / (ncol(values) - 1))
    scale[scale == 0] <- 1
  }
  list(X = sweep(sweep(X, 2L, rep(center, q), '-'), 2L, rep(scale, q), '/'),
       center = center, scale = scale)
}

smm_time_difference <- function(time) {
  if (!is.numeric(time) || !length(time) || any(!is.finite(time)) ||
      any(diff(time) <= 0)) {
    stop('time must be finite and strictly increasing.', call. = FALSE)
  }
  T <- length(time)
  D <- matrix(0, max(0L, T - 2L), T)
  if (T >= 3L) {
    h <- diff(time)
    for (j in seq_len(T - 2L)) {
      span <- h[j] + h[j + 1L]
      D[j, j + 0:2] <- sqrt(span / 2) * 2 / span *
        c(1 / h[j], -1 / h[j] - 1 / h[j + 1L], 1 / h[j + 1L])
    }
  }
  D
}

smm_fit_context <- function(X, y, standardize = TRUE) {
  packed <- .smm_pack(X)
  context <- new.env(parent = emptyenv())
  context$p <- packed$p; context$q <- packed$q
  context$labels <- .smm_labels(y, nrow(packed$X), FALSE)
  context$y <- y; context$standardize <- standardize
  context$preprocessing <- .smm_preprocess(packed$X, packed$p, packed$q, standardize)
  context$Z <- context$preprocessing$X
  context$center <- colMeans(context$Z)
  context$Zc <- sweep(context$Z, 2L, context$center, '-')
  context$signature <- digest::digest(list(context$Z, y, standardize), algo = 'sha256')
  context$cache <- new.env(parent = emptyenv())
  context$hits <- context$misses <- 0L
  class(context) <- 'smm_fit_context'
  context
}

.smm_cached <- function(context, key, build) {
  limit <- getOption('smm.solver_cache_limit', Inf)
  bounded <- startsWith(key, 'solver:') && is.finite(limit)
  if (bounded) {
    stopifnot(length(limit) == 1L, limit >= 1L)
    order <- context$solver_cache_order
    order <- order[order != key]
    if (!exists(key, context$cache, inherits = FALSE)) {
      while (length(order) >= limit) {
        rm(list = order[1L], envir = context$cache)
        order <- order[-1L]
      }
    }
    context$solver_cache_order <- c(order, key)
  }
  if (exists(key, context$cache, inherits = FALSE)) {
    context$hits <- context$hits + 1L
    return(get(key, context$cache, inherits = FALSE))
  }
  value <- build()
  assign(key, value, context$cache)
  context$misses <- context$misses + 1L
  value
}

.smm_banded_penalty <- function(q) {
  if (q < 3L) return(Matrix::Diagonal(q, 0))
  # D has three nonzero entries per row; D'D has bandwidth two.
  ii <- rep(seq_len(q - 2L), each = 3L)
  jj <- as.vector(t(outer(seq_len(q - 2L), 0:2, '+')))
  D <- Matrix::sparseMatrix(i = ii, j = jj, x = rep(c(1, -2, 1), q - 2L),
                            dims = c(q - 2L, q))
  Matrix::crossprod(D)
}
