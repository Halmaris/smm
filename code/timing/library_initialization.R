# Training-only LIBSVM starts; the existing ADMM solver still certifies every fit.
smm_library_initial <- function(context, parameters, initial = NULL) {
  stopifnot(requireNamespace('e1071', quietly = TRUE))
  p <- context$p; tt <- context$q; n <- length(context$y)
  gamma <- if (is.null(parameters$gamma)) 0 else parameters$gamma
  k <- length(context$labels$classes)
  Y <- outer(as.character(context$y), as.character(context$labels$classes),
    function(a, b) ifelse(a == b, 1, -1))
  costs <- apply(Y, 2L, function(z) parameters$Cbar /
    (2 * ifelse(z > 0, sum(z > 0), sum(z < 0))))
  h <- 1 / sqrt(max(1, mean(rowSums(context$Zc^2))))
  G <- crossprod(context$Z, Y * costs)
  zero <- vapply(seq_len(k), function(j)
    svd(matrix(G[, j], p, tt), nu = 0, nv = 0)$d[1L] <=
      parameters$tau * (1 - 1e-12), logical(1))
  zero_state <- list(context_id = context$signature,
    W = matrix(0, p * tt, k), b = numeric(k), S = matrix(0, p * tt, k),
    U = G, M = matrix(0, n, k), V = -costs / h, rho = 1)
  # For balanced OvR costs, alpha=C and nuclear dual G certify W=b=0
  # whenever ||G||op<=tau. Initializing its KKT state avoids a flat-loss stall.
  if (all(zero)) return(zero_state)
  key <- paste('libsvm_start', parameters$Cbar, gamma, sep = ':')
  start <- .smm_cached(context, key, function() {
    transform <- .smm_cached(context, paste('libsvm_transform', gamma, sep = ':'), function() {
      if (gamma == 0) return(list(X = context$Z, factor = NULL))
      Q <- Matrix::Diagonal(tt) + 2 * gamma * .smm_banded_penalty(tt)
      factor <- Matrix::Cholesky(Q, perm = FALSE, LDL = FALSE, super = FALSE)
      blocks <- matrix(aperm(array(t(context$Z), c(p, tt, n)), c(2L, 1L, 3L)), nrow = tt)
      solved <- as.matrix(Matrix::solve(factor, blocks, system = 'L'))
      X <- t(matrix(aperm(array(solved, c(tt, p, n)), c(2L, 1L, 3L)), nrow = p * tt))
      list(X = X, factor = factor)
    })
    k <- length(context$labels$classes)
    W <- matrix(0, p * tt, k); b <- numeric(k); alpha <- matrix(0, n, k)
    for (j in seq_len(k)) {
      signed <- factor(ifelse(context$y == context$labels$classes[j], 1, -1), levels = c(-1, 1))
      fit <- e1071::svm(transform$X, signed, kernel = 'linear',
        cost = parameters$Cbar / n, class.weights = n / (2 * table(signed)),
        scale = FALSE, tolerance = 1e-8, shrinking = FALSE, cachesize = 256)
      orientation <- if (fit$levels[fit$labels[1L]] == '1') 1 else -1
      W[, j] <- orientation * as.vector(crossprod(fit$SV, as.vector(fit$coefs)))
      b[j] <- -orientation * fit$rho
      alpha[fit$index, j] <- abs(as.vector(fit$coefs))
    }
    if (gamma > 0) {
      blocks <- matrix(aperm(array(W, c(p, tt, k)), c(2L, 1L, 3L)), nrow = tt)
      solved <- as.matrix(Matrix::solve(transform$factor, blocks, system = 'Lt'))
      W <- matrix(aperm(array(solved, c(tt, p, k)), c(2L, 1L, 3L)), nrow = p * tt)
    }
    Y <- outer(as.character(context$y), as.character(context$labels$classes),
      function(a, b) ifelse(a == b, 1, -1))
    h <- 1 / sqrt(max(1, mean(rowSums(context$Zc^2))))
    list(context_id = context$signature, W = W, b = b, S = W,
      U = matrix(0, p * tt, k), M = h * Y * sweep(context$Z %*% W, 2L, b, '+'),
      V = -alpha / h, rho = 1)
  })
  if (any(zero)) {
    for (name in c('W', 'S', 'U', 'M', 'V')) start[[name]][, zero] <- zero_state[[name]][, zero]
    start$b[zero] <- 0
  }
  if (is.null(initial)) return(start)
  objective <- function(state) {
    Y <- outer(as.character(context$y), as.character(context$labels$classes),
      function(a, b) ifelse(a == b, 1, -1))
    costs <- apply(Y, 2L, function(z) parameters$Cbar /
      (2 * ifelse(z > 0, sum(z > 0), sum(z < 0))))
    margin <- Y * sweep(context$Z %*% state$W, 2L, state$b, '+')
    colSums(state$W^2) / 2 + colSums(costs * pmax(0, 1 - margin)) +
      vapply(seq_len(ncol(state$W)), function(j) {
        w <- matrix(state$W[, j], p, tt)
        parameters$tau * sum(svd(w, nu = 0, nv = 0)$d) +
          gamma * sum((w %*% t(smm_time_difference(seq_len(tt))))^2)
      }, numeric(1))
  }
  use_initial <- !zero & objective(initial) < objective(start) - 1e-10
  if (any(use_initial)) {
    initial$U <- initial$U * initial$rho
    initial$V <- initial$V * initial$rho
    for (name in c('W', 'S', 'U', 'M', 'V'))
      start[[name]][, use_initial] <- initial[[name]][, use_initial]
    start$b[use_initial] <- initial$b[use_initial]
  }
  start
}
