# Multiclass simulation design ------------------------------------------

#' Enumerate the fixed class/rank/roughness/training-size factorial design.
#' @param cfg Generator settings from paper_design.
mc_plan <- function(cfg) {
  df <- expand.grid(K = cfg$classes, rank = cfg$ranks,
    roughness = cfg$roughness, train_per_class = cfg$train_per_class,
    stringsAsFactors = FALSE)
  df$scenario <- sprintf('s%02d', seq_len(nrow(df)))
  return(df)
}

#' Assign a deterministic RNG stream to one scenario, replication and generator step.
#' @param cfg Generator settings from paper_design.
#' @param scenario One-based numeric scenario index for the deterministic RNG streams.
#' @param replication One-based replication number.
#' @param stream Independent generator stream index.
mc_seed <- function(cfg, scenario, replication, stream) {
  # Disjoint deterministic streams; no dependence on execution order.
  return(as.integer(cfg$seed + scenario * 100000L + replication * 100L + stream))
}

#' Construct temporal second differences on the unit-spaced grid.
#' @param T Temporal grid length.
mc_D <- function(T) return(diff(diag(T), differences = 2L))

#' Draw orthonormal columns using the current RNG stream.
#' @param n Number of rows or observations.
#' @param k Number of orthonormal columns.
mc_orthogonal <- function(n, k) return(qr.Q(qr(matrix(rnorm(n * k), n, k)))[, seq_len(k), drop = FALSE])

#' Construct equal-Gram class coefficients with controlled rank and roughness.
#' @param cfg Generator settings from paper_design.
#' @param job One row of mc_plan.
#' @param replication One-based replication number.
#' @param scenario One-based numeric scenario index for the deterministic RNG streams.
mc_coefficients <- function(cfg, job, replication, scenario) {
  set.seed(mc_seed(cfg, scenario, replication, 1L))
  r <- job$rank; K <- job$K
  U <- mc_orthogonal(cfg$p, r)
  eig <- eigen(crossprod(mc_D(cfg$T)), symmetric = TRUE)
  ids <- if (job$roughness == 'smooth') order(eig$values)[seq_len(r)] else
    order(eig$values, decreasing = TRUE)[seq_len(r)]
  V <- eig$vectors[, ids, drop = FALSE]
  # A regular simplex in r^2 dimensions: same score covariance in every regime.
  H <- diag(K) - matrix(1 / K, K, K)
  E <- eigen(H, symmetric = TRUE)$vectors[, seq_len(K - 1L), drop = FALSE]
  A <- sqrt(K / (K - 1L)) * mc_orthogonal(r * r, K - 1L) %*% t(E)
  W <- vapply(seq_len(K), function(k)
    return(as.vector(U %*% matrix(A[, k], r, r) %*% t(V))), numeric(cfg$p * cfg$T))
  target <- (K / (K - 1L)) * H
  stopifnot(max(abs(crossprod(W) - target)) < 1e-10,
    max(abs(rowSums(W))) < 1e-10,
    all(vapply(seq_len(K), function(k)
      return(qr(matrix(W[, k], cfg$p, cfg$T), tol = 1e-9)$rank), 0L) == r))
  return(list(W = W, gram = crossprod(W),
    roughness = vapply(seq_len(K), function(k)
      return(sum((matrix(W[, k], cfg$p, cfg$T) %*% t(mc_D(cfg$T)))^2)), 0)))
}

#' Sample observations with fixed class counts by conditional rejection.
#' @param cfg Generator settings from paper_design.
#' @param job One row of mc_plan.
#' @param truth Generating coefficient matrices and their geometry.
#' @param per_class Required number of observations per class.
#' @param seed Integer RNG seed.
#' @param prefix Prefix for generated sample identifiers.
mc_sample <- function(cfg, job, truth, per_class, seed, prefix) {
  set.seed(seed)
  K <- job$K; d <- cfg$p * cfg$T
  accepted <- vector('list', K); counts <- integer(K); draws <- 0L
  while (any(counts < per_class)) {
    Z <- matrix(rnorm(512L * d), 512L, d)
    scores <- Z %*% truth$W + matrix(rnorm(512L * K, sd = cfg$noise_sd), 512L, K)
    label <- max.col(scores, ties.method = 'first')
    for (k in seq_len(K)) {
      ids <- head(which(label == k), per_class - counts[k])
      if (length(ids)) {
        accepted[[k]] <- rbind(accepted[[k]], Z[ids, , drop = FALSE])
        counts[k] <- counts[k] + length(ids)
      }
    }
    draws <- draws + 512L
    if (draws > 10000000L) stop('Conditional sampler exceeded its draw cap.')
  }
  Z <- do.call(rbind, accepted); y <- rep(seq_len(K), each = per_class)
  order <- sample.int(length(y)); Z <- Z[order, , drop = FALSE]; y <- y[order]
  return(list(X = lapply(seq_len(nrow(Z)), function(i) return(matrix(Z[i, ], cfg$p, cfg$T))),
    y = factor(y, levels = seq_len(K)), sample_id = paste0(prefix, seq_along(y)),
    candidate_draws = draws, sampling = 'fixed class counts by conditional rejection'))
}

#' Generate paired development/test samples and the fixed stratified folds.
#' @param cfg Generator settings from paper_design.
#' @param job One row of mc_plan.
#' @param replication One-based replication number.
#' @param scenario One-based numeric scenario index for the deterministic RNG streams.
mc_inputs <- function(cfg, job, replication, scenario) {
  truth <- mc_coefficients(cfg, job, replication, scenario)
  train <- mc_sample(cfg, job, truth, job$train_per_class,
    mc_seed(cfg, scenario, replication, 2L), 'train_')
  test <- mc_sample(cfg, job, truth, cfg$test_per_class,
    mc_seed(cfg, scenario, replication, 3L), 'test_')
  set.seed(mc_seed(cfg, scenario, replication, 4L))
  train$folds <- integer(length(train$y))
  for (k in levels(train$y)) {
    ids <- which(train$y == k)
    train$folds[ids] <- sample(rep(seq_len(cfg$folds), length.out = length(ids)))
  }
  stopifnot(all(table(train$y, train$folds) > 0),
    !any(train$sample_id %in% test$sample_id))
  return(list(train = train, test = test, truth = truth))
}

# Population and sample summaries ---------------------------------------

#' Compute pooled balanced accuracy and macro-F1 on all supplied predictions.
#' @param y Class labels; internal binary operations use -1/+1 and direct operations use 1..K.
#' @param pred Predicted class labels.
mc_metrics <- function(y, pred) {
  levels <- levels(y)
  stopifnot(length(pred) == length(y), !anyNA(pred), all(as.character(pred) %in% levels))
  cm <- table(factor(y, levels), factor(pred, levels)); tp <- diag(cm)
  f1 <- 2 * tp / (rowSums(cm) + colSums(cm))
  return(list(BA = mean(tp / rowSums(cm)), macro_F1 = mean(f1),
    recall = tp / rowSums(cm), confusion = cm))
}

#' Summarize coefficient effective ranks and scale-invariant temporal roughness.
#' @param fit Fitted direct model with matrix dimensions and coefficient columns.
mc_structure <- function(fit) {
  D <- mc_D(fit$q)
  one <- function(w) {
    W <- matrix(w, fit$p, fit$q); s <- svd(W, nu = 0, nv = 0)$d
    energy <- sum(s^2)
    return(c(effective_rank = if (energy > 0) sum(s)^2 / energy else 0,
      relative_roughness = if (energy > 0) sum((W %*% t(D))^2) / energy else 0))
  }
  pairs <- combn(ncol(fit$W), 2L)
  return(list(classes = apply(fit$W, 2L, one), contrasts = vapply(seq_len(ncol(pairs)),
    function(j) return(one(fit$W[, pairs[1L, j]] - fit$W[, pairs[2L, j]])), numeric(2)),
    coordinates = 'training-standardized; scale invariant summaries, not coefficient recovery'))
}
# Binary simulation design ----------------------------------------------

#' Generate the binary matrix-classification model used in the original study.
#' @param n_train Number of binary development observations.
#' @param n_test Number of binary test observations.
#' @param p Number of matrix rows.
#' @param T Temporal grid length.
#' @param rank Generating coefficient rank.
#' @param roughness Target relative squared second-difference norm.
#' @param correlation Separable AR(1) correlation parameter.
#' @param label_flip Independent label-flip probability.
#' @param signal_sd Standard deviation of the latent linear predictor.
#' @param seed Integer RNG seed.
simulate_smm <- function(n_train, n_test = 2000L, p = 8L, T = 32L,
                         rank = 2L, roughness = 0.5, correlation = 0,
                         label_flip = 0, signal_sd = 2, seed = 20261010L) {
  values <- c(n_train, n_test, p, T, rank)
  stopifnot(all(is.finite(values)), all(values == as.integer(values)),
    n_train >= 2L, n_test >= 1L, p >= 1L, T >= 3L, rank >= 1L,
    rank <= p, 2L * rank <= T, abs(correlation) < 1,
    is.finite(roughness), roughness >= 0,
    label_flip >= 0, label_flip < 0.5, signal_sd > 0,
    length(seed) == 1L, is.finite(seed))
  D <- diff(diag(T), differences = 2L)
  spectrum <- eigen(crossprod(D), symmetric = TRUE)
  order <- order(spectrum$values)
  eigenvalues <- pmax(spectrum$values[order], 0)
  Q <- spectrum$vectors[, order, drop = FALSE]
  low <- seq_len(rank); high <- (T - rank + 1L):T
  lower <- mean(eigenvalues[low]); upper <- mean(eigenvalues[high])
  if (roughness < lower - 1e-10 || roughness > upper + 1e-10)
    stop(sprintf('Roughness must lie in [%.6g, %.6g] for this rank.', lower, upper))
  mixing <- min(1, max(0, (roughness - lower) / (upper - lower)))
  V <- sqrt(1 - mixing) * Q[, low, drop = FALSE] +
    sqrt(mixing) * Q[, high, drop = FALSE]
  set.seed(seed)
  U <- qr.Q(qr(matrix(rnorm(p * rank), p, rank)))
  shape <- U %*% t(V) / sqrt(rank)
  row_covariance <- correlation^abs(outer(seq_len(p), seq_len(p), '-'))
  time_covariance <- correlation^abs(outer(seq_len(T), seq_len(T), '-'))
  signal_variance <- sum(shape * (row_covariance %*% shape %*% time_covariance))
  W <- signal_sd * shape / sqrt(signal_variance)
  row_root <- t(chol(row_covariance)); time_root <- t(chol(time_covariance))
  # Reset after constructing W so paired scenarios use the same input matrices.
  set.seed(seed + 1L)
  X <- lapply(seq_len(n_train + n_test), function(i)
    return(row_root %*% matrix(rnorm(p * T), p, T) %*% t(time_root)))
  eta <- vapply(X, function(x) return(sum(x * W)), numeric(1))
  probability <- plogis(eta)
  y <- ifelse(runif(length(X)) < probability, 1L, -1L)
  flip <- runif(length(X)) < label_flip
  y[flip] <- -y[flip]
  effective_probability <- label_flip + (1 - 2 * label_flip) * probability
  indices <- list(train = seq_len(n_train), test = n_train + seq_len(n_test))
  samples <- lapply(indices, function(i) return(list(X = X[i], y = factor(y[i], levels = c(-1, 1)),
    probability = effective_probability[i], linear_predictor = eta[i])))
  return(c(samples, list(W = W, coefficient_shape = shape, D = D,
    row_covariance = row_covariance, time_covariance = time_covariance,
    settings = list(n_train = n_train, n_test = n_test, p = p, T = T, rank = rank,
      roughness = roughness, correlation = correlation, label_flip = label_flip,
      signal_sd = signal_sd, seed = seed, mixing = mixing))))
}
