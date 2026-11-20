# Feature transformations and covariance estimates ----------------------

#' Construct orthonormal discrete cosine basis columns.
#' @param T Temporal grid length.
#' @param r Retained basis dimension.
.benchmark_basis <- function(T, r) {
  if (r > T || r < 1) stop('Invalid temporal basis dimension.')
  basis <- outer((seq_len(T) - 0.5) / T, 0:(r - 1L), function(t, k) cos(pi * t * k))
  basis[, 1L] <- 1 / sqrt(T)
  if (r > 1L) basis[, 2:r] <- basis[, 2:r, drop = FALSE] * sqrt(2 / T)
  return(basis)
}

#' Apply the mean/trend or truncated-cosine temporal representation.
#' @param X Matrix-valued observations or the vectorized design matrix for this operation.
#' @param method Comparator identifier.
#' @param parameters Named comparator hyperparameters.
.benchmark_transform <- function(X, method, parameters) {
  T <- ncol(X[[1L]])
  if (method == 'summary_svm') {
    t <- seq_len(T) - mean(seq_len(T))
    return(lapply(X, function(x) cbind(rowMeans(x), as.vector(x %*% t / sum(t^2)))))
  }
  if (method == 'functional_svm') {
    basis <- .benchmark_basis(T, parameters$components)
    return(lapply(X, function(x) x %*% basis))
  }
  return(X)
}

#' Shrink a covariance estimate toward a scaled identity.
#' @param S Empirical covariance matrix.
#' @param lambda Shrinkage proportion in (0, 1].
.benchmark_shrink <- function(S, lambda) {
  target <- mean(diag(S))
  if (target <= .Machine$double.eps) target <- 1
  return((1 - lambda) * S + lambda * target * diag(nrow(S)))
}

#' Accumulate covariance terms in observation order (R double precision).
#' @param E List of within-class residual matrices.
#' @param operation Function producing one covariance contribution.
covariance_sum <- function(E, operation) {
  value <- operation(E[[1L]])
  if (length(E) > 1L) for (i in 2:length(E)) value <- value + operation(E[[i]])
  return(value)
}

# Comparator fitting and prediction -------------------------------------

#' Fit a kernel, covariance or temporal-representation comparator on training data only.
#' @param X Matrix-valued observations or the vectorized design matrix for this operation.
#' @param y Class labels; internal binary operations use -1/+1 and direct operations use 1..K.
#' @param method Comparator identifier.
#' @param parameters Named comparator hyperparameters.
fit_comparator <- function(X, y, method, parameters) {
  context <- smm_context(X, y)
  p <- context$p; q <- context$q; preprocessing <- context$preprocessing
  labels <- context$labels
  Z <- preprocessing$X
  matrices <- if (method %in% c('summary_svm', 'functional_svm'))
    lapply(seq_len(nrow(Z)), function(i) matrix(Z[i, ], p, q)) else NULL
  cached <- function(key, build) return(if (is.null(context)) build() else .smm_cached(context, key, build))
  classes <- labels$classes
  group <- match(y, classes)
  n <- length(y); k <- length(classes)
  model <- list(method = method, parameters = parameters, classes = classes,
                p = p, q = q, center = preprocessing$center, scale = preprocessing$scale,
                factor_labels = labels$factor, converged = TRUE)
  if (method %in% c('summary_svm', 'functional_svm')) {
    inputs <- .benchmark_transform(matrices, method, parameters)
    model$fit <- fit_ovr(smm_context(inputs, y, standardize = FALSE),
      Cbar = parameters$Cbar, tau = 0, gamma = 0, relaxation = 1)
    model$converged <- model$fit$converged
  } else if (method == 'rbf_svm') {
    median_distance <- cached('median_distance', function() {
      distances <- as.numeric(dist(Z))^2
      value <- median(distances[distances > 0])
      if (is.finite(value)) value else 1
    })
    model$kernel_gamma <- parameters$bandwidth_multiplier / median_distance
    model$fits <- lapply(classes, function(label) {
      signed <- factor(ifelse(y == label, 1, -1), levels = c(-1, 1))
      weights <- n / (2 * table(signed))
      svm(Z, signed, kernel = 'radial', gamma = model$kernel_gamma,
                  cost = parameters$Cbar / n, class.weights = weights,
                  scale = FALSE, tolerance = 1e-4, shrinking = TRUE, cachesize = 256)
    })
    for (j in seq_along(model$fits)) {
      fit <- model$fits[[j]]
      signed <- ifelse(y == classes[j], 1, -1)
      costs <- parameters$Cbar / (2 * ifelse(signed > 0, sum(signed > 0), sum(signed < 0)))
      orientation <- if (fit$levels[fit$labels[1L]] == '1') 1 else -1
      pred <- predict(fit, Z, decision.values = TRUE)
      raw <- orientation * (as.vector(attr(pred, 'decision.values')) + fit$rho)
      model$fits[[j]]$rho <- -orientation * minimum_norm_intercept(raw, signed, costs)
    }
  } else if (method %in% c('shrinkage_lda', 'separable_lda')) {
    means <- cached('class_means', function()
      do.call(rbind, lapply(seq_len(k), function(j) colMeans(Z[group == j, , drop = FALSE]))))
    residual <- cached('class_residual', function() Z - means[group, , drop = FALSE])
    model$means <- means
    model$log_priors <- rep(-log(k), k)
    if (method == 'shrinkage_lda') {
      if (ncol(Z) > n) {
        target <- sum(residual^2) / ((n - k) * ncol(Z))
        if (target <= .Machine$double.eps) target <- 1
        diagonal <- parameters$lambda * target
        U <- residual * sqrt((1 - parameters$lambda) / (n - k))
        gram <- cached('lda_gram', function() tcrossprod(residual) / (n - k))
        R <- chol(diag(n) + (1 - parameters$lambda) * gram / diagonal)
        model$weights <- t(means) / diagonal - crossprod(U,
          .smm_chol_solve(R, U %*% t(means))) / diagonal^2
        model$linear_solver <- 'woodbury'
      } else {
        covariance <- cached('lda_covariance', function() crossprod(residual) / (n - k))
        covariance <- .benchmark_shrink(covariance, parameters$lambda)
        model$precision <- chol2inv(chol(covariance))
        model$weights <- model$precision %*% t(means)
        model$linear_solver <- 'primal'
      }
      model$offset <- -rowSums(means * t(model$weights)) / 2 + model$log_priors
    } else {
      covariances <- cached('separable_covariances', function() {
        E <- lapply(seq_len(n), function(i) matrix(residual[i, ], p, q))
        list(rows = covariance_sum(E, tcrossprod) / ((n - k) * q),
             times = covariance_sum(E, crossprod) / ((n - k) * p))
      })
      rows <- covariances$rows; times <- covariances$times
      time_scale <- mean(diag(times))
      if (time_scale <= .Machine$double.eps) times <- diag(q) else times <- times / time_scale
      model$row_covariance <- .benchmark_shrink(rows, parameters$lambda_row)
      model$time_covariance <- .benchmark_shrink(times, parameters$lambda_time)
      model$row_precision <- chol2inv(chol(model$row_covariance))
      model$time_precision <- chol2inv(chol(model$time_covariance))
      model$weights <- vapply(seq_len(k), function(j) {
        as.vector(model$row_precision %*% matrix(means[j, ], p, q) %*% model$time_precision)
      }, numeric(p * q))
      model$offset <- -rowSums(means * t(model$weights)) / 2 + model$log_priors
    }
  } else stop('Unknown comparator: ', method)
  class(model) <- 'smm_comparator'
  if (!model$converged) stop('Nonconverged comparator: ', method)
  return(model)
}

#' Predict labels or scores using training-estimated comparator preprocessing.
#' @param object Fitted comparator object.
#' @param newdata New observations with the training matrix shape.
#' @param type Return class labels or decision scores.
predict_comparator <- function(object, newdata, type = c('class', 'decision')) {
  type <- match.arg(type)
  packed <- .smm_pack(newdata)
  stopifnot(packed$p == object$p, packed$q == object$q)
  Z <- sweep(sweep(packed$X, 2L, rep(object$center, object$q), '-'),
             2L, rep(object$scale, object$q), '/')
  matrices <- lapply(seq_len(nrow(Z)), function(i) matrix(Z[i, ], object$p, object$q))
  if (object$method %in% c('summary_svm', 'functional_svm')) {
    inputs <- .benchmark_transform(matrices, object$method, object$parameters)
    scores <- predict_model(object$fit, inputs, type = 'decision')
  } else if (object$method == 'rbf_svm') {
    scores <- vapply(object$fits, function(fit) {
      prediction <- predict(fit, Z, decision.values = TRUE)
      orientation <- if (fit$levels[fit$labels[1L]] == '1') 1 else -1
      as.vector(attr(prediction, 'decision.values')) * orientation
    }, numeric(nrow(Z)))
    scores <- matrix(scores, nrow = nrow(Z), ncol = length(object$classes))
  } else scores <- sweep(Z %*% object$weights, 2L, object$offset, '+')
  colnames(scores) <- object$classes
  if (type == 'decision') return(scores)
  tolerance <- 64 * .Machine$double.eps * pmax(1, apply(abs(scores), 1L, max))
  labels <- object$classes[max.col(scores >= apply(scores, 1L, max) - tolerance, ties.method = 'first')]
  return(if (object$factor_labels) factor(labels, levels = object$classes) else labels)
}
