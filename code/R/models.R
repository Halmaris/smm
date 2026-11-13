#' Compile the study's FP32 kernel with the macOS Accelerate library.
#' @param source_file Path to admm_fp32.cpp.
#' @return Backend information; no data generation or fitting is performed.
compile_solver <- function(source_file) {
  if (Sys.info()[['sysname']] != 'Darwin') stop('This backend currently requires macOS Accelerate.')
  cache <- file.path(tempdir(), 'smm-native')
  dir.create(cache, showWarnings = FALSE)
  old <- Sys.getenv(c('PKG_LIBS', 'R_MAKEVARS_USER'), unset = NA_character_)
  on.exit(for (name in names(old)) {
    if (is.na(old[[name]])) Sys.unsetenv(name) else
      do.call(Sys.setenv, setNames(list(old[[name]]), name))
  })
  sdk <- system2('xcrun', '--show-sdk-path', stdout = TRUE)
  resource <- system2('/usr/bin/clang++', '-print-resource-dir', stdout = TRUE)
  makevars <- file.path(cache, 'Makevars')
  writeLines(c('CXX=/usr/bin/clang++ --no-default-config -std=gnu++17',
    paste('CXXFLAGS=-O2 -ffp-contract=off -nostdinc -nostdinc++',
      '-isystem', shQuote(file.path(sdk, 'usr/include/c++/v1')),
      '-isystem', shQuote(file.path(resource, 'include')),
      '-isystem', shQuote(file.path(sdk, 'usr/include'))),
    paste('CPPFLAGS=-iframework', shQuote(file.path(sdk, 'System/Library/Frameworks')))), makevars)
  Sys.setenv(PKG_LIBS = '-framework Accelerate', R_MAKEVARS_USER = makevars)
  Rcpp::sourceCpp(source_file, env = environment(compile_solver), cacheDir = cache, showOutput = FALSE)
  invisible(smm_fp32_info())
}

#' Fit the coupled multiclass objective with per-class nuclear and adjacency penalties.
#' @param context Partition-specific context from smm_context.
#' @param Cbar Multiplier of mean multiclass hinge loss.
#' @param tau Nonnegative nuclear penalty.
#' @param gamma Nonnegative adjacency penalty.
#' @param max_iter Maximum total ADMM updates, including the bounded FP64 correction.
#' @param relaxation ADMM overrelaxation, 1.5 for the new structured methods.
#' @return A fitted model with coefficients and achieved numerical diagnostics.
fit_direct <- function(context, Cbar = 1, tau = .1, gamma = .1, max_iter = 300000L, relaxation = 1.5) {
  stopifnot(Cbar > 0, tau >= 0, gamma >= 0, max_iter > 0)
  design <- context$Z; n <- nrow(design); d <- ncol(design)
  p <- context$p; q <- context$q; k <- length(context$labels$classes)
  yy <- match(context$y, context$labels$classes); costs <- rep(Cbar / n, n)
  h <- 1 / sqrt(max(1, mean(rowSums(design^2))))
  rho <- if (tau >= 16 || gamma >= 256) 2^ceiling(log2(max(16, tau, sqrt(gamma)))) else 1
  W <- S <- U <- matrix(0, d, k); M <- V <- matrix(0, n, k)
  inverse_H <- .smm_cached_inverse(context, gamma, 1)
  assess <- function(out) {
    raw <- design %*% out$W; scores <- h * raw
    pr <- sqrt(sum((out$W - out$S)^2) + sum((scores - out$M)^2))
    dr <- out$diagnostic_rho * sqrt(sum((out$S - out$S_old + h * crossprod(design, out$M - out$M_old))^2))
    pt <- sqrt((d + n) * k) * 1e-4 + 1e-3 * max(sqrt(sum(out$W^2) + sum(scores^2)), sqrt(sum(out$S^2) + sum(out$M^2)))
    dt <- sqrt(d * k) * 1e-4 + 1e-3 * out$rho * sqrt(sum((out$U + h * crossprod(design, out$V))^2))
    cert <- .msmm_structured_certificate(out$W, design, yy, costs, context, tau, gamma,
      out$rho, out$U, h * out$V, inverse_H, raw_scores = raw)
    list(pr = pr, dr = dr, pt = pt, dt = dt, certificate = cert,
      valid = isTRUE(pr <= pt && dr <= dt && cert$weak_duality_ok && cert$relative_gap <= 1e-3))
  }
  zero <- .structured_zero(design, yy, costs, p, q, tau)
  used <- 0L; chunks <- 0L; polished <- FALSE; polish_updates <- 0L
  state <- list(W = W, S = S, U = U, M = M, V = V, b = numeric(k))
  if (!is.null(zero)) {
    state$U <- zero$N / rho; state$V <- zero$B / (rho * h)
    out <- c(state, list(S_old = S, M_old = M, rho = rho, diagnostic_rho = rho))
    final <- assess(out)
  } else {
    packed <- .smm_cached(context, 'design:fp32', function() smm_pack_fp32(design))
    control <- list(family = 'direct', p = p, q = q, tau = tau, gamma = gamma, h = h,
      rho = rho, abs_tol = 1e-4, rel_tol = 1e-3, gap_tol = 1e-3, interval = 25L,
      center = numeric(d), relaxation = relaxation, refinement = 3L,
      accept_fp64 = function(out) assess(out)$valid)
    final <- NULL
    repeat {
      control$max_iter <- as.integer(min(5000L, max_iter - used))
      control$adaptive_rho <- chunks == 0L && tau < 16 && gamma < 256
      stalled <- !is.null(final) && !polished && .structured_precision_stall(
        c(final[c('pr', 'dr', 'pt', 'dt')], list(relative_gap = final$certificate$relative_gap)), 1e-3)
      if (stalled) {
        out <- .structured_polish(state, design, yy, costs, context, tau, gamma, h,
          control$rho, relaxation, min(100L, max_iter - used), assess)
        polished <- TRUE; polish_updates <- out$iterations
      } else out <- smm_admm_fp32(packed, packed, matrix(0, 0, 0), matrix(costs, n, 1),
        as.integer(yy), state, control,
        function(rho) .smm_fp32_factor(context, gamma, rho, h, FALSE),
        .smm_fp32_inverse(context, gamma, 1))
      final <- assess(out); used <- used + out$iterations; chunks <- chunks + 1L
      if (final$valid || used >= max_iter) break
      next_rho <- if (.structured_reduce_rho(final)) max(out$rho / 2, 1e-6) else out$rho
      state <- list(W = out$W, S = out$S, U = out$U * out$rho / next_rho,
        M = out$M, V = out$V * out$rho / next_rho, b = numeric(k))
      control$rho <- next_rho
    }
  }
  if (!final$valid) stop(sprintf('Nonconvergence after %d updates; gap %.3g.', used, final$certificate$relative_gap))
  structure(c(list(W = out$W, b = numeric(k), p = p, q = q,
    classes = context$labels$classes, center = context$preprocessing$center,
    scale = context$preprocessing$scale, converged = TRUE, iterations = used,
    primal_residual = final$pr, dual_residual = final$dr,
    primal_tolerance = final$pt, dual_tolerance = final$dt,
    precision = if (polished) 'FP32+FP64_polish' else 'FP32', polish_updates = polish_updates,
    Cbar = Cbar, tau = tau, gamma = gamma), final$certificate), class = 'smm_model')
}

#' Fit balanced one-versus-rest margins with a free intercept.
#' @param context Partition-specific context from smm_context.
#' @param Cbar Multiplier of each balanced binary mean hinge loss.
#' @param tau Nuclear penalty.
#' @param gamma Temporal or spatial penalty.
#' @param max_iter Maximum ADMM updates.
#' @param relaxation Overrelaxation; 1.5 for structured SMM, 1 for legacy margins.
#' @return Coefficients, intercepts and per-class achieved numerical diagnostics.
fit_ovr <- function(context, Cbar = 1, tau = .1, gamma = .1, max_iter = 300000L, relaxation = 1.5) {
  Z <- context$Z; n <- nrow(Z); d <- ncol(Z); p <- context$p; q <- context$q
  k <- length(context$labels$classes)
  Y <- outer(as.character(context$y), as.character(context$labels$classes), function(a, b) ifelse(a == b, 1, -1))
  costs <- apply(Y, 2L, function(y) Cbar / (2 * ifelse(y > 0, sum(y > 0), sum(y < 0))))
  h <- 1 / sqrt(max(1, mean(rowSums(context$Zc^2))))
  inverse <- .smm_cached_inverse(context, gamma, 1)
  assess <- function(out) {
    margin <- h * Y * sweep(Z %*% out$W, 2L, out$b, '+')
    pr <- sqrt(colSums((out$W - out$S)^2) + colSums((margin - out$M)^2))
    dz <- h * Y * (out$M - out$M_old)
    dr <- out$diagnostic_rho * sqrt(colSums((out$S - out$S_old + crossprod(Z, dz))^2) + colSums(dz)^2)
    pt <- sqrt(d + n) * 1e-4 + 1e-3 * pmax(sqrt(colSums(out$W^2) + colSums(margin^2)), sqrt(colSums(out$S^2) + colSums(out$M^2)))
    dt <- sqrt(d + 1) * 1e-4 + 1e-3 * out$rho * sqrt(colSums((out$U + h * crossprod(Z, Y * out$V))^2) + h^2 * colSums(Y * out$V)^2)
    certificates <- lapply(seq_len(k), function(j) .smm_certificate(out$W[, j], out$b[j], Z, Y[, j],
      costs[, j], p, q, tau, gamma, NULL, out$rho, out$U[, j], h * out$V[, j], inverse,
      details = FALSE, roughness_fn = context$roughness_fn))
    gaps <- vapply(certificates, function(x) x$relative_gap, 0)
    weak <- vapply(certificates, function(x) x$duality_gap >= -1e-10 * max(1, abs(x$objective), abs(x$dual_objective)), TRUE)
    list(valid = all(pr <= pt & dr <= dt & gaps <= 1e-3 & weak), certificates = certificates,
      diagnostics = data.frame(relative_gap = gaps, primal_residual = pr, dual_residual = dr,
        primal_tolerance = pt, dual_tolerance = dt))
  }
  initial <- list(W = matrix(0, d, k), S = matrix(0, d, k), U = matrix(0, d, k),
    M = matrix(0, n, k), V = matrix(0, n, k), b = numeric(k))
  control <- list(family = 'ovr', p = p, q = q, tau = tau, gamma = gamma, h = h,
    rho = 1, max_iter = as.integer(max_iter), abs_tol = 1e-4, rel_tol = 1e-3, gap_tol = 1e-3,
    adaptive_rho = TRUE, interval = 25L, center = context$center,
    relaxation = relaxation, refinement = 3L,
    accept_fp64 = function(out) assess(out)$valid)
  out <- smm_admm_fp32(Z, context$Zc, Y, costs, integer(), initial, control,
    function(rho) .smm_fp32_factor(context, gamma, rho, h, TRUE), .smm_fp32_inverse(context, gamma, 1))
  final <- assess(out)
  if (!final$valid) stop('OvR solver did not meet its numerical acceptance criteria.')
  # Canonicalize flat intercept minima and analytically certified zero columns.
  G <- crossprod(Z, Y * costs)
  zero <- vapply(seq_len(k), function(j)
    svd(matrix(G[, j], p, q), nu = 0, nv = 0)$d[1L] <= tau * (1 - 1e-12), TRUE)
  out$W[, zero] <- 0
  raw <- Z %*% out$W
  old_b <- out$b
  out$b <- vapply(seq_len(k), function(j)
    minimum_norm_intercept(raw[, j], Y[, j], costs[, j]), 0)
  for (j in seq_len(k)) {
    cert <- final$certificates[[j]]
    new_objective <- if (zero[j]) sum(costs[, j]) else cert$objective +
      sum(costs[, j] * (pmax(0, 1 - Y[, j] * (raw[, j] + out$b[j])) -
        pmax(0, 1 - Y[, j] * (raw[, j] + old_b[j]))))
    dual <- if (zero[j]) sum(costs[, j]) else cert$dual_objective
    gap <- abs(new_objective - dual) / max(1, abs(new_objective), abs(dual))
    stopifnot(new_objective <= cert$objective + 1e-9 * max(1, abs(cert$objective)), gap <= 1e-3)
    final$diagnostics$relative_gap[j] <- gap
  }
  structure(list(W = out$W, b = out$b, p = p, q = q, classes = context$labels$classes,
    center = context$preprocessing$center, scale = context$preprocessing$scale,
    Cbar = Cbar, tau = tau, gamma = gamma, converged = TRUE, iterations = out$iterations,
    diagnostics = final$diagnostics, precision = 'FP32', canonical_scores = TRUE), class = 'smm_model')
}

#' Predict with a fitted matrix classifier.
#' @param object Model returned by fit_direct or fit_ovr.
#' @param newdata New matrix-valued observations.
#' @param type Return class labels or the matrix of class scores.
predict_model <- function(object, newdata, type = 'class') {
  type <- match.arg(type, c('class', 'decision'))
  packed <- .smm_pack(newdata)
  stopifnot(packed$p == object$p, packed$q == object$q)
  Z <- sweep(sweep(packed$X, 2L, rep(object$center, object$q), '-'), 2L, rep(object$scale, object$q), '/')
  scores <- sweep(Z %*% object$W, 2L, object$b, '+')
  colnames(scores) <- object$classes
  if (type == 'decision') return(scores)
  if (isTRUE(object$canonical_scores)) scores <- scores >= apply(scores, 1L, max) -
    64 * .Machine$double.eps * pmax(1, apply(abs(scores), 1L, max))
  factor(object$classes[max.col(scores, ties.method = 'first')], levels = object$classes)
}

#' Select the minimum-absolute-value intercept among weighted hinge minimizers.
#' @param scores Training decision scores without an intercept.
#' @param signed Binary labels coded as -1 and 1.
#' @param costs Positive observation weights.
minimum_norm_intercept <- function(scores, signed, costs) {
  knots <- signed - scores; order <- order(knots)
  cumulative <- cumsum(costs[order] / sum(costs))
  target <- sum(costs[signed > 0]) / sum(costs)
  tolerance <- 64 * .Machine$double.eps
  j <- which(cumulative >= target - tolerance)[1L]
  lower <- knots[order[j]]
  upper <- if (abs(cumulative[j] - target) <= tolerance && j < length(knots))
    knots[order[j + 1L]] else lower
  max(lower, min(0, upper))
}

#' Fit one named study method through the common numerical backend.
#' @param X List of matrix-valued training observations.
#' @param y Training class labels.
#' @param method One of method_names().
#' @param parameters Named hyperparameters; Cbar, tau, gamma or comparator parameters.
#' @param channels Number of disconnected RGB planes for spatial methods.
#' @return A fitted classifier, usable with predict_method.
fit_method <- function(X, y, method, parameters, channels = 1L) {
  stopifnot(method %in% names(method_names()))
  if (method %in% c('rbf_svm', 'shrinkage_lda', 'separable_lda', 'summary_svm', 'functional_svm'))
    return(fit_comparator(X, y, method, parameters))
  spatial <- method %in% c('spatial_svm', 'spatial_smm', 'msmm_ss')
  context <- smm_context(X, y, if (spatial) 'spatial' else 'temporal',
    channels = if (spatial) channels else 1L)
  if (method %in% c('crammer_singer', 'msmm', 'msmm_ts', 'msmm_ss')) {
    tau <- if (method == 'crammer_singer') 0 else parameters$tau
    gamma <- if (method %in% c('msmm_ts', 'msmm_ss')) parameters$gamma else 0
    return(fit_direct(context, parameters$Cbar, tau, gamma,
      relaxation = if (method %in% c('msmm_ts', 'msmm_ss')) 1.5 else 1))
  }
  tau <- if (method %in% c('smm', 'temporal_smm', 'spatial_smm')) parameters$tau else 0
  gamma <- if (method %in% c('temporal_svm', 'temporal_smm', 'spatial_svm', 'spatial_smm')) parameters$gamma else 0
  fit_ovr(context, parameters$Cbar, tau, gamma,
    relaxation = if (method %in% c('temporal_smm', 'spatial_smm')) 1.5 else 1)
}

#' Predict using a model produced by any study method.
#' @param object Fitted model from fit_method.
#' @param X Matrix-valued test observations.
#' @param type Class labels or decision scores.
predict_method <- function(object, X, type = 'class') {
  if (inherits(object, 'smm_comparator')) predict_comparator(object, X, type) else
    predict_model(object, X, type)
}

#' List the scientific model families supported by the shared backend.
#' @return Named character vector of method identifiers and labels.
method_names <- function() {
  c(linear_svm = 'Linear SVM (OvR)', smm = 'SMM (OvR)',
    temporal_svm = 'SVM-TS (OvR)', temporal_smm = 'SMM-TS (OvR)',
    spatial_svm = 'SVM-SS (OvR)', spatial_smm = 'SMM-SS (OvR)',
    crammer_singer = 'Crammer-Singer', msmm = 'MSMM', msmm_ts = 'MSMM-TS',
    msmm_ss = 'MSMM-SS', rbf_svm = 'RBF SVM (OvR)', shrinkage_lda = 'Shrinkage LDA',
    separable_lda = 'Separable LDA', summary_svm = 'Mean-and-trend SVM',
    functional_svm = 'Temporal-basis SVM')
}
