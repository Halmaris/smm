# Numerical building blocks for the common FP32 solver and FP64 checks.
#' Validate a finite scalar numerical parameter.
#' @param x Numeric input or diagnostic record described by this function.
#' @param name Parameter name used in error messages.
#' @param lower Scalar lower bound.
#' @param strict Whether the lower bound is strict.
#' @param integer Whether an integral value is required.
.smm_scalar <- function(x, name, lower = 0, strict = FALSE, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      if (strict) x <= lower else x < lower) {
    stop(name, ' must be a finite scalar ', if (strict) '> ' else '>= ', lower,
         call. = FALSE)
  }
  if (integer && x != floor(x)) stop(name, ' must be an integer.', call. = FALSE)
  invisible(x)
}

#' Apply singular-value soft thresholding, the proximal map of the nuclear norm.
#' @param X Matrix-valued observations or the vectorized design matrix for this operation.
#' @param tau Nonnegative nuclear-norm penalty.
shrinkage <- function(X, tau) {
  .smm_scalar(tau, 'tau')
  if (!is.matrix(X) || !is.numeric(X) || !length(X) || any(!is.finite(X))) {
    stop('X must be a nonempty finite numeric matrix.', call. = FALSE)
  }
  dec <- svd(X)
  values <- pmax(dec$d - tau, 0)
  D <- sweep(dec$u, 2L, values, '*') %*% t(dec$v)
  list(D = D, nuc = sum(values), rk = sum(values > 0), singular_values = values)
}

#' Apply nuclear-norm proximal maps independently to coefficient columns.
#' @param A Stacked coefficient or score matrix.
#' @param p Number of matrix rows.
#' @param q Number of matrix columns.
#' @param threshold Singular-value threshold.
.smm_batch_shrinkage <- function(A, p, q, threshold) {
  if (threshold == 0) return(A)
  answer <- A
  for (j in seq_len(ncol(A)))
    answer[, j] <- as.vector(shrinkage(matrix(A[, j], p, q), threshold)$D)
  answer
}

#' Solve a positive-definite system from its upper Cholesky factor.
#' @param R Upper triangular Cholesky factor.
#' @param rhs Vector or matrix of right-hand sides.
.smm_chol_solve <- function(R, rhs) backsolve(R, forwardsolve(t(R), rhs))

# Apply (H_time kron I_p)^(-1) without constructing a (p*q)-square matrix.

#' Estimate row-wise scaling on this training partition and transform observations.
#' @param X Matrix-valued observations or the vectorized design matrix for this operation.
#' @param p Number of matrix rows.
#' @param q Number of matrix columns.
#' @param standardize Whether to estimate row-wise centering and scaling.
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

#' Construct a feasible binary dual point and evaluate the primal-dual gap.
#' @param w Vectorized coefficient matrix for a single margin.
#' @param b Unpenalized intercept.
#' @param X Matrix-valued observations or the vectorized design matrix for this operation.
#' @param y Class labels; internal binary operations use -1/+1 and direct operations use 1..K.
#' @param C Positive binary loss weights.
#' @param p Number of matrix rows.
#' @param q Number of matrix columns.
#' @param tau Nonnegative nuclear-norm penalty.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param D Temporal difference matrix, when supplied explicitly.
#' @param rho Positive ADMM penalty parameter.
#' @param u Coefficient-constraint scaled dual variable.
#' @param v Score-constraint scaled dual variable or vector to project.
#' @param inverse_H0 Function applying the inverse quadratic penalty.
#' @param details Whether to retain additional spectral diagnostics.
#' @param roughness_fn Function evaluating the structural quadratic penalty on a matrix.
.smm_certificate <- function(w, b, X, y, C, p, q, tau, gamma, D,
                             rho, u, v, inverse_H0, details = TRUE, roughness_fn = NULL) {
  W <- matrix(w, p, q)
  singular_values <- if (tau > 0 || details) svd(W, nu = 0L, nv = 0L)$d else numeric()
  margin <- y * as.vector(X %*% w + b)
  roughness <- if (gamma == 0 && !details) 0 else if (!is.null(roughness_fn)) roughness_fn(W) else if (!is.null(D)) sum((W %*% t(D))^2) else
    if (q < 3L) 0 else sum((W[, 1:(q-2L), drop = FALSE] -
      2 * W[, 2:(q-1L), drop = FALSE] + W[, 3:q, drop = FALSE])^2)
  objective <- sum(w^2) / 2 + tau * sum(singular_values) +
    gamma * roughness + sum(C * pmax(0, 1 - margin))
  alpha_raw <- -rho * v
  alpha <- pmin(C, pmax(0, alpha_raw))
  pos <- y > 0
  plus <- sum(alpha[pos])
  minus <- sum(alpha[!pos])
  # Scaling down the larger class total gives a feasible dual certificate.
  if (plus > minus) alpha[pos] <- alpha[pos] * minus / plus
  if (minus > plus) alpha[!pos] <- alpha[!pos] * plus / minus
  V <- matrix(rho * u, p, q)
  if (tau == 0) {
    V[] <- 0
  } else {
    spectral <- svd(V, nu = 0L, nv = 0L)$d[1L]
    if (spectral > tau) V <- V * tau / spectral
  }
  a <- as.vector(crossprod(X, y * alpha)) - as.vector(V)
  dual <- sum(alpha) - sum(a * inverse_H0(a)) / 2
  gap <- objective - dual
  relative_gap <- abs(gap) / max(1, abs(objective), abs(dual))
  list(objective = objective, dual_objective = dual, duality_gap = gap,
       relative_gap = relative_gap, singular_values = singular_values,
       roughness = roughness, alpha = alpha, nuclear_dual = V,
       dual_balance = abs(sum(y * alpha)),
       raw_dual_balance = abs(sum(y * alpha_raw)))
}

# Objective: 0.5*||W||_F^2 + tau*||W||_* + gamma*||W D'||_F^2
#            + sum_i C_i * max(0, 1 - y_i*(<W,X_i> + b)).
#' Project a vector onto the nonnegative simplex of a specified mass.
#' @param v Score-constraint scaled dual variable or vector to project.
#' @param radius Nonnegative simplex mass.
.msmm_simplex <- function(v, radius) {
  if (radius == 0) return(v * 0)
  u <- sort(v, decreasing = TRUE)
  css <- cumsum(u)
  active <- which(u > (css - radius) / seq_along(u))
  j <- max(active)
  pmax(v - (css[j] - radius) / j, 0)
}

# Proximal map of C * max(0, max_{k != y}(1 + z_k - z_y)).

#' Evaluate the coupled multiclass hinge proximal map row by row.
#' @param A Stacked coefficient or score matrix.
#' @param y Class labels; internal binary operations use -1/+1 and direct operations use 1..K.
#' @param cost Proximal loss weights.
.msmm_hinge_prox <- function(A, y, cost) {
  answer <- A
  for (i in seq_len(nrow(A))) {
    wrong <- seq_len(ncol(A)) != y[i]
    c <- 1 + A[i, wrong] - A[i, y[i]]
    u <- sort(c, decreasing = TRUE)
    cs <- cumsum(u)
    active <- which(u > cs / (seq_along(u) + 1))
    total <- if (length(active)) max(0, cs[max(active)] / (max(active) + 1)) else 0
    alpha <- if (total <= cost[i]) pmax(c - total, 0) else .msmm_simplex(c, cost[i])
    answer[i, wrong] <- A[i, wrong] - alpha
    answer[i, y[i]] <- A[i, y[i]] + sum(alpha)
  }
  answer
}

#' Evaluate the coupled multiclass hinge loss for each observation.
#' @param scores Observation-by-class decision scores.
#' @param y Class labels; internal binary operations use -1/+1 and direct operations use 1..K.
.msmm_loss <- function(scores, y) {
  wrong <- scores - scores[cbind(seq_along(y), y)] + 1
  wrong[cbind(seq_along(y), y)] <- 0
  apply(wrong, 1L, max)
}
#' Construct the sparse Gram matrix of temporal second differences.
#' @param q Number of matrix columns.
.smm_banded_penalty <- function(q) {
  if (q < 3L) return(Matrix::Diagonal(q, 0))
  # D has three nonzero entries per row; D'D has bandwidth two.
  ii <- rep(seq_len(q - 2L), each = 3L)
  jj <- as.vector(t(outer(seq_len(q - 2L), 0:2, '+')))
  D <- Matrix::sparseMatrix(i = ii, j = jj, x = rep(c(1, -2, 1), q - 2L),
                            dims = c(q - 2L, q))
  Matrix::crossprod(D)
}

#' Factor and apply the quadratic penalty inverse in temporal or spatial coordinates.
#' @param context Training-partition context containing design, adjacency and local caches.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param diagonal Positive diagonal term in the quadratic operator.
.smm_cached_inverse <- function(context, gamma, diagonal) {
  if (gamma == 0) return(function(rhs) rhs / diagonal)
  p <- context$p; q <- context$q
  key <- paste('inverse', gamma, diagonal, sep = ':')
  .smm_cached(context, key, function() {
    if (identical(context$structure, 'spatial'))
      return(.smm_spatial_inverse(context, gamma, diagonal))
    penalty <- .smm_cached(context, 'banded_penalty', function() .smm_banded_penalty(q))
    H <- Matrix::Diagonal(q, diagonal) + 2 * gamma * penalty
    factor <- Matrix::Cholesky(H, perm = FALSE, LDL = FALSE, super = FALSE)
    function(rhs) {
      vector_rhs <- is.null(dim(rhs))
      M <- matrix(rhs, nrow = p * q)
      blocks <- matrix(aperm(array(M, c(p, q, ncol(M))), c(2L, 1L, 3L)), nrow = q)
      solved <- as.matrix(Matrix::solve(factor, blocks, system = 'A'))
      answer <- matrix(aperm(array(solved, c(q, p, ncol(M))), c(2L, 1L, 3L)), nrow = p * q)
      if (vector_rhs) as.vector(answer) else answer
    }
  })
}

#' Build a primal or Woodbury solve for the ADMM coefficient update.
#' @param context Training-partition context containing design, adjacency and local caches.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param rho Positive ADMM penalty parameter.
#' @param h Score preconditioning scale.
#' @param centered Whether to use the centered design required for a free intercept.
.smm_cached_solver <- function(context, gamma, rho, h, centered) {
  key <- paste('solver', gamma, rho, format(h, digits = 17), centered, sep = ':')
  .smm_cached(context, key, function() {
    X <- if (centered) context$Zc else context$Z
    n <- nrow(X); d <- ncol(X); p <- context$p; q <- context$q
    diagonal <- 1 + rho
    gram <- .smm_cached(context, paste0('gram:', centered), function()
      if (d <= n) crossprod(X) else tcrossprod(X))
    if (d <= n) {
      H <- diagonal * diag(d)
      if (gamma > 0) {
        if (identical(context$structure, 'spatial')) {
          H <- H + 2 * gamma * (kronecker(diag(q), context$row_laplacian) +
            kronecker(context$column_laplacian, diag(p)))
        } else {
          penalty <- .smm_cached(context, 'banded_penalty', function() .smm_banded_penalty(q))
          H <- H + 2 * gamma * kronecker(as.matrix(penalty), diag(p))
        }
      }
      R <- chol(H + rho * h^2 * gram)
      return(function(rhs) .smm_chol_solve(R, rhs))
    }
    if (gamma == 0) {
      R <- chol(diag(n) / (rho * h^2) + gram / diagonal)
      return(function(rhs) rhs / diagonal - crossprod(X,
        .smm_chol_solve(R, X %*% rhs)) / diagonal^2)
    }
    inverse <- .smm_cached_inverse(context, gamma, diagonal)
    V <- inverse(t(X)); R <- chol(diag(n) / (rho * h^2) + X %*% V)
    function(rhs) {
      a <- inverse(rhs)
      a - V %*% .smm_chol_solve(R, X %*% a)
    }
  })
}
#' Apply a Kronecker operator through left and right matrix products.
#' @param rhs Vector or matrix of right-hand sides.
#' @param left Left matrix factor of the Kronecker action.
#' @param right Right matrix factor of the Kronecker action.
.smm_two_sided_action <- function(rhs, left, right) {
  p <- nrow(left); q <- nrow(right)
  M <- matrix(rhs, nrow = p * q); k <- ncol(M)
  B <- left %*% matrix(M, nrow = p)
  C <- right %*% matrix(aperm(array(B, c(p, q, k)), c(2L, 1L, 3L)), nrow = q)
  answer <- matrix(aperm(array(C, c(q, p, k)), c(2L, 1L, 3L)), nrow = p * q)
  if (is.null(dim(rhs))) as.vector(answer) else answer
}

#' Invert the spatial quadratic penalty in its separable eigenbasis.
#' @param context Training-partition context containing design, adjacency and local caches.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param diagonal Positive diagonal term in the quadratic operator.
.smm_spatial_inverse <- function(context, gamma, diagonal) {
  spectrum <- .smm_cached(context, 'spatial-spectrum', function() {
    r <- eigen(context$row_laplacian, symmetric = TRUE)
    c <- eigen(context$column_laplacian, symmetric = TRUE)
    list(row = r$vectors, column = c$vectors,
      values = as.vector(outer(pmax(r$values, 0), pmax(c$values, 0), '+')))
  })
  denominator <- diagonal + 2 * gamma * spectrum$values
  function(rhs) {
    transformed <- .smm_two_sided_action(rhs, t(spectrum$row), t(spectrum$column))
    .smm_two_sided_action(transformed / denominator, spectrum$row, spectrum$column)
  }
}
#' Project multiclass dual variables onto feasible sets and compute an FP64 gap.
#' @param W Stacked coefficient columns, one per class.
#' @param X Matrix-valued observations or the vectorized design matrix for this operation.
#' @param y Class labels; internal binary operations use -1/+1 and direct operations use 1..K.
#' @param costs Positive multiclass loss weights.
#' @param context Training-partition context containing design, adjacency and local caches.
#' @param tau Nonnegative nuclear-norm penalty.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param rho Positive ADMM penalty parameter.
#' @param U Coefficient-constraint scaled dual matrix.
#' @param V Score-constraint scaled dual matrix.
#' @param inverse_H Function applying the inverse quadratic penalty.
#' @param raw_scores Optional precomputed unscaled class scores.
.msmm_structured_certificate <- function(W, X, y, costs, context, tau, gamma,
                                         rho, U, V, inverse_H, raw_scores = NULL) {
  k <- ncol(W); p <- context$p; q <- context$q
  nuclear <- sum(vapply(seq_len(k), function(j)
    sum(svd(matrix(W[, j], p, q), nu = 0L, nv = 0L)$d), numeric(1)))
  roughness <- sum(vapply(seq_len(k), function(j)
    context$roughness_fn(matrix(W[, j], p, q)), numeric(1)))
  primal <- sum(W^2) / 2 + tau * nuclear + gamma * roughness +
    sum(costs * .msmm_loss(if (is.null(raw_scores)) X %*% W else raw_scores, y))
  B <- pmax(rho * V, 0)
  target <- cbind(seq_len(nrow(B)), y)
  B[target] <- 0
  mass <- rowSums(B)
  scaled <- which(mass > costs)
  # Preserve multiplication followed by division, rather than a rounded ratio.
  if (length(scaled)) B[scaled, ] <-
    B[scaled, , drop = FALSE] * costs[scaled] / mass[scaled]
  B[target] <- -rowSums(B)
  N <- rho * U
  if (tau == 0) N[] <- 0 else for (j in seq_len(k)) {
    spectral <- svd(matrix(N[, j], p, q), nu = 0L, nv = 0L)$d[1L]
    if (spectral > tau) N[, j] <- N[, j] * tau / spectral
  }
  linear <- crossprod(X, B) + N
  dual <- -sum(B[cbind(seq_along(y), y)]) - sum(linear * inverse_H(linear)) / 2
  gap <- primal - dual
  list(objective = primal, dual_objective = dual, duality_gap = gap,
    relative_gap = abs(gap) / max(1, abs(primal), abs(dual)),
    weak_duality_ok = is.finite(gap) && gap >= -1e-10 * max(1, abs(primal), abs(dual)),
    nuclear_norm = nuclear, roughness = roughness, score_dual = B, nuclear_dual = N)
}

#' Detect an exact zero-coefficient solution using the KKT spectral condition.
#' @param X Matrix-valued observations or the vectorized design matrix for this operation.
#' @param y Class labels; internal binary operations use -1/+1 and direct operations use 1..K.
#' @param costs Positive multiclass loss weights.
#' @param p Number of matrix rows.
#' @param q Number of matrix columns.
#' @param tau Nonnegative nuclear-norm penalty.
.structured_zero <- function(X, y, costs, p, q, tau) {
  k <- max(y); n <- nrow(X)
  B <- matrix(costs/(k-1L), n, k)
  B[cbind(seq_len(n), y)] <- -costs
  N <- -crossprod(X,B)
  norms <- vapply(seq_len(k),function(j) svd(matrix(N[,j],p,q),nu=0,nv=0)$d[1],numeric(1))
  if (!all(norms <= tau*(1-1e-12))) return(NULL)
  list(B=B,N=N)
}

#' Identify a dual-residual imbalance for the fixed penalty update rule.
#' @param check Residual and certificate assessment returned by the solver callback.
.structured_reduce_rho <- function(check) {
  isTRUE(!check$valid && all(is.finite(c(check$pr,check$dr,check$pt,check$dt,
    check$certificate$relative_gap))) && check$pr <= .1*check$pt && check$dr > check$dt)
}


# A small FP64 correction addresses a precision floor without relaxing acceptance.

#' Detect the bounded precision-floor condition for an FP64 correction.
#' @param x Numeric input or diagnostic record described by this function.
#' @param gap_tol Relative primal-dual gap tolerance.
.structured_precision_stall <- function(x, gap_tol) {
  !is.null(x) && isTRUE(all(is.finite(c(x$pr,x$dr,x$pt,x$dt,x$relative_gap))) &&
    x$pr <= .1*x$pt && x$dr <= .1*x$dt &&
    x$relative_gap > gap_tol && x$relative_gap <= 10*gap_tol)
}

#' Perform at most the requested number of FP64 correction updates with unchanged acceptance.
#' @param state Current ADMM primal and scaled dual variables.
#' @param design Training design matrix.
#' @param yy Class indices 1..K.
#' @param costs Positive multiclass loss weights.
#' @param context Training-partition context containing design, adjacency and local caches.
#' @param tau Nonnegative nuclear-norm penalty.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param h Score preconditioning scale.
#' @param rho Positive ADMM penalty parameter.
#' @param relaxation ADMM overrelaxation in (0, 2).
#' @param max_iter Maximum correction updates, included in the total budget.
#' @param assess Callback evaluating residuals and the FP64 certificate.
.structured_polish <- function(state, design, yy, costs, context, tau, gamma,
    h, rho, relaxation, max_iter, assess) {
  solve_w <- .smm_cached_solver(context, gamma, rho, h, FALSE)
  W<-state$W; S<-state$S; U<-state$U; M<-state$M; V<-state$V
  history<-list()
  for (iteration in seq_len(max_iter)) {
    W<-solve_w(rho*(S-U+h*crossprod(design,M-V)))
    S_old<-S; M_old<-M; scores<-h*(design%*%W)
    relaxed_W<-relaxation*W+(1-relaxation)*S_old
    relaxed_scores<-relaxation*scores+(1-relaxation)*M_old
    S<-.smm_batch_shrinkage(relaxed_W+U,context$p,context$q,tau/rho)
    M<-h*.msmm_hinge_prox((relaxed_scores+V)/h,yy,costs/(rho*h^2))
    U<-U+relaxed_W-S; V<-V+relaxed_scores-M
    out<-list(W=W,S=S,U=U,M=M,V=V,S_old=S_old,M_old=M_old,
      rho=rho,diagnostic_rho=rho,iterations=iteration)
    check<-assess(out)
    history[[iteration]]<-data.frame(iteration=iteration,objective=check$certificate$objective,
      dual_objective=check$certificate$dual_objective,relative_gap=check$certificate$relative_gap,
      primal_residual=check$pr,dual_residual=check$dr,rho=rho)
    if (check$valid) break
  }
  out$history<-do.call(rbind,history)
  out$certificate_checks<-iteration;out$fp64_acceptance_checks<-iteration
  out$fp64_rejected_stops<-0L
  out
}
#' Build or reuse the single-precision factor for an ADMM coefficient update.
#' @param context Training-partition context containing design, adjacency and local caches.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param rho Positive ADMM penalty parameter.
#' @param h Score preconditioning scale.
#' @param centered Whether to use the centered design required for a free intercept.
.smm_fp32_factor <- function(context, gamma, rho, h, centered) {
  key <- paste('solver:fp32', format(gamma, digits = 17), format(rho, digits = 17),
    format(h, digits = 17), centered, sep = ':')
  .smm_cached(context, key, function() {
    X <- if (centered) context$Zc else context$Z
    gram <- .smm_cached(context, paste0('gram:fp32:', centered),
      function() smm_gram_fp32(X))
    blocks <- if (is.null(context$row_blocks)) 1L else context$row_blocks
    if (blocks>1L) return(smm_factor_fp32(X,gram,context$p,context$q,gamma,rho,h,
      spatial=TRUE,row_blocks=blocks))
    smm_factor_fp32(X, gram, context$p, context$q, gamma, rho, h, spatial = identical(context$structure, 'spatial'))
  })
}

#' Build or reuse the single-precision quadratic penalty inverse.
#' @param context Training-partition context containing design, adjacency and local caches.
#' @param gamma Nonnegative quadratic adjacency penalty.
#' @param diagonal Positive diagonal term in the quadratic operator.
.smm_fp32_inverse <- function(context, gamma, diagonal) {
  .smm_cached(context, paste('native_inverse32', format(gamma, digits = 17),
    format(diagonal, digits = 17), sep = ':'), function() {
      blocks <- if (is.null(context$row_blocks)) 1L else context$row_blocks
      if (blocks>1L) return(smm_inverse_fp32(context$p,context$q,gamma,diagonal,
        spatial=TRUE,row_blocks=blocks))
      smm_inverse_fp32(context$p, context$q, gamma, diagonal, spatial = identical(context$structure, 'spatial'))
    })
}
