#' Tune one direct multiclass ablation on a fixed five-fold partition.
#' @param train List with matrix observations X, factor y and integer folds (1--5).
#' @param method M00, M10, M01 or M11: activation of nuclear and adjacency penalties.
#' @param seed TPE RNG seed, independent of the data-generator seed.
#' @param budget Total proposals including fixed, unpruned anchors.
#' @param gamma_max Upper endpoint of the positive adjacency-penalty domain.
#' @param structure Temporal or spatial adjacency.
#' @param channels Number of spatial colour planes (1 or 3).
#' @param fitter Function with the fit_direct interface, injectable for workflow tests.
#' @param predictor Function with the predict_model interface, injectable for tests.
#' @return Selected parameters, all trial scores/statuses and completed fold counts.
tune_direct <- function(train, method = 'M11', seed = 91001L, budget = 100L,
    gamma_max = 1e4, structure = 'temporal', channels = 1L,
    fitter = fit_direct, predictor = predict_model) {
  stopifnot(length(train$X) == length(train$y), length(train$folds) == length(train$y),
    is.factor(train$y), !anyNA(train$folds), setequal(train$folds, 1:5),
    all(table(train$y, train$folds) > 0))
  py <- import_from_path('tpe', file.path(.repo, 'python'))
  search <- py$Search(method, as.integer(seed), as.integer(budget), gamma_max)
  contexts <- lapply(1:5, function(fold) {
    ids <- which(train$folds != fold)
    return(smm_context(train$X[ids], train$y[ids], structure, channels = channels))
  })
  records <- vector('list', budget)
  for (i in seq_len(budget)) {
    request <- search$ask()
    prediction <- rep(NA_character_, length(train$y)); pruned <- FALSE
    for (fold in 1:5) {
      ids <- which(train$folds == fold)
      fit <- do.call(fitter, c(list(context = contexts[[fold]]), request$parameters))
      prediction[ids] <- as.character(predictor(fit, train$X[ids]))
      known <- !is.na(prediction)
      score <- mc_metrics(train$y[known], prediction[known])$BA
      pruned <- search$report(score, ba_upper_bound(train$y, prediction), as.integer(fold))
      if (pruned) break
    }
    score <- if (pruned) NA_real_ else mc_metrics(train$y, prediction)$BA
    search$tell(if (pruned) NULL else score, pruned = pruned)
    records[[i]] <- data.frame(trial = i, Cbar = request$parameters$Cbar,
      tau = request$parameters$tau, gamma = request$parameters$gamma,
      anchor = request$anchor, status = if (pruned) 'pruned' else 'complete',
      score = score, folds = fold)
  }
  return(list(selection = search$best(), trials = do.call(rbind, records), method = method,
    seed = seed, gamma_max = gamma_max))
}

#' Select a fixed-parameter model using pooled CV, without a TEST dependency.
#' @param X Development observations.
#' @param y Factor development labels.
#' @param folds Fixed fold assignments.
#' @param candidates Data frame of explicit candidate parameter values.
#' @param fit Function (X, y, parameters) that fits one candidate.
#' @param predict Function (model, X) returning class labels.
#' @return CV scores and the index of the best candidate; input order breaks ties.
select_grid <- function(X, y, folds, candidates, fit, predict) {
  stopifnot(length(X) == length(y), length(folds) == length(y), !anyNA(folds),
    all(table(y, folds) > 0), nrow(candidates) > 0)
  scores <- vapply(seq_len(nrow(candidates)), function(j) {
    pred <- rep(NA_character_, length(y))
    for (fold in sort(unique(folds))) {
      ids <- which(folds == fold)
      model <- fit(X[-ids], y[-ids], as.list(candidates[j, , drop = FALSE]))
      pred[ids] <- as.character(predict(model, X[ids]))
    }
    return(mc_metrics(y, pred)$BA)
  }, 0)
  return(list(index = which.max(scores), scores = scores))
}
