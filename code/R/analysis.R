#' Summarize paired effects with seeds averaged within each independent replication.
#' @param comparisons Data frame produced by run_additional: scenario, replication,
#'   tpe_seed, contrast, BA_difference, F1_difference; both prespecified seeds are required.
#' @return Replication means and exploratory t intervals/tests; Holm correction is
#'   applied across the provided scenarios/contrasts separately for each metric.
summarize_additional <- function(comparisons) {
  keys <- c('scenario', 'replication', 'tpe_seed', 'contrast')
  stopifnot(!anyDuplicated(comparisons[keys]), !anyNA(comparisons))
  blocks <- split(comparisons, interaction(comparisons[c('scenario', 'replication', 'contrast')], drop = TRUE))
  stopifnot(all(vapply(blocks, function(df) setequal(df$tpe_seed, c(91001L, 91002L)), TRUE)))
  means <- stats::aggregate(comparisons[c('BA_difference', 'F1_difference')],
    comparisons[c('scenario', 'replication', 'contrast')], mean)
  groups <- split(means, interaction(means[c('scenario', 'contrast')], drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(df) {
    stopifnot(setequal(df$replication, 1:10))
    do.call(rbind, lapply(c('BA_difference', 'F1_difference'), function(metric) {
      x <- df[[metric]]; n <- length(x); estimate <- mean(x); se <- sd(x) / sqrt(n)
      # Degenerate differences have no estimable sampling variance for a t test.
      p <- if (se > 0) 2 * pt(-abs(estimate / se), n - 1) else NA_real_
      data.frame(scenario = df$scenario[1], contrast = df$contrast[1], metric = metric,
        replications = n, estimate = estimate, lower = estimate - qt(.975, n - 1) * se,
        upper = estimate + qt(.975, n - 1) * se, p = p)
    }))
  }))
  summary$p_holm <- ave(summary$p, summary$metric, FUN = function(p) p.adjust(p, 'holm'))
  rownames(summary) <- NULL
  list(replication_means = means, summary = summary)
}
