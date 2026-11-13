# Summarize paired replications; x contains finite per-replication differences.
paired_interval <- function(x) {
  stopifnot(length(x) > 1L, all(is.finite(x)))
  se <- sd(x) / sqrt(length(x)); estimate <- mean(x)
  p <- if (se > 0) 2 * pt(-abs(estimate / se), length(x) - 1L) else NA_real_
  width <- qt(.975, length(x) - 1L) * se
  data.frame(n = length(x), mean = estimate, MC_SE = se,
    lower = estimate - width, upper = estimate + width, p = p)
}

# Score saved labels; truth and prediction are equally long class-label vectors.
classification_metrics <- function(truth, prediction) {
  stopifnot(length(truth) > 0L, length(truth) == length(prediction),
    !anyNA(truth), !anyNA(prediction))
  classes <- sort(unique(as.character(truth)))
  stopifnot(all(prediction %in% classes))
  tab <- table(factor(truth, classes), factor(prediction, classes))
  tp <- diag(tab)
  data.frame(n = length(truth), balanced_accuracy = mean(tp / rowSums(tab)),
    macro_f1 = mean(2 * tp / (rowSums(tab) + colSums(tab))))
}

# Summarize original four-variant ablations; df has scenario/replication/method/BA/macro_F1.
multiclass_contrasts <- function(df) {
  methods <- c('M00', 'M10', 'M01', 'M11')
  stopifnot(!anyNA(df[c('scenario', 'replication', 'method', 'BA', 'macro_F1')]),
    !anyDuplicated(df[c('scenario', 'replication', 'method')]),
    setequal(df$method, methods), all(df$BA >= 0 & df$BA <= 1),
    all(df$macro_F1 >= 0 & df$macro_F1 <= 1))
  rows <- list()
  for (scenario in sort(unique(df$scenario))) for (metric in c('BA', 'macro_F1')) {
    d <- df[df$scenario == scenario, ]
    reps <- sort(unique(d$replication))
    values <- vapply(methods, function(method) {
      z <- d[d$method == method, ]
      stopifnot(setequal(z$replication, reps))
      z[match(reps, z$replication), metric]
    }, numeric(length(reps)))
    effects <- list(nuclear_without_smoothing = values[, 'M10'] - values[, 'M00'],
      nuclear_given_smoothing = values[, 'M11'] - values[, 'M01'],
      smoothing_without_nuclear = values[, 'M01'] - values[, 'M00'],
      smoothing_given_nuclear = values[, 'M11'] - values[, 'M10'],
      interaction = values[, 'M11'] - values[, 'M10'] - values[, 'M01'] + values[, 'M00'])
    for (effect in names(effects)) rows[[length(rows) + 1L]] <- cbind(
      data.frame(scenario, metric, effect, primary = metric == 'BA' &&
        effect %in% c('nuclear_given_smoothing', 'smoothing_given_nuclear')),
      paired_interval(100 * effects[[effect]]))
  }
  result <- do.call(rbind, rows)
  result$p_holm <- NA_real_
  result$p_holm[result$primary] <- p.adjust(result$p[result$primary], 'holm')
  result
}
