# Summary statistics ----------------------------------------------------

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

# Paired factorial contrasts --------------------------------------------

# Summarize original four-variant ablations; df has scenario/replication/method/BA/macro_F1.
multiclass_contrasts <- function(df) {
  methods <- c('M00', 'M10', 'M01', 'M11')
  stopifnot(!anyNA(df[c('scenario', 'replication', 'method', 'BA', 'macro_F1')]),
    !anyDuplicated(df[c('scenario', 'replication', 'method')]),
    setequal(df$method, methods), all(df$BA >= 0 & df$BA <= 1),
    all(df$macro_F1 >= 0 & df$macro_F1 <= 1))
  # Align the four fitted variants by replication before taking differences.
  wide <- df |>
    tidyr::pivot_longer(c(BA, macro_F1), names_to = 'metric', values_to = 'score') |>
    tidyr::pivot_wider(id_cols = c(scenario, replication, metric),
                       names_from = method, values_from = score)
  stopifnot(!anyNA(wide[methods]))

  result <- wide |>
    dplyr::arrange(scenario, metric, replication) |>
    dplyr::transmute(scenario, replication, metric,
      nuclear_without_smoothing = M10 - M00,
      nuclear_given_smoothing = M11 - M01,
      smoothing_without_nuclear = M01 - M00,
      smoothing_given_nuclear = M11 - M10,
      interaction = M11 - M10 - M01 + M00) |>
    tidyr::pivot_longer(-c(scenario, replication, metric),
                        names_to = 'effect', values_to = 'difference') |>
    dplyr::summarise(paired_interval(100 * difference), .by = c(scenario, metric, effect)) |>
    dplyr::mutate(primary = metric == 'BA' &
      effect %in% c('nuclear_given_smoothing', 'smoothing_given_nuclear'), .after = effect) |>
    as.data.frame()

  # Adjust across the prespecified BA comparisons only.
  result$p_holm <- NA_real_
  result$p_holm[result$primary] <- p.adjust(result$p[result$primary], 'holm')
  result
}
