#' Summarize paired effects with seeds averaged within each independent replication.
#' @param comparisons Data frame produced by run_additional: scenario, replication,
#'   tpe_seed, contrast, BA_difference, F1_difference; both prespecified seeds are required.
#' @return Replication means and exploratory t intervals/tests; Holm correction is
#'   applied across the provided scenarios/contrasts separately for each metric.
summarize_additional <- function(comparisons) {
  keys <- c('scenario', 'replication', 'tpe_seed', 'contrast')
  stopifnot(!anyDuplicated(comparisons[keys]), !anyNA(comparisons))
  stopifnot(all(comparisons$tpe_seed %in% c(91001L, 91002L)))
  means <- comparisons %>%
    summarise(seeds = n(),
      across(c(BA_difference, F1_difference), mean),
      .by = c(scenario, replication, contrast)) %>%
    arrange(contrast, replication, scenario)
  stopifnot(all(means$seeds == 2L), all(means$replication %in% 1:10))
  means$seeds <- NULL

  # The independent unit is a replication, after averaging its two search seeds.
  summary <- means %>%
    pivot_longer(c(BA_difference, F1_difference), names_to = 'metric') %>%
    arrange(contrast, scenario, metric, replication) %>%
    summarise(paired_interval(value), .by = c(scenario, contrast, metric)) %>%
    mutate(p_holm = p.adjust(p, 'holm'), .by = metric) %>%
    select(scenario, contrast, metric, replications = n, estimate = mean,
           lower, upper, p, p_holm) %>%
    as.data.frame()
  stopifnot(all(summary$replications == 10L))
  return(list(replication_means = as.data.frame(means), summary = summary))
}
