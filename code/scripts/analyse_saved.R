# Analyse saved CSV files without loading or fitting a model.
# Arguments: mode (quality, multiclass, binary, timing, predictions), input CSV, output directory.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop('Usage: Rscript scripts/analyse_saved.R MODE INPUT.csv OUTPUT_DIR')
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(gsub('~+~', ' ', script, fixed = TRUE))))
source(file.path(root, 'R/libraries.R'))
source(file.path(root, 'R/saved_analysis.R'))
mode <- args[1]; df <- read.csv(args[2], check.names = FALSE)
out <- args[3]
if (dir.exists(out) && length(list.files(out))) stop('Output directory must be empty')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write <- function(x, name) return(write.csv(x, file.path(out, paste0(name, '.csv')), row.names = FALSE))
if (mode == 'quality') {
  stopifnot(names(df)[1] == 'dataset', !anyDuplicated(df$dataset))
  scores <- as.matrix(df[-1]); rownames(scores) <- df$dataset
  stopifnot(all(scores >= 0 & scores <= 1))
  source(file.path(root, 'R/rank_comparisons.R'))
  z <- compare_rank_matrix(scores)
  write(data.frame(method = colnames(scores), mean = colMeans(scores),
    mean_rank = z$mean_rank), 'summary')
  write(z$pairwise, 'pairwise'); write(z$omnibus, 'omnibus')
} else if (mode == 'multiclass') {
  stopifnot(setequal(df$scenario, sprintf('s%02d', 1:16)), nrow(df) == 3200L,
    setequal(df$replication, 1:50))
  write(multiclass_contrasts(df), 'contrasts')
} else if (mode == 'binary') {
  stopifnot(!anyDuplicated(df[c('scenario', 'replication', 'method')]),
    all(is.finite(df$balanced_accuracy)), all(is.finite(df$macro_f1)),
    all(df$balanced_accuracy >= 0 & df$balanced_accuracy <= 1),
    all(df$macro_f1 >= 0 & df$macro_f1 <= 1))
  counts <- count(df, scenario, method)
  stopifnot(!anyNA(df[c('scenario', 'replication', 'method')]),
    all(df$replication %in% 1:50), all(counts$n == 50L),
    nrow(counts) == n_distinct(df$scenario) * n_distinct(df$method))
  summary <- df %>%
    pivot_longer(c(balanced_accuracy, macro_f1), names_to = 'metric') %>%
    summarise(paired_interval(100 * value), .by = c(scenario, method, metric)) %>%
    select(-p)
  write(summary, 'summary')
} else if (mode == 'timing') {
  # Input consists of measured repeats only; warm-ups are kept in a separate file.
  stopifnot(!anyNA(df[c('dataset', 'method', 'repeat', 'seconds')]),
    !anyDuplicated(df[c('dataset', 'method', 'repeat')]),
    all(is.finite(df$seconds)), all(df$seconds > 0))
  write(aggregate(seconds ~ dataset + method, df, median), 'medians')
  write(aggregate(seconds ~ dataset + method, df, length), 'repeat_counts')
} else if (mode == 'predictions') {
  stopifnot(!anyNA(df[c('dataset', 'method', 'id', 'truth', 'prediction')]),
    !anyDuplicated(df[c('dataset', 'method', 'id')]))
  # Each method must score the same labelled observations within a dataset.
  pairs <- df %>%
    mutate(expected = n_distinct(method), .by = dataset) %>%
    summarise(valid = n() == first(expected) &
                       n_distinct(truth) == 1L, .by = c(dataset, id))
  stopifnot(all(pairs$valid))
  metrics <- df %>%
    summarise(classification_metrics(truth, prediction), .by = c(dataset, method))
  write(metrics, 'metrics')
} else stop('Unknown mode: ', mode)
