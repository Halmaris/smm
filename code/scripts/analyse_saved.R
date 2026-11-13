# Analyse saved CSV files without loading or fitting a model.
# Arguments: mode (quality, multiclass, binary, timing, predictions), input CSV, output directory.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop('Usage: Rscript scripts/analyse_saved.R MODE INPUT.csv OUTPUT_DIR')
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, 'R/saved_analysis.R'))
mode <- args[1]; df <- read.csv(args[2], check.names = FALSE)
out <- args[3]
if (dir.exists(out) && length(list.files(out))) stop('Output directory must be empty')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write <- function(x, name) write.csv(x, file.path(out, paste0(name, '.csv')), row.names = FALSE)
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
  rows <- list()
  for (s in unique(df$scenario)) for (method in unique(df$method)) {
    d <- df[df$scenario == s & df$method == method, ]
    stopifnot(setequal(d$replication, 1:50))
    for (metric in c('balanced_accuracy', 'macro_f1')) rows[[length(rows) + 1L]] <-
      cbind(data.frame(scenario = s, method, metric), paired_interval(100 * d[[metric]])[1:5])
  }
  write(do.call(rbind, rows), 'summary')
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
  rows <- list()
  for (dataset in unique(df$dataset)) {
    d <- df[df$dataset == dataset, ]; ref <- d[d$method == d$method[1], ]
    for (method in unique(d$method)) {
      z <- d[d$method == method, ]
      stopifnot(setequal(z$id, ref$id))
      z <- z[match(ref$id, z$id), ]; stopifnot(all(z$truth == ref$truth))
      rows[[length(rows) + 1L]] <- cbind(data.frame(dataset, method),
        classification_metrics(z$truth, z$prediction))
    }
  }
  write(do.call(rbind, rows), 'metrics')
} else stop('Unknown mode: ', mode)
