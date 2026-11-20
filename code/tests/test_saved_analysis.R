# Run from code/: Rscript --vanilla tests/test_saved_analysis.R. No model fitting.
source('R/libraries.R')
source('R/saved_analysis.R')
x <- classification_metrics(c('a', 'a', 'b', 'b'), c('a', 'a', 'a', 'a'))
stopifnot(abs(x$balanced_accuracy - .5) < 1e-12, abs(x$macro_f1 - 1/3) < 1e-12)
df <- expand.grid(scenario = 's01', replication = 1:5,
  method = c('M00', 'M10', 'M01', 'M11'), stringsAsFactors = FALSE)
df$BA <- .5 + match(df$method, c('M00', 'M10', 'M01', 'M11')) * df$replication / 100
df$macro_F1 <- df$BA
z <- multiclass_contrasts(df)
set.seed(1)
stopifnot(isTRUE(all.equal(z, multiclass_contrasts(df[sample(nrow(df)), ]))))
stopifnot(abs(z$mean[z$metric == 'BA' & z$effect == 'nuclear_given_smoothing'] - 3) < 1e-12)
stopifnot(inherits(try(multiclass_contrasts(df[-1, ]), silent = TRUE), 'try-error'))
stopifnot(inherits(try(multiclass_contrasts(rbind(df, df[1, ])), silent = TRUE), 'try-error'))
stopifnot(is.na(paired_interval(rep(1, 5))$p))
cat('Saved-result analysis checks passed.\n')
