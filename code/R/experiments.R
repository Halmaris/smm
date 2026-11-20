# Study design ----------------------------------------------------------

#' Return the fixed design of the completed multiclass ablation.
#' @return Generator settings; calling this function generates no data and fits no model.
paper_design <- function() {
  return(list(seed = 20261107L, classes = c(3L, 5L), ranks = c(2L, 8L),
    roughness = c('smooth', 'rough'), train_per_class = c(25L, 100L),
    p = 8L, T = 32L, noise_sd = .5, test_per_class = 1000L,
    replications = 50L, folds = 5L))
}

#' Enumerate the minimal additional experiments without running them.
#' @param stage Core equal-budget comparison or boundary sensitivity analysis.
#' @return One row per tuning; fixed budgets include anchors.
additional_plan <- function(stage = 'core') {
  stopifnot(stage %in% c('core', 'boundary'))
  methods <- if (stage == 'core') c('M00', 'M10', 'M01', 'M11') else c('M01', 'M11')
  df <- expand.grid(scenario = if (stage == 'core') c('s09', 's13', 's15') else 's09',
    replication = 1:10, tpe_seed = c(91001L, 91002L), method = methods,
    stringsAsFactors = FALSE)
  df$budget <- unname(c(M00 = 34L, M10 = 33L, M01 = 33L, M11 = 100L)[df$method])
  df$gamma_max <- if (stage == 'core') 1e4 else 1e5
  df$stage <- stage
  return(df)
}

# Model selection and held-out assessment -------------------------------

#' Choose the best single-penalty family using development CV alone.
#' @param tuned Named tuning objects for M00, M10 and M01.
#' @return Method name; exact score ties retain the order M00, M10, M01.
select_portfolio <- function(tuned) {
  methods <- c('M00', 'M10', 'M01')
  scores <- vapply(tuned[methods], function(x) x$selection$score, 0)
  return(methods[which.max(scores)])
}

#' Refit a frozen development choice and evaluate the held-out test observations.
#' @param inputs List with separate train and test lists containing X and factor y.
#' @param choice Tuning selection with a parameters list.
#' @param fitter Function with the fit_direct interface.
#' @param predictor Function with the predict_model interface.
#' @return TEST BA, macro-F1 and predictions; no TEST-dependent model selection.
assess_choice <- function(inputs, choice, fitter = fit_direct, predictor = predict_model) {
  model <- do.call(fitter, c(list(context = smm_context(inputs$train$X, inputs$train$y)),
    choice$parameters))
  pred <- predictor(model, inputs$test$X)
  scores <- mc_metrics(inputs$test$y, pred)
  return(list(BA = scores$BA, macro_F1 = scores$macro_F1, prediction = pred))
}

# Experiment blocks -----------------------------------------------------

#' Execute one explicitly requested additional-experiment block, sequentially.
#' @param inputs Frozen train/test inputs including original folds; attribute identity
#'   contains the source input-file checksum, as set by scripts/run_additional.R.
#' @param scenario One of s09, s13, s15; boundary stage uses only s09.
#' @param replication Original replication number, 1--10.
#' @param seed TPE seed, 91001 or 91002.
#' @param stage Core or boundary.
#' @param core Completed matching core-stage object, required for boundary comparisons.
#' @param tuner Function with the tune_direct interface, injectable for unit tests.
#' @param assessor Function with the assess_choice interface, injectable for unit tests.
#' @return Tuning traces, frozen choices, TEST comparisons and the exact refit count.
run_additional <- function(inputs, scenario, replication, seed, stage = 'core', core = NULL,
    tuner = tune_direct, assessor = assess_choice) {
  if (is.null(attr(inputs, 'identity'))) stop('Supply the identity of the frozen input file as an attribute.')
  plan <- additional_plan(stage)
  rows <- plan[plan$scenario == scenario & plan$replication == replication & plan$tpe_seed == seed, ]
  stopifnot(nrow(rows) > 0)
  if (stage == 'boundary') stopifnot(!is.null(core), core$stage == 'core',
    core$scenario == scenario, core$replication == replication, core$tpe_seed == seed,
    identical(core$input_identity, attr(inputs, 'identity')))
  tuned <- setNames(lapply(seq_len(nrow(rows)), function(i)
    tuner(inputs$train, rows$method[i], seed, rows$budget[i], rows$gamma_max[i])), rows$method)
  evaluations <- list(); refits <- 0L
  if (stage == 'core') {
    selected <- select_portfolio(tuned)
    for (method in c('M11', selected)) {
      evaluations[[method]] <- assessor(inputs, tuned[[method]]$selection)
      refits <- refits + 1L
    }
    comparison <- data.frame(scenario = scenario, replication = replication, tpe_seed = seed,
      contrast = 'M11 - CV-selected portfolio', selected = selected,
      BA_difference = evaluations$M11$BA - evaluations[[selected]]$BA,
      F1_difference = evaluations$M11$macro_F1 - evaluations[[selected]]$macro_F1)
  } else {
    comparison <- list()
    for (method in names(tuned)) {
      baseline <- core$evaluations[[method]]
      if (is.null(baseline)) {
        baseline <- assessor(inputs, core$tuned[[method]]$selection)
        refits <- refits + 1L
      }
      expanded <- assessor(inputs, tuned[[method]]$selection)
      refits <- refits + 1L
      evaluations[[method]] <- list(original = baseline, expanded = expanded)
      comparison[[method]] <- data.frame(scenario = scenario, replication = replication, tpe_seed = seed,
        contrast = paste(method, 'expanded - original gamma domain'), selected = method,
        BA_difference = expanded$BA - baseline$BA,
        F1_difference = expanded$macro_F1 - baseline$macro_F1)
    }
    comparison <- bind_rows(comparison)
  }
  return(list(stage = stage, scenario = scenario, replication = replication, tpe_seed = seed,
    input_identity = attr(inputs, 'identity'), tuned = tuned, evaluations = evaluations,
    comparison = comparison, refits = refits))
}

#' Run one block of the original four-way direct multiclass ablation.
#' @param inputs Fixed development/test observations and development fold assignments.
#' @param seed TPE seed; supply the archived method-specific seeds to reproduce a block.
#' @return Four tuning traces and held-out assessments; no execution on sourcing this file.
run_ablation <- function(inputs, seed) {
  methods <- c('M00', 'M10', 'M01', 'M11')
  stopifnot(length(seed) %in% c(1L, 4L))
  seed <- rep(seed, length.out = 4L)
  tuned <- setNames(lapply(seq_along(methods), function(i)
    tune_direct(inputs$train, methods[i], seed[i], budget = 100L)), methods)
  return(list(tuned = tuned, evaluations = lapply(tuned, function(x) assess_choice(inputs, x$selection))))
}
