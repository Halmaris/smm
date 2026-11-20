# Across-dataset inference from a finite task-by-method score matrix.
# p: pairwise p-values in combn(k, 2) order; k: number of methods.
bergmann_hommel <- local({
  engine <- NULL
  adjust <- function(p, k) {
    stopifnot(k >= 2L, k <= 12L, length(p) == choose(k, 2),
      all(is.finite(p)), all(p >= 0 & p <= 1))
    if (is.null(engine)) {
      env <- new.env()
      if (Sys.info()[['sysname']] == 'Darwin') {
        old_makevars <- Sys.getenv('R_MAKEVARS_USER', unset = NA_character_)
        on.exit(if (is.na(old_makevars)) Sys.unsetenv('R_MAKEVARS_USER') else
          Sys.setenv(R_MAKEVARS_USER = old_makevars), add = TRUE)
        sdk <- system2('xcrun', '--show-sdk-path', stdout = TRUE)
        resource <- system2('/usr/bin/clang++', '-print-resource-dir', stdout = TRUE)
        makevars <- tempfile('smm-rank-Makevars-')
        on.exit(unlink(makevars), add = TRUE)
        writeLines(c('CXX=/usr/bin/clang++ -std=gnu++17',
          paste('CXXFLAGS=-O2 -nostdinc -nostdinc++',
            '-isystem', shQuote(file.path(sdk, 'usr/include/c++/v1')),
            '-isystem', shQuote(file.path(resource, 'include')),
            '-isystem', shQuote(file.path(sdk, 'usr/include'))),
          paste('CPPFLAGS=-iframework', shQuote(file.path(sdk, 'System/Library/Frameworks')))), makevars)
        Sys.setenv(R_MAKEVARS_USER = makevars)
      }
      cppFunction(code = '
        Rcpp::List bh_partitions(Rcpp::NumericVector p, int k) {
          int m = p.size();
          std::vector<int> a(k, 0), left, right;
          for (int j = 0; j<k-1; ++j) for (int l = j+1; l<k; ++l) {
            left.push_back(j); right.push_back(l);
          }
          Rcpp::NumericVector adjusted(m);
          int count = 0;
          std::function<void(int,int)> visit = [&](int pos, int maximum) {
            if (pos == k) {
              std::vector<int> members;
              double smallest = 1.0;
              for (int h = 0; h<m; ++h) if (a[left[h]] == a[right[h]]) {
                members.push_back(h); smallest = std::min(smallest, p[h]);
              }
              if (members.empty()) return;
              ++count;
              double local = std::min(1.0, members.size()*smallest);
              for (int h : members) adjusted[h] = std::max(adjusted[h], local);
              return;
            }
            for (int block = 0; block<=maximum+1; ++block) {
              a[pos] = block; visit(pos+1, std::max(maximum, block));
            }
          };
          visit(1,0);
          return Rcpp::List::create(Rcpp::Named("adjusted") = adjusted,
            Rcpp::Named("exhaustive_sets") = count);
        }', includes = '#include <functional>', env = env)
      engine <<- env$bh_partitions
    }
    result <- engine(p, as.integer(k))
    result$closed_p <- result$adjusted
    # Conservative monotone envelope used in the published example and scmamp.
    levels <- sort(unique(p))
    envelope <- cummax(vapply(levels, function(x) max(result$closed_p[p == x]), numeric(1)))
    result$adjusted <- envelope[match(p, levels)]
    return(result)
  }
  adjust
})

# Transitive near-equality classes avoid breaking ties through CSV round-off.
# x: scores (larger is better); tolerance: numerical tie threshold.
metric_ranks <- function(x, tolerance = 1e-12) {
  ord <- order(x)
  groups <- cumsum(c(TRUE, diff(x[ord]) > tolerance))
  z <- numeric(length(x)); z[ord] <- groups
  return(rank(-z, ties.method = 'average'))
}

# k: method count; pairs: index pairs; rejected: pairwise rejection flags.
nonrejection_cliques <- function(k, pairs, rejected) {
  adjacent <- matrix(TRUE, k, k)
  adjacent[cbind(pairs[1L, rejected], pairs[2L, rejected])] <- FALSE
  adjacent <- adjacent & t(adjacent)
  candidates <- lapply(seq_len(2^k - 1L), function(mask)
    which(as.logical(intToBits(mask)[seq_len(k)])))
  candidates <- Filter(function(s) length(s) > 1L && all(adjacent[s, s]), candidates)
  is_maximal <- function(s) {
    contains_s <- function(t) return(length(t) > length(s) && all(s %in% t))
    return(!any(vapply(candidates, contains_s, logical(1))))
  }
  return(Filter(is_maximal, candidates))
}

# scores: task-by-method matrix; alpha: familywise significance level.
compare_rank_matrix <- function(scores, alpha = .05) {
  stopifnot(is.matrix(scores), all(is.finite(scores)), nrow(scores) > 1L,
    ncol(scores) > 1L, !is.null(colnames(scores)))
  ranks <- t(apply(scores, 1L, metric_ranks))
  dimnames(ranks) <- dimnames(scores)
  N <- nrow(ranks); k <- ncol(ranks); mean_rank <- colMeans(ranks)
  Q0 <- 12 * N / (k * (k + 1)) * sum((mean_rank - (k + 1) / 2)^2)
  tie_sum <- sum(apply(ranks, 1L, function(x) sum(table(x)^3 - table(x))))
  tie_factor <- 1 - tie_sum / (N * (k^3 - k))
  Q <- if (tie_factor > 0) Q0 / tie_factor else 0
  F <- if (Q >= N * (k - 1)) Inf else (N - 1) * Q / (N * (k - 1) - Q)
  p_friedman <- pchisq(Q, k - 1, lower.tail = FALSE)
  p_iman <- pf(F, k - 1, (k - 1) * (N - 1), lower.tail = FALSE)
  if (tie_factor > 0) stopifnot(abs(unname(friedman.test(ranks)$statistic) - Q) < 1e-8)
  gate <- p_friedman < alpha && p_iman < alpha
  pairs <- combn(k, 2L)
  differences <- mean_rank[pairs[1L, ]] - mean_rank[pairs[2L, ]]
  # Demsar's untied variance is retained for conservative handling of ties.
  se <- sqrt(k * (k + 1) / (6 * N))
  raw_p <- 2 * pnorm(-abs(differences) / se)
  bh <- bergmann_hommel(raw_p, k)
  rejected <- gate & bh$adjusted < alpha
  ci_width <- qnorm(1 - alpha / (2 * ncol(pairs))) * se
  pairwise <- data.frame(first = colnames(scores)[pairs[1L, ]],
    second = colnames(scores)[pairs[2L, ]], rank_difference = differences,
    standard_error = se, raw_p = raw_p, exhaustive_set_p = bh$closed_p,
    bergmann_hommel_p = bh$adjusted,
    rejected = rejected, bonferroni_lower_95 = differences - ci_width,
    bonferroni_upper_95 = differences + ci_width, row.names = NULL)
  return(list(ranks = ranks, mean_rank = mean_rank, pairwise = pairwise,
    cliques = nonrejection_cliques(k, pairs, rejected),
    omnibus = data.frame(N = N, k = k, tie_factor = tie_factor,
      friedman_Q = Q, friedman_df = k - 1, friedman_p = p_friedman,
      iman_F = F, iman_df1 = k - 1, iman_df2 = (k - 1) * (N - 1),
      iman_p = p_iman, alpha = alpha, gate = gate,
      exhaustive_sets = bh$exhaustive_sets, rejected_pairs = sum(rejected))))
}
