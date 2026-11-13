# Models and numerical conventions

## Coupled multiclass models

For training observations `X_i` and class matrices `W_k`, the fitted objective is

$$
\frac12\sum_k\|W_k\|_F^2 + \tau\sum_k\|W_k\|_* +
\gamma\sum_k Q(W_k) + \frac{\bar C}{n}\sum_i
\max\{0,\max_{k\ne y_i}(1+\langle W_k-W_{y_i},X_i\rangle)\}.
$$

There is no intercept and no class reweighting in this coupled objective.
`Cbar` is the multiplier of mean loss, not the per-observation cost.

| Variant | Nuclear penalty | Adjacency penalty |
|---|---|---|
| M00 | Fixed at zero | Fixed at zero |
| M10 | Tuned | Fixed at zero |
| M01 | Fixed at zero | Tuned |
| M11 | Tuned | Tuned |

“Tuned” includes an explicit zero branch. M00 is the coupled Crammer–Singer loss;
the variants are not OvR classifiers. Temporal `Q(W)` is the squared Frobenius norm
of second differences along columns. Spatial `Q(W)` sums squared first differences
on the two-dimensional grid. RGB uses three disconnected grids stacked along rows.

## Balanced one-versus-rest models

For each class, `fit_ovr` minimizes

$$
\frac12\|W\|_F^2 + \tau\|W\|_* + \gamma Q(W)
+\sum_i C_i\max\{0,1-y_i(\langle W,X_i\rangle+b)\},
\qquad C_i=\frac{\bar C}{2n_{y_i}}.
$$

The intercept is unpenalized. Flat weighted-hinge intercept minima use the value
of smallest absolute magnitude. Analytically certified zero coefficients are
canonicalized. Numerical prediction ties use the first class, with a
`64 * machine-epsilon` score tolerance for these canonicalized classifiers.
The coupled direct variants retain their original first-exact-maximum convention.

## Preprocessing and comparators

For each matrix row, center and scale across development observations and matrix
columns, using the sample standard deviation; constant rows use scale one.
Each CV training partition estimates its own transformation. TEST uses only the
transformation fitted on the full development sample.

`fit_method(X, y, method, parameters, channels)` uses these parameters:

| Methods | Parameters |
|---|---|
| `linear_svm`, `crammer_singer` | `Cbar` |
| `smm`, `msmm` | `Cbar`, `tau` |
| `temporal_svm`, `spatial_svm` | `Cbar`, `gamma` |
| `temporal_smm`, `spatial_smm`, `msmm_ts`, `msmm_ss` | `Cbar`, `tau`, `gamma` |
| `rbf_svm` | `Cbar`, `bandwidth_multiplier` |
| `shrinkage_lda` | `lambda` |
| `separable_lda` | `lambda_row`, `lambda_time` |
| `summary_svm` | `Cbar` |
| `functional_svm` | `Cbar`, `components` |

RBF uses balanced OvR e1071 models, `gamma = bandwidth_multiplier / median squared
nonzero training distance`, tolerance `1e-4`, shrinking enabled and 256 MB cache.
It evaluates kernels through e1071 rather than maintaining the old shared-kernel
cache. This is the quality-study objective/settings, not the historical timing
engine. LDA uses uniform class priors and covariance shrinkage toward a scaled
identity. Separable LDA estimates row and column covariance factors. Temporal
comparators use channel means/slopes or a truncated orthonormal cosine basis,
after the common training transformation, without a second standardization.

## Numerical acceptance

The single C++ backend performs ADMM linear algebra and proximal steps in FP32.
Periodic diagnostics run at iteration 1, every 25 updates and the end of a chunk.
FP32 certificates are screening values. At a candidate stop, the R callback
recomputes residuals and a feasible primal–dual certificate in FP64; the returned
solution must pass that check. This is not an FP32-only acceptance rule.

Defaults are absolute residual tolerance `1e-4`, relative residual tolerance
`1e-3`, relative gap `1e-3`, and at most 300,000 updates. Structured SMM variants
use relaxation 1.5; the named legacy margins use 1. Direct structured fits retain
the bounded, at-most-100-update FP64 precision-floor correction, included in the
total update limit. This is an internal correction, not another selectable solver.
A failure to satisfy acceptance raises an error; it is never counted as pruning
or silently accepted. There is no automatic retry or tolerance relaxation.

The compact OvR wrapper starts from zero; it does not recreate every legacy
library initialization or checkpoint. The native acceptance thresholds and
statistical objective are retained, but historical solver trajectories and timings
must be taken from the original frozen evidence, not inferred from this package.
