# Small synthetic solver tests, not paper experiments.
args <- commandArgs(FALSE)
script <- sub('^--file=', '', args[startsWith(args, '--file=')][1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root, 'R', 'load.R'))
m <- load_smm(root)
set.seed(31)
X <- replicate(30, matrix(rnorm(12),3,4),simplify=FALSE)
y <- factor(rep(1:3,10))
ctx <- m$smm_context(X,y)
for (tau in c(0,.02)) for(gamma in c(0,.1)) {
  fit <- m$fit_direct(ctx,10,tau,gamma)
  stopifnot(fit$relative_gap<=1e-3,fit$primal_residual<=fit$primal_tolerance,
    fit$dual_residual<=fit$dual_tolerance,length(m$predict_model(fit,X))==30)
}
expanded <- m$fit_direct(ctx,10,.02,1e5)
stopifnot(expanded$relative_gap<=1e-3)
fit <- m$fit_direct(ctx,1,100,0)
stopifnot(fit$iterations==0,all(fit$W==0),fit$relative_gap<1e-12)
for (structure in c('temporal','spatial')) {
  rgb <- replicate(30, matrix(rnorm(36),9,4),simplify=FALSE)
  cxt <- m$smm_context(rgb,y,structure,channels=if(structure=='spatial')3 else 1)
  fit <- m$fit_direct(cxt,10,.02,.1)
  stopifnot(fit$relative_gap<=1e-3)
  ovr <- m$fit_ovr(cxt,10,.02,.1)
  stopifnot(all(ovr$diagnostics$relative_gap<=1e-3))
}
zero <- m$fit_ovr(ctx,1,100,0)
stopifnot(all(zero$W==0),all(zero$b==0))
parameters <- list(Cbar=10,tau=.02,gamma=.1,components=2L,bandwidth_multiplier=1,
  lambda=.1,lambda_row=.1,lambda_time=.1)
for (method in names(m$method_names())) {
  model <- m$fit_method(X,y,method,parameters)
  pred <- m$predict_method(model,X)
  stopifnot(length(pred)==length(y),!anyNA(pred))
}
# Native inverse/linear factors agree with an independently assembled FP64 system.
rgb <- m$smm_context(replicate(10,matrix(rnorm(18),6,3),simplify=FALSE),
  factor(rep(1:2,5)),'spatial',channels=3)
rhs <- matrix(rnorm(36),18,2)
expected <- m$.smm_cached_inverse(rgb,.3,2)(rhs)
actual <- m$smm_apply_inverse_fp32(rhs,6L,3L,.3,2,TRUE,3L)
stopifnot(max(abs(expected-actual))<1e-5)
factor <- m$.smm_fp32_factor(rgb,.3,1,.5,FALSE)
expected <- m$.smm_cached_solver(rgb,.3,1,.5,FALSE)(rhs)
actual <- m$smm_apply_factor_fp32(rgb$Z,rhs,factor,6L,3L)
stopifnot(max(abs(expected-actual))<1e-5)
cat('PASS: four ablations, exact zeros, temporal/RGB fits, all named methods and native algebra.\n')
