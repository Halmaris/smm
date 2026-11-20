# Algebraic and workflow tests; no research campaign is executed.
args <- commandArgs(FALSE)
script <- sub('^--file=', '', args[startsWith(args, '--file=')][1])
root <- dirname(dirname(normalizePath(gsub('~+~', ' ', script, fixed = TRUE))))
source(file.path(root, 'R', 'load.R'))
m <- load_smm(root, compile = FALSE)
stopifnot(identical(m$.smm_pack(array(1:12, c(1,3,4)))$p, 1L))
set.seed(41)
X <- replicate(30, matrix(rnorm(12), 3, 4), simplify = FALSE)
y <- factor(rep(1:3, 10))
ctx <- m$smm_context(X, y)
D <- diff(diag(4), differences = 2)
H <- diag(12) * 2 + .6 * kronecker(crossprod(D), diag(3))
rhs <- matrix(rnorm(24), 12, 2)
stopifnot(max(abs(m$.smm_cached_inverse(ctx, .3, 2)(rhs) - solve(H, rhs))) < 1e-12)
for (centered in c(FALSE, TRUE)) {
  Z <- if (centered) ctx$Zc else ctx$Z
  expected <- solve(diag(12) * 2 + .6 * kronecker(crossprod(D), diag(3)) + .25 * crossprod(Z), rhs)
  stopifnot(max(abs(m$.smm_cached_solver(ctx, .3, 1, .5, centered)(rhs) - expected)) < 1e-12)
}
# Exercise the observation-space (Woodbury) branch separately.
ctx_small <- m$smm_context(X[1:6], y[1:6])
expected <- solve(H + .25 * crossprod(ctx_small$Z), rhs)
stopifnot(max(abs(m$.smm_cached_solver(ctx_small, .3, 1, .5, FALSE)(rhs) - expected)) < 1e-12)
RGB <- replicate(30, matrix(rnorm(36), 9, 4), simplify = FALSE)
rgb <- m$smm_context(RGB, y, 'spatial', channels = 3)
W <- matrix(rep(c(rep(1,3),rep(7,3),rep(-2,3)), 4),9,4)
stopifnot(rgb$roughness_fn(W) == 0, rgb$row_laplacian[3,4] == 0, rgb$row_laplacian[6,7] == 0)
Hrgb <- diag(36) * 2 + .6 * (kronecker(diag(4), rgb$row_laplacian) + kronecker(rgb$column_laplacian, diag(9)))
z <- rnorm(36)
stopifnot(max(abs(m$.smm_cached_inverse(rgb,.3,2)(z) - solve(Hrgb,z))) < 1e-12)
stopifnot(all.equal(m$.msmm_simplex(c(-1,3,2),2),c(0,1.5,.5)))
stopifnot(max(abs(m$shrinkage(diag(c(3,1)),2)$D - diag(c(1,0)))) < 1e-12)
truth <- factor(c('a','a','b','b'))
stopifnot(m$ba_upper_bound(truth,c('a','b',NA,NA)) == .75,
  m$minimum_norm_intercept(c(0,0),c(-1,1),c(.5,.5)) == 0)
plan <- m$additional_plan()
stopifnot(nrow(plan) == 240, sum(plan$budget) == 12000,
  nrow(m$additional_plan('boundary')) == 40, sum(m$additional_plan('boundary')$budget) == 2660)
# Mock fitting validates orchestration without spending a research budget.
train <- list(X=X,y=y,folds=m$stratified_folds(y,7))
inputs <- list(train=train,test=list(X=X,y=y));attr(inputs,'identity') <- 'synthetic-test'
tuner <- function(train,method,seed,budget,gamma_max) {
  list(selection=list(score=if(method=='M00') .8 else .7, parameters=list(method=method,gamma_max=gamma_max)),
    budget=budget, method=method)
}
assessor <- function(inputs,choice) list(BA=.75,macro_F1=.7)
core <- m$run_additional(inputs,'s09',1,91001,tuner=tuner,assessor=assessor)
boundary <- m$run_additional(inputs,'s09',1,91001,'boundary',core,tuner,assessor)
stopifnot(core$refits==2,core$comparison$selected=='M00',boundary$refits==3)
core$evaluations$M01 <- core$evaluations$M00
stopifnot(m$run_additional(inputs,'s09',1,91001,'boundary',core,tuner,assessor)$refits==2)
# A paired analysis has ten independent replication means, not twenty seed units.
df <- expand.grid(scenario=c('s09','s13','s15'),replication=1:10,tpe_seed=c(91001L,91002L))
df$contrast <- 'M11 - CV-selected portfolio';df$BA_difference <- sin(df$replication)/100
df$F1_difference <- cos(df$replication)/100
summary <- m$summarize_additional(df)
stopifnot(nrow(summary$summary)==6,nrow(summary$replication_means)==30,
  all(summary$summary$replications==10))
cat('PASS: scaling, temporal/spatial algebra, RGB boundaries, proximal maps, CV bound, plans, paired workflow and inference.\n')
