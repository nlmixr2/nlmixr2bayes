# Benchmark: the SAME 1-cmt oral ODE population model on theo_sd, written
# (A) natively in Stan with Stan's own ODE solver (ode_rk45_tol, tolerances
# matched to rxode2's 1e-8), and (B) as an nlmixr2 model through
# est="stan" (the linked rxode2/nlmixr2est likelihood with analytic
# forward-sensitivity gradients).  Identical priors, identical non-centred
# eta structure (3 diagonal blocks, half-Cauchy SDs), identical sampler
# settings.  Reports wall time (sampling only, warm compiles), worst bulk
# ESS, ESS/s, and a per-gradient-evaluation microbenchmark.
#
# FAIR SAMPLER SETTINGS: both sides run Stan defaults (adapt_delta=0.8,
# max_treedepth=10; these are also stanControl's defaults, pinned here
# explicitly so the benchmark stays matched if either side's defaults
# move).
# FAIR CORE BUDGET: both methods get the same number of cores
# (BENCH_CORES), each spent the way the method naturally parallelizes --
# native rstan runs its chains in parallel processes (its only built-in
# parallelism without reduce_sum), nlmixr2bayes runs chains sequentially
# with BENCH_CORES subject-parallel threads inside every gradient.
library(nlmixr2bayes)
BENCH_CORES <- max(2L, min(4L, rxode2::getRxThreads()))
cat(sprintf("core budget per method: %d\n", BENCH_CORES))
d <- nlmixr2data::theo_sd
ids <- unique(d$ID)
N <- length(ids)
obs <- d[d$EVID == 0, ]
dose <- vapply(ids, function(i) d$AMT[d$ID == i & d$EVID != 0][1], numeric(1))
s <- integer(N); e <- integer(N); k <- 0L
tm <- numeric(0); dv <- numeric(0)
for (i in seq_len(N)) {
  di <- obs[obs$ID == ids[i], ]
  s[i] <- k + 1L; k <- k + nrow(di); e[i] <- k
  tm <- c(tm, di$TIME); dv <- c(dv, di$DV)
}
stanData <- list(N = N, nObs = length(dv), time = tm, dv = dv,
                 dose = dose, s = s, e = e)

nativeCode <- "
functions {
  vector onecmt(real t, vector y, real ka, real ke) {
    vector[2] dy;
    dy[1] = -ka * y[1];
    dy[2] = ka * y[1] - ke * y[2];
    return dy;
  }
}
data {
  int<lower=1> N;
  int<lower=1> nObs;
  vector[nObs] time;
  vector[nObs] dv;
  vector[N] dose;
  array[N] int s;
  array[N] int e;
}
parameters {
  real tka;
  real tcl;
  real tv;
  real<lower=0> add_sd;
  real<lower=0> sd_ka;
  real<lower=0> sd_cl;
  real<lower=0> sd_v;
  matrix[N, 3] z;
}
model {
  tka ~ normal(0.45, 1);
  tcl ~ normal(1, 1);
  tv ~ normal(3.45, 1);
  add_sd ~ cauchy(0, 2.5);
  sd_ka ~ cauchy(0, 1.936);  // 2.5*sqrt(0.6)
  sd_cl ~ cauchy(0, 1.369);  // 2.5*sqrt(0.3)
  sd_v  ~ cauchy(0, 0.791);  // 2.5*sqrt(0.1)
  to_vector(z) ~ std_normal();
  for (i in 1:N) {
    real ka = exp(tka + sd_ka * z[i, 1]);
    real cl = exp(tcl + sd_cl * z[i, 2]);
    real v  = exp(tv  + sd_v  * z[i, 3]);
    vector[2] y0 = [dose[i], 0]';
    array[e[i] - s[i] + 1] vector[2] sol =
      ode_rk45_tol(onecmt, y0, -1e-6, to_array_1d(time[s[i]:e[i]]),
                   1e-8, 1e-8, 100000, ka, cl / v);
    for (j in s[i]:e[i]) {
      dv[j] ~ normal(sol[j - s[i] + 1][2] / v, add_sd);
    }
  }
}
"

odeMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- c(0, 0.7)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    prior(tka) ~ dnorm(0.45, 1)
    prior(tcl) ~ dnorm(1, 1)
    prior(tv) ~ dnorm(3.45, 1)
    prior(add.sd) ~ dcauchy(0, 2.5)
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

report <- function(tag, sf, wall) {
  su <- rstan::summary(sf)$summary
  keep <- grepl("^(tka|tcl|tv|add_sd|sd_)", rownames(su))
  worst <- min(su[keep, "n_eff"])
  cat(sprintf("%-14s wall %7.1f s | worst bulk ESS %6.0f | ESS/s %7.3f\n",
              tag, wall, worst, worst / wall))
}

## ---- native Stan ---------------------------------------------------------
cat("compiling native model...\n")
# rstan::stan_model leaks PKG_CPPFLAGS/PKG_LIBS/USE_CXX17 into the session,
# which breaks later rxode2 C compiles; snapshot and restore (nlmixr2bayes's
# own stanCompile() does this internally)
.leak <- c("PKG_CPPFLAGS", "PKG_LIBS", "USE_CXX17", "PKG_CXXFLAGS")
.old <- Sys.getenv(.leak, unset = NA_character_)
smN <- rstan::stan_model(model_code = nativeCode)
for (.v in .leak) {
  if (is.na(.old[[.v]])) Sys.unsetenv(.v) else
    do.call(Sys.setenv, stats::setNames(list(.old[[.v]]), .v))
}
iniN <- function() {
  list(tka = 0.45, tcl = 1, tv = 3.45, add_sd = 0.7,
       sd_ka = sqrt(0.6), sd_cl = sqrt(0.3), sd_v = sqrt(0.1),
       z = matrix(0, N, 3))
}
t0 <- proc.time()[["elapsed"]]
sfN <- rstan::sampling(smN, data = stanData, chains = 2, iter = 1000,
                       warmup = 500, seed = 42, refresh = 0,
                       cores = min(2L, BENCH_CORES),
                       init = list(iniN(), iniN()))
wallN <- proc.time()[["elapsed"]] - t0
report("native ode_rk45", sfN, wallN)

## ---- nlmixr2bayes (warm compile: run twice, report the second) -----------
cat("nlmixr2bayes cold run (compiles + caches)...\n")
fit1 <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  odeMod, d, est = "stan",
  control = stanControl(chains = 2L, iter = 1000L, warmup = 500L,
                        seed = 42L, cores = BENCH_CORES, calcTables = FALSE,
                        adapt_delta = 0.8, max_treedepth = 10L,
                        print = 0L,
                        # matched solver FAMILY + mode: Stan's ode_rk45
                        # is non-stiff Dormand-Prince WITH dense output;
                        # rxode2's twin is dop853 + dense=TRUE (large
                        # internal steps, interpolation at observations)
                        rxControl = rxode2::rxControl(method = "dop853",
                                                      dense = TRUE,
                                                      atol = 1e-8,
                                                      rtol = 1e-8),
                        ofv = "none", onDiagnostic = "none"))))
t0 <- proc.time()[["elapsed"]]
fit2 <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  odeMod, d, est = "stan",
  control = stanControl(chains = 2L, iter = 1000L, warmup = 500L,
                        seed = 42L, cores = BENCH_CORES, calcTables = FALSE,
                        adapt_delta = 0.8, max_treedepth = 10L,
                        print = 0L,
                        # matched solver FAMILY + mode: Stan's ode_rk45
                        # is non-stiff Dormand-Prince WITH dense output;
                        # rxode2's twin is dop853 + dense=TRUE (large
                        # internal steps, interpolation at observations)
                        rxControl = rxode2::rxControl(method = "dop853",
                                                      dense = TRUE,
                                                      atol = 1e-8,
                                                      rtol = 1e-8),
                        ofv = "none", onDiagnostic = "none"))))
wallL <- proc.time()[["elapsed"]] - t0
report("nlmixr2bayes", fit2$env$stanfit, wallL)

## ---- linked, retry ladder + FD fallback disabled -------------------------
# isolates the failure-cascade cost: failed solves reject immediately
# with -Inf instead of relaxing tolerances / retrying with FD gradients
t0 <- proc.time()[["elapsed"]]
fit3 <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  odeMod, d, est = "stan",
  control = stanControl(chains = 2L, iter = 1000L, warmup = 500L,
                        seed = 42L, cores = BENCH_CORES, calcTables = FALSE,
                        adapt_delta = 0.8, max_treedepth = 10L,
                        print = 0L,
                        maxOdeRecalc = 0L, fallbackFD = FALSE,
                        rxControl = rxode2::rxControl(method = "dop853",
                                                      dense = TRUE,
                                                      atol = 1e-8,
                                                      rtol = 1e-8),
                        ofv = "none", onDiagnostic = "none"))))
wallN0 <- proc.time()[["elapsed"]] - t0
report("linked no-retry", fit3$env$stanfit, wallN0)

## ---- per-gradient-evaluation microbenchmark ------------------------------
upN <- rstan::unconstrain_pars(sfN, list(
  tka = 0.45, tcl = 1, tv = 3.45, add_sd = 0.7, sd_ka = sqrt(0.6),
  sd_cl = sqrt(0.3), sd_v = sqrt(0.1), z = matrix(0, N, 3)))
tN <- system.time(for (i in 1:50) rstan::grad_log_prob(sfN, upN))[["elapsed"]]
cat(sprintf("native  grad eval: %6.2f ms\n", 1000 * tN / 50))
sfL <- fit2$env$stanfit
npar <- rstan::get_num_upars(sfL)
h <- stanLinkSetup(odeMod, d, thetaSens = TRUE, cores = BENCH_CORES,
                   rxControl = rxode2::rxControl(method = "dop853",
                                                 dense = TRUE,
                                                 atol = 1e-8, rtol = 1e-8))
.map <- nlmixr2bayes:::.stanMap(rxode2::rxode2(odeMod))
.Call(nlmixr2bayes:::`_nlmixr2bayes_setThetaBase`, as.double(h$initPar))
.Call(nlmixr2bayes:::`_nlmixr2bayes_setMuRef`, as.integer(.map$muRefIdx))
upL <- rep(0.1, npar)
tL <- system.time(for (i in 1:50) rstan::grad_log_prob(sfL, upL))[["elapsed"]]
cat(sprintf("linked  grad eval: %6.2f ms\n", 1000 * tL / 50))
.Call(nlmixr2bayes:::`_nlmixr2bayes_clearThetaBase`)
stanLinkFree()

## ---- posterior sanity ----------------------------------------------------
sN <- rstan::summary(sfN, pars = c("tka", "tcl", "tv", "add_sd"))$summary
sL <- rstan::summary(sfL, pars = c("tka", "tcl", "tv", "add_sd"))$summary
cat("posterior means (native vs linked):\n")
for (p in c("tka", "tcl", "tv", "add_sd")) {
  cat(sprintf("  %-7s %8.4f  vs %8.4f\n", p, sN[p, "mean"], sL[p, "mean"]))
}
