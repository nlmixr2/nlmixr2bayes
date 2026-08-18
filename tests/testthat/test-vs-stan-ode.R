## The strongest independent check available: solve the same model with Stan's
## own ODE integrator and autodiff sensitivities, and require the log density
## and gradient to agree.  Nothing is shared between the two routes except the
## data, so agreement here rules out a whole class of wiring mistakes that
## finite differencing against ourselves cannot see.
library(testthat)

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
"
truth <- c(lka = log(1.1), lcl = log(4), lv = log(30))

test_that("the bridge agrees with Stan's own ODE solver", {
  skipUnlessStan(nlmixr2 = FALSE)

  nObs <- 20L
  ts <- seq(0.5, 24, length.out = nObs)
  relTol <- 1e-10
  absTol <- 1e-10

  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, ts)
  h <- rxsRegister(pkModel, events = ev, sens = names(truth),
                   output = "center", atol = absTol, rtol = relTol)
  on.exit(rxsRelease(h))

  cp <- rxsSolve(h, truth)[, 1L] / exp(truth[["lv"]])
  set.seed(9)
  cpObs <- exp(log(cp) + stats::rnorm(nObs, 0, 0.1))

  smBridge <- rxsStanModel(system.file("stan", "pk_1cmt_oral.stan",
                                       package = "nlmixr2bayes"),
                           modelName = "xcheck_bridge")
  fitBridge <- rstan::sampling(smBridge, chains = 0,
                               data = list(nObs = nObs,
                                           handle = as.integer(unclass(h)),
                                           cpObs = cpObs))

  smStan <- rstan::stan_model(system.file("stan", "pk_1cmt_oral_stanode.stan",
                                          package = "nlmixr2bayes"),
                              model_name = "xcheck_stanode")
  fitStan <- rstan::sampling(smStan, chains = 0,
                             data = list(nObs = nObs, ts = ts, cpObs = cpObs,
                                         amt = 100, relTol = relTol,
                                         absTol = absTol, maxSteps = 100000L,
                                         stiff = 1L))

  ## Several points, not just the truth: a constant offset would hide here.
  for (delta in list(c(0, 0, 0, 0), c(0.3, -0.2, 0.15, 0.4),
                     c(-0.5, 0.4, -0.3, -0.2))) {
    u <- c(truth[["lka"]], truth[["lcl"]], truth[["lv"]], log(0.1)) + delta

    lpB <- rstan::log_prob(fitBridge, u, adjust_transform = FALSE)
    lpS <- rstan::log_prob(fitStan, u, adjust_transform = FALSE)
    expect_equal(lpB, lpS, tolerance = 1e-6)

    gB <- as.numeric(rstan::grad_log_prob(fitBridge, u, adjust_transform = FALSE))
    gS <- as.numeric(rstan::grad_log_prob(fitStan, u, adjust_transform = FALSE))
    expect_equal(gB, gS, tolerance = 1e-5)
  }
})
