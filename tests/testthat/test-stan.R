## End-to-end: a Stan program whose ODE solve is done by rxode2, with rxode2's
## analytic sensitivities carried onto Stan's reverse tape.
##
## Slow (compiles a Stan model), so it is skipped unless RXSTAN_STAN_TESTS is set.
library(testthat)

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
cp <- center / v
"

truth <- c(lka = log(1.1), lcl = log(4), lv = log(30))

makeFixture <- function(tag = NULL) {
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by = 1.5))

  h <- rxsRegister(pkModel, events = ev,
                   sens = c("lka", "lcl", "lv"), output = "center",
                   atol = 1e-10, rtol = 1e-10)

  ## Seed immediately before the draw: rxSolve() advances R's RNG but the fast
  ## path does not, so seeding earlier would make the data depend on which
  ## solve route ran.
  cp <- rxsSolve(h, truth)[, 1L] / exp(truth[["lv"]])
  set.seed(42)
  cpObs <- exp(log(cp) + stats::rnorm(length(cp), 0, 0.1))

  name <- if (is.null(tag)) "rxstan_pk" else "rxstan_pk_fd"
  sm <- stanModelFor(system.file("stan", "pk_1cmt_oral.stan", package = "nlmixr2bayes"),
                     name, tag = tag)

  fit <- rstan::sampling(sm, chains = 0,
                         data = list(nObs = length(cpObs),
                                     handle = as.integer(unclass(h)),
                                     cpObs = cpObs))
  list(handle = h, model = sm, fit = fit, cpObs = cpObs)
}

test_that("Stan gradients through rxode2 match finite differences", {
  skipUnlessStan(nlmixr2 = FALSE)

  f <- makeFixture()
  on.exit(rxsRelease(f$handle))

  ## unconstrained: lka, lcl, lv, log(sigma)
  upars <- c(truth[["lka"]], truth[["lcl"]], truth[["lv"]], log(0.1))
  chk <- rxsCheckGradient(f$fit, upars)

  print(chk)
  expect_true(all(chk$relDiff < 1e-6),
              info = paste(utils::capture.output(print(chk)), collapse = "\n"))
})

test_that("the analytic policy and the finite-difference policy agree", {
  skipUnlessStan(nlmixr2 = FALSE)

  fa <- makeFixture(tag = NULL)
  ff <- makeFixture(tag = "::rxstan::finite_diff_tag")
  on.exit({ rxsRelease(fa$handle); rxsRelease(ff$handle) })

  upars <- c(truth[["lka"]], truth[["lcl"]], truth[["lv"]], log(0.1))

  ga <- rstan::grad_log_prob(fa$fit, upars, adjust_transform = FALSE)
  gf <- rstan::grad_log_prob(ff$fit, upars, adjust_transform = FALSE)

  ## Same log density either way; only the derivative route differs.
  expect_equal(rstan::log_prob(fa$fit, upars, adjust_transform = FALSE),
               rstan::log_prob(ff$fit, upars, adjust_transform = FALSE),
               tolerance = 1e-10)
  expect_equal(as.numeric(ga), as.numeric(gf), tolerance = 1e-5)
})

test_that("the bridge samples and recovers the simulated truth", {
  skipUnlessStan(nlmixr2 = FALSE)

  f <- makeFixture()
  on.exit(rxsRelease(f$handle))

  ## Started at the truth on purpose.  A one-compartment oral profile is
  ## genuinely bimodal -- the flip-flop mode swaps ka and cl/v and fits just as
  ## well -- so a diffuse start tests PK identifiability, not this bridge.
  init <- list(list(lka = truth[["lka"]], lcl = truth[["lcl"]],
                    lv = truth[["lv"]], sigma = 0.1))

  s <- rstan::sampling(f$model, chains = 1, iter = 1000, warmup = 500,
                       seed = 11, refresh = 0, init = init,
                       data = list(nObs = length(f$cpObs),
                                   handle = as.integer(unclass(f$handle)),
                                   cpObs = f$cpObs))
  post <- as.matrix(s, pars = c("lka", "lcl", "lv"))
  ci <- apply(post, 2, stats::quantile, probs = c(0.025, 0.975))
  print(ci)

  ## The chain has to actually move, or "covers the truth" means nothing.
  expect_gt(length(unique(post[, "lka"])), nrow(post) / 2)

  for (nm in c("lka", "lcl", "lv")) {
    expect_gte(truth[[nm]], ci[1, nm])
    expect_lte(truth[[nm]], ci[2, nm])
  }
})
