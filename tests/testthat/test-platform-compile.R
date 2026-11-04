## A real Stan compile + link + sample on every CI platform.
##
## The macOS/Windows CI jobs run only the CRAN-visible subset
## (NLMIXR2BAYES_TEST_CRAN_ONLY forces NOT_CRAN=false), so every other test
## that compiles Stan is skipped there.  That hid, for weeks, that no Stan
## compile had ever been shown to work on those platforms -- it surfaced only
## as a vignette failure.  These two tests are the smallest proof that each
## subsystem compiles, links against rxode2/nlmixr2est and runs: they are
## gated on NLMIXR2BAYES_TEST_PLATFORM_COMPILE (set by the macOS/Windows jobs
## in R-CMD-check.yaml) or NOT_CRAN, and never run on CRAN itself.
library(testthat)

skipUnlessPlatformCompile <- function() {
  skip_if_not_installed("rstan")
  if (
    !isTRUE(as.logical(Sys.getenv("NLMIXR2BAYES_TEST_PLATFORM_COMPILE", "false"))) &&
      !identical(Sys.getenv("NOT_CRAN"), "true")
  ) {
    skip("set NLMIXR2BAYES_TEST_PLATFORM_COMPILE=true (or NOT_CRAN=true) to compile Stan")
  }
}

test_that("likelihood-level linking: est='nuts' compiles, links and samples", {
  skipUnlessPlatformCompile()
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(
      .estMod,
      .linkData(),
      est = "nuts",
      control = stanControl(
        chains = 1L,
        iter = 100L,
        warmup = 50L,
        seed = 42L,
        cores = 1L,
        print = 0L,
        calcTables = FALSE,
        ofv = "none",
        onDiagnostic = "none"
      )
    )
  ))
  expect_true(inherits(.fit, "nlmixr2FitCore"))
  expect_equal(names(.fit$ui$theta), c("tcl", "tv", "add.sd"))
  expect_true(all(is.finite(.fit$ui$theta)))
})

test_that("rxstan bridge: a hand-coded .stan program compiles, links and samples", {
  skipUnlessPlatformCompile()
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by = 4))
  h <- rxsRegister(
    "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
",
    events = ev,
    sens = c("lka", "lcl", "lv"),
    output = "center"
  )
  on.exit(rxsRelease(h))
  .truth <- c(lka = log(1.1), lcl = log(4), lv = log(30))
  .cp <- rxsSolve(h, unname(.truth))[, 1] / exp(.truth[["lv"]])
  sm <- rxsStanModel(
    system.file("stan", "pk_1cmt_oral.stan", package = "nlmixr2bayes"),
    modelName = "platformcheck"
  )
  .fit <- suppressWarnings(rstan::sampling(
    sm,
    data = list(nObs = length(.cp), handle = h, cpObs = .cp),
    chains = 1,
    iter = 100,
    warmup = 50,
    seed = 42,
    refresh = 0,
    init = function(chain_id = 1) list(lka = 0.1, lcl = 1.4, lv = 3.4, sigma = 0.2)
  ))
  .draws <- as.matrix(.fit, pars = c("lka", "lcl", "lv"))
  expect_equal(dim(.draws), c(50L, 3L))
  expect_true(all(is.finite(.draws)))
})
