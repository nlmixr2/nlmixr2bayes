# G9 (cross-session): bitwise-identical draws for a fixed seed across FRESH
# R sessions.  Catches reproducibility leaks the in-session gate cannot:
# package-load-order state, lazy-initialized globals, cache-restored model
# objects behaving differently from freshly compiled ones.
# NLMIXR2STAN_SLOW: two callr sessions, each links + samples (the compile
# itself is warm through the on-disk cache).

test_that("fixed seed: bitwise-identical draws across fresh sessions (G9)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("callr")
  skip_if_not(nzchar(Sys.getenv("NLMIXR2STAN_SLOW")),
              "set NLMIXR2STAN_SLOW=TRUE for the cross-session gate")
  .run <- function() {
    callr::r(function() {
      .mod <- function() {
        ini({
          tcl <- 1
          tv <- 3
          add.sd <- c(0, 0.5)
          eta.cl ~ 0.1
          prior(tcl) ~ dnorm(1, 2)
          prior(add.sd) ~ dcauchy(0, 2.5)
        })
        model({
          cl <- exp(tcl + eta.cl)
          v <- exp(tv)
          cp <- 100 / v * exp(-cl / v * time)
          cp ~ add(add.sd)
        })
      }
      set.seed(42)
      .tt <- c(0.5, 1, 2, 4, 8)
      .d <- do.call(rbind, lapply(1:4, function(id) {
        data.frame(ID = id, TIME = .tt,
                   DV = 5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5),
                   AMT = 0, EVID = 0)
      }))
      .fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
        .mod, .d, est = "stan",
        control = nlmixr2stan::stanControl(chains = 1L, iter = 400L,
                                           warmup = 200L, seed = 42L,
                                           likCores = 1L,
                                           calcTables = FALSE,
                                           onDiagnostic = "none"))))
      rstan::extract(.fit$env$stanfit, permuted = FALSE)
    })
  }
  expect_identical(.run(), .run())
})
