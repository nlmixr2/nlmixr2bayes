# Simulate-and-recover (the "is it right?" gate that exercises the whole
# pipeline against a KNOWN truth): data generated from known theta/Omega with
# a correlated 2x2 omega block, fit with est="stan", and the generative
# values must land inside the posterior.  This is also the only end-to-end
# exercise of the LKJ + SD omega path (the default fixture has a single eta).
# Heavy (a fresh program shape compiles + a real HMC run): NLMIXR2STAN_SLOW.

test_that("simulate-and-recover: correlated omega rho inside the 90% CrI", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("mvtnorm")
  skip_if_not(nzchar(Sys.getenv("NLMIXR2STAN_SLOW")),
              "set NLMIXR2STAN_SLOW=TRUE for the simulate-and-recover gate")
  .mod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      eta.cl + eta.v ~ c(0.1,
                         0.06, 0.1)
      prior(tcl) ~ dnorm(1, 2)
      prior(tv) ~ dnorm(3, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  # ---- simulate from the truth: rho = 0.6 ---------------------------------
  .nid <- 40L
  .tt <- c(0.25, 0.5, 1, 2, 4, 8)
  .omTrue <- matrix(c(0.1, 0.06, 0.06, 0.1), 2, 2)
  .truth <- list(tcl = 1, tv = 3, addSd = 0.3, rho = 0.6)
  .d <- rxode2::rxWithSeed(1234, {
    .eta <- mvtnorm::rmvnorm(.nid, sigma = .omTrue)
    do.call(rbind, lapply(seq_len(.nid), function(.i) {
      .cl <- exp(.truth$tcl + .eta[.i, 1])
      .v <- exp(.truth$tv + .eta[.i, 2])
      .cp <- 100 / .v * exp(-.cl / .v * .tt)
      data.frame(ID = .i, TIME = .tt,
                 DV = .cp + stats::rnorm(length(.tt), 0, .truth$addSd),
                 AMT = 0, EVID = 0)
    }))
  })
  # ---- fit ----------------------------------------------------------------
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(
      .mod, .d, est = "stan",
      control = stanControl(chains = 2L, iter = 1500L, warmup = 500L,
                            seed = 77L, adapt_delta = 0.95, cores = 2L,
                            onDiagnostic = "none"))))
  .sf <- .fit$env$stanfit
  # ---- recovery -----------------------------------------------------------
  .om <- rstan::extract(.sf, pars = "omegaOut")$omegaOut  # draws x 2 x 2
  .rho <- .om[, 1, 2] / sqrt(.om[, 1, 1] * .om[, 2, 2])
  .ci <- stats::quantile(.rho, c(0.05, 0.95))
  expect_true(.ci[1] < .truth$rho && .truth$rho < .ci[2],
              label = paste0("rho 90% CrI [", signif(.ci[1], 3), ", ",
                             signif(.ci[2], 3), "] covers 0.6"))
  # thetas recovered within 3 posterior SDs
  .sum <- rstan::summary(.sf, pars = c("tcl", "tv", "add_sd"))$summary
  expect_lt(abs(.sum["tcl", "mean"] - .truth$tcl), 3 * .sum["tcl", "sd"])
  expect_lt(abs(.sum["tv", "mean"] - .truth$tv), 3 * .sum["tv", "sd"])
  expect_lt(abs(.sum["add_sd", "mean"] - .truth$addSd),
            3 * .sum["add_sd", "sd"] + 0.05)
  # the fit's omega is the posterior-mean covariance with the right names
  expect_equal(dimnames(.fit$omega)[[1]], c("eta.cl", "eta.v"))
  expect_true(all(eigen(.fit$omega)$values > 0))
})
