# G13: centred and non-centred eta parameterizations target the same
# posterior by construction; sampling the same model both ways must agree.
# A failure here almost always means a missing Jacobian in the centred path
# or a wrong Cholesky factor in the non-centred one.
#
# The data are simulated FROM the model so omega is identified: with
# model-mismatched data (the .linkData fixture) the omega posterior is
# prior-dominated with a half-Cauchy tail, and tail-driven MOMENTS disagree
# by Monte Carlo alone (observed: posterior sd 96 vs 3 at agreeing bulk).
# Omega is therefore compared on quantiles, thetas on mean/sd.

test_that("centred == non-centred posterior (G13)", {
  skip_if_not_installed("rstan")
  skip_on_cran()
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
  .tt <- c(0.5, 1, 2, 4, 8)
  .d <- rxode2::rxWithSeed(11, {
    do.call(rbind, lapply(1:8, function(id) {
      .cl <- exp(1 + stats::rnorm(1, 0, 0.3))
      .f <- 100 / exp(3) * exp(-.cl / exp(3) * .tt)
      data.frame(ID = id, TIME = .tt,
                 DV = .f + stats::rnorm(length(.tt), 0, 0.5),
                 AMT = 0, EVID = 0)
    }))
  })
  .fitWith <- function(etaParam) {
    suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      .mod, .d, est = "stan",
      control = stanControl(chains = 2L, iter = 4000L, warmup = 1000L,
                            seed = 73L, etaParam = etaParam,
                            adapt_delta = 0.95, cores = 1L,
                            calcTables = FALSE, onDiagnostic = "none"))))
  }
  .nc <- .fitWith("noncentered")
  .ce <- .fitWith("centered")
  # the centred program really is centred (samples eta directly)
  expect_true(any(grepl("etaP_eta_cl", .ce$env$stanCode, fixed = TRUE)))
  expect_false(any(grepl("etaP_eta_cl", .nc$env$stanCode, fixed = TRUE)))
  .sN <- rstan::summary(.nc$env$stanfit,
                        pars = c("tcl", "tv", "add_sd"))$summary
  .sC <- rstan::summary(.ce$env$stanfit,
                        pars = c("tcl", "tv", "add_sd"))$summary
  for (.p in c("tcl", "tv", "add_sd")) {
    .tol <- 3 * sqrt(.sN[.p, "se_mean"]^2 + .sC[.p, "se_mean"]^2)
    expect_lt(abs(.sN[.p, "mean"] - .sC[.p, "mean"]), .tol,
              label = paste0(.p, " mean |", signif(.sN[.p, "mean"], 4), " - ",
                             signif(.sC[.p, "mean"], 4), "|"))
    expect_lt(abs(.sN[.p, "sd"] - .sC[.p, "sd"]) / .sN[.p, "sd"], 0.1,
              label = paste0(.p, " sd rel diff"))
  }
  # omega: quantile agreement on the SD scale (robust to the prior tail)
  .omN <- sqrt(rstan::extract(.nc$env$stanfit, pars = "omegaOut")$omegaOut[, 1, 1])
  .omC <- sqrt(rstan::extract(.ce$env$stanfit, pars = "omegaOut")$omegaOut[, 1, 1])
  .qN <- stats::quantile(.omN, c(0.25, 0.5, 0.75))
  .qC <- stats::quantile(.omC, c(0.25, 0.5, 0.75))
  expect_true(all(abs(.qN - .qC) / .qN < 0.15),
              label = paste0("omega sd quartiles nc(",
                             paste(signif(.qN, 3), collapse = ", "),
                             ") vs c(",
                             paste(signif(.qC, 3), collapse = ", "), ")"))
})
