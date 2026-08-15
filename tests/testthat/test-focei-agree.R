# G6: agreement with est="focei" point estimates on nlmixr2data::theo_sd
# under weak priors.  Not a correctness proof (the two need not agree
# exactly even when both are right) -- a cheap catch for gross parameter
# mis-mapping on a real dataset.  Heavy (a full HMC run over an ODE-family
# model), so it runs only under NLMIXR2STAN_SLOW=TRUE; linCmt() keeps each
# gradient evaluation analytic.

.theoStan <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- c(0, 0.7)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    prior(tka) ~ dnorm(0, 10)
    prior(tcl) ~ dnorm(0, 10)
    prior(tv) ~ dnorm(0, 10)
    prior(add.sd) ~ dcauchy(0, 5)
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.theoFocei <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- c(0, 0.7)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

test_that("G6: posterior means agree with est='focei' on theo_sd", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("nlmixr2data")
  skip_if_not(nzchar(Sys.getenv("NLMIXR2STAN_SLOW")),
              "set NLMIXR2STAN_SLOW=TRUE for the G6 agreement gate")
  .d <- nlmixr2data::theo_sd
  .focei <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.theoFocei, .d, est = "focei",
                        control = nlmixr2est::foceiControl(print = 0L,
                                                           covMethod = ""))))
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.theoStan, .d, est = "stan",
                        control = stanControl(chains = 2L, iter = 1500L,
                                              warmup = 500L, seed = 7L,
                                              adapt_delta = 0.95,
                                              onDiagnostic = "none"))))
  .thF <- .focei$ui$theta
  .thS <- .fit$ui$theta
  .sd <- sqrt(diag(.fit$cov))
  for (.p in c("tka", "tcl", "tv")) {
    # weak priors: the posterior should cover the MLE
    expect_lt(abs(.thS[[.p]] - .thF[[.p]]), max(2 * .sd[[.p]], 0.1),
              label = paste0(.p, " (stan=", signif(.thS[[.p]], 4),
                             " focei=", signif(.thF[[.p]], 4), ")"))
    # log-scale thetas within ~15% on the natural scale
    expect_lt(abs(.thS[[.p]] - .thF[[.p]]), 0.15, label = paste0(.p, " log"))
  }
  # omega SDs within 40% (a weakly-informative half-Cauchy still shrinks a
  # 12-subject variance; pretending otherwise would be a fake gate)
  .omF <- sqrt(diag(.focei$omega))
  .omS <- sqrt(diag(.fit$omega))
  for (.e in names(.omF)) {
    expect_lt(abs(.omS[[.e]] - .omF[[.e]]) / .omF[[.e]], 0.4,
              label = paste0("omega SD ", .e, " (stan=", signif(.omS[[.e]], 4),
                             " focei=", signif(.omF[[.e]], 4), ")"))
  }
})
