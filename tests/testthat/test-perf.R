# G12: performance tripwire on nlmixr2data::theo_sd with linCmt().  The
# plan's budget: wall clock <= 15 min single-core for 4x2000, worst-case
# bulk ESS/s >= 1.  This gate runs a scaled-down version (2 x 1000),
# REPORTS the numbers (they are the recorded baseline), and asserts only
# generous tripwire bounds -- machine variance is real, so a hard budget
# would fail on slow CI for reasons that are not regressions (the plan
# classifies G12 as a tripwire, not a build break).  NLMIXR2STAN_SLOW.

test_that("theo_sd linCmt() performance tripwire (G12)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("nlmixr2data")
  skip_if_not(nzchar(Sys.getenv("NLMIXR2STAN_SLOW")),
              "set NLMIXR2STAN_SLOW=TRUE for the performance gate")
  .mod <- function() {
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
      linCmt() ~ add(add.sd)
    })
  }
  .t0 <- proc.time()[["elapsed"]]
  .fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .mod, nlmixr2data::theo_sd, est = "stan",
    control = stanControl(chains = 2L, iter = 1000L, warmup = 500L,
                          seed = 42L, likCores = 2L, calcTables = FALSE,
                          onDiagnostic = "none"))))
  .wall <- proc.time()[["elapsed"]] - .t0
  .sum <- rstan::summary(.fit$env$stanfit,
                         pars = c("tka", "tcl", "tv", "add_sd"))$summary
  .worstEss <- min(.sum[, "n_eff"])
  .essPerSec <- .worstEss / .wall
  cat(sprintf("\nG12 baseline: wall %.1f s (2x1000, warm compile), worst bulk ESS %.0f, ESS/s %.3f\n",
              .wall, .worstEss, .essPerSec))
  # tripwires (generous: 3x the plan budget scaled to this run size)
  expect_lt(.wall, 45 * 60)
  expect_gt(.essPerSec, 0.05)
  expect_true(inherits(.fit, "nlmixr2FitData"))
})
