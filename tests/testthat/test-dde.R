## The case that justifies the bridge: a delay differential equation.  Stan has
## no DDE solver at all, so there is no pure-Stan oracle here -- correctness
## rests on rxode2's delayed sensitivities agreeing with central differences.
library(testthat)

ddeModel <- "
kin  <- exp(lkin)
kout <- exp(lkout)
tau  <- exp(ltau)
R(0) <- 1
d/dt(R) <- kin - kout * delay(R, tau)
"

test_that("a delay differential equation solves with delayed sensitivities", {
  m <- rxode2::rxode2(ddeModel, calcSens = c("lkin", "lkout"))
  expect_equal(unname(rxode2::rxModelVars(m)$flags[["hasDelay"]]), 1L)
  expect_true(all(
    c("rx__sens_R_BY_lkin__", "rx__sens_R_BY_lkout__") %in%
      rxode2::rxState(m)
  ))
})

test_that("the bridge returns analytic DDE sensitivities", {
  tau <- 1.5
  tms <- seq(0.5, 20, by = 1)
  ev <- rxode2::et(tms)
  h <- rxsRegister(
    ddeModel,
    events = ev,
    sens = c("lkin", "lkout"),
    output = "R",
    params = c(ltau = log(tau)),
    atol = 1e-10,
    rtol = 1e-10,
    method = "dop853"
  )
  on.exit(rxsRelease(h))

  kin <- 0.7
  kout <- 0.4
  p <- c(log(kin), log(kout))
  got <- rxsSolve(h, p)

  ## Closed form while t <= tau: the constant pre-history R(s <= 0) = 1 makes
  ## the equation R' = kin - kout, so R = 1 + (kin - kout) t and the
  ## sensitivities are exactly kin*t and -kout*t.  This is an oracle that owes
  ## nothing to finite differencing.
  early <- which(tms <= tau)
  expect_equal(got[early, 2L], kin * tms[early], tolerance = 1e-9)
  expect_equal(got[early, 3L], -kout * tms[early], tolerance = 1e-9)

  ## Past the first delay interval there is no closed form, so fall back to
  ## central differences.  The step has to stay well above the solver noise
  ## floor: differencing a solve accurate to ~1e-9 with a 1e-6 step amplifies
  ## that noise to ~1e-3.
  fd <- vapply(
    seq_along(p),
    function(j) {
      step <- 1e-4
      pp <- p
      pp[j] <- pp[j] + step
      pm <- p
      pm[j] <- pm[j] - step
      (rxsSolve(h, pp)[, 1L] - rxsSolve(h, pm)[, 1L]) / (2 * step)
    },
    numeric(attr(h, "ny"))
  )

  ## Not identically zero, i.e. the delay really is being differentiated.
  expect_gt(max(abs(got[, -1L])), 1e-3)
  expect_equal(got[, -1L], fd, tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("an estimated delay() lag is detected and warned about (nlmixr2est#1122)", {
  skip_on_cran()
  # nlmixr2est's theta sensitivity is identically zero for a parameter that
  # reaches the model only as a delay() lag, so est="nuts" cannot sample it;
  # the fit must say so rather than return a silently unmixed posterior
  .lagMod <- function(fixTau) {
    if (fixTau) {
      function() {
        ini({
          lkin <- log(0.7)
          ltau <- fix(log(1.5))
          add.sd <- 0.3
          eta.kin ~ 0.05
        })
        model({
          kin <- exp(lkin + eta.kin)
          tau <- exp(ltau)
          R(0) <- 1
          d/dt(R) <- kin - 0.4 * delay(R, tau)
          R ~ add(add.sd)
        })
      }
    } else {
      function() {
        ini({
          lkin <- log(0.7)
          ltau <- log(1.5)
          add.sd <- 0.3
          eta.kin ~ 0.05
        })
        model({
          kin <- exp(lkin + eta.kin)
          tau <- exp(ltau)
          R(0) <- 1
          d/dt(R) <- kin - 0.4 * delay(R, tau)
          R ~ add(add.sd)
        })
      }
    }
  }
  # resolved through the model's own assignments: delay(R, tau), tau <- exp(ltau)
  expect_equal(.stanDelayLagPars(rxode2::rxode2(.lagMod(FALSE))), "ltau")
  # fix()ed: nothing to sample, nothing to warn about
  expect_equal(.stanDelayLagPars(rxode2::rxode2(.lagMod(TRUE))), character(0))
  # a literal lag, and a model with no delay at all
  .litMod <- function() {
    ini({
      lkin <- log(0.7)
      add.sd <- 0.3
      eta.kin ~ 0.05
    })
    model({
      kin <- exp(lkin + eta.kin)
      R(0) <- 1
      d/dt(R) <- kin - 0.4 * delay(R, 1.5)
      R ~ add(add.sd)
    })
  }
  expect_equal(.stanDelayLagPars(rxode2::rxode2(.litMod)), character(0))
  expect_equal(.stanDelayLagPars(rxode2::rxode2(.linkMod)), character(0))
})
