# G2: the target must be a pure function of state.  NUTS assumes
# reversibility: if the density at a point depends on what was evaluated
# before it (ODE tolerance stickiness, retry state, per-subject caches), the
# sampler is silently wrong.  The ODE retry machinery is ENABLED (a hard
# point re-solves with relaxed tolerance instead of rejecting -- better than
# truncating the posterior to -Inf); purity comes from the C batch entries
# resetting all retry/sticky/tolerance state per batch and per subject, so
# the same point reproduces the same retry ladder.  This gate proves that
# bitwise, with the retry machinery live, at the exact C entries the
# compiled Stan model calls -- no Stan needed.
#
# The fixture is a DOSED ODE model on purpose: the history-dependent state
# this gate exists to catch (LSODA solver work arrays, tolerance stickiness,
# per-subject solve caches, the event-handling path) is only in play when
# the solver actually runs -- the closed-form .linkMod would pass trivially.

.odeDetMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- 0.5
    eta.cl ~ 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d / dt(central) <- -cl / v * central
    cp <- central / v
    cp ~ add(add.sd)
  })
}

.odeDetData <- function() {
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 100, EVID = 1),
          data.frame(ID = id, TIME = .tt,
                     DV = 5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5),
                     AMT = 0, EVID = 0))
  }))
}

test_that("cond batch: 500 interleaved evaluations are bitwise identical (G2)", {
  skip_on_cran()
  h <- stanLinkSetup(.odeDetMod, .odeDetData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  .linkSetTheta(h$initPar)
  set.seed(17)
  .etaA <- matrix(stats::rnorm(4, 0, 0.2), 4, 1)
  .etaB <- matrix(stats::rnorm(4, 0, 0.4), 4, 1)
  .refA <- .condBatch(.etaA)
  .refB <- .condBatch(.etaB)
  .ok <- TRUE
  for (.i in 1:250) {
    .a <- .condBatch(.etaA)
    .b <- .condBatch(.etaB)
    if (!identical(.a$value, .refA$value) ||
          !identical(.a$grad, .refA$grad) ||
          !identical(.b$value, .refB$value) ||
          !identical(.b$grad, .refB$grad)) {
      .ok <- FALSE
      break
    }
  }
  expect_true(.ok)
  # thread-count invariance: the same batch at cores=2 is bitwise
  # identical to the serial result (retry/tolerance state is subject-scoped
  # and reset per batch, so the parallel schedule cannot leak between
  # subjects)
  stanSetCores(2L)
  .a2 <- .condBatch(.etaA)
  stanSetCores(1L)
  expect_identical(.a2$value, .refA$value)
  expect_identical(.a2$grad, .refA$grad)
})

test_that("tier-2 batch: interleaved theta+eta evaluations are bitwise identical (G2)", {
  skip_on_cran()
  h <- stanLinkSetup(.odeDetMod, .odeDetData(), thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2bayes:::`_nlmixr2bayes_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setMuRef`, 1L)
  .bt <- function(theta, e) {
    .Call(nlmixr2bayes:::`_nlmixr2bayes_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  set.seed(23)
  .eta <- matrix(stats::rnorm(4, 0, 0.2), 4, 1)
  .thA <- c(1.05, 2.95, 0.55)
  .thB <- c(0.9, 3.1, 0.45)
  .refA <- .bt(.thA, .eta)
  .refB <- .bt(.thB, .eta)
  .ok <- TRUE
  for (.i in 1:100) {
    .a <- .bt(.thA, .eta)
    .b <- .bt(.thB, .eta)
    if (!identical(.a$value, .refA$value) ||
          !identical(.a$gradEta, .refA$gradEta) ||
          !identical(.a$gradTheta, .refA$gradTheta) ||
          !identical(.b$value, .refB$value) ||
          !identical(.b$gradTheta, .refB$gradTheta)) {
      .ok <- FALSE
      break
    }
  }
  expect_true(.ok)
})
