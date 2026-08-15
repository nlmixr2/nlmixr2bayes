# G2: the target must be a pure function of state.  NUTS assumes
# reversibility: if the density at a point depends on what was evaluated
# before it (ODE tolerance stickiness, retry state, per-subject caches), the
# sampler is silently wrong.  stanControl pins the retry machinery off; this
# gate proves the batch entries are bitwise identical under heavy
# interleaved evaluation -- at the exact C entries the compiled Stan model
# calls, no Stan needed.

test_that("cond batch: 500 interleaved evaluations are bitwise identical (G2)", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
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
})

test_that("tier-2 batch: interleaved theta+eta evaluations are bitwise identical (G2)", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
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
