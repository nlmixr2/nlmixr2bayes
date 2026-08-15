# The tier-2 C path (theta sampled too): value + d/d(eta) + d/d(theta)
# against central differences.  This validates OUR gradient assembly -- in
# particular the mu-reference identity (d/dtheta_p = d/deta_k for a theta
# mu-referencing eta k) that the C shim adds on top of nlmixr2est's forward
# sensitivities.

test_that("cond_batch_theta: value + both gradients FD-agree", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  expect_true(h$thetaSens)
  # tcl mu-references eta.cl (theta index 1); tv + add.sd carry sensitivities
  expect_equal(h$thetaSensIdx, c(2L, 3L))
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)

  set.seed(31)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.2), h$nid, h$neta)
  th <- c(1.05, 2.95, 0.55)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  got <- .bt(th, eta)
  expect_equal(got$nBad, 0L)
  .h <- 1e-5
  # d/d(eta)
  fdE <- matrix(0, h$nid, h$neta)
  for (k in seq_len(h$neta)) {
    up <- eta; up[, k] <- up[, k] + .h
    dn <- eta; dn[, k] <- dn[, k] - .h
    fdE[, k] <- (.bt(th, up)$value - .bt(th, dn)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradEta), as.numeric(fdE), tolerance = 1e-4)
  # d/d(theta): mu-referenced (tcl), sensitivity (tv), residual (add.sd)
  fdT <- matrix(0, h$nid, length(th))
  for (t in seq_along(th)) {
    up <- th; up[t] <- up[t] + .h
    dn <- th; dn[t] <- dn[t] - .h
    fdT[, t] <- (.bt(up, eta)$value - .bt(dn, eta)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradTheta), as.numeric(fdT), tolerance = 1e-3)
  # and the mu-referenced column really is the eta gradient
  expect_equal(got$gradTheta[, 1], got$gradEta[, 1], tolerance = 1e-10)
})

test_that("cond_batch_theta refuses without the tier-2 state installed", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
  expect_error(
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`,
          as.double(c(1, 3, 0.5)), matrix(0, h$nid, h$neta)),
    "status -101")
})
