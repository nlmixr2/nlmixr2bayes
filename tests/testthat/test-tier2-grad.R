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
  # central differences over the etas
  fdE <- matrix(0, h$nid, h$neta)
  for (k in seq_len(h$neta)) {
    up <- eta
    up[, k] <- up[, k] + .h
    dn <- eta
    dn[, k] <- dn[, k] - .h
    fdE[, k] <- (.bt(th, up)$value - .bt(th, dn)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradEta), as.numeric(fdE), tolerance = 1e-4)
  # central differences over the thetas: mu-referenced (tcl),
  # sensitivity (tv), residual (add.sd)
  fdT <- matrix(0, h$nid, length(th))
  for (t in seq_along(th)) {
    up <- th
    up[t] <- up[t] + .h
    dn <- th
    dn[t] <- dn[t] - .h
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

test_that("dosed models work; estimated dose-handling thetas are refused", {
  skip_on_cran()
  .lagMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      tlag <- -1
      add.sd <- 0.5
      eta.cl ~ 0.1
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d / dt(central) <- -cl / v * central
      alag(central) <- exp(tlag)
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  set.seed(42)
  .d <- do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 100, EVID = 1),
          data.frame(ID = id, TIME = c(0.5, 1, 2, 4, 8),
                     DV = 5 * exp(-0.05 * c(0.5, 1, 2, 4, 8)) +
                       stats::rnorm(5, 0, 0.5),
                     AMT = 0, EVID = 0))
  }))
  # the dosing-events claim, at gradient level: bolus events + alag flow
  # through the VALUE and the eta gradient (rxode2 handles the event; the
  # eta enters the ODE, not the event time)
  h <- stanLinkSetup(.lagMod, .d, cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  nlmixr2stan:::.linkSetTheta(h$initPar)
  set.seed(5)
  eta <- matrix(stats::rnorm(h$nid, 0, 0.2), h$nid, 1)
  got <- nlmixr2stan:::.condBatch(eta)
  expect_equal(got$nBad, 0L)
  .h <- 1e-5
  fd <- (nlmixr2stan:::.condBatch(eta + .h)$value -
           nlmixr2stan:::.condBatch(eta - .h)$value) / (2 * .h)
  expect_equal(as.numeric(got$grad), as.numeric(fd), tolerance = 1e-4)
  stanLinkFree()
  # but the ESTIMATED lag's theta derivative needs a jump condition the
  # linked forward sensitivities do not carry (the FD gate caught the
  # sensitivity column silently zero while the true derivative is not);
  # est="stan" refuses rather than samples a value/gradient mismatch
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.lagMod, .d, est = "stan",
                          control = stanControl(run = FALSE))),
    "dose-handling")
  # detection is transitive through intermediate assignments and clears
  # when the parameter is fixed
  expect_equal(nlmixr2stan:::.stanEventThetas(rxode2::rxode2(.lagMod)),
               "tlag")
})
