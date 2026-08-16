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

test_that("dosed models work; dose-handling theta gradients carry the jump", {
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
  # detection of dose-handling thetas is transitive through intermediate
  # assignments
  expect_equal(nlmixr2stan:::.stanEventThetas(rxode2::rxode2(.lagMod)),
               "tlag")
  if (!nlmixr2stan:::.stanHasEventThetaSens()) {
    # an nlmixr2est without nlmixr2/nlmixr2est#946 advertises the ESTIMATED
    # lag theta in the sensitivity index but leaves its column silently
    # zero; est="stan" refuses rather than samples a value/gradient mismatch
    expect_error(
      suppressMessages(
        nlmixr2est::nlmixr2(.lagMod, .d, est = "stan",
                            control = stanControl(run = FALSE))),
      "dose-handling")
    skip("nlmixr2est lacks the #946 event-jump theta sensitivities")
  }
  # with the #946 fix the derivative through the event is real: the full
  # tier-2 comparison -- value + d/d(eta) + d/d(theta) vs central
  # differences, INCLUDING the alag theta whose dependence runs through the
  # event time, not the ODE right-hand side
  h2 <- stanLinkSetup(.lagMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  # tcl mu-references eta.cl; tv, tlag, add.sd all carry sensitivities
  expect_equal(h2$thetaSensIdx, c(2L, 3L, 4L))
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h2$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  th <- c(1.05, 2.95, -1.05, 0.55)
  got <- .bt(th, eta)
  expect_equal(got$nBad, 0L)
  # the alag column is real, not silently zero ...
  expect_true(all(abs(got$gradTheta[, 3]) > 1e-3))
  .h <- 1e-5
  # ... and every theta column FD-agrees, the jump-carried one included
  fdT <- matrix(0, h2$nid, length(th))
  for (t in seq_along(th)) {
    up <- th
    up[t] <- up[t] + .h
    dn <- th
    dn[t] <- dn[t] - .h
    fdT[, t] <- (.bt(up, eta)$value - .bt(dn, eta)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradTheta), as.numeric(fdT), tolerance = 1e-3)
  # the eta gradient still FD-agrees under the tier-2 path
  fdE <- (.bt(th, eta + .h)$value - .bt(th, eta - .h)$value) / (2 * .h)
  expect_equal(as.numeric(got$gradEta), as.numeric(fdE), tolerance = 1e-4)
  # est="stan" now accepts the model (no dose-handling refusal)
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.lagMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
})

test_that("mu-referenced covariate coefficient gradients FD-agree", {
  skip_on_cran()
  .covMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      wt.cl <- 0.75
      add.sd <- 0.5
      eta.cl ~ 0.1
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl + wt.cl * WT + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  set.seed(7)
  .wt <- c(0.15, -0.1, 0.05, -0.2) # already centered/log-scaled
  .d <- do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = c(0.5, 1, 2, 4, 8),
               DV = 5 * exp(-0.05 * c(0.5, 1, 2, 4, 8)) +
                 stats::rnorm(5, 0, 0.5),
               WT = .wt[id], AMT = 0, EVID = 0)
  }))
  # the map resolves the coefficient to its (theta, eta) pair, and the
  # per-subject values come out in id order (the EVID=9 NA rows ignored)
  .ui <- rxode2::rxode2(.covMod)
  .map <- nlmixr2stan:::.stanMap(.ui)
  expect_equal(.map$muRefCov$thetaIdx, 3L)
  expect_equal(.map$muRefCov$etaIdx, 1L)
  .env <- new.env(parent = emptyenv())
  .env$table <- nlmixr2est::tableControl()
  nlmixr2est::.foceiPreProcessData(.d, .env, .ui, rxode2::rxControl())
  .cv <- nlmixr2stan:::.stanMuRefCovValues(.map, .env$dataSav)
  expect_equal(as.numeric(.cv$val), .wt)
  expect_false(any(.cv$timeVarying))
  # a time-varying covariate cannot factor through the scatter identity --
  # it is FLAGGED (not refused): the coefficient rides the
  # forward-sensitivity model like any other structural theta, the same
  # way the other nlmixr2est methods treat a time-varying regressor
  .dTv <- .d
  .dTv$WT[.dTv$ID == 1][3] <- 0.4
  .envTv <- new.env(parent = emptyenv())
  .envTv$table <- nlmixr2est::tableControl()
  nlmixr2est::.foceiPreProcessData(.dTv, .envTv, .ui, rxode2::rxControl())
  .cvTv <- nlmixr2stan:::.stanMuRefCovValues(.map, .envTv$dataSav)
  expect_true(.cvTv$timeVarying[1])
  # --- regime 1: the theta-sensitivity model carries the coefficient -------
  # upstream classifies wt.cl as a plain structural theta, so it gets an
  # exact forward sensitivity; est="stan" must NOT also scatter (2x bug)
  h <- stanLinkSetup(.covMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  expect_true(3L %in% h$thetaSensIdx)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, as.integer(.map$muRefIdx))
  # the est-side rule: no scatter for sensitivity-covered coefficients
  expect_length(which(!(.map$muRefCov$thetaIdx %in% h$thetaSensIdx)), 0L)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  set.seed(11)
  eta <- matrix(stats::rnorm(4, 0, 0.2), 4, 1)
  th <- c(1.05, 2.95, 0.7, 0.55)
  got <- .bt(th, eta)
  expect_equal(got$nBad, 0L)
  # the chain-rule identity d/d(wt.cl) = WT_i * d/d(eta.cl) holds for the
  # forward-sensitivity column too (same derivative, independent route)
  expect_equal(got$gradTheta[, 3], .wt * got$gradEta[, 1], tolerance = 1e-4)
  # and every theta column FD-agrees, the covariate coefficient included
  .h <- 1e-5
  fdT <- matrix(0, 4, length(th))
  for (t in seq_along(th)) {
    up <- th
    up[t] <- up[t] + .h
    dn <- th
    dn[t] <- dn[t] - .h
    fdT[, t] <- (.bt(up, eta)$value - .bt(dn, eta)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradTheta), as.numeric(fdT), tolerance = 1e-3)
  stanLinkFree()
  .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
  # --- regime 2: no theta-sensitivity model; the scatter is the source -----
  h2 <- stanLinkSetup(.covMod, .d, thetaSens = FALSE, cores = 1L)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h2$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, as.integer(.map$muRefIdx))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRefCov`,
        as.integer(.map$muRefCov$thetaIdx),
        as.integer(.map$muRefCov$etaIdx - 1L), .cv$val)
  got2 <- .bt(th, eta)
  expect_equal(got2$nBad, 0L)
  # tcl (mu-ref) and wt.cl (scatter) columns match regime 1's; tv/add.sd
  # are zero here (no sensitivity model), which is why est="stan" loads one
  # whenever such thetas are estimated
  expect_equal(got2$gradTheta[, 1], got$gradTheta[, 1], tolerance = 1e-8)
  expect_equal(got2$gradTheta[, 3], .wt * got2$gradEta[, 1],
               tolerance = 1e-12)
  expect_equal(as.numeric(got2$gradTheta[, 3]), as.numeric(fdT[, 3]),
               tolerance = 1e-3)
  stanLinkFree()
  .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
  # --- time-varying covariate: the forward-sensitivity route --------------
  h3 <- stanLinkSetup(.covMod, .dTv, thetaSens = TRUE, cores = 1L)
  expect_true(3L %in% h3$thetaSensIdx)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h3$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, as.integer(.map$muRefIdx))
  got3 <- .bt(th, eta)
  expect_equal(got3$nBad, 0L)
  .h <- 1e-5
  up <- th
  up[3] <- up[3] + .h
  dn <- th
  dn[3] <- dn[3] - .h
  fdTv <- (.bt(up, eta)$value - .bt(dn, eta)$value) / (2 * .h)
  expect_equal(as.numeric(got3$gradTheta[, 3]), as.numeric(fdTv),
               tolerance = 1e-3)
  # est="stan" accepts the model end-to-end (codegen path)
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.covMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
})
