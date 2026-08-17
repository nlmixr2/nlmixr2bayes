# Smoke stage 0 (issue #1, Spec 5): the cross-package C path -- our
# R_RegisterCCallable'd nlmixr2bayes_cond_batch through nlmixr2est's pointer
# table -- validated against foceiLikRun() and central differences, with NO
# Stan in the picture.  This is the stage that catches a sign error in the
# linkage before Stan ever runs.

test_that("stanLinkSetup loads, guards re-entry, reports the handle", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  expect_equal(h$nid, 4L)
  expect_equal(h$neta, 1L)
  expect_equal(h$scale, "natural")
  expect_true(is.character(h$setupHash))
  # re-entrancy lock: a second link while one is active errors
  expect_error(stanLinkSetup(.linkMod, .linkData()), "already linked")
  expect_true(stanLinkFree())
  expect_false(stanLinkFree())
})

test_that("stage 0: value matches foceiLikRun(type='cond') exactly", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  set.seed(7)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.2), h$nid, h$neta)
  ref <- nlmixr2est::foceiLikRun(h$initPar, eta, type = "cond")
  got <- nlmixr2bayes:::.condBatch(eta)
  expect_equal(got$nBad, 0L)
  expect_equal(as.numeric(got$value), as.numeric(ref), tolerance = 1e-12)
})

test_that("stage 0: gradient matches central differences (the sign test)", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  set.seed(11)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.25), h$nid, h$neta)
  nlmixr2bayes:::.linkSetTheta(h$initPar)
  got <- nlmixr2bayes:::.condBatch(eta)
  .h <- 1e-5
  fd <- matrix(0, h$nid, h$neta)
  for (k in seq_len(h$neta)) {
    up <- eta
    up[, k] <- up[, k] + .h
    dn <- eta
    dn[, k] <- dn[, k] - .h
    fd[, k] <- (nlmixr2bayes:::.condBatch(up)$value -
                  nlmixr2bayes:::.condBatch(dn)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$grad), as.numeric(fd), tolerance = 1e-4)
  expect_false(isTRUE(all.equal(as.numeric(got$grad), as.numeric(-fd),
                                tolerance = 1e-2)))
})

test_that("stage 0: natural scale means theta is iniDf's est", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  # scale="natural": the theta part of initPar IS the ini() estimates
  expect_equal(h$initPar[seq_len(h$ntheta)], c(1, 3, 0.5))
})

test_that("fused single-solve tier-2 entry (#958) is used and consistent", {
  skip_on_cran()
  skip_if_not(nlmixr2bayes:::.stanHasCombSens(),
              "nlmixr2est lacks the combined-sensitivity build (#958)")
  .mod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      kout <- 0.2
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center - kout * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .d <- nlmixr2data::theo_sd
  .h <- stanLinkSetup(.mod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2bayes:::`_nlmixr2bayes_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  # the combined build is loaded: dims flag 0x80
  .dm <- .Call(nlmixr2bayes:::`_nlmixr2bayes_dims`)
  expect_true(bitwAnd(.dm[["flags"]], 0x80L) != 0L)
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setThetaBase`, as.double(.h$initPar))
  .map <- nlmixr2bayes:::.stanMap(rxode2::rxode2(.mod))
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setMuRef`, as.integer(.map$muRefIdx))
  set.seed(3)
  .eta <- matrix(stats::rnorm(.h$nid * .h$neta, 0, 0.2), .h$nid, .h$neta)
  .th <- .h$initPar[seq_len(.h$ntheta)]
  .g <- .Call(nlmixr2bayes:::`_nlmixr2bayes_condBatchTheta`, as.double(.th),
              .eta)
  # value + eta gradient from the fused entry are BITWISE identical to the
  # separate condBatch entry on the same combined load (the fused theta
  # pass reads the very solve the value pass produced)
  .sep <- nlmixr2bayes:::.condBatch(.eta)
  expect_identical(.g$value, .sep$value)
  expect_identical(.g$gradEta, .sep$grad)
  # assembled theta gradient (fused sens columns + mu-ref scatter) FD-agrees
  for (.jj in seq_len(.h$ntheta)) {
    .hs <- 1e-5 * max(1, abs(.th[.jj]))
    .tp <- .th; .tp[.jj] <- .tp[.jj] + .hs
    .tm <- .th; .tm[.jj] <- .tm[.jj] - .hs
    .vp <- .Call(nlmixr2bayes:::`_nlmixr2bayes_condBatchTheta`,
                 as.double(.tp), .eta)$value
    .vm <- .Call(nlmixr2bayes:::`_nlmixr2bayes_condBatchTheta`,
                 as.double(.tm), .eta)$value
    .fd <- (.vp - .vm) / (2 * .hs)
    expect_lt(max(abs((.g$gradTheta[, .jj] - .fd) / (abs(.fd) + 1e-8))),
              1e-4)
  }
})
