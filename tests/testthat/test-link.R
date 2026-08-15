# Smoke stage 0 (issue #1, Spec 5): the cross-package C path -- our
# R_RegisterCCallable'd nlmixr2stan_cond_batch through nlmixr2est's pointer
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
  got <- nlmixr2stan:::.condBatch(eta)
  expect_equal(got$nBad, 0L)
  expect_equal(as.numeric(got$value), as.numeric(ref), tolerance = 1e-12)
})

test_that("stage 0: gradient matches central differences (the sign test)", {
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  set.seed(11)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.25), h$nid, h$neta)
  nlmixr2stan:::.linkSetTheta(h$initPar)
  got <- nlmixr2stan:::.condBatch(eta)
  .h <- 1e-5
  fd <- matrix(0, h$nid, h$neta)
  for (k in seq_len(h$neta)) {
    up <- eta; up[, k] <- up[, k] + .h
    dn <- eta; dn[, k] <- dn[, k] - .h
    fd[, k] <- (nlmixr2stan:::.condBatch(up)$value -
                  nlmixr2stan:::.condBatch(dn)$value) / (2 * .h)
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
