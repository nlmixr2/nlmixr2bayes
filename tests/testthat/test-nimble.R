# NIMBLE linked-likelihood backend (phase 1, R/nimbleGen.R, R/nimbleShim.R):
# the custom distribution dFoceiCondLik wraps nlmixr2bayes_cond_batch_theta
# through a compiled shim rather than re-deriving the FOCEi conditional in
# BUGS syntax.  Mirrors test-link.R/test-stan-smoke.R's structure: a
# ground-truth value comparison first (catches a sign/layout bug before any
# MCMC noise can hide it), then an end-to-end smoke test.  Uses the same
# small analytic .linkMod/.linkData fixture (helper-link.R) so nimble's
# compile step (the slow part) stays fast.

test_that(".nimbleAssertSupported refuses unsupported model shapes", {
  # mirrors .stanMap()'s actual column structure (R/stanMap.R) so a shape
  # .stanMap() never emits can't accidentally pass/fail this test for the
  # wrong reason -- theta needs >= 2 rows since a single theta is itself
  # one of the refused shapes (the length-1 distribution-call-site collapse).
  .baseMap <- list(
    nMix = 1L,
    theta = data.frame(name = c("tcl", "tv"), par = c("tcl", "tv"),
                       est = c(1, 3), lower = c(-Inf, -Inf),
                       upper = c(Inf, Inf), fix = c(FALSE, FALSE),
                       ntheta = c(1L, 2L), stringsAsFactors = FALSE),
    muRefCov = data.frame(thetaIdx = integer(0), etaIdx = integer(0),
                          covariate = character(0), name = character(0),
                          stringsAsFactors = FALSE),
    blocks = list(list(members = "eta.cl", k = 1L, fix = FALSE)))
  expect_true(nlmixr2bayes:::.nimbleAssertSupported(.baseMap))

  .mix <- .baseMap
  .mix$nMix <- 2L
  expect_error(nlmixr2bayes:::.nimbleAssertSupported(.mix), "mixtures")

  .cov <- .baseMap
  .cov$muRefCov <- data.frame(thetaIdx = 1L, etaIdx = 1L, covariate = "wt",
                              name = "tcl_wt", stringsAsFactors = FALSE)
  expect_error(nlmixr2bayes:::.nimbleAssertSupported(.cov), "covariates")

  .thetaFix <- .baseMap
  .thetaFix$theta$fix[1] <- TRUE
  expect_error(nlmixr2bayes:::.nimbleAssertSupported(.thetaFix), "fix\\(\\)ed")

  .oneTheta <- .baseMap
  .oneTheta$theta <- .baseMap$theta[1, , drop = FALSE]
  expect_error(nlmixr2bayes:::.nimbleAssertSupported(.oneTheta),
              "at least 2 population parameters")

  .corr <- .baseMap
  .corr$blocks <- list(list(members = c("eta.cl", "eta.v"), k = 2L,
                            fix = c(FALSE, FALSE)))
  expect_error(nlmixr2bayes:::.nimbleAssertSupported(.corr), "diagonal omega")

  .fix <- .baseMap
  .fix$blocks <- list(list(members = "eta.cl", k = 1L, fix = TRUE))
  expect_error(nlmixr2bayes:::.nimbleAssertSupported(.fix), "fixed-variance")
})

test_that("the NIMBLE shim compiles and caches", {
  skip_if_not_installed("nimble")
  skip_on_cran()
  .o1 <- nlmixr2bayes:::.nimbleShimCompile()
  expect_true(file.exists(.o1))
  .o2 <- nlmixr2bayes:::.nimbleShimCompile()
  expect_identical(.o1, .o2) # cache hit, not a fresh compile
})

test_that("dFoceiCondLik's compiled log-density matches the ground truth", {
  skip_if_not_installed("nimble")
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), thetaSens = FALSE, cores = 1L)
  on.exit({
    .Call(nlmixr2bayes:::`_nlmixr2bayes_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .map <- nlmixr2bayes:::.stanMap(rxode2::rxode2(.linkMod))
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setMuRef`, as.integer(.map$muRefIdx))

  set.seed(99)
  .eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.15), h$nid, h$neta)
  .thetaCases <- list(
    init = h$initPar[seq_len(h$ntheta)],
    perturbed = { .t <- h$initPar[seq_len(h$ntheta)]; .t[1] <- .t[1] + 0.2; .t }
  )

  nlmixr2bayes:::.nimbleCondLikSetup()
  # same flattening pattern as .nimbleBuildCode(): only a 1D range crosses
  # the distribution-call boundary (see the eta=double(1) design note in
  # R/nimbleGen.R -- a whole-array 2D slice collapses rank when nid or neta
  # is exactly 1, which a single-eta model like .linkMod hits).
  .code <- nimble::nimbleCode({
    for (i in 1:nid) {
      for (j in 1:neta) {
        etaFlat[(i - 1) * neta + j] <- eta[i, j]
      }
    }
    zero ~ dFoceiCondLik(theta[1:ntheta], etaFlat[1:(nid * neta)], ntheta, nid, neta)
  })
  .m <- nimble::nimbleModel(.code,
    constants = list(ntheta = h$ntheta, nid = h$nid, neta = h$neta),
    data = list(zero = 0),
    inits = list(theta = .thetaCases$init, eta = .eta),
    dimensions = list(theta = h$ntheta, eta = c(h$nid, h$neta)),
    calculate = FALSE)
  .cm <- nimble::compileNimble(.m)
  .cm$eta <- .eta
  for (.nm in names(.thetaCases)) {
    .cm$theta <- .thetaCases[[.nm]]
    .cm$calculate()
    .ref <- .Call(nlmixr2bayes:::`_nlmixr2bayes_condBatchTheta`,
                  as.double(.thetaCases[[.nm]]), .eta)
    expect_equal(.cm$logProb_zero, sum(.ref$value), tolerance = 1e-8,
                info = .nm)
  }
})

test_that("nimbleLinkedSample() runs end to end and returns finite samples", {
  skip_if_not_installed("nimble")
  skip_on_cran()
  # .linkMod: 4 subjects, 1 eta (eta.cl), 3 free thetas (tcl, tv, add.sd)
  .samples <- nimbleLinkedSample(.linkMod, .linkData(), niter = 20L,
                                 nburnin = 10L, cores = 1L, seed = 1)
  expect_equal(dim(.samples), c(10L, 7L))
  expect_true(all(c("theta[1]", "theta[2]", "theta[3]",
                    "eta[1, 1]", "eta[4, 1]") %in% colnames(.samples)))
  expect_true(all(is.finite(.samples)))
})

test_that("nimbleLinkedSample() refuses a single-subject/single-eta model", {
  skip_if_not_installed("nimble")
  skip_on_cran()
  # nid * neta == 1: the same length-1 distribution-call-site collapse as
  # the theta check, this time for etaFlat[1:(nid*neta)]. Caught before any
  # NIMBLE compilation (right after stanLinkSetup() reports h$nid), so this
  # is cheap even though it links the real problem.
  .oneSubject <- .linkData()[.linkData()$ID == 1, ]
  expect_error(
    nimbleLinkedSample(.linkMod, .oneSubject, niter = 5L, nburnin = 1L, cores = 1L),
    "at least 2 total")
})
