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

# ---- phase 2: ini({}) priors + bounds (R/nimblePriors.R) -------------------

test_that(".nimbleThetaPriorText translates ini() priors with correct bounds/density", {
  skip_if_not_installed("nimble")
  skip_on_cran()
  # tcl: bounded [0,5], no explicit prior() -> default normal, T(0,5)
  # tv: unbounded (logistic's own support), explicit prior() -> no truncation
  # add.sd: explicit gamma prior, positive support -> T(,0,) (a no-op given
  #   gamma's own boundary, but still exercises the support-promotion path)
  .priorMod <- function() {
    ini({
      tcl <- c(0, 1, 5)
      tv <- 3
      prior(tv) ~ dlogis(3, 1)
      add.sd <- 0.5
      prior(add.sd) ~ dgamma(2, 1)
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxode2(.priorMod)
  .map <- nlmixr2bayes:::.stanMap(.ui)
  .pri <- stanPriors(.ui)
  .thetaSdVec <- pmax(3 * abs(.map$theta$est), 1)

  .priRow <- function(nm) {
    .w <- which(.pri$pop$name == nm & .pri$pop$kind != "multivariate")
    if (length(.w) == 1L) .pri$pop[.w, ] else NULL
  }
  .idx <- function(nm) match(nm, .map$theta$name)
  .textFor <- function(nm) {
    .i <- .idx(nm)
    nlmixr2bayes:::.nimbleThetaPriorText(
      nm, .map$theta$est[.i], .thetaSdVec[.i], .priRow(nm),
      .map$theta$lower[.i], .map$theta$upper[.i], .map$theta$name)
  }
  .txtTcl <- .textFor("tcl")
  .txtTv <- .textFor("tv")
  .txtSd <- .textFor("add.sd")

  expect_identical(.txtTcl, "T(dnorm(1, sd = 3), 0, 5)")
  expect_identical(.txtTv, "dlogis(location = 3, scale = 1)")
  expect_match(.txtSd, "^T\\(dgamma\\(shape = 2, rate = 1\\), 0(\\.0+)?, \\)$")

  .code <- eval(parse(text = paste0(
    "nimble::nimbleCode({\n",
    "  theta1 ~ ", .txtTcl, "\n",
    "  theta2 ~ ", .txtTv, "\n",
    "  theta3 ~ ", .txtSd, "\n",
    "})")))
  # no nimbleExternalCall involved here, so calculate=TRUE at nimbleModel()
  # build time works fine (unlike the dFoceiCondLik-containing models above,
  # which need compileNimble() first) -- no C++ compile needed for this test.
  .m <- nimble::nimbleModel(.code, inits = list(theta1 = 2, theta2 = 3, theta3 = 1),
                            calculate = TRUE)
  .refTcl <- stats::dnorm(2, 1, 3, log = TRUE) -
    log(stats::pnorm(5, 1, 3) - stats::pnorm(0, 1, 3))
  .refTv <- stats::dlogis(3, 3, 1, log = TRUE)
  .refSd <- stats::dgamma(1, shape = 2, rate = 1, log = TRUE)
  expect_equal(.m$logProb_theta1, .refTcl, tolerance = 1e-8)
  expect_equal(.m$logProb_theta2, .refTv, tolerance = 1e-8)
  expect_equal(.m$logProb_theta3, .refSd, tolerance = 1e-8)

  # out-of-bounds tcl is rejected (truncation actually enforced, not just
  # cosmetic in the generated text)
  .m$theta1 <- -1
  .m$calculate()
  expect_true(is.infinite(.m$logProb_theta1) && .m$logProb_theta1 < 0)
})

test_that("a prior argument that references another parameter translates to its theta[] slot", {
  skip_if_not_installed("nimble")
  skip_on_cran()
  # .nimbleArgText()'s is.name(arg) branch (a bare parameter-name reference,
  # as opposed to a numeric literal) had no test coverage: the mixed-prior
  # test above only exercises literal arguments.
  .hierMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      prior(tv) ~ normal(tcl, 1)
      add.sd <- 0.5
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxode2(.hierMod)
  .map <- nlmixr2bayes:::.stanMap(.ui)
  .pri <- stanPriors(.ui)
  .w <- which(.pri$pop$name == "tv" & .pri$pop$kind != "multivariate")
  .tvIdx <- match("tv", .map$theta$name)
  .tclIdx <- match("tcl", .map$theta$name)
  .txt <- nlmixr2bayes:::.nimbleThetaPriorText(
    "tv", .map$theta$est[.tvIdx], 3, .pri$pop[.w, ],
    .map$theta$lower[.tvIdx], .map$theta$upper[.tvIdx], .map$theta$name)
  expect_identical(.txt, sprintf("dnorm(theta[%d], sd = 1)", .tclIdx))

  .code <- eval(parse(text = paste0(
    "nimble::nimbleCode({\n",
    "  theta[", .tclIdx, "] <- 2\n",
    "  theta[", .tvIdx, "] ~ ", .txt, "\n",
    "})")))
  .m <- nimble::nimbleModel(.code, dimensions = list(theta = 2),
                            inits = list(theta = c(0, 0)), calculate = FALSE)
  .m$theta[.tvIdx] <- 1.5
  .m$calculate()
  expect_equal(.m$logProb_theta[.tvIdx], stats::dnorm(1.5, 2, 1, log = TRUE),
              tolerance = 1e-8)
})

test_that(".nimbleThetaPriorText refuses a prior distribution outside the catalog", {
  .priRow <- data.frame(prior = "dweibull(1, 2)", lower = -Inf, upper = Inf,
                        stringsAsFactors = FALSE)
  expect_error(
    nlmixr2bayes:::.nimbleThetaPriorText("x", 1, 1, .priRow, -Inf, Inf, "x"),
    "not yet supported")
})

test_that(".nimbleAssertSupported refuses priors it cannot honor", {
  .baseMap <- list(
    nMix = 1L,
    theta = data.frame(name = c("tcl", "tv"), fix = c(FALSE, FALSE),
                       stringsAsFactors = FALSE),
    muRefCov = data.frame(thetaIdx = integer(0)),
    blocks = list(list(members = "eta.cl", k = 1L, fix = FALSE)))

  .mvPrior <- data.frame(name = "tcl", kind = "multivariate",
                         stringsAsFactors = FALSE)
  expect_error(
    nlmixr2bayes:::.nimbleAssertSupported(.baseMap, priPop = .mvPrior),
    "multivariate")

  .omegaPrior <- data.frame(name = "eta.cl", stringsAsFactors = FALSE)
  expect_error(
    nlmixr2bayes:::.nimbleAssertSupported(.baseMap, priOmega = .omegaPrior),
    "omega")
})

test_that(".nimbleCondLikSetup() restores .GlobalEnv bindings wiped independently of its own cache", {
  skip_if_not_installed("nimble")
  skip_on_cran()
  # dFoceiCondLik/rFoceiCondLik/condBatchExternal live in TWO places: the
  # package-internal .nimbleGenEnv cache (what gates re-doing the expensive
  # shim-compile/registerDistributions() work) and .GlobalEnv (what
  # nimbleModel()'s plain-name lookup actually needs). Something as ordinary
  # as the user's own rm(list=ls()) clears the second without touching the
  # first -- remove just the three target objects (not a blanket
  # rm(list=ls()), to avoid disturbing anything else this test run has in
  # .GlobalEnv) and confirm a second setup call restores them rather than
  # short-circuiting on the still-populated cache.
  nlmixr2bayes:::.nimbleCondLikSetup()
  expect_true(exists("dFoceiCondLik", envir = .GlobalEnv, inherits = FALSE))
  rm(list = c("dFoceiCondLik", "rFoceiCondLik", "condBatchExternal"),
     envir = .GlobalEnv)
  expect_false(exists("dFoceiCondLik", envir = .GlobalEnv, inherits = FALSE))

  nlmixr2bayes:::.nimbleCondLikSetup()
  expect_true(exists("dFoceiCondLik", envir = .GlobalEnv, inherits = FALSE))
  expect_true(exists("rFoceiCondLik", envir = .GlobalEnv, inherits = FALSE))
  expect_true(exists("condBatchExternal", envir = .GlobalEnv, inherits = FALSE))
})
