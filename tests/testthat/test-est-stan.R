# The full est="stan" pipeline: compile, link, sample, and the nlmixr2 fit
# contract.  Needs the Stan toolchain; the first run compiles (~1-2 min).

test_that("nlmixr2(est='stan') returns a first-class nlmixr2 fit", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .fit <- suppressWarnings(suppressMessages(
                                            nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                                                                control = stanControl(chains = 2L, iter = 400L,
                                                                                      warmup = 200L, seed = 42L,
                                                                                      cores = 1L,
                                                                                      onDiagnostic = "none"))))
  expect_true(inherits(.fit, "nlmixr2FitData"))
  .env <- .fit$env
  # point estimates pushed back into the ui (the env's fullTheta/theta are
  # consumed and dropped by the output machinery) + posterior covariance
  expect_equal(names(.fit$ui$theta), c("tcl", "tv", "add.sd"))
  # the fixture's data decay 0.05 = cl/v with v ~= exp(3) implies cl ~= 1,
  # i.e. tcl ~= 0 (the prior at 1 pulls up); loose sanity bound
  expect_true(abs(.fit$ui$theta[["tcl"]]) < 1.5)
  expect_equal(dimnames(.fit$cov)[[1]], c("tcl", "tv", "add.sd"))
  expect_equal(.fit$covMethod, "stan.posterior")
  # omega is the posterior mean of L L' -- symmetric PD
  expect_true(isSymmetric(.fit$omega))
  expect_true(all(eigen(.fit$omega)$values > 0))
  # posterior-mean etas survive the automatic FOCEi objective row (H3): the
  # FOCEi EBEs are kept separately
  expect_true(is.data.frame(.fit$etaObf))
  expect_equal(names(.fit$etaObf)[1], "ID")
  expect_true(is.numeric(.fit$etaObf$OBJI))
  if (exists("etaObfFocei", .env)) {
    expect_false(isTRUE(all.equal(.fit$etaObf, get("etaObfFocei", .env))))
  }
  # Bayesian accessors ride on the fit environment
  expect_s4_class(.env$stanfit, "stanfit")
  # the documented retrieval path: fit$stanfit reaches the rstan object
  expect_s4_class(.fit$stanfit, "stanfit")
  expect_true(is.character(.fit$stanCode))
  expect_true(is.character(.env$stanCode))
  expect_true(is.data.frame(.env$posteriorSummary))
  # rownames read like the model, not its Stan translation (#1): theta/error
  # parameters use their original (dotted) names, not the Stan-mangled
  # identifiers (add_sd) or the redundant theta[] array duplicate; the solo
  # eta.cl omega block shows up as om.eta.cl (natural variance scale), not
  # the internal sd_bN/Lcorr_bN/omega_bN sampling parameterization
  .psRn <- rownames(.env$posteriorSummary)
  expect_true(all(c("tcl", "tv", "add.sd", "om.eta.cl") %in% .psRn))
  expect_false(any(grepl("^(theta\\[|add_sd$|sd_b|Lcorr_b|omega_b|L_b)",
                        .psRn)))
  expect_true(is.list(.env$stanDiagnostics))
  # the FOCEi comparability row was made current (ofv="focei" default)
  expect_true("FOCEi" %in% rownames(.fit$objDf))
  # tables were built from the posterior point estimates
  expect_true("IPRED" %in% names(.fit))
  # posterior means should sit near the ini() estimates for this
  # well-behaved fixture (loose sanity bound, not a calibration claim)
  expect_true(abs(.fit$ui$theta[["tv"]] - 3) < 1)
  # ---- G11 contract breadth ----------------------------------------------
  # print/summary run clean (the fit is presentable, not just structured)
  expect_no_error(invisible(utils::capture.output(print(.fit))))
  expect_no_error(invisible(utils::capture.output(summary(.fit$parFixedDf))))
  # $theta follows the standard nlmixr2 fit contract (named numeric; the
  # finalize step replaces the pre-finalize bounds data.frame, same as focei)
  expect_true(is.numeric(.fit$env$theta))
  expect_equal(names(.fit$env$theta), c("tcl", "tv", "add.sd"))
  # G11d: the posterior-mean etaObf survives a setOfv round-trip -- the
  # FOCEi ofv step overwrites $etaObf with EBEs unless stashed/restored
  .etaBefore <- .fit$etaObf
  .obj <- .fit$objDf["FOCEi", "OBJF"]
  suppressWarnings(suppressMessages(nlmixr2est::setOfv(.fit, "FOCEi")))
  expect_equal(.fit$etaObf, .etaBefore)
  expect_equal(.fit$objDf["FOCEi", "OBJF"], .obj, tolerance = 1e-8)
  # D7: subject-level WAIC/LOO rows and the loo objects
  if (requireNamespace("loo", quietly = TRUE)) {
    expect_s3_class(.fit$env$loo, "loo")
    expect_s3_class(.fit$env$waic, "waic")
    expect_true(any(grepl("^WAIC", rownames(.fit$objDf))))
    expect_true(any(grepl("^LOO", rownames(.fit$objDf))))
    # deviance-scale OBJF = -2 * elpd
    expect_equal(.fit$objDf[grep("^LOO", rownames(.fit$objDf)), "OBJF"],
                 -2 * .fit$env$loo$estimates["elpd_loo", "Estimate"],
                 tolerance = 1e-8)
  }
})

test_that("iteration print + parameter history over the sampler", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(nlmixr2bayes:::.stanHasIterPrint())
  # the scale.h table prints via Rprintf (stdout); nlmixr2's own progress
  # goes to the message stream -- capture both
  .out <- utils::capture.output(type = "output", {
    .msg <- utils::capture.output(type = "message", {
      .fit <- suppressWarnings(
        nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                            control = stanControl(chains = 1L, iter = 200L,
                                                  warmup = 100L, seed = 42L,
                                                  cores = 1L, print = 100L,
                                                  calcTables = FALSE,
                                                  ofv = "none",
                                                  onDiagnostic = "none")))
    })
  })
  # the familiar scale.h iteration table printed during sampling
  expect_true(any(grepl("Function Val", c(.out, .msg), fixed = TRUE)))
  # ... and the printed rows were recorded as parHistData with natural-scale
  # thetas + the ACTUAL omega variance columns (om.<eta>), objf = -2*sum(ll)
  .ph <- .fit$parHistData
  expect_true(is.data.frame(.ph))
  expect_true(nrow(.ph) > 0L)
  expect_true(all(c("iter", "type", "objf", "tcl", "tv", "add.sd",
                    "om.eta.cl", "nEval") %in% names(.ph)))
  expect_true(all(is.finite(.ph$objf[.ph$type == "Scaled"])))
  # #1: nlmixr2est's own "iter"/"#" column counts PRINTED rows, not raw
  # evaluations, so it always steps by 1 regardless of the print cadence --
  # nEval is the true per-evaluation count (the C tick shim overwrites it on
  # every raw call, ungated), and should step by roughly the requested
  # cadence (100) between consecutive printed rows, not by 1
  .nEval <- .ph$nEval[.ph$type == "Scaled"]
  expect_true(length(.nEval) > 1L)
  expect_true(median(diff(.nEval)) > 10)
  # print = 0 turns the machinery off entirely
  .fit0 <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                        control = stanControl(chains = 1L, iter = 200L,
                                              warmup = 100L, seed = 42L,
                                              cores = 1L, print = 0L,
                                              calcTables = FALSE,
                                              ofv = "none",
                                              onDiagnostic = "none"))))
  expect_true(is.null(.fit0$parHistData) ||
                nrow(.fit0$parHistData) == 0L)
})

test_that("Windows PSOCK chain parallelization runs", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(identical(.Platform$OS.type, "windows"))
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                        control = stanControl(chains = 2L, chainCores = 2L,
                                              cores = 2L, iter = 100L,
                                              warmup = 50L, seed = 42L,
                                              print = 0L, calcTables = FALSE,
                                              ofv = "none",
                                              onDiagnostic = "none"))))
  expect_s4_class(.fit$stanfit, "stanfit")
  expect_equal(.fit$stanfit@sim$chains, 2L)
  expect_equal(length(.fit$stanfit@sim$samples), 2L)
})
