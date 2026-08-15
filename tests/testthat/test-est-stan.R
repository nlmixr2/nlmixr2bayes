# The full est="stan" pipeline: compile, link, sample, and the nlmixr2 fit
# contract.  Needs the Stan toolchain; the first run compiles (~1-2 min).

test_that("nlmixr2(est='stan') returns a first-class nlmixr2 fit", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .fit <- suppressWarnings(suppressMessages(
                                            nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                                                                control = stanControl(chains = 2L, iter = 400L,
                                                                                      warmup = 200L, seed = 42L,
                                                                                      likCores = 1L,
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
  expect_true(is.character(.env$stanCode))
  expect_true(is.data.frame(.env$posteriorSummary))
  expect_true(is.list(.env$stanDiagnostics))
  # the FOCEi comparability row was made current (ofv="focei" default)
  expect_true("FOCEi" %in% rownames(.fit$objDf))
  # tables were built from the posterior point estimates
  expect_true("IPRED" %in% names(.fit))
  # posterior means should sit near the ini() estimates for this
  # well-behaved fixture (loose sanity bound, not a calibration claim)
  expect_true(abs(.fit$ui$theta[["tv"]] - 3) < 1)
})
