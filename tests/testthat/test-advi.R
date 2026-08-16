# ADVI (rstan::vb) as an alternative inference algorithm on the SAME
# generated program and linked likelihood: stanControl(algorithm=
# "meanfield"/"fullrank").  The whole nlmixr2 fit contract holds; the
# Pareto-k diagnostic replaces Rhat/ESS/divergences.

test_that("stanControl algorithm argument validates", {
  expect_error(stanControl(algorithm = "bogus"))
  expect_equal(stanControl(algorithm = "meanfield")$algorithm, "meanfield")
  expect_equal(stanControl()$algorithm, "NUTS")
})

test_that("ADVI meanfield: a complete nlmixr2 fit through rstan::vb", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                        control = stanControl(algorithm = "meanfield",
                                              seed = 42L, cores = 1L,
                                              print = 0L,
                                              onDiagnostic = "none"))))
  expect_true(inherits(.fit, "nlmixr2FitData"))
  .env <- .fit$env
  expect_equal(.env$method, "Stan (ADVI meanfield)")
  # the vb result flows through the same posterior->fit machinery
  expect_equal(names(.fit$ui$theta), c("tcl", "tv", "add.sd"))
  expect_equal(dimnames(.fit$cov)[[1]], c("tcl", "tv", "add.sd"))
  expect_equal(.fit$covMethod, "stan.posterior")
  expect_true(isSymmetric(.fit$omega))
  expect_true(all(eigen(.fit$omega)$values > 0))
  expect_true(is.data.frame(.fit$etaObf))
  expect_equal(names(.fit$etaObf)[1], "ID")
  # loose sanity: the approximation lands in the same region as NUTS
  expect_true(abs(.fit$ui$theta[["tcl"]]) < 1.5)
  expect_true(abs(.fit$ui$theta[["tv"]] - 3) < 1)
  # ADVI diagnostics: khat replaces Rhat/ESS
  .dx <- .env$stanDiagnostics
  expect_true(is.na(.dx$maxRhat))
  expect_true(is.na(.dx$minEss))
  expect_true(is.numeric(.dx$khat))
  expect_true(grepl("ADVI meanfield", .env$extra, fixed = TRUE))
  expect_true(grepl("Pareto khat", .env$extra, fixed = TRUE))
  # WAIC/LOO rows still built from the (approximate) posterior draws
  if (requireNamespace("loo", quietly = TRUE)) {
    expect_true(any(grepl("^WAIC", rownames(.fit$objDf))))
  }
})

test_that("ADVI fullrank runs and labels itself", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                        control = stanControl(algorithm = "fullrank",
                                              seed = 42L, cores = 1L,
                                              print = 0L,
                                              calcTables = FALSE,
                                              ofv = "none",
                                              onDiagnostic = "none"))))
  # calcTables=FALSE -> core fit (no table build), still the full contract
  expect_true(inherits(.fit, "nlmixr2FitCore"))
  expect_equal(.fit$env$method, "Stan (ADVI fullrank)")
})

test_that("a high khat is surfaced through onDiagnostic", {
  skip_on_cran()
  # pure-R check on the diagnostic gate (no Stan run needed): a synthetic
  # vb-shaped object with a stored PSIS khat above 0.7 must warn
  .sf <- structure(list(), class = "fakeStanfit")
  .sim <- list(diagnostics = list(psis = list(pareto_k = 1.2)))
  .fake <- methods::setClass("fakeVb", representation(sim = "list"))(sim = .sim)
  .ctl <- list(algorithm = "meanfield", onDiagnostic = "warn")
  expect_warning(nlmixr2stan:::.stanDiagnosticsVb(.fake, .ctl),
                 "khat")
  .ctl$onDiagnostic <- "none"
  .dx <- nlmixr2stan:::.stanDiagnosticsVb(.fake, .ctl)
  expect_equal(.dx$khat, 1.2)
})
