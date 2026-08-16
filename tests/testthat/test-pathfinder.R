# Pathfinder: the pure-R core against analytic targets (no Stan needed),
# and the est="pathfinder" / est="advi" sugar methods end to end.

test_that("pathfinder core recovers an analytic correlated Gaussian", {
  skip_on_cran()
  set.seed(7)
  .s <- matrix(c(2, 1.2, 1.2, 1.5), 2, 2)
  .si <- solve(.s)
  .mu <- c(3, -1)
  .fn <- function(x) -0.5 * as.numeric(t(x - .mu) %*% .si %*% (x - .mu))
  .gr <- function(x) as.numeric(-.si %*% (x - .mu))
  .r <- nlmixr2stan:::.pathfinderMulti(.fn, .gr, c(-5, 5), paths = 4L,
                                       jitterSd = 1, drawsPerPath = 2000L,
                                       nDraws = 4000L)
  expect_equal(.r$nPathsOk, 4L)
  # a Gaussian target is exactly representable: khat must be small
  expect_true(is.finite(.r$khat) && .r$khat < 0.5)
  expect_equal(colMeans(.r$draws), .mu, tolerance = 0.1)
  expect_equal(stats::cov(.r$draws), .s, tolerance = 0.15)
})

test_that("pathfinder core survives a -Inf region (half support)", {
  skip_on_cran()
  set.seed(11)
  # log-density of Exp(1) on x > 0: -Inf below zero exercises the
  # line-search rejection path
  .fn <- function(x) if (x[1] <= 0) -Inf else -x[1]
  .gr <- function(x) -1
  .r <- nlmixr2stan:::.pathfinderMulti(.fn, .gr, 2, paths = 2L,
                                       jitterSd = 0.5,
                                       drawsPerPath = 500L, nDraws = 500L)
  # a boundary-running target has no interior optimum; either outcome
  # (NULL or finite draws) is acceptable -- the point is no crash
  expect_true(is.null(.r) || all(is.finite(.r$draws)))
})

test_that("est='advi' and est='pathfinder' sugar dispatch + naming", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  # algorithm contradictions error instead of silently overriding
  .env <- new.env()
  .env$control <- stanControl(algorithm = "pathfinder")
  expect_error(nlmixr2Est.advi(.env), "advi")
  .env2 <- new.env()
  .env2$control <- stanControl(algorithm = "meanfield")
  expect_error(nlmixr2Est.pathfinder(.env2), "pathfinder")
})

test_that("est='pathfinder': a complete nlmixr2 fit", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("loo")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.estMod, .linkData(), est = "pathfinder",
                        control = stanControl(seed = 42L, cores = 1L,
                                              print = 0L,
                                              pathfinderPaths = 2L,
                                              vbOutputSamples = 300L,
                                              onDiagnostic = "none"))))
  expect_true(inherits(.fit, "nlmixr2FitCore"))
  .env <- .fit$env
  expect_equal(.env$est, "pathfinder")
  expect_equal(.env$method, "Stan (Pathfinder)")
  expect_true(is.numeric(.env$stanDiagnostics$khat))
  expect_true(grepl("Pathfinder", .env$extra, fixed = TRUE))
  expect_equal(names(.fit$ui$theta), c("tcl", "tv", "add.sd"))
  expect_true(isSymmetric(.fit$omega))
  expect_true(all(eigen(.fit$omega)$values > 0))
  expect_true(is.data.frame(.fit$etaObf))
  expect_true(abs(.fit$ui$theta[["tv"]] - 3) < 1)
})

test_that("est='advi': sugar records its own est name", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.estMod, .linkData(), est = "advi",
                        control = stanControl(seed = 42L, cores = 1L,
                                              print = 0L,
                                              calcTables = FALSE,
                                              ofv = "none",
                                              onDiagnostic = "none"))))
  expect_equal(.fit$env$est, "advi")
  expect_equal(.fit$env$method, "Stan (ADVI meanfield)")
})
