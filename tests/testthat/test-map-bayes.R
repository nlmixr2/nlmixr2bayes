# The individual (MAP-Bayes) posterior p(eta | y, theta, Omega, sigma) with
# every population parameter fix()ed (#22): no prior() lines are needed, a
# wholly fixed correlated omega block is a constant Cholesky factor, and the
# program samples only the etas against the eta-only linked conditional.

.mapBayesMod <- function() {
  ini({
    tcl <- fix(log(0.5))
    tv <- fix(log(5))
    eta.cl + eta.v ~ fix(0.09,
                         0.03, 0.04)
    add.sd <- fix(0.2)
    prop.sd <- fix(0.1)
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(centr) <- -cl / v * centr
    cp <- centr / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

# the reprex in #22 (and rxode2#1593)
.mapBayesData <- function() {
  data.frame(
    id = 1,
    time = c(0, 2, 11.9, 12, 24, 26),
    amt = c(100, NA, NA, 100, 100, NA),
    dv = c(NA, 17.5, 6.1, NA, NA, 13.9),
    evid = c(1, 0, 0, 1, 1, 0)
  )
}

test_that("only an all-population-fixed model with etas skips the prior requirement", {
  .ui <- function(f) rxode2::rxode2(f)
  expect_true(.stanAllPopFixed(.ui(.mapBayesMod)))
  # one estimated theta: not the individual posterior, priors still required
  .oneFree <- function() {
    ini({
      tcl <- log(0.5)
      tv <- fix(log(5))
      eta.cl ~ fix(0.09)
      add.sd <- fix(0.2)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_false(.stanAllPopFixed(.ui(.oneFree)))
  # an estimated omega: the same
  .omFree <- function() {
    ini({
      tcl <- fix(log(0.5))
      tv <- fix(log(5))
      eta.cl ~ 0.09
      add.sd <- fix(0.2)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_false(.stanAllPopFixed(.ui(.omFree)))
  # no etas at all: nothing to sample
  .noEta <- function() {
    ini({
      tcl <- fix(log(0.5))
      add.sd <- fix(0.2)
    })
    model({
      cp <- 100 * exp(-exp(tcl) * time)
      cp ~ add(add.sd)
    })
  }
  expect_false(.stanAllPopFixed(.ui(.noEta)))
})

test_that("all-population-fixed codegen: constant L, eta-only conditional, no priors", {
  skip_on_cran()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.mapBayesMod, .mapBayesData(), est = "nuts", control = stanControl(run = FALSE))
  )
  expect_s3_class(.code, "nlmixr2bayesCode")
  .lines <- strsplit(.code$code, "\n")[[1]]
  # every theta (structural and residual) dissolved into the model
  expect_true(any(grepl("vector[0] theta;", .lines, fixed = TRUE)))
  expect_true(any(grepl("vector nlmixr2_cond_all(matrix etaMat);", .lines, fixed = TRUE)))
  expect_false(any(grepl("nlmixr2_cond_all2", .lines, fixed = TRUE)))
  expect_true(any(grepl(
    "matrix[2,2] L_eta_cl_eta_v = [[0.3, 0], [0.1, 0.173205080756888]];",
    .lines,
    fixed = TRUE
  )))
  # the only sampled parameters are the standardized etas
  .par <- .lines[seq(which(.lines == "parameters {") + 1L, which(.lines == "transformed parameters {") - 2L)]
  expect_equal(trimws(.par), "matrix[N, 2] z_eta_cl_eta_v;")
  expect_true(any(grepl("logLikSubj = nlmixr2_cond_all(eta);", .lines, fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code, allow_undefined = TRUE))
  }
})

test_that("all-population-fixed est=\"nuts\" recovers the exact individual posterior (#22)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .mapBayesMod,
    .mapBayesData(),
    est = "nuts",
    control = stanControl(
      chains = 2L,
      iter = 3000L,
      warmup = 1000L,
      seed = 22L,
      chainCores = 1L,
      cores = 1L,
      print = 0L,
      onDiagnostic = "none"
    )
  )))
  expect_true(inherits(.fit, "nlmixr2FitData"))
  # the fixed population parameters come back as given
  expect_equal(
    .fit$theta,
    c(tcl = log(0.5), tv = log(5), add.sd = 0.2, prop.sd = 0.1)
  )
  expect_equal(unname(.fit$omega), matrix(c(0.09, 0.03, 0.03, 0.04), 2))
  expect_null(.fit$cov)
  # the reference posterior (an rxode2 MCMC and babelmixr2 fmeMcmc, the
  # rxode2#1593 reply): eta.cl 0.070/0.072, eta.v 0.046/0.085 (mean/sd);
  # 4000 draws put the Monte Carlo error of each mean near 0.002
  .eta <- rstan::extract(.fit$stanfit, pars = "eta")$eta
  expect_equal(mean(.eta[, 1, 1]), 0.070, tolerance = 0.01 / 0.070)
  expect_equal(stats::sd(.eta[, 1, 1]), 0.072, tolerance = 0.1)
  expect_equal(mean(.eta[, 1, 2]), 0.046, tolerance = 0.01 / 0.046)
  expect_equal(stats::sd(.eta[, 1, 2]), 0.085, tolerance = 0.1)
  # the posterior-mean etas are what the fit reports
  expect_equal(.fit$etaObf$eta.cl, mean(.eta[, 1, 1]))
  expect_equal(.fit$etaObf$eta.v, mean(.eta[, 1, 2]))
  # with no population parameter sampled, the etas are what gets diagnosed
  .dx <- .fit$env$stanDiagnostics
  expect_true(is.finite(.dx$maxRhat))
  expect_true(.dx$minEss > 100)
})
