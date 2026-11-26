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
  .assert <- function(f) .stanAssert(rxode2::rxode2(f))
  expect_true(.assert(.mapBayesMod))
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
  expect_error(.assert(.oneFree), "prior(tcl)", fixed = TRUE)
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
  expect_error(.assert(.omFree), "prior(eta.cl)", fixed = TRUE)
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
  skip_if_not(.stanHasNlmApi(), "nlmixr2est lacks the nlm C API (#953)")
  expect_error(.assert(.noEta), "nothing to sample")
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
  # the centered parameterization samples the etas against the same
  # constant 2 x 2 factor
  .codeC <- suppressMessages(nlmixr2est::nlmixr2(
    .mapBayesMod,
    .mapBayesData(),
    est = "nuts",
    control = stanControl(run = FALSE, etaParam = "centered")
  ))
  .linesC <- strsplit(.codeC$code, "\n")[[1]]
  expect_true(any(grepl(
    "etaP_eta_cl_eta_v[i] ~ multi_normal_cholesky(rep_vector(0, 2), L_eta_cl_eta_v);",
    .linesC,
    fixed = TRUE
  )))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .codeC$code, allow_undefined = TRUE))
  }
})

test_that("a correlated block with only some elements fixed is refused, off-diagonal included", {
  skip_on_cran()
  .run <- function(f) {
    suppressMessages(nlmixr2est::nlmixr2(f, .mapBayesData(), est = "nuts", control = stanControl(run = FALSE)))
  }
  # a fixed covariance between free variances: the diagonal fix flags alone
  # would not see it, and the block would be sampled with the fix() ignored
  # (rxode2 refuses a prior on such a block, so the block takes the default
  # priors behind a free, prior-carrying theta)
  .fixCov <- function() {
    ini({
      tcl <- log(0.5)
      tv <- fix(log(5))
      eta.cl + eta.v ~ c(0.09, fix(0.03), 0.04)
      add.sd <- fix(0.2)
      prior(tcl) ~ dnorm(0, 2)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_error(.run(.fixCov), "mixes fixed and estimated")
  # fixed variances around a free covariance: no prior() line can be
  # suggested, so the refusal is the block's, not an empty "for example:"
  .fixVar <- function() {
    ini({
      tcl <- fix(log(0.5))
      tv <- fix(log(5))
      eta.cl + eta.v ~ c(fix(0.09), 0.03, fix(0.04))
      add.sd <- fix(0.2)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_error(.run(.fixVar), "mixes fixed and estimated")
  # one free variance: the prior advice cannot help (rxode2 refuses a prior
  # on such a block), so the block's refusal comes first
  .oneVar <- function() {
    ini({
      tcl <- fix(log(0.5))
      tv <- fix(log(5))
      eta.cl + eta.v ~ c(fix(0.09), 0.03, 0.04)
      add.sd <- fix(0.2)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_error(.run(.oneVar), "mixes fixed and estimated")
})

test_that("all-population-fixed est=\"pathfinder\" runs (no theta columns to rebuild)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .mapBayesMod,
    .mapBayesData(),
    est = "pathfinder",
    control = stanControl(seed = 22L, cores = 1L, print = 0L, vbOutputSamples = 500L, onDiagnostic = "none")
  )))
  expect_true(inherits(.fit, "nlmixr2FitData"))
  # an approximation, so only a loose agreement with the exact posterior
  expect_equal(.fit$etaObf$eta.cl, 0.070, tolerance = 0.05 / 0.070)
  expect_equal(.fit$etaObf$eta.v, 0.046, tolerance = 0.05 / 0.046)
  # the fixed omega cells are summarized (from the generated quantities),
  # and no eta or theta cell leaks in
  expect_equal(rownames(.fit$posteriorSummary), c("om.eta.cl", "cov.eta.v.eta.cl", "om.eta.v"))
  expect_equal(.fit$posteriorSummary[, "mean"], c(0.09, 0.03, 0.04))
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
  .etaSum <- rstan::summary(.fit$stanfit)$summary[c("eta[1,1]", "eta[1,2]"), ]
  expect_equal(.dx$maxRhat, max(.etaSum[, "Rhat"]))
  expect_equal(.dx$minEss, min(.etaSum[, "n_eff"]))
  expect_true(.dx$minEss > 100)
})
