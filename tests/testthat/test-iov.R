# First-class IOV: nlmixr2est's preprocessing hook (.uiApplyIov, enabled by
# attr(nlmixr2Est.stan, "iov")) expands `iov.x ~ v | OCC` into per-occasion
# FIXED unit-variance etas scaled by an estimated sd theta in the model
# code.  The generator renders a fixed eta as a CONSTANT Cholesky factor
# (z still sampled -- a random effect with known variance), which also
# gives plain fixed-variance etas for free.  Fixed etas with modeled
# correlations are refused.

.iovMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- c(0, 0.5)
    prior(tcl) ~ dnorm(1, 2)
    prior(add.sd) ~ dcauchy(0, 2.5)
    eta.cl ~ 0.1
    iov.cl ~ 0.02 | OCC
  })
  model({
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.iovData <- function() {
  do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = c(1, 4, 8), DV = c(2, 1.4, 0.8),
                     AMT = 0, EVID = 0, OCC = 1),
          data.frame(ID = id, TIME = c(25, 28, 32), DV = c(2.1, 1.5, 0.9),
                     AMT = 0, EVID = 0, OCC = 2))
  }))
}

test_that("IOV expands to fixed unit-variance etas + an sd theta", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iovMod, .iovData(), est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2bayesCode")
  .lines <- strsplit(.code$code, "\n")[[1]]
  # the IOV magnitude is an ordinary theta ...
  expect_true(any(grepl("theta[4] = iov_cl;", .lines, fixed = TRUE)))
  # ... and the occasion etas are fixed unit-variance blocks: constant L,
  # no omega parameter, z still sampled (named from the occasion eta, not
  # an opaque "_b2"/"_b3")
  expect_true(any(grepl("matrix[1,1] L_rx_iov_cl_1 = [[1]];", .lines,
                        fixed = TRUE)))
  expect_true(any(grepl("matrix[1,1] L_rx_iov_cl_2 = [[1]];", .lines,
                        fixed = TRUE)))
  expect_true(any(grepl("to_vector(z_rx_iov_cl_1) ~ std_normal();", .lines,
                        fixed = TRUE)))
  expect_false(any(grepl("omega_rx_iov_cl_1|sd_rx_iov_cl_1", .lines)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("the assembled IOV target's gradient FD-agrees end to end", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  skip_if_not_installed("rstan")
  skip_if_not_installed("numDeriv")
  .d <- .iovData()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iovMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .sm <- stanCompile(.code$code)
  # link the REWRITTEN ui (the hook's expansion) exactly as est="stan" does
  h <- stanLinkSetup(.code$ui, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(`_nlmixr2bayes_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(`_nlmixr2bayes_setThetaBase`, as.double(h$initPar))
  .Call(`_nlmixr2bayes_setMuRef`,
        as.integer(.code$map$muRefIdx))
  .sf <- rstan::sampling(.sm, data = .code$data, chains = 1, iter = 2,
                         warmup = 1, refresh = 0, cores = 1,
                         show_messages = FALSE,
                         init = list(list(tcl = 1, tv = 3, add_sd = 0.5,
                                          iov_cl = sqrt(0.02),
                                          sd_eta_cl = sqrt(0.1),
                                          z_eta_cl = matrix(0, 4, 1),
                                          z_rx_iov_cl_1 = matrix(0, 4, 1),
                                          z_rx_iov_cl_2 = matrix(0, 4, 1))))
  set.seed(9)
  .pt <- list(tcl = 1.05, tv = 2.95, add_sd = 0.55, iov_cl = 0.2,
              sd_eta_cl = 0.3,
              z_eta_cl = matrix(stats::rnorm(4, 0, 0.4), 4, 1),
              z_rx_iov_cl_1 = matrix(stats::rnorm(4, 0, 0.4), 4, 1),
              z_rx_iov_cl_2 = matrix(stats::rnorm(4, 0, 0.4), 4, 1))
  .up <- rstan::unconstrain_pars(.sf, .pt)
  .gA <- rstan::grad_log_prob(.sf, .up)
  attributes(.gA) <- NULL
  .gN <- numDeriv::grad(function(u) rstan::log_prob(.sf, u), .up,
                        method = "Richardson")
  # every block: thetas (incl. the IOV sd through the forward
  # sensitivities), the id-level eta, and BOTH occasion-eta blocks
  expect_equal(.gA, .gN, tolerance = 1e-4)
})

test_that("plain fixed-variance etas work; fixed correlated blocks refuse", {
  skip_on_cran()
  .fixEta <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      eta.cl ~ fix(0.04)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .d <- .iovData()[, c("ID", "TIME", "DV", "AMT", "EVID")]
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.fixEta, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("matrix[1,1] L_eta_cl = [[0.2]];", .lines,
                        fixed = TRUE)))
  expect_false(any(grepl("sd_eta_cl", .lines)))
  # the centered parameterization composes with a fixed block: eta sampled
  # directly against the CONSTANT L
  .codeC <- suppressMessages(
    nlmixr2est::nlmixr2(.fixEta, .d, est = "stan",
                        control = stanControl(run = FALSE,
                                              etaParam = "centered")))
  .linesC <- strsplit(.codeC$code, "\n")[[1]]
  expect_true(any(grepl("matrix[1,1] L_eta_cl = [[0.2]];", .linesC,
                        fixed = TRUE)))
  expect_true(any(grepl("etaP_eta_cl[i] ~ multi_normal_cholesky", .linesC,
                        fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .codeC$code,
                               allow_undefined = TRUE))
  }
  # fixed etas with modeled correlations: out of scope, must error
  .fixBlk <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      eta.cl + eta.v ~ fix(0.04,
                           0.01, 0.04)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.fixBlk, .d, est = "stan",
                          control = stanControl(run = FALSE))),
    "correlations")
})

test_that("IOV models are refused while the upstream sens column is wrong (#952)", {
  skip_on_cran()
  skip_if(.stanHasIovSens())
  # with the iov attribute gated off, the | occ etas survive to the
  # id-only assertion and are refused rather than sampled with a
  # ~14% wrong gradient on the IOV magnitude theta
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.iovMod, .iovData(), est = "stan",
                          control = stanControl(run = FALSE))))
})

# issue #15: nlmixr2est's IOV rewrite builds the magnitude theta by copying
# the model's FIRST theta row, so without the R/stanIov.R repair the
# magnitude tracks theta #1's prior and the declared prior(iov.x) -- which
# lives on the deleted eta row -- never reaches the program.

.iovPriorMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- c(0, 0.5)
    prior(tcl) ~ dnorm(-7.25, 0.125)  # deliberately distinctive
    prior(add.sd) ~ dcauchy(0, 2.5)
    prior(iov.cl) ~ dcauchy(0, 1)
    eta.cl ~ 0.1
    iov.cl ~ 0.02 | OCC
  })
  model({
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

test_that("a declared prior(iov.x) lands on the IOV magnitude theta", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iovPriorMod, .iovData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("iov_cl ~ cauchy(0, 1);", .lines, fixed = TRUE)))
  # the first theta's prior stays on the first theta and nowhere else
  expect_true(any(grepl("tcl ~ normal(-7.25, 0.125);", .lines, fixed = TRUE)))
  expect_false(any(grepl("iov_cl ~ normal", .lines, fixed = TRUE)))
  # the magnitude IS an SD: lower-bounded, so cauchy() is the half-Cauchy
  # idiom (literal arguments -> constraint declared, T[,] omitted)
  expect_true(any(grepl("real<lower=0> iov_cl;", .lines, fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("an undeclared IOV magnitude gets the omega-SD default, not a leak", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  # .iovMod declares prior(tcl) ~ dnorm(1, 2) but nothing on iov.cl
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iovMod, .iovData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_false(any(grepl("iov_cl ~ normal", .lines, fixed = TRUE)))
  # diagOmegaSdPrior at 2.5 * sd0, exactly as an omega diagonal gets
  expect_true(any(grepl("^\\s*iov_cl ~ cauchy\\(0, 0\\.3535533905932",
                        .lines)))
  # ... announced, and NOT reported as flat
  expect_true(any(grepl("default prior iov_cl ~ cauchy", .code$notes)))
  expect_false(any(grepl("flat[^\n]*iov", .code$notes)))
  # the template stays honest for a non-default diagOmegaSdPrior too
  .code2 <- suppressMessages(
    nlmixr2est::nlmixr2(.iovMod, .iovData(), est = "stan",
                        control = stanControl(run = FALSE,
                                              diagOmegaSdPrior = "normal(0, %s)")))
  expect_true(any(grepl("iov_cl ~ normal(0, ", strsplit(.code2$code, "\n")[[1]],
                        fixed = TRUE)))
})

test_that("a matrix prior on an IOV magnitude is refused, not mis-emitted", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  .mvIov <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      prior(iov.cl) ~ invWishart(4)
      eta.cl ~ 0.1
      iov.cl ~ 0.02 | OCC
    })
    model({
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  # which side refuses depends on the nlmixr2est: once the IOV rewrite
  # carries the declared prior onto the magnitude theta itself, rxode2's own
  # validation gets there first ("applies to a covariance matrix, but
  # 'iov.cl' is a population estimate"); before that, the repair here is
  # what sees it.  Either way it is refused, never mis-emitted.
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.mvIov, .iovData(), est = "stan",
                          control = stanControl(run = FALSE))),
    "must be univariate|covariance matrix")
})

test_that("prior(iov.x) alone satisfies the est=\"stan\" prior requirement", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  # with the leak, the declared prior vanished and theta #1 had none, so
  # this model was refused as "declares none"
  .onlyIov <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      prior(iov.cl) ~ dcauchy(0, 1)
      eta.cl ~ 0.1
      iov.cl ~ 0.02 | OCC
    })
    model({
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.onlyIov, .iovData(), est = "stan",
                        control = stanControl(run = FALSE)))
  expect_true(any(grepl("iov_cl ~ cauchy(0, 1);",
                        strsplit(.code$code, "\n")[[1]], fixed = TRUE)))
})

test_that("the IOV prior capture is checked, not trusted", {
  skip_on_cran()
  # a rewritten ui whose magnitude thetas the capture never saw must error
  # rather than fall back to whatever the rewrite left behind
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iovPriorMod, .iovData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .old <- .stanIovEnv$priors
  on.exit(.stanIovEnv$priors <- .old, add = TRUE)
  .stanIovEnv$priors <- character(0)
  expect_error(stanPriors(.code$ui), "did not see this model")
})

test_that("the capture hook is registered ahead of nlmixr2est's IOV rewrite", {
  skip_on_cran()
  .hooks <- nlmixr2est::preProcessHooks()
  expect_true(".stanCaptureIovPriors" %in% .hooks)
  expect_true(which(.hooks == ".stanCaptureIovPriors") <
                which(.hooks == ".uiApplyIov"))
})

test_that("several IOV parameters over several occasion variables", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  .twoIov <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(-7.25, 0.125)
      prior(add.sd) ~ dcauchy(0, 2.5)
      prior(iov.v) ~ dcauchy(0, 3)   # declared on the SECOND one only
      eta.cl ~ 0.1
      iov.cl ~ 0.02 | OCC
      iov.v ~ 0.03 | VIS
    })
    model({
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + iov.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .d <- .iovData()
  .d$VIS <- .d$OCC
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.twoIov, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  # each magnitude gets ITS OWN prior: the declared one, the default one
  expect_true(any(grepl("iov_v ~ cauchy(0, 3);", .lines, fixed = TRUE)))
  expect_true(any(grepl("^\\s*iov_cl ~ cauchy\\(0, 0\\.3535533905932", .lines)))
  expect_false(any(grepl("iov_cl ~ normal|iov_v ~ normal", .lines)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("a parameter-dependent IOV prior is truncated at the SD's zero", {
  skip_on_cran()
  skip_if_not(.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  # the magnitude theta is lower-bounded at 0, so a prior whose arguments
  # reference another parameter needs the T[0, ] normalizer (rule 2 of
  # R/stanPriors.R) rather than the half-Cauchy literal idiom
  .depIov <- function() {
    ini({
      tcl <- 1
      tv <- 3
      tauScale <- c(0, 1)
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      prior(tauScale) ~ dcauchy(0, 1)
      prior(iov.cl) ~ dcauchy(0, tauScale)
      eta.cl ~ 0.1
      iov.cl ~ 0.02 | OCC
    })
    model({
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv) * tauScale / tauScale
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.depIov, .iovData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("iov_cl ~ cauchy(0, tauScale) T[0, ];", .lines,
                        fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("a non-sd iovXform expansion is refused, not read as an SD", {
  skip_on_cran()
  # est="stan" leaves iovXform at nlmixr2est's "sd" default, which is what
  # makes the magnitude theta an SD; the backTransform records the xform
  # that was actually used, so a "logvar" expansion must refuse rather than
  # emit an SD-scale prior on a log-variance
  .iniDf <- data.frame(ntheta = c(1, 2), neta1 = c(NA_real_, NA_real_),
                       name = c("tcl", "iov.cl"),
                       backTransform = c(NA_character_, "nlmixr2iovLogvarCv"),
                       stringsAsFactors = FALSE)
  expect_error(.stanIovMagnitude(.iniDf), "standard deviation")
  .iniDf$backTransform[2] <- "nlmixr2iovSdCv"
  expect_equal(.stanIovMagnitude(.iniDf), "iov.cl")
})
