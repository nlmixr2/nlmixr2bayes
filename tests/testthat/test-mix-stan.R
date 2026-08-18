# Finite mixtures (nlmixr2/nlmixr2est#955 + the Stan Users Guide
# log-sum-exp marginalization).  Component-specific etas ride the blessed
# component-major batch layout; the eta priors factor OUT of the
# marginalization (each keeps its ordinary z ~ std_normal), the mixing
# weight sits INSIDE as a pure Stan autodiff parameter (the
# component-conditional is p-free, FD-verified upstream), and the
# membership posteriors land in generated quantities.

.mixMod <- function() {
  ini({
    tcl1 <- 1
    tcl2 <- 2
    tv <- 3
    p1 <- 0.3
    add.sd <- 0.5
    eta.cl ~ 0.1
    prior(tcl1) ~ dnorm(1, 2)
    prior(tcl2) ~ dnorm(2, 2)
    prior(p1) ~ dbeta(2, 2)
    prior(add.sd) ~ dcauchy(0, 2.5)
  })
  model({
    cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.mixData <- function() {
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = .tt,
               DV = 5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5),
               AMT = 0, EVID = 0)
  }))
}

.stanHasMixApi <- function() !identical(.Call(nlmixr2bayes:::`_nlmixr2bayes_nMix`), -2L)

test_that("mixture codegen: component-major etas + log_sum_exp + membership", {
  skip_on_cran()
  skip_if_not(.stanHasMixApi(),
              "nlmixr2est lacks the blessed mixture layout (#955)")
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.mixMod, .mixData(), est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2bayesCode")
  .lines <- strsplit(.code$code, "\n")[[1]]
  # the mixing probability is constrained to (0,1) regardless of iniDf
  expect_true(any(grepl("real<lower=0,upper=1> p1;", .lines, fixed = TRUE)))
  # component-specific etas: 2*N rows through the same block
  expect_true(any(grepl("matrix[2*N, 1] z_eta_cl;", .lines, fixed = TRUE)))
  expect_true(any(grepl("matrix[2*N, 1] eta;", .lines, fixed = TRUE)))
  # the marginalization and the membership posteriors
  expect_true(any(grepl("log_sum_exp(log(p1) + llCond[i], log1p(-p1) + llCond[N + i])",
                        .lines, fixed = TRUE)))
  expect_true(any(grepl("matrix[N, 2] mixProbOut;", .lines, fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
  # K > 2 refuses with a clear message
  .mix3 <- function() {
    ini({
      tcl1 <- 1
      tcl2 <- 2
      tcl3 <- 3
      tv <- 3
      p1 <- 0.3
      p2 <- 0.3
      add.sd <- 0.5
      eta.cl ~ 0.1
      prior(tcl1) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl), p2,
                exp(tcl3 + eta.cl))
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.mix3, .mixData(), est = "stan",
                          control = stanControl(run = FALSE))),
    "2-component")
})

test_that("mixture tier-2 shim: component-conditional value + gradients FD-agree", {
  skip_on_cran()
  skip_if_not(.stanHasMixApi(),
              "nlmixr2est lacks the blessed mixture layout (#955)")
  .d <- .mixData()
  h <- stanLinkSetup(.mixMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2bayes:::`_nlmixr2bayes_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  expect_identical(.Call(nlmixr2bayes:::`_nlmixr2bayes_nMix`), 2L)
  .map <- nlmixr2bayes:::.stanMap(rxode2::rxode2(.mixMod))
  expect_equal(.map$nMix, 2L)
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2bayes:::`_nlmixr2bayes_setMuRef`, as.integer(.map$muRefIdx))
  .bt <- function(theta, e) {
    .Call(nlmixr2bayes:::`_nlmixr2bayes_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  set.seed(7)
  eta <- matrix(stats::rnorm(8, 0, 0.2), 8, 1) # component-major 2 x 4
  .th <- h$initPar[seq_len(5)]
  got <- .bt(.th, eta)
  expect_equal(got$nBad, 0L)
  .h <- 1e-5
  # eta gradient per expanded row
  fdE <- (.bt(.th, eta + .h)$value - .bt(.th, eta - .h)$value) / (2 * .h)
  expect_equal(as.numeric(got$gradEta), as.numeric(fdE), tolerance = 1e-4)
  # theta columns FD-agree; tcl1 moves only component-1 rows and the mixing
  # probability moves NOTHING (p-free conditional)
  for (t in seq_len(5)) {
    up <- .th
    up[t] <- up[t] + .h
    dn <- .th
    dn[t] <- dn[t] - .h
    fd <- (.bt(up, eta)$value - .bt(dn, eta)$value) / (2 * .h)
    expect_equal(as.numeric(got$gradTheta[, t]), as.numeric(fd),
                 tolerance = 1e-3, info = paste0("theta ", t))
  }
  .p1Idx <- .map$mixProbIdx[1]
  expect_true(all(got$gradTheta[, .p1Idx] == 0))
  # positive controls first (a bug zeroing all columns must NOT pass), then
  # the per-component selectivity
  expect_true(any(abs(got$gradTheta[1:4, 1]) > 1e-3)) # tcl1 x component 1
  expect_true(any(abs(got$gradTheta[5:8, 2]) > 1e-3)) # tcl2 x component 2
  expect_true(all(abs(got$gradTheta[5:8, 1]) < 1e-10)) # tcl1 x component 2
  expect_true(all(abs(got$gradTheta[1:4, 2]) < 1e-10)) # tcl2 x component 1
})

test_that("a fix()ed mixing probability inlines as a literal", {
  skip_on_cran()
  skip_if_not(.stanHasMixApi(),
              "nlmixr2est lacks the blessed mixture layout (#955)")
  .fixP <- function() {
    ini({
      tcl1 <- 1
      tcl2 <- 2
      tv <- 3
      p1 <- fix(0.3)
      add.sd <- 0.5
      eta.cl ~ 0.1
      prior(tcl1) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  # literalFix would dissolve it entirely; test the kept-in-vector path
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.fixP, .mixData(), est = "stan",
                        control = stanControl(run = FALSE,
                                              literalFix = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("log_sum_exp(log(0.3)", .lines, fixed = TRUE)))
  expect_false(any(grepl("real<lower=0,upper=1> p1;", .lines, fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("mixture end to end: assembled gradient + membership + fit contract", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(.stanHasMixApi(),
              "nlmixr2est lacks the blessed mixture layout (#955)")
  .d <- .mixData()
  .fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .mixMod, .d, est = "stan",
    control = stanControl(chains = 1L, iter = 500L, warmup = 250L,
                          seed = 42L, cores = 1L, calcTables = FALSE,
                          ofv = "none", onDiagnostic = "none"))))
  expect_true(inherits(.fit, "nlmixr2FitCore"))
  # membership posteriors: one row per subject, rows sum to 1
  .mp <- .fit$env$mixProb
  expect_true(is.data.frame(.mp))
  expect_equal(nrow(.mp), 4L)
  expect_equal(as.numeric(.mp$mix1 + .mp$mix2), rep(1, 4), tolerance = 1e-6)
  # the etaObf is the membership-weighted eta (one row per PHYSICAL subject)
  expect_equal(nrow(.fit$env$etaObf), 4L)
  # assembled gradient FD through the whole mixture target, INCLUDING the
  # mixing probability (whose gradient is pure Stan autodiff of the
  # log_sum_exp weights)
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    h <- stanLinkSetup(.mixMod, .d, thetaSens = TRUE, cores = 1L)
    on.exit({
      .Call(nlmixr2bayes:::`_nlmixr2bayes_clearThetaBase`)
      stanLinkFree()
    }, add = TRUE)
    .map <- nlmixr2bayes:::.stanMap(rxode2::rxode2(.mixMod))
    .Call(nlmixr2bayes:::`_nlmixr2bayes_setThetaBase`, as.double(h$initPar))
    .Call(nlmixr2bayes:::`_nlmixr2bayes_setMuRef`, as.integer(.map$muRefIdx))
    .sf <- .fit$env$stanfit
    set.seed(3)
    .pt <- list(tcl1 = 1.05, tcl2 = 1.95, tv = 2.95, p1 = 0.35,
                add_sd = 0.55, sd_eta_cl = 0.3,
                z_eta_cl = matrix(stats::rnorm(8, 0, 0.4), 8, 1))
    .up <- rstan::unconstrain_pars(.sf, .pt)
    .gA <- rstan::grad_log_prob(.sf, .up)
    attributes(.gA) <- NULL
    .gN <- numDeriv::grad(function(u) rstan::log_prob(.sf, u), .up,
                          method = "Richardson")
    expect_equal(.gA, .gN, tolerance = 1e-4)
  }
})
