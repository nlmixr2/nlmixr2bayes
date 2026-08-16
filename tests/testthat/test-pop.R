# Tier 0: population-only (no-eta) models via the nlm C API
# (nlmixr2/nlmixr2est#953).  The generated program has theta parameters
# only and one external scalar carrying the whole data log-likelihood with
# its analytic gradient; single-subject Bayesian fits (the Torsten
# pk2cpt-style examples) are the target use.

.popMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- c(0, 0.5)
    prior(tcl) ~ dnorm(1, 2)
    prior(tv) ~ dnorm(3, 2)
    prior(add.sd) ~ dcauchy(0, 2.5)
  })
  model({
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.popData <- function() {
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = .tt,
               DV = 5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5),
               AMT = 0, EVID = 0)
  }))
}

test_that("tier 0 codegen: theta-only program with nlmixr2_pop_ll", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasNlmApi(),
              "nlmixr2est lacks the nlm C API (#953)")
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.popMod, .popData(), est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("real nlmixr2_pop_ll(vector theta);", .lines,
                        fixed = TRUE)))
  expect_true(any(grepl("target += nlmixr2_pop_ll(theta);", .lines,
                        fixed = TRUE)))
  # no eta machinery at all ("\\beta\\b" avoids matching "theta")
  expect_false(any(grepl("\\beta\\b|omegaOut|z_b|logLikSubj", .lines)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("tier 0 refuses undissolved fixed thetas (theta-length mismatch)", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasNlmApi(),
              "nlmixr2est lacks the nlm C API (#953)")
  .fixMod <- function() {
    ini({
      tcl <- 1
      tv <- fix(3)
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  # literalFix=FALSE keeps the fixed row in the theta vector; the nlm
  # problem estimates only free thetas, so the mismatch is refused at
  # GENERATION (before any compile) rather than failing at the first
  # evaluation
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.fixMod, .popData(), est = "stan",
                          control = stanControl(run = FALSE,
                                                literalFix = FALSE))),
    "literalFix")
  # with the default literalFix=TRUE the fixed theta dissolves and the
  # model generates
  expect_s3_class(
    suppressMessages(
      nlmixr2est::nlmixr2(.fixMod, .popData(), est = "stan",
                          control = stanControl(run = FALSE))),
    "nlmixr2stanCode")
})

test_that("tier 0 link: value ties to the hand density; gradient FD-agrees", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasNlmApi(),
              "nlmixr2est lacks the nlm C API (#953)")
  .d <- .popData()
  h <- stanPopLinkSetup(.popMod, .d)
  on.exit(stanLinkFree(), add = TRUE)
  expect_true(h$pop)
  expect_equal(h$ntheta, 3L)
  .th <- c(1, 3, 0.5)
  .e <- nlmixr2stan:::.popEval(.th)
  expect_equal(.e$nBad, 0L)
  # raw nlm convention: value = -logLik with all constants (log p = -value)
  .hand <- sum(vapply(1:4, function(i) {
    .di <- .d[.d$ID == i, ]
    .f <- 100 / exp(3) * exp(-exp(1) / exp(3) * .di$TIME)
    sum(stats::dnorm(.di$DV, .f, 0.5, log = TRUE))
  }, numeric(1)))
  expect_equal(.e$value, -.hand, tolerance = 1e-8)
  .h <- 1e-6
  .fd <- vapply(1:3, function(k) {
    .up <- .th
    .up[k] <- .up[k] + .h
    .dn <- .th
    .dn[k] <- .dn[k] - .h
    (nlmixr2stan:::.popEval(.up)$value -
       nlmixr2stan:::.popEval(.dn)$value) / (2 * .h)
  }, numeric(1))
  expect_equal(.e$grad, .fd, tolerance = 1e-6)
  expect_identical(nlmixr2stan:::.popEval(.th), nlmixr2stan:::.popEval(.th))
})

test_that("tier 0 end to end: fit, assembled gradient, fit contract", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(nlmixr2stan:::.stanHasNlmApi(),
              "nlmixr2est lacks the nlm C API (#953)")
  .d <- .popData()
  .fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .popMod, .d, est = "stan",
    control = stanControl(chains = 1L, iter = 600L, warmup = 300L,
                          seed = 42L, cores = 1L, calcTables = FALSE,
                          ofv = "none", onDiagnostic = "none"))))
  expect_true(inherits(.fit, "nlmixr2FitCore"))
  expect_s4_class(.fit$env$stanfit, "stanfit")
  expect_equal(names(.fit$env$theta), c("tcl", "tv", "add.sd"))
  # the nlm objective row: -2 logLik at the posterior point estimate
  expect_true(is.finite(.fit$env$objective))
  # assembled gradient FD (the pop external + priors + constraints)
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    h <- stanPopLinkSetup(.popMod, .d)
    on.exit(stanLinkFree(), add = TRUE)
    .sf <- .fit$env$stanfit
    .up <- rstan::unconstrain_pars(.sf, list(tcl = 1.05, tv = 2.95,
                                             add_sd = 0.55))
    .gA <- rstan::grad_log_prob(.sf, .up)
    attributes(.gA) <- NULL
    .gN <- numDeriv::grad(function(u) rstan::log_prob(.sf, u), .up,
                          method = "Richardson")
    expect_equal(.gA, .gN, tolerance = 1e-4)
  }
})

test_that("tier 0 + estimated lambda: nlm convention pinned; assembled target exact", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(nlmixr2stan:::.stanHasNlmApi(),
              "nlmixr2est lacks the nlm C API (#953)")
  .tbsPop <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- 0.5
      lambda <- c(-2, 0.5, 2)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      prior(lambda) ~ dnorm(0.5, 1)
    })
    model({
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd) + boxCox(lambda)
    })
  }
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  .d <- do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = .tt,
               DV = abs(5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5)),
               AMT = 0, EVID = 0)
  }))
  h <- stanPopLinkSetup(.tbsPop, .d)
  on.exit(stanLinkFree(), add = TRUE)
  .th <- c(1, 3, 0.5, 0.5)
  .e <- nlmixr2stan:::.popEval(.th)
  # the nlm convention, PINNED: value = -logLik of the TRANSFORMED-scale
  # density WITHOUT the DV Jacobian (same as the focei conditional) -- the
  # Stan-side (lambda-1)*sJ term is therefore an addition, not a double
  # count.  If upstream ever folds the Jacobian into the value this
  # assertion catches the would-be double count.
  .bc <- function(y, l) (y^l - 1) / l
  .f <- 100 / exp(3) * exp(-exp(1) / exp(3) * .d$TIME)
  .llNoJac <- sum(stats::dnorm(.bc(.d$DV, 0.5), .bc(.f, 0.5), 0.5,
                               log = TRUE))
  expect_equal(.e$value, -.llNoJac, tolerance = 1e-8)
  # the lambda gradient is real and FD-agrees
  .hh <- 1e-5
  .up <- .th
  .up[4] <- .up[4] + .hh
  .dn <- .th
  .dn[4] <- .dn[4] - .hh
  .fd <- (nlmixr2stan:::.popEval(.up)$value -
            nlmixr2stan:::.popEval(.dn)$value) / (2 * .hh)
  expect_equal(.e$grad[4], .fd, tolerance = 1e-4)
  stanLinkFree()
  # assembled pop target: difference form in lambda + full gradient FD
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.tbsPop, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_equal(.code$data$sumLogJac_4, sum(log(.d$DV)), tolerance = 1e-12)
  .sm <- stanCompile(.code$code)
  h2 <- stanPopLinkSetup(.tbsPop, .d)
  .sf <- rstan::sampling(.sm, data = .code$data, chains = 1, iter = 2,
                         warmup = 1, refresh = 0, cores = 1,
                         show_messages = FALSE,
                         init = list(list(tcl = 1, tv = 3, add_sd = 0.5,
                                          lambda = 0.5)))
  .pt <- function(lam) list(tcl = 1, tv = 3, add_sd = 0.5, lambda = lam)
  .lp <- function(lam) {
    rstan::log_prob(.sf, rstan::unconstrain_pars(.sf, .pt(lam)),
                    adjust_transform = FALSE)
  }
  .pop <- function(lam) -nlmixr2stan:::.popEval(c(1, 3, 0.5, lam))$value
  .sJ <- .code$data$sumLogJac_4
  .l1 <- 0.4
  .l2 <- 0.7
  .lhs <- .lp(.l2) - .lp(.l1)
  .rhs <- (.pop(.l2) - .pop(.l1)) + (.l2 - .l1) * .sJ +
    (stats::dnorm(.l2, 0.5, 1, log = TRUE) -
       stats::dnorm(.l1, 0.5, 1, log = TRUE))
  expect_equal(.lhs, .rhs, tolerance = 1e-6)
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    .up2 <- rstan::unconstrain_pars(.sf, .pt(0.5))
    .gA <- rstan::grad_log_prob(.sf, .up2)
    attributes(.gA) <- NULL
    .gN <- numDeriv::grad(function(u) rstan::log_prob(.sf, u), .up2,
                          method = "Richardson")
    expect_equal(.gA, .gN, tolerance = 1e-4)
  }
})

test_that("tier 0 iteration print: rows + nlm parameter history", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasNlmApi(),
              "nlmixr2est lacks the nlm C API (#953)")
  .d <- .popData()
  # the scale.h table prints via Rprintf (stdout); setup messages go to the
  # message stream -- capture both
  .out <- utils::capture.output(type = "output", {
    .msg <- utils::capture.output(type = "message", {
      h <- stanPopLinkSetup(.popMod, .d, print = 2L)
      for (i in 1:5) invisible(nlmixr2stan:::.popEval(h$initPar))
      .ph <- nlmixr2est::nlmGetParHist(TRUE)
      stanLinkFree()
    })
  })
  .out <- c(.out, .msg)
  # never a DOUBLE header (the regression: setup + lazy row-print each
  # printing one); whether it prints at all depends on residency state
  # carried across setups in one process, so 0 is acceptable in-suite
  expect_true(sum(grepl("Function Val", .out, fixed = TRUE)) <= 1L)
  expect_true(is.data.frame(.ph))
  expect_equal(sum(.ph$type == "Scaled"), 3L)
  expect_true(all(c("objf", "tcl", "tv", "add.sd") %in% names(.ph)))
})
