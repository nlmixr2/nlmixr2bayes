# Estimated transform-both-sides lambda (unlocked by nlmixr2/nlmixr2est#949
# + the Stan-side DV-transform Jacobian).  The linked conditional evaluates
# the TRANSFORMED-scale density with an exact d/dlambda column; the
# generator adds target += (lambda - 1) * sumLogJac (a pure data statistic)
# so the assembled target is the correct untransformed-scale density with
# an exact lambda gradient end to end.

.tbsMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- 0.5
    lambda <- c(-2, 0.5, 2)
    eta.cl ~ 0.1
    prior(tcl) ~ dnorm(1, 2)
    prior(add.sd) ~ dcauchy(0, 2.5)
    prior(lambda) ~ dnorm(0.5, 1)
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd) + boxCox(lambda)
  })
}

.tbsData <- function() {
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  # includes a dose row (EVID=1, NA DV) per subject: the Jacobian statistic
  # must count OBSERVATIONS only
  do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 0, EVID = 2),
          data.frame(ID = id, TIME = .tt,
                     DV = abs(5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5)),
                     AMT = 0, EVID = 0))
  }))
}

# observation DVs only (the Jacobian statistic's domain)
.tbsObsDv <- function(d) d$DV[d$EVID == 0 & !is.na(d$DV)]

test_that("estimated lambda: the linked d/dlambda column FD-agrees (#949)", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasTbsLambdaSens(),
              "nlmixr2est lacks the #949 lambda sensitivities")
  h <- stanLinkSetup(.tbsMod, .tbsData(), thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  expect_true(4L %in% h$thetaSensIdx)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  .eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  .th <- c(1, 3, 0.5, 0.5)
  got <- .bt(.th, .eta)
  expect_equal(got$nBad, 0L)
  expect_true(all(abs(got$gradTheta[, 4]) > 1e-3)) # real, not silently zero
  .h <- 1e-5
  for (t in 2:4) {
    up <- .th
    up[t] <- up[t] + .h
    dn <- .th
    dn[t] <- dn[t] - .h
    fd <- (.bt(up, .eta)$value - .bt(dn, .eta)$value) / (2 * .h)
    expect_equal(as.numeric(got$gradTheta[, t]), as.numeric(fd),
                 tolerance = 1e-3, info = paste0("theta ", t))
  }
})

test_that("estimated lambda: the generator emits the Jacobian statistic", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasTbsLambdaSens(),
              "nlmixr2est lacks the #949 lambda sensitivities")
  .d <- .tbsData()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.tbsMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("real sumLogJac_4;", .lines, fixed = TRUE)))
  expect_true(any(grepl("target += (lambda - 1) * sumLogJac_4;", .lines,
                        fixed = TRUE)))
  # the statistic is exactly sum(log DV) for Box-Cox
  expect_equal(.code$data$sumLogJac_4, sum(log(.tbsObsDv(.d))),
               tolerance = 1e-12)
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
  # Yeo-Johnson variant: signed statistic
  .yjMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- 0.5
      lambda <- c(-2, 0.5, 2)
      eta.cl ~ 0.1
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      prior(lambda) ~ dnorm(0.5, 1)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd) + yeoJohnson(lambda)
    })
  }
  .codeY <- suppressMessages(
    nlmixr2est::nlmixr2(.yjMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .dvObs <- .tbsObsDv(.d)
  expect_equal(.codeY$data$sumLogJac_4,
               sum(log1p(.dvObs[.dvObs >= 0])) -
                 sum(log1p(-.dvObs[.dvObs < 0])), tolerance = 1e-12)
  # a FIXED lambda emits no Jacobian machinery (constant shift)
  .fixMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- 0.5
      lambda <- fix(0.5)
      eta.cl ~ 0.1
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd) + boxCox(lambda)
    })
  }
  .codeF <- suppressMessages(
    nlmixr2est::nlmixr2(.fixMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_false(any(grepl("sumLogJac", strsplit(.codeF$code, "\n")[[1]],
                         fixed = TRUE)))
  # Box-Cox with non-positive DV is refused with the endpoint named
  .dBad <- .d
  .dBad$DV[3] <- -1
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.tbsMod, .dBad, est = "stan",
                          control = stanControl(run = FALSE))),
    "strictly positive")
})

test_that("estimated lambda: the assembled target is the untransformed-scale density", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(nlmixr2stan:::.stanHasTbsLambdaSens(),
              "nlmixr2est lacks the #949 lambda sensitivities")
  .d <- .tbsData()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.tbsMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .sm <- stanCompile(.code$code)
  h <- stanLinkSetup(.tbsMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .sf <- rstan::sampling(.sm, data = .code$data, chains = 1, iter = 2,
                         warmup = 1, refresh = 0, cores = 1,
                         show_messages = FALSE,
                         init = list(list(tcl = 1, tv = 3, add_sd = 0.5,
                                          lambda = 0.5, sd_b1 = sqrt(0.1),
                                          z_b1 = matrix(0, 4, 1))))
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  .eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  .sJ <- .code$data$sumLogJac_4
  # difference form in lambda ONLY (everything else cancels): the target
  # difference must equal the conditional difference + the Jacobian
  # difference + the lambda-prior difference, exactly
  .pt <- function(lam) {
    list(tcl = 1, tv = 3, add_sd = 0.5, lambda = lam, sd_b1 = sqrt(0.1),
         z_b1 = .eta / sqrt(0.1))
  }
  .l1 <- 0.4
  .l2 <- 0.7
  .lp <- function(lam) {
    rstan::log_prob(.sf, rstan::unconstrain_pars(.sf, .pt(lam)),
                    adjust_transform = FALSE)
  }
  .cond <- function(lam) sum(.bt(c(1, 3, 0.5, lam), .eta)$value)
  .lhs <- .lp(.l2) - .lp(.l1)
  .rhs <- (.cond(.l2) - .cond(.l1)) +
    (.l2 - .l1) * .sJ +
    (stats::dnorm(.l2, 0.5, 1, log = TRUE) -
       stats::dnorm(.l1, 0.5, 1, log = TRUE))
  expect_equal(.lhs, .rhs, tolerance = 1e-6)
  # and Stan's assembled gradient (incl. the Jacobian's d/dlambda = sJ)
  # agrees with finite differences of its own log density
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    .up <- rstan::unconstrain_pars(.sf, .pt(0.5))
    .gA <- rstan::grad_log_prob(.sf, .up)
    attributes(.gA) <- NULL
    .gN <- numDeriv::grad(function(u) rstan::log_prob(.sf, u), .up,
                          method = "Richardson")
    expect_equal(.gA, .gN, tolerance = 1e-4)
  }
})

test_that("estimated lambda on a logit base: statistic on the transformed DV; target exact", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(nlmixr2stan:::.stanHasTbsLambdaSens(),
              "nlmixr2est lacks the #949 lambda sensitivities")
  # logitNorm + yeoJohnson: the lambda acts on t = logit((y-lo)/(hi-lo)),
  # so the Jacobian statistic must be the YJ statistic ON t, not on the raw
  # DV (the review question that prompted this gate); the linked d/dlambda
  # column is exact for the combined transform (FD 9e-11)
  .lgMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- 0.5
      lambda <- c(-2, 0.5, 2)
      eta.cl ~ 0.1
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      prior(lambda) ~ dnorm(0.5, 1)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ logitNorm(add.sd, 0, 20) + yeoJohnson(lambda)
    })
  }
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  .d <- do.call(rbind, lapply(1:4, function(id) {
    data.frame(ID = id, TIME = .tt,
               DV = pmin(19.5, pmax(0.5, 5 * exp(-0.05 * .tt) +
                                      stats::rnorm(5, 0, 0.5))),
               AMT = 0, EVID = 0)
  }))
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.lgMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  # the statistic is the YJ statistic on the LOGIT-transformed DV
  .t <- stats::qlogis(.d$DV / 20)
  expect_equal(.code$data$sumLogJac_4,
               sum(log1p(.t[.t >= 0])) - sum(log1p(-.t[.t < 0])),
               tolerance = 1e-12)
  # assembled-target exactness in difference form (lambda only)
  .sm <- stanCompile(.code$code)
  h <- stanLinkSetup(.lgMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .sf <- rstan::sampling(.sm, data = .code$data, chains = 1, iter = 2,
                         warmup = 1, refresh = 0, cores = 1,
                         show_messages = FALSE,
                         init = list(list(tcl = 1, tv = 3, add_sd = 0.5,
                                          lambda = 0.5, sd_b1 = sqrt(0.1),
                                          z_b1 = matrix(0, 4, 1))))
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  .eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  .sJ <- .code$data$sumLogJac_4
  .pt <- function(lam) {
    list(tcl = 1, tv = 3, add_sd = 0.5, lambda = lam, sd_b1 = sqrt(0.1),
         z_b1 = .eta / sqrt(0.1))
  }
  .lp <- function(lam) {
    rstan::log_prob(.sf, rstan::unconstrain_pars(.sf, .pt(lam)),
                    adjust_transform = FALSE)
  }
  .cond <- function(lam) sum(.bt(c(1, 3, 0.5, lam), .eta)$value)
  .l1 <- 0.4
  .l2 <- 0.7
  .lhs <- .lp(.l2) - .lp(.l1)
  .rhs <- (.cond(.l2) - .cond(.l1)) + (.l2 - .l1) * .sJ +
    (stats::dnorm(.l2, 0.5, 1, log = TRUE) -
       stats::dnorm(.l1, 0.5, 1, log = TRUE))
  expect_equal(.lhs, .rhs, tolerance = 1e-6)
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    .up <- rstan::unconstrain_pars(.sf, .pt(0.5))
    .gA <- rstan::grad_log_prob(.sf, .up)
    attributes(.gA) <- NULL
    .gN <- numDeriv::grad(function(u) rstan::log_prob(.sf, u), .up,
                          method = "Richardson")
    expect_equal(.gA, .gN, tolerance = 1e-4)
  }
})
