# G15: the go/no-go for the committed censoring + ll() scope (D9).  The
# linked conditional follows nlmixr2est's NONMEM (adjLik) convention: a
# normal observation contributes -0.5*err^2/r - 0.5*log(r) (no -0.5*log(2pi))
# and a censored observation contributes log Phi(...) - 0.5*log(2pi)
# (censEst.h doCensNormal1).  Relative to the textbook density the value is
# therefore shifted by (nUncensored - nCensored) * 0.5*log(2pi) per subject
# -- a DATA-dependent constant, constant in every parameter Stan samples, so
# the target is a valid log density up to an additive constant.  These tests
# pin (a) that constant exactly (if it ever became parameter-dependent the
# equality breaks), and (b) the non-constant parts: the censored CDF term
# depends on the residual SD and the prediction, so eta AND theta gradients
# must carry it (FD-verified).

.censMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- 0.5
    eta.cl ~ 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.censData <- function(limit = NA_real_) {
  set.seed(42)
  .tt <- c(0.5, 1, 2, 4, 8)
  .d <- do.call(rbind, lapply(1:4, function(id) {
    .dv <- 5 * exp(-0.05 * .tt) + stats::rnorm(5, 0, 0.5)
    .cens <- as.integer(.dv < 4.5)
    .dv[.cens == 1] <- 4.5
    data.frame(ID = id, TIME = .tt, DV = .dv, CENS = .cens, AMT = 0, EVID = 0)
  }))
  if (!is.na(limit)) .d$LIMIT <- limit
  .d
}

# hand-computed textbook conditional for the fixture at theta=(1,3,.5).
# With a LIMIT column, UNCENSORED rows are M2: the truncated-normal density
# dnorm/P(X > limit) (Beal 2001) -- not the plain density.
.censHand <- function(d, eta, limit = NA_real_) {
  vapply(1:4, function(i) {
    .di <- d[d$ID == i, ]
    .f <- 100 / exp(3) * exp(-exp(1 + eta[i, 1]) / exp(3) * .di$TIME)
    .obs <- stats::dnorm(.di$DV, .f, 0.5, log = TRUE)
    if (!is.na(limit)) {
      .obs <- .obs - stats::pnorm(limit, .f, 0.5, log.p = TRUE,
                                  lower.tail = FALSE)
    }
    .cen <- if (is.na(limit)) {
      stats::pnorm(.di$DV, .f, 0.5, log.p = TRUE)      # M3: log Phi
    } else {                                            # M4: interval/tail
      log(stats::pnorm(.di$DV, .f, 0.5) - stats::pnorm(limit, .f, 0.5)) -
        stats::pnorm(limit, .f, 0.5, log.p = TRUE, lower.tail = FALSE)
    }
    sum(ifelse(.di$CENS == 1, .cen, .obs))
  }, numeric(1))
}

test_that("M3 censoring: value = textbook + the documented constant; gradients carry the CDF terms", {
  skip_on_cran()
  .d <- .censData()
  h <- stanLinkSetup(.censMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  .th <- c(1, 3, 0.5)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  got <- .bt(.th, .eta)
  expect_equal(got$nBad, 0L)
  # the adjLik constant, exactly: (nUncens - nCens) * 0.5*log(2pi) per
  # subject.  Parameter-dependence of anything dropped would break this.
  .nc <- tapply(.d$CENS, .d$ID, sum)
  .shift <- (5 - 2 * as.numeric(.nc)) * 0.5 * log(2 * pi)
  expect_equal(got$value, .censHand(.d, .eta) + .shift, tolerance = 1e-8)
  # eta gradient carries the censored CDF term
  .h <- 1e-5
  fdE <- (.bt(.th, .eta + .h)$value - .bt(.th, .eta - .h)$value) / (2 * .h)
  expect_equal(as.numeric(got$gradEta), as.numeric(fdE), tolerance = 1e-4)
  # theta gradients: tv (structural) and add.sd -- the censored
  # log Phi((loq-f)/sd) depends on BOTH, so neither column is constant
  fdT <- matrix(0, 4, 3)
  for (t in 1:3) {
    up <- .th
    up[t] <- up[t] + .h
    dn <- .th
    dn[t] <- dn[t] - .h
    fdT[, t] <- (.bt(up, .eta)$value - .bt(dn, .eta)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradTheta), as.numeric(fdT), tolerance = 1e-3)
  expect_true(all(abs(fdT[, 3]) > 1e-3)) # add.sd really enters the CDF
})

test_that("M4 censoring (LIMIT column): value constant + FD gradients", {
  skip_on_cran()
  .d <- .censData(limit = 0)
  h <- stanLinkSetup(.censMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  .th <- c(1, 3, 0.5)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  got <- .bt(.th, .eta)
  expect_equal(got$nBad, 0L)
  # under M2 (uncensored + LIMIT) the adjLik constant cancels against the
  # missing -0.5*log(2pi), so uncensored rows are the EXACT truncated
  # density; only the censored M4 rows carry the -0.5*log(2pi) constant
  .nc <- tapply(.d$CENS, .d$ID, sum)
  .shift <- -as.numeric(.nc) * 0.5 * log(2 * pi)
  expect_equal(got$value, .censHand(.d, .eta, limit = 0) + .shift,
               tolerance = 1e-8)
  .h <- 1e-5
  fdE <- (.bt(.th, .eta + .h)$value - .bt(.th, .eta - .h)$value) / (2 * .h)
  expect_equal(as.numeric(got$gradEta), as.numeric(fdE), tolerance = 1e-4)
  fdT <- matrix(0, 4, 3)
  for (t in 1:3) {
    up <- .th
    up[t] <- up[t] + .h
    dn <- .th
    dn[t] <- dn[t] - .h
    fdT[, t] <- (.bt(up, .eta)$value - .bt(dn, .eta)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradTheta), as.numeric(fdT), tolerance = 1e-3)
})

test_that("ll() endpoint: twin of add() up to the per-obs constant; FD gradients", {
  skip_on_cran()
  # the conditional path sums llikObs directly (npEvalCondLik), so an ll()
  # observation contributes the user's expression EXACTLY -- the nNonNormal
  # adjLik adjustment applies only to the fInd->llik finalization, not here.
  # The add() twin omits -0.5*log(2pi) per obs, so the twin difference is
  # exactly nobs * 0.5*log(2pi) per subject -- constant in every parameter
  .llMod <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- 0.5
      eta.cl ~ 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) -
        0.5 * ((DV - cp) / add.sd)^2
    })
  }
  .d <- .linkData()
  .eta <- matrix(c(-0.1, 0.05, 0.2, -0.15), 4, 1)
  .th <- c(1, 3, 0.5)
  .bt <- function(theta, e) {
    .Call(nlmixr2stan:::`_nlmixr2stan_condBatchTheta`, as.double(theta),
          as.matrix(e))
  }
  h <- stanLinkSetup(.censMod, .d, thetaSens = TRUE, cores = 1L)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .addV <- .bt(.th, .eta)
  stanLinkFree()
  .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
  h2 <- stanLinkSetup(.llMod, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h2$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  got <- .bt(.th, .eta)
  expect_equal(got$nBad, 0L)
  # ll() value = the full textbook density; add() twin = same + 5*0.5*log(2pi)
  .full <- vapply(1:4, function(i) {
    .di <- .d[.d$ID == i, ]
    .f <- 100 / exp(3) * exp(-exp(1 + .eta[i, 1]) / exp(3) * .di$TIME)
    sum(stats::dnorm(.di$DV, .f, 0.5, log = TRUE))
  }, numeric(1))
  expect_equal(got$value, .full, tolerance = 1e-8)
  expect_equal(.addV$value, .full + 5 * 0.5 * log(2 * pi), tolerance = 1e-8)
  # gradients agree with the twin exactly (the constant differentiates away)
  expect_equal(got$gradEta, .addV$gradEta, tolerance = 1e-6)
  # and FD-agree in their own right
  .h <- 1e-5
  fdE <- (.bt(.th, .eta + .h)$value - .bt(.th, .eta - .h)$value) / (2 * .h)
  expect_equal(as.numeric(got$gradEta), as.numeric(fdE), tolerance = 1e-4)
  fdT <- matrix(0, 4, 3)
  for (t in 1:3) {
    up <- .th
    up[t] <- up[t] + .h
    dn <- .th
    dn[t] <- dn[t] - .h
    fdT[, t] <- (.bt(up, .eta)$value - .bt(dn, .eta)$value) / (2 * .h)
  }
  expect_equal(as.numeric(got$gradTheta), as.numeric(fdT), tolerance = 1e-3)
})
