## Covariates.  rxode2 already reads these straight from the event table, so
## the bridge needs no change; what has to be right is that codegen notices
## them, tells baseline from time-varying, and declares them in Stan at the
## right length.  Getting the length wrong would silently collapse a moving
## covariate to its first value per subject.
library(testthat)

wtModel <- function() {
  ini({
    tka <- 0.5
    tcl <- 1.0
    tv <- 3.4
    wtcl <- 0.75
    add.sd <- 0.2
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + wtcl * log(wt / 70))
    v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

## A circadian input: the covariate moves within every subject.
circModel <- function() {
  ini({
    tkin <- 0
    tkout <- 0
    amp <- 0.1
    add.sd <- 0.05
  })
  model({
    kin <- exp(tkin) + amp * cos(ctime / 24)
    kout <- exp(tkout)
    R(0) <- 1
    d/dt(R) <- kin - kout * R
    resp <- R
    resp ~ add(add.sd)
  })
}

wtData <- function(nsub = 6L, wts = NULL) {
  if (is.null(wts)) wts <- seq(50, 100, length.out = nsub)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = "depot")
    e <- rxode2::et(e, c(0.5, 1, 2, 4, 8, 12, 24))
    x <- as.data.frame(e)
    x$id <- i
    x$dv <- 1
    x$wt <- wts[i]
    x
  }))
}

circData <- function(nsub = 3L) {
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(seq(0, 48, by = 4))
    x <- as.data.frame(e)
    x$id <- i
    x$dv <- 1
    x$ctime <- (x$time + 4 * i) %% 24
    x
  }))
}

test_that("a baseline covariate is declared once per subject", {
  skipUnlessStan(nlmixr2 = FALSE)
  d <- wtData()
  g <- rxsStanFromUi(wtModel, d)
  on.exit(rxsRelease(g$handle), add = TRUE)

  expect_match(g$code, "vector[nSub] cov_wt;", fixed = TRUE)
  expect_match(g$code, "real wt = cov_wt[s];", fixed = TRUE)
  ## One value per subject, in subject order -- not one per row.
  expect_length(g$standata$cov_wt, length(unique(d$id)))
  expect_equal(g$standata$cov_wt, seq(50, 100, length.out = 6L))
})

test_that("a time-varying covariate is declared once per observation", {
  skipUnlessStan(nlmixr2 = FALSE)
  d <- circData()
  g <- rxsStanFromUi(circModel, d)
  on.exit(rxsRelease(g$handle), add = TRUE)

  expect_match(g$code, "vector[nObs] cov_ctime;", fixed = TRUE)
  expect_match(g$code, "real ctime = cov_ctime[i];", fixed = TRUE)
  obs <- d[d$evid == 0, ]
  expect_length(g$standata$cov_ctime, nrow(obs))
  expect_equal(g$standata$cov_ctime, obs$ctime)
})

test_that("a symbol that is neither parameter nor column is named", {
  skipUnlessStan(nlmixr2 = FALSE)
  d <- wtData()
  d$wt <- NULL
  ## Must not surface as rxode2's expanded-model dump.
  expect_error(rxsStanFromUi(wtModel, d), "uses wt", fixed = TRUE)
})

test_that("covariates leave the fast path armed", {
  skipUnlessStan(nlmixr2 = FALSE)
  for (case in list(list(wtModel, wtData()), list(circModel, circData()))) {
    g <- rxsStanFromUi(case[[1]], case[[2]])
    expect_true(rxsFastAvailable(g$handle))
    rxsRelease(g$handle)
  }
})

test_that("gradients are right with a covariate in the transform", {
  skipUnlessStan(nlmixr2 = FALSE)
  d <- wtData()
  set.seed(11)
  g <- rxsStanFromUi(wtModel, d, atol = 1e-10, rtol = 1e-10)
  on.exit(rxsRelease(g$handle), add = TRUE)
  sm <- stanModelFor(g$code, "cov_wt_grad")
  fit <- rstan::sampling(sm, data = g$standata, chains = 1, iter = 1,
                         refresh = 0, seed = 1)
  err <- rxsCheckGradient(fit, rep(0.1, rstan::get_num_upars(fit)))
  expect_lt(max(err$relDiff), 1e-6)
})
