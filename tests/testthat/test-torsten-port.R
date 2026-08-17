# The Torsten-comparison vignette's ported models must be REAL: each port
# goes through est="stan" codegen (run=FALSE) and, with rstan, the complete
# generated program passes stanc.  (Sampling-level agreement with Torsten
# itself would need cmdstan + Torsten installed; codegen acceptance is what
# keeps the vignette honest.)

test_that("Torsten poppk2cpt port generates and parses", {
  skip_on_cran()
  .poppk2cpt <- function() {
    ini({
      lcl <- log(10)
      lq <- log(15)
      lv1 <- log(35)
      lv2 <- log(105)
      lka <- log(2.5)
      prior(lcl) ~ dnorm(log(10), 0.25)
      prior(lq) ~ dnorm(log(15), 0.5)
      prior(lv1) ~ dnorm(log(35), 0.25)
      prior(lv2) ~ dnorm(log(105), 0.5)
      prior(lka) ~ dnorm(log(2.5), 1)
      lnorm.sd <- c(0, 0.2)
      prior(lnorm.sd) ~ dcauchy(0, 1)
      eta.cl + eta.q + eta.v1 + eta.v2 + eta.ka ~
        c(0.06,
          0.01, 0.06,
          0.01, 0.01, 0.06,
          0.01, 0.01, 0.01, 0.06,
          0.01, 0.01, 0.01, 0.01, 0.06)
    })
    model({
      cl <- exp(lcl + eta.cl)
      q <- exp(lq + eta.q)
      v <- exp(lv1 + eta.v1)
      vp <- exp(lv2 + eta.v2)
      ka <- exp(lka + eta.ka)
      cp <- linCmt()
      cp ~ lnorm(lnorm.sd)
    })
  }
  .d <- do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 1000, EVID = 1),
          data.frame(ID = id, TIME = c(1, 2, 4, 8, 12, 24),
                     DV = c(20, 25, 18, 10, 6, 2), AMT = 0, EVID = 0))
  }))
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.poppk2cpt, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  # the announced default omega prior mirrors Torsten's hand-coded choice
  expect_true(any(grepl("lkj_corr_cholesky", strsplit(.code$code, "\n")[[1]],
                        fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("Torsten Friberg-Karlsson port (two endpoints) generates and parses", {
  skip_on_cran()
  .fk <- function() {
    ini({
      lcl <- log(10)
      lq <- log(15)
      lv1 <- log(35)
      lv2 <- log(105)
      lka <- log(2)
      lmtt <- log(125)
      lcirc0 <- log(5)
      lalpha <- log(3e-4)
      lgamma <- log(0.17)
      prior(lcl) ~ dnorm(log(10), 0.25)
      prior(lmtt) ~ dnorm(log(125), 0.25)
      prior(lcirc0) ~ dnorm(log(5), 0.25)
      prior(lalpha) ~ dnorm(log(3e-4), 1)
      prior(lgamma) ~ dnorm(log(0.17), 0.5)
      lnorm.sd <- c(0, 0.2)
      lnormNeut.sd <- c(0, 0.3)
      prior(lnorm.sd) ~ dcauchy(0, 1)
      prior(lnormNeut.sd) ~ dcauchy(0, 1)
      eta.cl ~ 0.06
      eta.mtt ~ 0.04
      eta.circ0 ~ 0.04
    })
    model({
      cl <- exp(lcl + eta.cl)
      q <- exp(lq)
      v1 <- exp(lv1)
      v2 <- exp(lv2)
      ka <- exp(lka)
      mtt <- exp(lmtt + eta.mtt)
      circ0 <- exp(lcirc0 + eta.circ0)
      alpha <- exp(lalpha)
      gamma <- exp(lgamma)
      k10 <- cl / v1
      k12 <- q / v1
      k21 <- q / v2
      ktr <- 4 / mtt
      d / dt(depot) <- -ka * depot
      d / dt(cent) <- ka * depot - (k10 + k12) * cent + k21 * peri
      d / dt(peri) <- k12 * cent - k21 * peri
      conc <- cent / v1
      edrug <- alpha * conc
      d / dt(prol) <- ktr * prol * (1 - edrug) * (circ0 / circ)^gamma -
        ktr * prol
      d / dt(tr1) <- ktr * (prol - tr1)
      d / dt(tr2) <- ktr * (tr1 - tr2)
      d / dt(tr3) <- ktr * (tr2 - tr3)
      d / dt(circ) <- ktr * tr3 - ktr * circ
      prol(0) <- circ0
      tr1(0) <- circ0
      tr2(0) <- circ0
      tr3(0) <- circ0
      circ(0) <- circ0
      conc ~ lnorm(lnorm.sd)
      circ ~ lnorm(lnormNeut.sd)
    })
  }
  .d <- do.call(rbind, lapply(1:3, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 80, EVID = 1,
                     CMT = "depot"),
          data.frame(ID = id, TIME = c(1, 4, 12, 24),
                     DV = c(1.5, 1.1, 0.5, 0.2), AMT = 0, EVID = 0,
                     CMT = "conc"),
          data.frame(ID = id, TIME = c(24, 96, 168, 336),
                     DV = c(4.9, 3.2, 2.1, 4.4), AMT = 0, EVID = 0,
                     CMT = "circ"))
  }))
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.fk, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("Torsten pk2cpt port (single patient, no etas) generates and parses", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasNlmApi(),
              "nlmixr2est lacks the nlm C API (#953)")
  # the single-patient two-compartment example, Torsten's informative
  # lognormal priors verbatim (normal on the log-scale parameters)
  .pk2cpt <- function() {
    ini({
      lcl <- log(10)
      lq <- log(15)
      lv1 <- log(35)
      lv2 <- log(105)
      lka <- log(2.5)
      prior(lcl) ~ dnorm(log(10), 0.25)
      prior(lq) ~ dnorm(log(15), 0.5)
      prior(lv1) ~ dnorm(log(35), 0.25)
      prior(lv2) ~ dnorm(log(105), 0.5)
      prior(lka) ~ dnorm(log(2.5), 1)
      lnorm.sd <- c(0, 0.2)
      prior(lnorm.sd) ~ dcauchy(0, 1)
    })
    model({
      cl <- exp(lcl)
      q <- exp(lq)
      v <- exp(lv1)
      vp <- exp(lv2)
      ka <- exp(lka)
      cp <- linCmt()
      cp ~ lnorm(lnorm.sd)
    })
  }
  .d <- rbind(data.frame(ID = 1, TIME = 0, DV = NA_real_, AMT = 1000,
                         EVID = 1),
              data.frame(ID = 1, TIME = c(0.5, 1, 2, 4, 8, 12, 24),
                         DV = c(15, 22, 24, 18, 12, 8, 3), AMT = 0,
                         EVID = 0))
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.pk2cpt, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("nlmixr2_pop_ll", .lines, fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("IOV occasion-indicator pattern generates and parses", {
  skip_on_cran()
  .iov <- function() {
    ini({
      lcl <- log(10)
      lv <- log(35)
      prior(lcl) ~ dnorm(log(10), 0.5)
      lnorm.sd <- c(0, 0.2)
      prior(lnorm.sd) ~ dcauchy(0, 1)
      eta.cl ~ 0.06
      eta.occ1 ~ 0.02
      eta.occ2 ~ 0.02
    })
    model({
      kappa <- eta.occ1 * (OCC == 1) + eta.occ2 * (OCC == 2)
      cl <- exp(lcl + eta.cl + kappa)
      v <- exp(lv)
      cp <- linCmt()
      cp ~ lnorm(lnorm.sd)
    })
  }
  .d <- do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 100, EVID = 1,
                     OCC = 1),
          data.frame(ID = id, TIME = c(1, 4, 8), DV = c(2, 1.4, 0.8),
                     AMT = 0, EVID = 0, OCC = 1),
          data.frame(ID = id, TIME = 24, DV = NA_real_, AMT = 100, EVID = 1,
                     OCC = 2),
          data.frame(ID = id, TIME = c(25, 28, 32), DV = c(2.1, 1.5, 0.9),
                     AMT = 0, EVID = 0, OCC = 2))
  }))
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iov, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("Torsten effCpt port (effect-compartment population PK/PD) generates and parses", {
  skip_on_cran()
  .effCpt <- function() {
    ini({
      lcl <- log(10)
      lq <- log(15)
      lv1 <- log(35)
      lv2 <- log(105)
      lka <- log(2)
      lke0 <- log(1)
      lec50 <- log(100)
      prior(lcl) ~ dnorm(log(10), 0.25)
      prior(lke0) ~ dnorm(log(1), 0.25)
      prior(lec50) ~ dnorm(log(100), 0.5)
      lnorm.sd <- c(0, 0.2)
      resp.sd <- c(0, 5)
      prior(lnorm.sd) ~ dcauchy(0, 1)
      prior(resp.sd) ~ dcauchy(0, 5)
      eta.cl ~ 0.06
      eta.v1 ~ 0.06
      eta.ke0 ~ 0.04
      eta.ec50 ~ 0.04
    })
    model({
      cl <- exp(lcl + eta.cl)
      q <- exp(lq)
      v1 <- exp(lv1 + eta.v1)
      v2 <- exp(lv2)
      ka <- exp(lka)
      ke0 <- exp(lke0 + eta.ke0)
      ec50 <- exp(lec50 + eta.ec50)
      k10 <- cl / v1
      k12 <- q / v1
      k21 <- q / v2
      d / dt(depot) <- -ka * depot
      d / dt(cent) <- ka * depot - (k10 + k12) * cent + k21 * peri
      d / dt(peri) <- k12 * cent - k21 * peri
      cp <- 1000 * cent / v1
      d / dt(ce) <- ke0 * (cp - ce)
      resp <- 100 * ce / (ec50 + ce)
      cp ~ lnorm(lnorm.sd)
      resp ~ add(resp.sd)
    })
  }
  .d <- do.call(rbind, lapply(1:3, function(id) {
    rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 100, EVID = 1,
                     CMT = "depot"),
          data.frame(ID = id, TIME = c(0.5, 2, 8, 24),
                     DV = c(400, 900, 500, 150), AMT = 0, EVID = 0,
                     CMT = "cp"),
          data.frame(ID = id, TIME = c(1, 4, 12, 24),
                     DV = c(40, 70, 60, 30), AMT = 0, EVID = 0,
                     CMT = "resp"))
  }))
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.effCpt, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})
