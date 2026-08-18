## Multiple endpoints.  The bridge already returns every state, so one solve
## serves them all; what has to be right is that each observation row is
## scored against ITS endpoint's prediction and error model, and that log_lik
## branches the same way the model block does.
library(testthat)

pkpd <- function() {
  ini({
    tka <- 0.5
    tcl <- 1.0
    tv <- 3.4
    tec50 <- 1.0
    add.pk <- 0.2
    add.pd <- 0.1
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    ec50 <- exp(tec50)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    d/dt(eff) <- 1 - eff * (1 + (center / v) / (ec50 + center / v))
    eff(0) <- 1
    cp <- center / v
    resp <- eff
    cp ~ add(add.pk)
    resp ~ add(add.pd)
  })
}

## Parent and metabolite: two endpoints reading different states.
parentMetab <- function() {
  ini({
    tka <- 0.5
    tcl <- 1.0
    tv <- 3.4
    tclm <- 0.5
    add.p <- 0.2
    add.m <- 0.1
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    clm <- exp(tclm)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    d/dt(metab) <- cl / v * center - clm * metab
    parent <- center / v
    meta <- metab
    parent ~ add(add.p)
    meta ~ add(add.m)
  })
}

epData <- function(eps, nsub = 4L, times = c(1, 4, 8, 12, 24), seed = 3) {
  d <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = "depot")
    for (ep in eps) e <- rxode2::et(e, times, cmt = ep)
    x <- as.data.frame(e)
    x$id <- i
    x$dv <- 1
    x
  }))
  set.seed(seed)
  d
}

test_that("two endpoints each get their own prediction and error model", {
  skip_if_not_installed("nlmixr2")
  d <- epData(c("cp", "resp"))
  g <- rxsStanFromUi(pkpd, d)
  on.exit(rxsRelease(g$handle))

  ## Every observation row is assigned an endpoint, alternating here.
  expect_equal(length(g$standata$dvid), g$standata$nObs)
  expect_setequal(unique(g$standata$dvid), 1:2)
  expect_true(grepl("array[nObs] int<lower=1, upper=2> dvid;", g$code,
                    fixed = TRUE))

  ## Each endpoint scored with its own sd.
  expect_true(grepl("if (dvid[i] == 1) pred[i] = cp;", g$code, fixed = TRUE))
  expect_true(grepl("else if (dvid[i] == 2) pred[i] = resp;", g$code,
                    fixed = TRUE))
  expect_true(grepl("normal_lpdf(dv[i] | pred[i], add_pk)", g$code,
                    fixed = TRUE))
  expect_true(grepl("normal_lpdf(dv[i] | pred[i], add_pd)", g$code,
                    fixed = TRUE))

  ## log_lik must branch identically to the model block: each density appears
  ## exactly twice, once in each.
  for (sd in c("add_pk", "add_pd")) {
    hits <- gregexpr(sprintf("normal_lpdf(dv[i] | pred[i], %s)", sd), g$code,
                     fixed = TRUE)[[1]]
    expect_length(hits, 2L)
  }
})

test_that("a parent/metabolite model reads the states it should", {
  skip_if_not_installed("nlmixr2")
  d <- epData(c("parent", "meta"))
  g <- rxsStanFromUi(parentMetab, d)
  on.exit(rxsRelease(g$handle))

  expect_equal(g$states, c("depot", "center", "metab"))
  expect_true(grepl("if (dvid[i] == 1) pred[i] = parent;", g$code,
                    fixed = TRUE))
  expect_true(grepl("else if (dvid[i] == 2) pred[i] = meta;", g$code,
                    fixed = TRUE))
})

test_that("a single endpoint still emits the unbranched program", {
  skip_if_not_installed("nlmixr2")
  one <- function() {
    ini({ tka <- 0.5; tcl <- 1.0; tv <- 3.4; add.sd <- 0.2 })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  e <- rxode2::et(amt = 100, cmt = "depot")
  e <- rxode2::et(e, c(1, 4, 8, 12))
  d <- as.data.frame(e); d$id <- 1L; d$dv <- 1

  g <- rxsStanFromUi(one, d)
  on.exit(rxsRelease(g$handle))

  ## No dvid machinery at all when there is nothing to choose between.
  expect_false(grepl("dvid", g$code, fixed = TRUE))
  expect_null(g$standata$dvid)
  expect_true(grepl("        pred[i] = cp;", g$code, fixed = TRUE))
})

test_that("an endpoint the data does not identify is refused", {
  skip_if_not_installed("nlmixr2")

  ## cmt naming something the model has no endpoint for.
  d <- epData(c("cp", "resp"))
  d$cmt[d$cmt == "resp"] <- "effect"
  expect_error(rxsStanFromUi(pkpd, d), "endpoint the model does not have")

  ## No way at all to tell the endpoints apart.
  d2 <- epData(c("cp", "resp"))
  d2$cmt <- NULL
  expect_error(rxsStanFromUi(pkpd, d2), "dvid.*cmt|cmt.*dvid")
})

test_that("two endpoints recover their own residual error", {
  skipUnlessStan()

  d <- epData(c("cp", "resp"), nsub = 6L)

  ## Simulate from the truth so each endpoint has its own noise level.
  sim <- rxode2::rxode2("
ka <- exp(tka)
cl <- exp(tcl)
v  <- exp(tv)
ec50 <- exp(tec50)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - cl / v * center
d/dt(eff)    <-  1 - eff * (1 + (center / v) / (ec50 + center / v))
eff(0) <- 1
cp <- center / v
resp <- eff
")
  s <- rxode2::rxSolve(sim, params = c(tka = 0.5, tcl = 1.0, tv = 3.4,
                                       tec50 = 1.0),
                       events = d, returnType = "data.frame", cores = 1L)
  isObs <- d$evid == 0
  isPk <- isObs & d$cmt == "cp"
  set.seed(9)
  d$dv[isPk] <- s$cp[d$cmt[isObs] == "cp"] + stats::rnorm(sum(isPk), 0, 0.2)
  d$dv[isObs & d$cmt == "resp"] <-
    s$resp[d$cmt[isObs] == "resp"] +
    stats::rnorm(sum(isObs & d$cmt == "resp"), 0, 0.1)

  ## atol/rtol tighter than the 1e-8 default: the effect compartment's Emax
  ## feedback makes this system harder to integrate than a plain PK model, and
  ## at the default the gradient check sits at 1.4e-5.  Tightening to 1e-10
  ## brings it to 3.1e-8, so it is solver accuracy rather than a wrong
  ## derivative -- 1e-12 is no better (1.4e-7), which is the fd noise floor.
  g <- rxsStanFromUi(pkpd, d, atol = 1e-10, rtol = 1e-10)
  on.exit(rxsRelease(g$handle))
  sm <- stanModelFor(g$code, "rxstan_pkpd")

  fit <- rstan::sampling(sm, data = g$standata, chains = 0)
  set.seed(4)
  chk <- rxsCheckGradient(fit, stats::rnorm(rstan::get_num_upars(fit), 0, 0.2))
  expect_true(all(chk$relDiff < 1e-6),
              info = paste(utils::capture.output(print(chk)), collapse = "\n"))

  s2 <- rstan::sampling(sm, data = g$standata, chains = 2, iter = 700,
                        warmup = 350, seed = 11, refresh = 0,
                        init = rxsInit(g, jitter = 0.1))
  post <- as.matrix(s2, pars = c("add_pk", "add_pd"))
  m <- colMeans(post)
  cat("\n  add_pk (0.2):", round(m[["add_pk"]], 3),
      " add_pd (0.1):", round(m[["add_pd"]], 3), "\n")

  ## The point of separate endpoints: the two error terms are distinguished
  ## rather than pooled into one.
  expect_lt(abs(m[["add_pk"]] - 0.2), 0.1)
  expect_lt(abs(m[["add_pd"]] - 0.1), 0.05)
  expect_gt(m[["add_pk"]], m[["add_pd"]])
})
