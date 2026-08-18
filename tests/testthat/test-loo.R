## log_lik is only worth emitting if it agrees with what was actually fitted,
## so these tests check the two against each other rather than just checking
## that loo() runs.
library(testthat)

oneCmt <- function() {
  ini({
    tka <- 0.5
    tcl <- 1.0
    tv <- 3.4
    eta.cl ~ 0.09
    add.sd <- 0.2
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

simData <- function(nsub = 6L, seed = 42) {
  times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = "depot")
    e <- rxode2::et(e, times)
    d <- as.data.frame(e)
    d$id <- i
    d$dv <- 0
    d
  }))
  set.seed(seed)
  etaCl <- stats::rnorm(nsub, 0, 0.3)
  m <- rxode2::rxode2("
ka <- exp(tka)
cl <- exp(tcl + eta_cl)
v  <- exp(tv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - cl / v * center
cp <- center / v
")
  pm <- cbind(tka = 0.5, tcl = 1.0, tv = 3.4, eta_cl = etaCl)
  s <- rxode2::rxSolve(m, params = pm, events = dat,
                       returnType = "data.frame", cores = 1L)
  dat$dv[dat$evid == 0] <- s$cp + stats::rnorm(nrow(s), 0, 0.2)
  dat
}

test_that("the generated program always emits log_lik", {
  skip_if_not_installed("nlmixr2")
  gen <- rxsStanFromUi(oneCmt, simData(3L))
  on.exit(rxsRelease(gen$handle))

  expect_match(gen$code, "generated quantities \\{")
  expect_match(gen$code, "vector\\[nObs\\] log_lik;")

  ## pred must live in transformed parameters, or generated quantities would
  ## have to solve the system a second time.
  tp <- sub(".*transformed parameters \\{", "", gen$code)
  tp <- sub("\\nmodel \\{.*", "", tp)
  expect_match(tp, "vector\\[nObs\\] pred;")
  expect_match(tp, "rx_solve\\(handle, p\\)")

  ## The model block and log_lik must use the SAME density expression.
  expect_equal(
    length(gregexpr("normal_lpdf(dv[i] | pred[i], add_sd)", gen$code,
                    fixed = TRUE)[[1]]),
    2L)
})

test_that("rxsInit starts every chain from the ini() values", {
  skip_if_not_installed("nlmixr2")
  gen <- rxsStanFromUi(oneCmt, simData(3L))
  on.exit(rxsRelease(gen$handle))

  expect_equal(gen$inits$tka, 0.5)
  expect_equal(gen$inits$tcl, 1.0)
  expect_equal(gen$inits$omega_eta_cl, sqrt(0.09))
  expect_equal(gen$inits$add_sd, 0.2)
  expect_equal(dim(gen$inits$z), c(1L, 3L))
  expect_true(all(gen$inits$z == 0))

  f <- rxsInit(gen, jitter = 0)
  expect_equal(f(1L)$tka, 0.5)

  set.seed(1)
  g <- rxsInit(gen, jitter = 0.5)
  expect_false(identical(g(1L)$tka, g(2L)$tka))
  ## jitter must not disturb the positively constrained parameters
  expect_equal(g(1L)$add_sd, 0.2)
  expect_equal(g(1L)$omega_eta_cl, sqrt(0.09))
})

test_that("log_lik matches the density the model actually used", {
  skipUnlessStan()

  dat <- simData(6L)
  gen <- rxsStanFromUi(oneCmt, dat)
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_loglik")

  s <- rstan::sampling(sm, data = gen$standata, chains = 1, iter = 400,
                       warmup = 200, seed = 3, refresh = 0,
                       init = rxsInit(gen, jitter = 0.1))

  ll <- rstan::extract(s, "log_lik")$log_lik
  pred <- rstan::extract(s, "pred")$pred
  sd <- rstan::extract(s, "add_sd")$add_sd

  expect_equal(dim(ll), c(nrow(pred), gen$standata$nObs))
  expect_true(all(is.finite(ll)))

  ## Recompute the density in R from the saved pred and sd.
  manual <- t(vapply(seq_len(nrow(ll)), function(d) {
    stats::dnorm(gen$standata$dv, pred[d, ], sd[d], log = TRUE)
  }, numeric(ncol(ll))))
  expect_equal(ll, manual, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("loo works on the generated program", {
  skipUnlessStan()
  skip_if_not_installed("loo")

  dat <- simData(6L)
  gen <- rxsStanFromUi(oneCmt, dat)
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_loo")

  s <- rstan::sampling(sm, data = gen$standata, chains = 2, iter = 600,
                       warmup = 300, seed = 4, refresh = 0,
                       init = rxsInit(gen, jitter = 0.1))

  l <- loo::loo(loo::extract_log_lik(s, "log_lik"))
  expect_true(is.finite(l$estimates["elpd_loo", "Estimate"]))
})

test_that("optimizing() and ADVI run through the bridge", {
  skipUnlessStan()

  dat <- simData(6L)
  gen <- rxsStanFromUi(oneCmt, dat)
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_opt")

  ## Both are first-order, so precomputed_gradients is enough for them.
  o <- rstan::optimizing(sm, data = gen$standata, seed = 5,
                         init = rxsInit(gen, jitter = 0)(1L), verbose = FALSE)
  expect_equal(o$return_code, 0L)
  expect_true(is.finite(o$value))
  expect_equal(unname(o$par["tcl"]), 1.0, tolerance = 0.5)

  v <- rstan::vb(sm, data = gen$standata, seed = 5, iter = 2000,
                 output_samples = 200, init = rxsInit(gen, jitter = 0)(1L))
  expect_s4_class(v, "stanfit")
  expect_true(is.finite(mean(rstan::extract(v, "tcl")$tcl)))
})
