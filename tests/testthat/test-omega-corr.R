## Correlated between-subject variability: nlmixr2 records off-diagonal omega
## entries as their own iniDf rows, and the generated Stan uses an SD vector
## plus a Cholesky correlation factor.
library(testthat)

corModel <- function() {
  ini({
    tka <- 0.0953
    tcl <- 1.386
    tv <- 3.401
    eta.ka + eta.cl ~ c(0.09,
                        0.045, 0.0625)
    eta.v ~ 0.04
    add.sd <- 0.3
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

test_that("omega blocks are the connected components of the off-diagonals", {
  om <- matrix(0, 4, 4)
  diag(om) <- 1
  om[1, 2] <- om[2, 1] <- 0.3
  om[3, 4] <- om[4, 3] <- 0.2
  expect_equal(nlmixr2bayes:::.rxsOmegaBlocks(om), list(c(1L, 2L), c(3L, 4L)))

  diag4 <- diag(4)
  expect_equal(nlmixr2bayes:::.rxsOmegaBlocks(diag4), list(1L, 2L, 3L, 4L))

  ## A chain 1-2, 2-3 is one block of three, not two of two.
  chain <- diag(3)
  chain[1, 2] <- chain[2, 1] <- 0.1
  chain[2, 3] <- chain[3, 2] <- 0.1
  expect_equal(nlmixr2bayes:::.rxsOmegaBlocks(chain), list(c(1L, 2L, 3L)))
})

simCor <- function(nsub = 8L, seed = 5) {
  times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = "depot")
    e <- rxode2::et(e, times)
    d <- as.data.frame(e)
    d$id <- i
    d$dv <- 0
    d
  }))

  theta <- c(0.0953, 1.386, 3.401)
  omega <- matrix(c(0.09, 0.045, 0, 0.045, 0.0625, 0, 0, 0, 0.04), 3, 3)

  set.seed(seed)
  L <- chol(omega)
  eta <- matrix(stats::rnorm(3L * nsub), nrow = nsub) %*% L
  phi <- sweep(eta, 2, theta, "+")

  m <- rxode2::rxode2("
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - cl / v * center
cp <- center / v
")
  pm <- cbind(lka = phi[, 1], lcl = phi[, 2], lv = phi[, 3])
  s <- rxode2::rxSolve(m, params = pm, events = dat, returnType = "data.frame",
                       cores = 1L, atol = 1e-10, rtol = 1e-10)
  dat$dv[dat$evid == 0] <- s$cp + stats::rnorm(nrow(s), 0, 0.3)
  list(data = dat, theta = theta, omega = omega, eta = eta, nsub = nsub)
}

test_that("a correlated block becomes a Cholesky factor with an LKJ prior", {
  skip_if_not_installed("nlmixr2")
  sim <- simCor(3L)
  gen <- rxsStanFromUi(corModel, sim$data, lkjEta = 3)
  on.exit(rxsRelease(gen$handle))

  expect_match(gen$code, "vector<lower=0>[2] omega_blk1;", fixed = TRUE)
  expect_match(gen$code, "cholesky_factor_corr[2] L_blk1;", fixed = TRUE)
  expect_match(gen$code, "L_blk1 ~ lkj_corr_cholesky(3);", fixed = TRUE)
  expect_match(gen$code, "diag_pre_multiply(omega_blk1, L_blk1)", fixed = TRUE)

  ## The uncorrelated eta keeps the plain scalar treatment.
  expect_match(gen$code, "real<lower=0> omega_eta_v;", fixed = TRUE)
  expect_match(gen$code, "eta[3, s] = omega_eta_v * z[3, s];", fixed = TRUE)

  ## The correlation is reported, otherwise it cannot be read off the fit.
  expect_match(gen$code, "corr_blk1 = multiply_lower_tri_self_transpose(L_blk1)",
               fixed = TRUE)

  ## Off-diagonal rows must not be mistaken for etas.
  expect_equal(gen$etaNames, c("eta.ka", "eta.cl", "eta.v"))
  expect_equal(attr(gen$handle, "nBlock"), 6L)
})

test_that("a correlated model has correct gradients", {
  skipUnlessStan()

  sim <- simCor(4L)
  gen <- rxsStanFromUi(corModel, sim$data)
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_corr")
  fit <- rstan::sampling(sm, data = gen$standata, chains = 0)

  set.seed(21)
  for (trial in 1:3) {
    u <- stats::runif(rstan::get_num_upars(fit), -0.5, 0.5)
    chk <- rxsCheckGradient(fit, u)
    expect_true(all(chk$relDiff < 1e-4),
                info = paste(utils::capture.output(
                  print(chk[order(-chk$relDiff), ][1:3, ])), collapse = "\n"))
  }
})

test_that("the correlation itself is recovered, not just the marginals", {
  skipUnlessStan()

  sim <- simCor(30L, seed = 12)
  gen <- rxsStanFromUi(corModel, sim$data)
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_corr_fit")

  s <- rstan::sampling(sm, data = gen$standata, chains = 1, iter = 700,
                       warmup = 350, seed = 3, refresh = 0)

  rho <- as.matrix(s, pars = "corr_blk1[1,2]")
  ci <- stats::quantile(rho, c(0.05, 0.95))
  truth <- sim$omega[1, 2] / sqrt(sim$omega[1, 1] * sim$omega[2, 2])
  cat("\ntrue rho:", truth, " 90% CI:", ci, "\n")

  expect_gte(truth, ci[[1]])
  expect_lte(truth, ci[[2]])
})
