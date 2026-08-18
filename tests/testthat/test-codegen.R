## Codegen is only worth having if what it emits is trustworthy, so the
## generated program is checked the same way the hand-written ones are:
## gradients against finite differences, and simulate-and-recover.
library(testthat)

oneCmt <- function() {
  ini({
    tka <- 0.0953
    tcl <- 1.386
    tv <- 3.401
    eta.ka ~ 0.09
    eta.cl ~ 0.0625
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

simData <- function(nsub = 5L, seed = 77) {
  times <- c(0.5, 1, 2, 4, 8, 12, 24)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = "depot")
    e <- rxode2::et(e, times)
    d <- as.data.frame(e)
    d$id <- i
    d$dv <- 0
    d
  }))

  set.seed(seed)
  theta <- c(tka = 0.0953, tcl = 1.386, tv = 3.401)
  omega <- c(0.3, 0.25, 0.2)
  eta <- matrix(stats::rnorm(3L * nsub), nrow = nsub, byrow = TRUE)
  eta <- sweep(eta, 2, omega, "*")

  m <- rxode2::rxode2("
ka  <- exp(tka + eta_ka)
cl  <- exp(tcl + eta_cl)
v   <- exp(tv + eta_v)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - cl / v * center
cp <- center / v
")
  pm <- cbind(tka = theta[["tka"]], tcl = theta[["tcl"]], tv = theta[["tv"]],
              eta_ka = eta[, 1], eta_cl = eta[, 2], eta_v = eta[, 3])
  s <- rxode2::rxSolve(m, params = pm, events = dat, returnType = "data.frame",
                       cores = 1L, atol = 1e-10, rtol = 1e-10)
  dat$dv[dat$evid == 0] <- s$cp + stats::rnorm(nrow(s), 0, 0.3)
  list(data = dat, theta = theta, omega = omega, eta = eta)
}

test_that("expressions translate to Stan or are refused outright", {
  expect_equal(nlmixr2bayes:::.rxsExprToStan(quote(a + b)), "(a + b)")
  expect_equal(nlmixr2bayes:::.rxsExprToStan(quote(exp(a))), "exp(a)")
  expect_equal(nlmixr2bayes:::.rxsExprToStan(quote(a^2)), "pow(a, 2)")
  expect_equal(nlmixr2bayes:::.rxsExprToStan(quote(eta.ka)), "eta_ka")
  expect_error(nlmixr2bayes:::.rxsExprToStan(quote(besselJ(a, b))), "unsupported")
})

test_that("codegen produces the layout the model implies", {
  skip_if_not_installed("nlmixr2")
  sim <- simData(3L)
  gen <- rxsStanFromUi(oneCmt, sim$data)
  on.exit(rxsRelease(gen$handle))

  expect_equal(gen$thetaNames, c("tka", "tcl", "tv"))
  expect_equal(gen$etaNames, c("eta.ka", "eta.cl", "eta.v"))
  expect_equal(gen$errNames, "add.sd")
  expect_equal(gen$stanNames$err, "add_sd")
  expect_equal(gen$states, c("depot", "center"))

  ## 3 thetas + 3 etas per subject
  expect_equal(attr(gen$handle, "nBlock"), 6L)
  expect_equal(attr(gen$handle, "nBlocks"), 3L)
  expect_equal(gen$standata$nSub, 3L)
  expect_equal(gen$standata$nObs, 21L)
})

## --- ini() bounds, fix() and prior overrides --------------------------------

boundedModel <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- fix(3.45)
    eta.ka ~ fix(0.6)
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - (cl / v) * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

boundedData <- function() {
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, c(0.5, 2, 6, 12))
  d <- as.data.frame(ev)
  d$id <- 1L
  d$dv <- 1
  d
}

test_that("ini() bounds become Stan constraints", {
  skip_if_not_installed("nlmixr2")
  gen <- rxsStanFromUi(boundedModel, boundedData())
  on.exit(rxsRelease(gen$handle))

  expect_match(gen$code, "real<upper=4\\.60517", fixed = FALSE)
  expect_match(gen$code, "real<lower=0> add_sd;", fixed = TRUE)
})

test_that("fix() removes the parameter from sampling, like a constant() prior", {
  skip_if_not_installed("nlmixr2")
  gen <- rxsStanFromUi(boundedModel, boundedData())
  on.exit(rxsRelease(gen$handle))

  lines <- strsplit(gen$code, "\n")[[1]]
  parBlock <- lines[which(lines == "parameters {"):which(lines == "transformed parameters {")]

  ## The fixed theta is a constant, never a tightly-constrained parameter.
  expect_false(any(grepl("\\btv\\b", parBlock)))
  expect_match(gen$code, "real tv = 3.45", fixed = TRUE)

  ## A fixed eta fixes its between-subject SD at sqrt(variance); z stays free.
  expect_false(any(grepl("omega_eta_ka", parBlock)))
  expect_match(gen$code, sprintf("real omega_eta_ka = %.17g", sqrt(0.6)),
               fixed = TRUE)
  expect_match(gen$code, "eta[1, s] = omega_eta_ka * z[1, s];", fixed = TRUE)

  ## And it drops out of the sensitivity system: tka, tcl + 2 etas, not 5.
  expect_equal(attr(gen$handle, "nBlock"), 4L)
  expect_equal(gen$thetaNames, c("tka", "tcl"))
})

test_that("priors can be overridden by name, and typos are refused", {
  skip_if_not_installed("nlmixr2")
  gen <- rxsStanFromUi(boundedModel, boundedData(),
                       priors = list(tka = "normal(0, 2)",
                                     add.sd = "exponential(1)",
                                     eta.v = "cauchy(0, 1)"))
  on.exit(rxsRelease(gen$handle))

  expect_match(gen$code, "tka ~ normal(0, 2);", fixed = TRUE)
  expect_match(gen$code, "add_sd ~ exponential(1);", fixed = TRUE)
  expect_match(gen$code, "omega_eta_v ~ cauchy(0, 1);", fixed = TRUE)

  expect_error(rxsStanFromUi(boundedModel, boundedData(),
                             priors = list(nosuchpar = "normal(0, 1)")),
               "unknown parameter")
})

test_that("priorSd controls the default theta prior", {
  skip_if_not_installed("nlmixr2")
  gen <- rxsStanFromUi(boundedModel, boundedData(), priorSd = 2.5)
  on.exit(rxsRelease(gen$handle))
  expect_match(gen$code, "tka ~ normal(0.45", fixed = TRUE)
  expect_match(gen$code, ", 2.5);", fixed = TRUE)
})

test_that("a model with a fixed parameter still has correct gradients", {
  skipUnlessStan()

  gen <- rxsStanFromUi(boundedModel, boundedData())
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_fixed")
  fit <- rstan::sampling(sm, data = gen$standata, chains = 0)

  set.seed(11)
  u <- stats::runif(rstan::get_num_upars(fit), -0.5, 0.5)
  chk <- rxsCheckGradient(fit, u)
  expect_true(all(chk$relDiff < 1e-4),
              info = paste(utils::capture.output(print(chk)), collapse = "\n"))
})

test_that("the generated program has correct gradients and recovers the truth", {
  skipUnlessStan()

  sim <- simData(5L)
  gen <- rxsStanFromUi(oneCmt, sim$data)
  on.exit(rxsRelease(gen$handle))

  sm <- stanModelFor(gen$code, "rxstan_gen")

  fit <- rstan::sampling(sm, data = gen$standata, chains = 0)

  set.seed(4)
  npar <- rstan::get_num_upars(fit)
  for (trial in 1:3) {
    u <- stats::runif(npar, -0.6, 0.6)
    chk <- rxsCheckGradient(fit, u)
    expect_true(all(chk$relDiff < 1e-4),
                info = paste(utils::capture.output(
                  print(chk[order(-chk$relDiff), ][1:3, ])), collapse = "\n"))
  }

  init <- list(list(theta = as.numeric(sim$theta), omega = sim$omega,
                    add_sd = 0.3,
                    z = t(sweep(sim$eta, 2, sim$omega, "/"))))
  s <- rstan::sampling(sm, data = gen$standata, chains = 1, iter = 600,
                       warmup = 300, seed = 8, refresh = 0, init = init)

  post <- as.matrix(s, pars = "theta")
  ci <- apply(post, 2, stats::quantile, probs = c(0.025, 0.975))
  print(ci)

  expect_gt(length(unique(post[, 1])), nrow(post) / 2)
  for (k in 1:3) {
    expect_gte(sim$theta[[k]], ci[1, k])
    expect_lte(sim$theta[[k]], ci[2, k])
  }
})
