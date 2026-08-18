## End-to-end hierarchical model: population parameters plus per-subject random
## effects, with rxode2 solving every subject in one call and supplying the
## block-diagonal sensitivities.
library(testthat)

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
"
sens <- c("lka", "lcl", "lv")
theta <- c(lka = log(1.1), lcl = log(4), lv = log(30))
omega <- c(0.3, 0.25, 0.2)

popFixture <- function(nsub = 6L) {
  times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
  ev <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = "depot")
    e <- rxode2::et(e, times)
    d <- as.data.frame(e)
    d$id <- i
    d
  }))

  h <- rxsRegister(pkModel, events = ev, sens = sens, output = "center",
                   perSubject = TRUE, atol = 1e-10, rtol = 1e-10)

  set.seed(2024)
  eta <- matrix(stats::rnorm(3L * nsub), nrow = nsub, byrow = TRUE)
  phi <- sweep(sweep(eta, 2, omega, "*"), 2, theta, "+")
  p <- as.numeric(t(phi))

  cp <- rxsSolve(h, p)[, 1L] / rep(exp(phi[, 3]), each = length(times))
  cpObs <- exp(log(cp) + stats::rnorm(length(cp), 0, 0.1))

  list(handle = h, nsub = nsub, phi = phi, p = p, cpObs = cpObs,
       subj = rep(seq_len(nsub), each = length(times)))
}

test_that("a population model samples and recovers the population parameters", {
  skipUnlessStan(nlmixr2 = FALSE)

  f <- popFixture()
  on.exit(rxsRelease(f$handle))

  sm <- stanModelFor(system.file("stan", "pk_1cmt_oral_pop.stan",
                                 package = "nlmixr2bayes"),
                     "rxstan_pk_pop")

  dat <- list(nSub = f$nsub, nObs = length(f$cpObs),
              handle = as.integer(unclass(f$handle)),
              cpObs = f$cpObs, subj = f$subj)

  fit0 <- rstan::sampling(sm, data = dat, chains = 0)

  ## Gradients first: if the block-diagonal wiring were wrong, the random
  ## effects would get each other's partials and this is where it shows.
  upars <- c(theta, omega, as.numeric(t(sweep(f$phi, 2, theta) / rep(omega, each = f$nsub))),
             log(0.1))
  upars[4:6] <- log(omega)
  chk <- rxsCheckGradient(fit0, upars)
  print(utils::head(chk[order(-chk$relDiff), ], 5))
  expect_true(all(chk$relDiff < 1e-5))

  ## Started at the truth: the one-compartment oral profile is bimodal
  ## (flip-flop), which is a property of the PK model, not of this bridge.
  init <- list(list(theta = as.numeric(theta), omega = omega, sigma = 0.1,
                    z = t(sweep(f$phi, 2, theta) / rep(omega, each = f$nsub))))

  s <- rstan::sampling(sm, data = dat, chains = 1, iter = 600, warmup = 300,
                       seed = 5, refresh = 0, init = init)

  post <- as.matrix(s, pars = "theta")
  ci <- apply(post, 2, stats::quantile, probs = c(0.025, 0.975))
  print(ci)

  expect_gt(length(unique(post[, 1])), nrow(post) / 2)
  for (k in 1:3) {
    expect_gte(theta[[k]], ci[1, k])
    expect_lte(theta[[k]], ci[2, k])
  }
})
