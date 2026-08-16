# G14: prior-only recovery.  With the likelihood term stripped from the
# generated program, every sampled parameter's marginal must equal its
# declared prior (one-sample KS at alpha = 0.01/nPar).  This isolates the
# prior codegen -- support promotion, truncation branches, the default
# LKJ + half-Cauchy omega path, name mangling -- completely from the
# likelihood, so when a posterior-level gate (G4/G6) fails it is immediately
# clear which half is broken.

test_that("prior-only sampling recovers every declared prior (G14)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .mod <- function() {
    ini({
      tcl <- 1
      fbio <- c(0, 0.4, 1)
      add.sd <- c(0, 0.5)
      eta.cl + eta.v ~ c(0.1,
                         0.02, 0.1)
      prior(tcl) ~ dnorm(1, 2)
      prior(fbio) ~ dbeta(2, 3)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(3 + eta.v)
      cp <- fbio * 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .d <- .linkData()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.mod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  # neutralize the ONE likelihood evaluation (the model block computes
  # llCond once and adds sum(llCond)); everything else stays
  # byte-identical -- the iteration tick's external call is a no-op when
  # printing is not armed, so it can stay
  .w <- grep("] llCond = nlmixr2_cond_all2(eta, theta);", .lines,
             fixed = TRUE)
  expect_length(.w, 1L)
  .lines[.w] <- '    vector[N] llCond = rep_vector(0, N);' 
  # generated quantities would call the (unlinked) external function
  .w2 <- grep("logLikSubj = nlmixr2_cond_all2", .lines, fixed = TRUE)
  if (length(.w2) == 1L) .lines[.w2] <- "  logLikSubj = rep_vector(0, N);"
  .sm <- stanCompile(paste(.lines, collapse = "\n"))
  .sf <- rxode2::rxWithSeed(99, {
    rstan::sampling(.sm, data = .code$data, chains = 2L, iter = 6000L,
                    warmup = 1000L, thin = 5L, seed = 99L, refresh = 0L)
  })
  .alpha <- 0.01 / 6
  .ks <- function(x, F) suppressWarnings(stats::ks.test(x, F))$p.value
  .ex <- rstan::extract(.sf)
  # thetas: plain normal, unit-support beta, half-Cauchy (truncation branch
  # (b): constraint declared, no T[] -- the dropped normalizer is constant)
  expect_gt(.ks(.ex$tcl, function(q) stats::pnorm(q, 1, 2)), .alpha)
  expect_gt(.ks(.ex$fbio, function(q) stats::pbeta(q, 2, 3)), .alpha)
  expect_gt(.ks(.ex$add_sd,
                function(q) 2 * (stats::pcauchy(q, 0, 2.5) - 0.5)), .alpha)
  # omega block (default LKJ(2) + half-Cauchy SDs, scale 2.5*sqrt(diag init))
  .om <- .ex$omegaOut
  .s <- 2.5 * sqrt(0.1)
  .hc <- function(q) 2 * (stats::pcauchy(q, 0, .s) - 0.5)
  expect_gt(.ks(sqrt(.om[, 1, 1]), .hc), .alpha)
  expect_gt(.ks(sqrt(.om[, 2, 2]), .hc), .alpha)
  # LKJ(2), k=2: the correlation is 2*Beta(2,2) - 1
  .r <- .om[, 1, 2] / sqrt(.om[, 1, 1] * .om[, 2, 2])
  expect_gt(.ks(.r, function(q) stats::pbeta((q + 1) / 2, 2, 2)), .alpha)
})
