# G5: simulation-based calibration.  For each replicate, draw
# (tcl, omega, add.sd) FROM THE PRIORS, simulate data from that draw, fit
# with est="stan", and record the rank of the prior draw among L thinned
# posterior draws.  If (and only if) the sampler targets the correct
# posterior, the ranks are uniform on {0..L} -- this is the whole-pipeline
# correctness gate that catches a mis-normalized prior, a wrong Jacobian,
# or a biased target that every fixed-truth test can miss.
#
# NLMIXR2STAN_SLOW runs a smoke-scale N (default 40, ~5-10 min after the
# compile); the nightly-scale run sets NLMIXR2STAN_SBC_N=1000.  Replicates
# whose fit diverges are redrawn (a discard rate above 20% fails the gate
# -- at nightly scale the plan's threshold is 5%).

test_that("simulation-based calibration: ranks are uniform (G5)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not(nzchar(Sys.getenv("NLMIXR2STAN_SLOW")),
              "set NLMIXR2STAN_SLOW=TRUE for the SBC gate")
  .n <- as.integer(Sys.getenv("NLMIXR2STAN_SBC_N", "40"))
  .mod <- function() {
    ini({
      tcl <- 1
      add.sd <- c(0, 0.4)
      eta.cl ~ 0.01
      prior(tcl) ~ dnorm(1, 0.3)
      prior(add.sd) ~ dlnorm(-1, 0.3)
    })
    model({
      cl <- exp(tcl + eta.cl)
      cp <- 100 / exp(3) * exp(-cl / exp(3) * time)
      cp ~ add(add.sd)
    })
  }
  # the default omega prior the generator announces for this model:
  # sd ~ half-Cauchy(0, 2.5 * sqrt(0.01))
  .sdScale <- 2.5 * sqrt(0.01)
  .tt <- c(0.5, 1, 2, 4, 8)
  .nid <- 6L
  .oneRep <- function(repSeed) {
    set.seed(repSeed)
    .tcl0 <- stats::rnorm(1, 1, 0.3)
    .om0 <- abs(stats::rcauchy(1, 0, .sdScale))
    .add0 <- stats::rlnorm(1, -1, 0.3)
    .eta <- stats::rnorm(.nid, 0, .om0)
    .d <- do.call(rbind, lapply(seq_len(.nid), function(i) {
      .f <- 100 / exp(3) * exp(-exp(.tcl0 + .eta[i]) / exp(3) * .tt)
      data.frame(ID = i, TIME = .tt,
                 DV = .f + stats::rnorm(length(.tt), 0, .add0),
                 AMT = 0, EVID = 0)
    }))
    .fit <- try(suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
      .mod, .d, est = "stan",
      control = stanControl(chains = 1L, iter = 1200L, warmup = 400L,
                            thin = 8L, seed = repSeed, likCores = 1L,
                            adapt_delta = 0.95, calcTables = FALSE,
                            ofv = "none", onDiagnostic = "none")))),
      silent = TRUE)
    if (inherits(.fit, "try-error")) return(NULL)
    .sf <- .fit$env$stanfit
    .sp <- rstan::get_sampler_params(.sf, inc_warmup = FALSE)
    if (sum(vapply(.sp, function(x) sum(x[, "divergent__"]),
                   numeric(1))) > 0) {
      return(NULL) # divergent replicate: discard and redraw
    }
    .ex <- rstan::extract(.sf, pars = c("tcl", "add_sd", "omegaOut"))
    c(tcl = sum(.ex$tcl < .tcl0),
      add.sd = sum(.ex$add_sd < .add0),
      omSd = sum(sqrt(.ex$omegaOut[, 1, 1]) < .om0),
      L = length(.ex$tcl))
  }
  .ranks <- list()
  .discard <- 0L
  .seed <- 5000L
  while (length(.ranks) < .n) {
    .seed <- .seed + 1L
    .r <- .oneRep(.seed)
    if (is.null(.r)) {
      .discard <- .discard + 1L
      next
    }
    .ranks[[length(.ranks) + 1L]] <- .r
  }
  .rk <- do.call(rbind, .ranks)
  .L <- .rk[1, "L"]
  cat(sprintf("\nSBC: %d replicates, %d discarded (%.0f%%), L = %d\n",
              .n, .discard, 100 * .discard / (.n + .discard), .L))
  expect_lt(.discard / (.n + .discard), 0.20)
  # uniformity per parameter: randomized-PIT KS against U(0,1)
  set.seed(1)
  for (.p in c("tcl", "add.sd", "omSd")) {
    .u <- (.rk[, .p] + stats::runif(.n)) / (.L + 1)
    .ks <- suppressWarnings(stats::ks.test(.u, "punif"))
    cat(sprintf("SBC %s: KS p = %.4f\n", .p, .ks$p.value))
    expect_gt(.ks$p.value, 0.001)
  }
})
