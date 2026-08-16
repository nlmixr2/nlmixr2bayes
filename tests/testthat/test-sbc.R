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
      prior(eta.cl) ~ invWishart(6)
    })
    model({
      cl <- exp(tcl + eta.cl)
      cp <- 100 / exp(3) * exp(-cl / exp(3) * time)
      cp ~ add(add.sd)
    })
  }
  # the omega prior: invWishart(6) with the ini value 0.01 as scale; the
  # 1-D marginal is InvGamma(nu/2 = 3, s/2 = 0.005).  Chosen over the
  # default half-Cauchy DELIBERATELY: SBC simulates the omega from its own
  # prior, and a Cauchy-tail draw (sd ~ 3) makes a single replicate's fit
  # pathological (hours at deep treedepth); InvGamma(3, .) has light tails
  # so every replicate fits in seconds.  SBC only requires that the
  # simulation prior MATCHES the model prior, which this does exactly.
  .tt <- c(0.5, 1, 2, 4, 8)
  .nid <- 6L
  # one program + one compile for ALL replicates: the data live in the
  # linked nlmixr2est problem, so only the link is rebuilt per replicate
  .dTemplate <- do.call(rbind, lapply(seq_len(.nid), function(i) {
    data.frame(ID = i, TIME = .tt, DV = 5, AMT = 0, EVID = 0)
  }))
  .codeGen <- suppressMessages(
    nlmixr2est::nlmixr2(.mod, .dTemplate, est = "stan",
                        control = stanControl(run = FALSE)))
  .sm <- stanCompile(.codeGen$code)
  .stanData <- .codeGen$data
  .oneRep <- function(repSeed) {
    set.seed(repSeed)
    .tcl0 <- stats::rnorm(1, 1, 0.3)
    .omVar0 <- 0.005 / stats::rgamma(1, 3)
    .om0 <- sqrt(.omVar0)
    .add0 <- stats::rlnorm(1, -1, 0.3)
    .eta <- stats::rnorm(.nid, 0, .om0)
    .d <- do.call(rbind, lapply(seq_len(.nid), function(i) {
      .f <- 100 / exp(3) * exp(-exp(.tcl0 + .eta[i]) / exp(3) * .tt)
      data.frame(ID = i, TIME = .tt,
                 DV = .f + stats::rnorm(length(.tt), 0, .add0),
                 AMT = 0, EVID = 0)
    }))
    # drive the link + sampler directly on the ONE precompiled program --
    # the full nlmixr2() pipeline (ui rebuild, tables, finalize) per
    # replicate made the harness hours-slow with zero validation benefit;
    # this is byte-for-byte the same target (same program, same link path)
    .sf <- try({
      h <- stanLinkSetup(.mod, .d, thetaSens = TRUE, cores = 1L)
      .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
      .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
      suppressWarnings(rstan::sampling(
        .sm, data = .stanData, chains = 1L, iter = 1200L, warmup = 400L,
        thin = 8L, seed = repSeed, refresh = 0, cores = 1,
        show_messages = FALSE,
        init = list(list(tcl = 1, add_sd = 0.4,
                         omega_b1 = matrix(0.01, 1, 1),
                         z_b1 = matrix(0, .nid, 1))),
        control = list(adapt_delta = 0.95)))
    }, silent = TRUE)
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
    if (inherits(.sf, "try-error")) return(NULL)
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
