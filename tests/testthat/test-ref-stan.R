# G4: the strongest single correctness test.  The same model hand-written
# NATIVELY in Stan -- analytic solution, no external function -- with the
# same priors and the same non-centred parameterization.  The linked fit's
# posterior must match the reference within Monte Carlo error on every
# shared parameter.  This validates the whole chain at once: codegen, prior
# emission, eta-prior ownership, the external likelihood and its gradients.
#
# The two programs deliberately share parameter names (tcl, tv, add_sd,
# sd_eta_cl), so the comparison reads straight off rstan::summary.  The linked
# conditional omits the per-observation -0.5*log(2*pi) constant; constants
# do not move a posterior, which is exactly what this test demonstrates.

.refStanCode <- function() {
  paste(
        "data {",
        "  int<lower=1> N;",
        "  int<lower=1> Nobs;",
        "  array[Nobs] int<lower=1, upper=N> id;",
        "  vector[Nobs] time;",
        "  vector[Nobs] dv;",
        "}",
        "parameters {",
        "  real tcl;",
        "  real tv;",
        "  real<lower=0> add_sd;",
        "  real<lower=0> sd_eta_cl;",
        "  matrix[N, 1] z_eta_cl;",
        "}",
        "model {",
        "  tcl ~ normal(1, 2);",
        "  tv ~ normal(3, 2);",
        "  add_sd ~ cauchy(0, 2.5);",
        # the exact default the generator announces: cauchy(0, 2.5*sqrt(0.1))
        paste0("  sd_eta_cl ~ cauchy(0, ",
               format(2.5 * sqrt(0.1), digits = 15, trim = TRUE,
                      scientific = FALSE), ");"),
        "  to_vector(z_eta_cl) ~ std_normal();",
        "  {",
        "    real v = exp(tv);",
        "    for (o in 1:Nobs) {",
        "      real cl = exp(tcl + sd_eta_cl * z_eta_cl[id[o], 1]);",
        "      real cp = 100 / v * exp(-cl / v * time[o]);",
        "      dv[o] ~ normal(cp, add_sd);",
        "    }",
        "  }",
        "}",
        sep = "\n")
}

test_that("G4: the linked posterior matches a hand-written pure-Stan reference", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  d <- .linkData()
  .iter <- if (nzchar(Sys.getenv("NLMIXR2STAN_SLOW"))) 10000L else 3000L
  .warm <- if (nzchar(Sys.getenv("NLMIXR2STAN_SLOW"))) 2000L else 1000L
  # ---- the linked fit ------------------------------------------------------
  .fit <- suppressWarnings(suppressMessages(
                                            nlmixr2est::nlmixr2(.estMod, d, est = "stan",
                                                                control = stanControl(chains = 2L, iter = .iter,
                                                                                      warmup = .warm, seed = 101L,
                                                                                      adapt_delta = 0.95,
                                                                                      cores = 1L,
                                                                                      onDiagnostic = "none"))))
  .sfL <- .fit$env$stanfit
  # ---- the reference fit ---------------------------------------------------
  .smR <- stanCompile(.refStanCode())
  .dataR <- list(N = 4L, Nobs = nrow(d), id = as.integer(d$ID),
                 time = d$TIME, dv = d$DV)
  .sfR <- rstan::sampling(.smR, data = .dataR, chains = 2, iter = .iter,
                          warmup = .warm, seed = 202, refresh = 0,
                          control = list(adapt_delta = 0.95))
  # ---- compare -------------------------------------------------------------
  .pars <- c("tcl", "tv", "add_sd", "sd_eta_cl")
  .sL <- rstan::summary(.sfL, pars = .pars)$summary
  .sR <- rstan::summary(.sfR, pars = .pars)$summary
  for (.p in .pars) {
    .tol <- max(3 * sqrt(.sL[.p, "se_mean"]^2 + .sR[.p, "se_mean"]^2), 0.03)
    expect_lt(abs(.sL[.p, "mean"] - .sR[.p, "mean"]), .tol,
              label = paste0("posterior mean of ", .p,
                             " (linked=", signif(.sL[.p, "mean"], 4),
                             " ref=", signif(.sR[.p, "mean"], 4), ")"))
    expect_lt(abs(.sL[.p, "sd"] - .sR[.p, "sd"]) / .sR[.p, "sd"], 0.2,
              label = paste0("posterior sd of ", .p))
  }
})
