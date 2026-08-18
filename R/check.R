#' Compile a Stan program against the rxstan bridge
#'
#' Thin wrapper over [rstan::stanc()] + [rstan::stan_model()] that arranges the
#' namespace surgery `includes=` needs and keeps the model name unobfuscated so
#' the namespace is predictable.
#'
#' @param file path to a `.stan` file declaring `vector rx_solve(int, vector);`
#' @param modelName C++-safe model name
#' @param tag,pre,post passed to [rxsStanIncludes()]
#' @param ... passed to [rstan::stan_model()]
#' @return a `stanmodel`
#' @author Lukas A. Widmer
#' @export
rxsStanModel <- function(file, modelName = "rxstan_model",
                         tag = NULL, pre = NULL, post = NULL, ...) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("rxsStanModel() needs rstan", call. = FALSE)
  }

  ## rstan::stan_model() leaves PKG_CPPFLAGS and PKG_LIBS pointing at the Stan
  ## headers for the rest of the session, which makes every later rxode2 model
  ## build fail (its generated C picks up Stan's C++ Eigen header).  Put them
  ## back so registering another handle after compiling still works.
  vars <- c("PKG_CPPFLAGS", "PKG_CXXFLAGS", "PKG_LIBS")
  saved <- Sys.getenv(vars, names = TRUE, unset = NA)
  on.exit({
    set <- saved[!is.na(saved)]
    if (length(set)) do.call(Sys.setenv, as.list(set))
    unset <- names(saved)[is.na(saved)]
    if (length(unset)) Sys.unsetenv(unset)
  }, add = TRUE)

  sc <- rstan::stanc(file = file, model_name = modelName,
                     allow_undefined = TRUE, obfuscate_model_name = FALSE)
  rstan::stan_model(stanc_ret = sc,
                    includes = rxsStanIncludes(modelName, tag = tag,
                                               pre = pre, post = post),
                    obfuscate_model_name = FALSE, ...)
}

#' Compare a Stan model's gradient against finite differences of its log density
#'
#' The core correctness check for the bridge: if rxode2's analytic
#' sensitivities are wired onto the tape correctly, `grad_log_prob()` and a
#' central difference of `log_prob()` agree.
#'
#' It is also the only guard against a subtler failure.  The value and the
#' gradient come from two *different* integrations -- the state equations and
#' the sensitivity equations -- so they agree only as far as the solver
#' tolerances make them.  Stan's own ODE solvers differentiate the actual
#' numerical trajectory, so there value and gradient are consistent by
#' construction even when both are inaccurate; here they are not.  An
#' inconsistency acts like a non-conservative force in the leapfrog integrator:
#' the Hamiltonian drifts, acceptance falls and the posterior can be biased
#' with nothing visibly wrong.  This function catches exactly that, because
#' `log_prob()` walks the value path and `grad_log_prob()` the sensitivity one.
#'
#' So re-run it whenever `atol`/`rtol` are loosened from the [rxsRegister()]
#' defaults, and treat a relative difference much worse than about `1e-6` as a
#' reason to tighten them again rather than as noise.
#'
#' @param fit a `stanfit` (one from `sampling(..., chains = 0)` is enough)
#' @param upars unconstrained parameter vector to test at
#' @param h finite-difference step.  Differencing a solve amplifies solver
#'   noise by `1/h`, so very small steps make the check worse, not better.
#' @return a data frame with the analytic gradient, the numeric gradient and
#'   their relative difference
#' @author Lukas A. Widmer
#' @export
rxsCheckGradient <- function(fit, upars, h = 1e-5) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("rxsCheckGradient() needs rstan", call. = FALSE)
  }
  lp <- function(u) rstan::log_prob(fit, u, adjust_transform = FALSE)

  analytic <- as.numeric(rstan::grad_log_prob(fit, upars,
                                              adjust_transform = FALSE))
  numeric <- vapply(seq_along(upars), function(i) {
    step <- h * max(1, abs(upars[i]))
    up <- upars; up[i] <- up[i] + step
    um <- upars; um[i] <- um[i] - step
    (lp(up) - lp(um)) / (2 * step)
  }, numeric(1))

  scale <- pmax(1, abs(analytic), abs(numeric))
  data.frame(par = seq_along(upars), analytic = analytic, numeric = numeric,
             relDiff = abs(analytic - numeric) / scale)
}
