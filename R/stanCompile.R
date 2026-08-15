# Compile the generated Stan program against the linked likelihood via
# rstan::stan_model(allow_undefined=TRUE, includes=), cached by content hash.

#' The phase-0 (eta-only) Stan program: theta and Omega fixed, Stan samples
#' the etas with the full, normalized eta prior on the Stan side.  The
#' observations never enter Stan's data block -- they live in nlmixr2est's
#' loaded problem.
#' @noRd
.stanPhase0Code <- function() {
  paste(
        "functions {",
        "  // external (allow_undefined): value + analytic d/d(eta) supplied by",
        "  // the linked rxode2/nlmixr2est likelihood via nlmixr2stan_lp.hpp",
        "  vector nlmixr2_cond_all(matrix etaMat);",
        "}",
        "data {",
        "  int<lower=1> N;              // subjects, in the link handle's idLvl order",
        "  int<lower=1> P;              // etas per subject",
        "  cholesky_factor_cov[P] L;    // Omega = L L', fixed in phase 0",
        "}",
        "parameters {",
        "  matrix[N, P] eta;",
        "}",
        "model {",
        "  // Stan owns the FULL, normalized eta prior; nlmixr2_cond_all returns",
        "  // the CONDITIONAL log p(y_i | eta_i) only, so nothing is double-counted",
        "  for (i in 1:N) {",
        "    eta[i]' ~ multi_normal_cholesky(rep_vector(0, P), L);",
        "  }",
        "  target += sum(nlmixr2_cond_all(eta));",
        "}",
        sep = "\n")
}

#' Path to the injected external-function header
#' @noRd
.stanLpHeader <- function() {
  .h <- system.file("include", "nlmixr2stan_lp.hpp", package = "nlmixr2stan")
  if (!nzchar(.h)) stop("nlmixr2stan_lp.hpp not found", call. = FALSE) # nocov
  # forward slashes: rstan uses includes= as a sub() replacement, where
  # backslashes are escape characters
  normalizePath(.h, winslash = "/")
}

#' Compile a Stan program against the linked likelihood
#'
#' Wraps [rstan::stan_model()] with `allow_undefined=TRUE` and the
#' `nlmixr2stan_lp.hpp` header injected via `includes=`.  Compiled models are
#' cached by a hash of the program, the header, and the rstan/nlmixr2stan
#' versions.
#'
#' @param code Stan program text (default: the phase-0 eta-only program)
#' @param cache use (and populate) the on-disk cache
#' @param cacheDir cache directory (default
#'   `tools::R_user_dir("nlmixr2stan", "cache")`)
#' @param verbose show the compiler output
#' @return an `rstan::stanmodel`
#' @export
#' @author Matthew L. Fidler
stanCompile <- function(code = .stanPhase0Code(), cache = TRUE,
                        cacheDir = NULL, verbose = FALSE) {
  rxode2::rxReq("rstan")
  .stanAssertBuildOk()
  .hpp <- .stanLpHeader()
  .key <- digest::digest(list(code, readLines(.hpp),
                              as.character(utils::packageVersion("rstan")),
                              as.character(utils::packageVersion("nlmixr2stan")),
                              R.version.string))
  .dir <- if (is.null(cacheDir)) tools::R_user_dir("nlmixr2stan", "cache") else cacheDir
  .rds <- file.path(.dir, paste0("stanmodel-", .key, ".rds"))
  if (cache && file.exists(.rds)) {
    .m <- tryCatch(readRDS(.rds), error = function(e) NULL)
    if (!is.null(.m)) return(.m)
  }
  .m <- rstan::stan_model(model_code = code, model_name = "nlmixr2stan",
                          allow_undefined = TRUE,
                          includes = paste0("\n#include \"", .hpp, "\"\n"),
                          auto_write = FALSE, verbose = verbose)
  if (cache) {
    dir.create(.dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(.m, .rds)
  }
  .m
}
