# Compile the generated Stan program against the linked likelihood via
# rstan::stan_model(allow_undefined=TRUE, includes=), cached by content hash.

# In-session memo of already-compiled stanmodel objects, keyed by the same
# content hash as the on-disk cache.  A stanmodel round-tripped through
# saveRDS()/readRDS() -- even within the SAME session -- carries a
# `.CXXDSOMISC$module` that is a distinct, serialization-broken copy; when
# the model's DLL happens to already be loaded in-process (true whenever
# THIS session compiled it), rstan's own `is_dso_loaded()` guard reports
# TRUE and skips reloading, so that stale copy's dead external pointer
# reaches `Module(mustStart = TRUE)` and crashes with "NULL value passed
# for DllInfo".  Returning the original live object -- never a
# freshly-read copy -- for anything this session already compiled avoids
# the bug outright, rather than merely detecting and recompiling around it.
.stanCompileEnv <- new.env(parent = emptyenv())

#' The phase-0 (eta-only) Stan program: theta and Omega fixed, Stan samples
#' the etas with the full, normalized eta prior on the Stan side.  The
#' observations never enter Stan's data block -- they live in nlmixr2est's
#' loaded problem.
#' @noRd
.stanPhase0Code <- function() {
  paste(
        "functions {",
        "  // external (allow_undefined): value + analytic d/d(eta) supplied by",
        "  // the linked rxode2/nlmixr2est likelihood via nlmixr2bayes_lp.hpp",
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

#' Whether a cached `stanmodel` object's compiled module can actually be
#' loaded.  A `stanmodel` read back via [readRDS()] carries a `dso` slot
#' whose C++ module pointer is only meaningful in the session that compiled
#' it (external pointers do not survive (de)serialization); rstan is
#' normally supposed to transparently reload the module from the embedded
#' DSO bytes on first use, but that reload can itself fail (e.g. a session
#' near R's loaded-DLL limit), surfacing as `Error in Module(module,
#' mustStart = TRUE): ... NULL value passed for DllInfo` deep inside
#' [rstan::sampling()].  Probing with the same call rstan uses internally
#' catches a stale/corrupt cache entry before it reaches the user.
#' @noRd
.stanModelUsable <- function(m) {
  isTRUE(methods::is(m, "stanmodel")) &&
    isTRUE(tryCatch({
      rstan:::mk_cppmodule(m)
      TRUE
    }, error = function(e) FALSE))
}

#' Path to the injected external-function header
#' @noRd
.stanLpHeader <- function() {
  .h <- system.file("include", "nlmixr2bayes_lp.hpp", package = "nlmixr2bayes")
  if (!nzchar(.h)) stop("nlmixr2bayes_lp.hpp not found", call. = FALSE) # nocov
  # forward slashes: rstan uses includes= as a sub() replacement, where
  # backslashes are escape characters
  normalizePath(.h, winslash = "/")
}

#' Compile a Stan program against the linked likelihood
#'
#' Wraps [rstan::stan_model()] with `allow_undefined=TRUE` and the
#' `nlmixr2bayes_lp.hpp` header injected via `includes=`.  Compiled models are
#' cached by a hash of the program, the header, and the rstan/nlmixr2bayes
#' versions.
#'
#' @param code Stan program text (default: the phase-0 eta-only program)
#' @param cache use (and populate) the on-disk cache
#' @param cacheDir cache directory (default
#'   `tools::R_user_dir("nlmixr2bayes", "cache")`)
#' @param verbose show the compiler output
#' @return an `rstan::stanmodel`
#' @export
#' @author Matthew L Fidler
stanCompile <- function(code = .stanPhase0Code(), cache = TRUE,
                        cacheDir = NULL, verbose = FALSE) {
  rxode2::rxReq("rstan")
  .stanAssertBuildOk()
  .hpp <- .stanLpHeader()
  .key <- digest::digest(list(code, readLines(.hpp),
                              as.character(utils::packageVersion("rstan")),
                              as.character(utils::packageVersion("nlmixr2bayes")),
                              R.version.string))
  if (cache && exists(.key, envir = .stanCompileEnv, inherits = FALSE)) {
    return(get(.key, envir = .stanCompileEnv, inherits = FALSE))
  }
  .dir <- if (is.null(cacheDir)) tools::R_user_dir("nlmixr2bayes", "cache") else cacheDir
  .rds <- file.path(.dir, paste0("stanmodel-", .key, ".rds"))
  if (cache && file.exists(.rds)) {
    .m <- tryCatch(readRDS(.rds), error = function(e) NULL)
    if (!is.null(.m) && .stanModelUsable(.m)) {
      assign(.key, .m, envir = .stanCompileEnv)
      return(.m)
    }
    if (!is.null(.m)) {
      # a stale/corrupt cache entry: fall through and recompile rather than
      # handing back a stanmodel that will crash sampling() with the
      # DllInfo error above
      unlink(.rds)
    }
  }
  # rstan::stan_model LEAKS compiler environment variables into the session
  # (PKG_CPPFLAGS with a forced `-include .../Eigen.hpp`, PKG_LIBS with the
  # Stan libraries, USE_CXX17) and never restores them, which then breaks
  # every subsequent rxode2 C model compile in the session ("compilation
  # terminated": gcc compiling C with a force-included C++ header).  est=
  # "stan" compiles rxode2 models AFTER this point (the link), so snapshot
  # and restore around the Stan compile.
  .leak <- c("PKG_CPPFLAGS", "PKG_LIBS", "USE_CXX17", "PKG_CXXFLAGS")
  .old <- Sys.getenv(.leak, unset = NA_character_)
  on.exit({
    for (.v in .leak) {
      if (is.na(.old[[.v]])) {
        Sys.unsetenv(.v)
      } else {
        do.call(Sys.setenv, stats::setNames(list(.old[[.v]]), .v))
      }
    }
  }, add = TRUE)
  # a fresh program shape takes 1-2 minutes to compile; say so, or the
  # user reasonably concludes the session froze
  cli::cli_inform("compiling Stan model (typically 1-2 minutes)...")
  .m <- rstan::stan_model(model_code = code, model_name = "nlmixr2bayes",
                          allow_undefined = TRUE,
                          includes = paste0("\n#include \"", .hpp, "\"\n"),
                          auto_write = FALSE, verbose = verbose)
  cli::cli_inform("done")
  if (cache) {
    assign(.key, .m, envir = .stanCompileEnv)
    dir.create(.dir, recursive = TRUE, showWarnings = FALSE)
    # write-then-rename: a second session compiling the same model shape
    # concurrently must never observe a partially-written .rds at .rds's
    # final path
    .tmp <- paste0(.rds, ".tmp", Sys.getpid())
    saveRDS(.m, .tmp)
    file.rename(.tmp, .rds)
  }
  .m
}
