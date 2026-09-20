library(testthat)
library(nlmixr2bayes)

# -------------------------------------------------------------------------
# Test policy
#
# Compiling dominates this suite's wall time: every est="stan"/est="nuts" test
# does a real Stan compile (~1-2 min) the first time a given generated-program
# shape is seen, and the hand-coded .stan bridge tests compile a full Stan
# program each.  The gates below exist so CI can split that cost up instead of
# losing a job to a hosted runner reclaiming it mid-run.
#
#   * NLMIXR2BAYES_TEST_CRAN_ONLY=true -- run only the CRAN-visible subset by
#     forcing NOT_CRAN=false.  Most files here are behind skip_on_cran(), so
#     this is quick; it is what the macOS/Windows/devel/oldrel jobs use, where
#     the point is proving the package BUILDS and checks (install, examples,
#     Rd, vignettes, compiled code), not re-running the sampler.
#   * NLMIXR2BAYES_TEST_SHARD="i/n" -- split the file list across n jobs.
#   * testthat workers: SERIAL on CI or CRAN.  Nothing here is parallel-safe to
#     assume otherwise, and rxode2/nlmixr2est already parallelize inside a
#     single likelihood evaluation.
#   * rxode2 within-solve threads: capped to 2 on CRAN, per CRAN's two-core
#     policy; on CI and locally rxode2 manages its own threads.
#
# The two gates that stay OFF unless a workflow sets them (see
# .github/workflows/slow-tests.yaml) are RXSTAN_STAN_TESTS (the hand-coded
# .stan compile suite, via skipUnlessStan() in helper-rxstan.R) and
# NLMIXR2STAN_SLOW (G5 simulation-based calibration, G6 FOCEi agreement).
# -------------------------------------------------------------------------
if (isTRUE(as.logical(Sys.getenv("NLMIXR2BAYES_TEST_CRAN_ONLY", "false")))) {
  Sys.setenv("NOT_CRAN" = "false") # nolint
}

.onCran <- !identical(Sys.getenv("NOT_CRAN"), "true")
.onCI <- isTRUE(as.logical(Sys.getenv("CI", "false")))

if (.onCI || .onCran) {
  options(Ncpus = 1L)
  Sys.setenv(TESTTHAT_CPUS = "1")
  Sys.setenv(TESTTHAT_PARALLEL = "FALSE")
}
if (.onCran) {
  rxode2::setRxThreads(2L)
}

.filter <- NULL
.shard <- Sys.getenv("NLMIXR2BAYES_TEST_SHARD")
if (nzchar(.shard)) {
  .p <- strsplit(.shard, "/", fixed = TRUE)[[1]]
  .i <- suppressWarnings(as.integer(.p[1]))
  .n <- suppressWarnings(as.integer(.p[2]))
  if (length(.p) != 2L || is.na(.i) || is.na(.n) || .n < 1L || .i < 1L || .i > .n) {
    stop(sprintf("NLMIXR2BAYES_TEST_SHARD=%s must be \"i/n\" with 1 <= i <= n", .shard))
  }
  # Deal the files out round-robin over a SORTED list so the split is stable
  # across jobs and platforms, and every file lands in exactly one shard.
  .all <- sort(sub(
    "^test-",
    "",
    sub(
      "\\.R$",
      "",
      basename(Sys.glob(file.path("testthat", "test-*.R")))
    )
  ))
  # An empty list would make every shard match nothing and report a green run
  # having tested nothing at all -- fail loudly instead.
  if (length(.all) == 0L) {
    stop("NLMIXR2BAYES_TEST_SHARD set but no test files found in ./testthat")
  }
  .mine <- .all[seq_along(.all) %% .n == (.i %% .n)]
  .filter <- if (length(.mine) == 0L) "^$" else paste0("^(", paste(.mine, collapse = "|"), ")$")
}
# Locally (and on CRAN) .filter stays NULL -> run everything.

# stop_on_failure MUST stay TRUE (the default): with it FALSE, R CMD check
# never errors on a test failure, so `checking tests` passes no matter what.
test_check("nlmixr2bayes", filter = .filter, perl = TRUE)
