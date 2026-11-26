library(testthat)
library(nlmixr2bayes)

# -------------------------------------------------------------------------
# Test policy
#
# Compiling dominates this suite's wall time: every est="stan"/est="nuts" test
# does a real Stan compile (~1-2 min) the first time a given generated-program
# shape is seen, and the hand-coded .stan bridge tests compile a full Stan
# program each.  Most files are therefore behind skip_on_cran().
#
#   * testthat workers: SERIAL on CI or CRAN.  Nothing here is parallel-safe to
#     assume otherwise, and rxode2/nlmixr2est already parallelize inside a
#     single likelihood evaluation.
#   * rxode2 within-solve threads: capped to 2 on CRAN, per CRAN's two-core
#     policy; on CI and locally rxode2 manages its own threads.
#
# CI test selection lives in .github/workflows/ci.yaml, not here: the fast
# tier runs tests/ci-core.txt and the nightly full tier runs every file, both
# sharded by .github/scripts/ci-shards.R and run through testthat::test_dir()
# (.github/workflows/_tests.yaml).  R CMD check (NOT_CRAN=false) runs the
# CRAN-visible subset below.  Coverage sets NLMIXR2BAYES_TEST_FILTER to the
# files it can afford.
#
# The opt-in gates, all off unless a workflow sets them:
#   * RXSTAN_STAN_TESTS -- the hand-coded .stan compile suite (skipUnlessStan()
#     in helper-rxstan.R); every _tests.yaml shard sets it.
#   * NLMIXR2STAN_SLOW -- G5 simulation-based calibration, G6 FOCEi agreement;
#     only the weekly slow-tests.yaml sets it.
#   * NLMIXR2BAYES_TEST_PLATFORM_COMPILE -- set by the macOS/Windows
#     CRAN-mirror checks in the full tier, opting test-platform-compile.R back
#     in under NOT_CRAN=false: one real Stan compile per subsystem, so those
#     platforms prove Stan builds there at all.
# -------------------------------------------------------------------------
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
if (nzchar(Sys.getenv("NLMIXR2BAYES_TEST_FILTER"))) {
  .filter <- Sys.getenv("NLMIXR2BAYES_TEST_FILTER")
}
# Locally (and on CRAN) .filter stays NULL -> run everything.

# stop_on_failure MUST stay TRUE (the default): with it FALSE, R CMD check
# never errors on a test failure, so `checking tests` passes no matter what.
test_check("nlmixr2bayes", filter = .filter, perl = TRUE)
