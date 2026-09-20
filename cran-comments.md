## Submission

This is a new submission.

nlmixr2bayes adds a Stan-based Bayesian estimation method to 'nlmixr2'.
Rather than re-expressing the model in the Stan language, it links Stan
directly to the compiled 'rxode2'/'nlmixr2est' log-likelihood and its
analytic gradients, so ODE solving, the residual-error models and censoring
stay in one place.

## Test environments

* local: Ubuntu 24.04, R 4.6.1
* GitHub Actions: ubuntu-latest (R release, R devel, R oldrel-1),
  macos-latest (R release), windows-latest (R release)

## R CMD check results

<!-- Fill in from the final `R CMD check --as-cran` run before submitting. -->

## Notes the maintainer expects, and why

* **`Unexported object imported by a ':::' call: 'rstan:::mk_cppmodule'`**
  (R/stanCompile.R).  This package caches compiled Stan models on disk.  A
  cached `stanmodel` whose compiled module can no longer be instantiated
  otherwise fails much later, inside `rstan::sampling()`, as
  `"NULL value passed for DllInfo"` -- with no indication that the cache is
  the cause.  Probing with the same internal call rstan itself uses is the
  only way we have found to detect that before it reaches the user.  The
  call is wrapped in `tryCatch()`, so if the internal is removed or changes
  signature the probe reports "not usable" and the model is recompiled
  rather than erroring; nothing breaks.  We would gladly switch to a
  documented rstan entry point if one is added.

* **`Suggests` contains 'StanEstimators', which is not on CRAN.**  It is
  declared in `Additional_repositories`
  (https://andrjohns.r-universe.dev).  It is the Pathfinder backend: rstan
  does not expose Stan's Pathfinder service, and CmdStan cannot work here
  at all, because a separate executable cannot reach the in-process linked
  likelihood.  It is reached only through
  `requireNamespace("StanEstimators", quietly = TRUE)`; when it is absent,
  `est = "pathfinder"` reports that cleanly and every other algorithm is
  unaffected.

* **Package size.**  `inst/cache` holds the fitted models the vignettes and
  README display.  Every number shown in the documentation is produced by
  really running the code, and caching the fits is what keeps the vignettes
  to seconds instead of hours of sampling.  The fits are already stripped of
  rstan's compiled shared object before being saved (roughly 15 MB per fit
  down to 1.5 MB) while keeping every read-only accessor working.

* **Compiled code.**  The package calls 'nlmixr2est''s exported C
  function-pointer tables (installed at load time via `R_GetCCallable`)
  rather than reimplementing any likelihood or gradient math.  A runtime
  API-version check refuses to run against a mismatched 'nlmixr2est' build
  instead of misbehaving silently.
