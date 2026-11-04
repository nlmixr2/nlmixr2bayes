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
* win-builder: R-devel (pending)

## R CMD check results

Local `R CMD check --as-cran`: 0 errors | 0 warnings | 1 notes

* `New submission` (CRAN incoming feasibility).

## Notes the maintainer expects, and why

* **Precomputed vignettes.**  The vignettes demonstrate Stan fits, which
  take a Stan compile and minutes of sampling each.  They are knitted ahead
  of time from `vignettes/*.Rmd.orig` (not shipped), so the shipped
  `vignettes/*.Rmd` contain static output and rebuilding them compiles and
  samples nothing.  Every number shown was produced by really running the
  code.

* **Compiled code.**  The package calls 'nlmixr2est''s exported C
  function-pointer tables (installed at load time via `R_GetCCallable`)
  rather than reimplementing any likelihood or gradient math.  A runtime
  API-version check refuses to run against a mismatched 'nlmixr2est' build
  instead of misbehaving silently.
