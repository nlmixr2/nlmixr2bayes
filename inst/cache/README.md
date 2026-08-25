# Documentation fit cache

These files are the fits the nlmixr2bayes README and vignettes display.  Each
document runs its code for real and caches the expensive results here with
[nlmixr2save](https://github.com/nlmixr2/nlmixr2save)'s `:=` operator, so the
rendered pages show live output without re-running an hour of Stan sampling on
every build.

The cache lives under `inst/`, so it is **installed with the package**: if you
re-render a vignette from an installed nlmixr2bayes you get these fits back
instead of refitting.

- `../../vignettes/cacheSetup.R` points `:=` here (`bayesCache("<prefix>-")`).
- `../../vignettes/precompute.R` (re)populates this directory by rendering the
  documents:

  ```sh
  Rscript vignettes/precompute.R           # fit only what is missing
  Rscript vignettes/precompute.R --clean   # refit everything
  ```

A cached `stanfit` is stored without rstan's compiled shared object (see
`saveFitItem.stanfit` in `R/stanSave.R`), which is ~99% of its serialized size
and is unusable in another session anyway.  Everything that reads a finished
posterior -- `rstan::summary()`, `extract()`, `traceplot()`,
`get_sampler_params()`, `posterior::as_draws_df()` -- still works on a loaded
fit; only re-entering the compiled model (`log_prob()`, more `sampling()`)
needs the session that compiled it.

Refresh the cache whenever a displayed model, its data, or something that
changes the numbers is updated.
