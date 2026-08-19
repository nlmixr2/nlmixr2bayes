# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

`nlmixr2bayes` gives nlmixr2 a Stan-based Bayesian estimation method, via **two
independent, complementary subsystems** that share a `.Rbuildignore`d C++
runtime but have no code paths in common otherwise:

1. **Likelihood-level linking** (`est = "stan"` / `"nuts"` / `"advi"` /
   `"pathfinder"`) — the primary, actively developed engine. Stan is linked
   directly to the *compiled* rxode2/nlmixr2est likelihood (ODE solving,
   residual-error models, censoring all stay in rxode2/nlmixr2est) rather than
   the model being re-expressed in the Stan language. Stan gets the
   log-likelihood and analytic gradients through `precomputed_gradients()`
   and samples theta/Omega/eta jointly.
2. **ODE-level hand-coded Stan models** (the "rxstan" bridge —
   `rxsRegister()`/`rxsStanModel()`/`rxsStanFromUi()`): write your own `.stan`
   program and declare rxode2 as the ODE solver backend (`rx_solve()`), for
   when you need custom Stan code, DDEs, or full control the automatic path
   doesn't expose.

Read `README.md` for the full "why two approaches" pitch, the speed/gradient
benchmark numbers, and the correctness-oracle table — it's kept up to date
and isn't duplicated here.

## Development commands

Load the package for interactive development (this repo has no
`NEWS.md`/install workflow beyond standard R package tooling):

```r
devtools::load_all(".")
```

Run tests. `tests/testthat.R` uses `test_check()`, which needs the package
*installed*; for iteration, load + run individual files instead (much
faster, and is how this repo's own test suite is normally driven):

```r
devtools::load_all(quiet = TRUE)
testthat::test_file("tests/testthat/test-est-stan.R")
```

`skip_on_cran()` gates most of the suite. `devtools::test()` sets
`NOT_CRAN` automatically; a plain `Rscript -e 'testthat::test_file(...)'`
does not, so set it explicitly:

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-est-stan.R")'
```

Other test-gating env vars:
- `RXSTAN_STAN_TESTS=1` — the rxstan-bridge hand-coded `.stan` compilation
  tests (`test-omega-corr.R`, `test-lincmt.R`, etc.) are skipped without it.
- `NLMIXR2STAN_SLOW=TRUE` — nightly-scale gates (G5 simulation-based
  calibration in `test-sbc.R`, G6 FOCEi-agreement in `test-focei-agree.R`).
  `NLMIXR2STAN_SBC_N` overrides the SBC replicate count (default 40).

Most `est="stan"`/`est="nuts"` tests do a **real Stan compile** the first
time a given generated-program shape is seen (~1-2 min), then hit an
on-disk + in-session cache (`stanControl(cache=)`, keyed by content hash) —
`tests/testthat/helper-rxstan.R` also memoizes compiled `.stan` models
across files for the rxstan-bridge suite. Compiling dominates suite
wall-time, so prefer running the one or two files relevant to a change over
the full suite while iterating; check `uptime`/`ps` before kicking off a
full run or a benchmark, since a concurrent heavy compile will both slow
down and invalidate the other.

Lint (`.lintr`: camelCase, internal functions `^\.[a-z][a-zA-Z0-9]*$`, S3
methods `name.class`, 120-col lines, cyclomatic complexity ≤ 50):

```r
lintr::lint_package()
```

Regenerate docs after adding/changing a roxygen-tagged export:

```r
devtools::document()
```

**Gotcha:** this repo's `man/*.Rd` and `NAMESPACE` have drifted from what a
clean `roxygen2::roxygenise()` run produces even though
`Config/roxygen2/version` matches the installed roxygen2 — a full
regeneration will reformat many unrelated `Rd` files (adding `\author{}`
blocks, etc.) and generate fresh docs for previously-undocumented `rxs*`
exports, because `R/complete.R`'s S3 methods are missing `@export`/
`@exportS3Method` tags. When adding a single new export, it's safer to
hand-edit `NAMESPACE` and add one matching `man/*.Rd` file in the existing
style than to run `document()` and sort through an unrelated diff.

## Architecture

### Likelihood-level linking (`est="stan"` family)

Pipeline, roughly `stanControl.R` (options + validation) → `stanEst.R`
(`nlmixr2Est.stan`, the orchestrator) → `stanMap.R` (`.stanMap()`: builds
the theta/eta/omega-block parameter map from the ui's `iniDf`) →
`stanPriors.R` (`.stanParName()` mangles nlmixr2 names to Stan identifiers,
`.stanBlockId()` derives an omega block's Stan suffix from its member eta
names, prior emission) → `stanGen.R` (`.stanGenerate()`: assembles the
actual `.stan` program text) → `stanCompile.R` (cached `rstan::stan_model()`)
→ `stanLink.R` (links the compiled rxode2/nlmixr2est likelihood; the C side)
→ `rstan::sampling()`/`vb()`/Pathfinder → posterior mapped back into an
nlmixr2 fit (theta/omega/etaObf, `$posteriorSummary`, WAIC/LOO, diagnostics).

- **Two generated-program shapes.** A model *with* etas ("tier 2") gets a
  `functions{}` block declaring two external functions the compiled Stan
  model calls into: `nlmixr2_cond_all2(matrix etaMat, vector theta)` (value +
  analytic gradient, via nlmixr2est) and `nlmixr2_iter_tick(vector par, real
  objf)` (an **ungated** print hook — fires on every raw log-density
  evaluation; nlmixr2est's C side decides internally whether to actually
  print/record, at the `stanControl(print=)` cadence). A model with *no*
  etas ("tier 0" — population-only/single-subject) instead uses
  nlmixr2est's separate NLM C API (`nlmixr2_pop_ll(theta)`) with its own,
  unrelated print/history mechanism (`nlmixr2est::nlmGetParHist`).
- **Declared Stan parameter names read like the model.** Thetas (including
  error-model and covariate-coefficient parameters — ordinary theta rows)
  get Stan identifiers straight from `.stanParName()` (dots → underscores).
  Each omega block's internal sampling parameters (`sd<id>`/`Lcorr<id>`/
  `omega<id>`/`L<id>`/`z<id>`/`etaP<id>`) are suffixed by `.stanBlockId()`
  from the block's member eta names, not a raw block index — so
  `fit$stanCode`/`fit$stanfit` are directly legible, not just
  `fit$posteriorSummary` (which additionally relabels the derived
  `omegaOut` cells to the `om.<eta>`/`cov.<eta1>.<eta2>` convention and
  drops the internal per-block parameterization as redundant with it).
- **Sugar estimation methods** (`stanSugar.R`): `est="nuts"`/`"advi"`/
  `"pathfinder"` are thin wrappers that force the matching
  `stanControl(algorithm=)` and dispatch to `nlmixr2Est.stan`, the way
  `"foce"`/`"focei"` are members of nlmixr2est's own focei family — letting
  a user pick the *algorithm* instead of the *tool*. Each needs a matching
  `getValidNlmixrCtl.<name>`/`nmObjGetControl.<name>` shim plus a
  `nlmixr2Est.<name>` S3 method registered in `NAMESPACE`; errors (not
  silent overrides) when the control already names a conflicting algorithm.
- **Parallelism has two axes**, not one: chains (`chainCores`, forked
  processes on unix, PSOCK workers on Windows) and subject-parallel OpenMP
  *within* a single evaluation (`cores`). They're meant to be mutually
  exclusive — forked chains default `cores` to 1 (OpenMP-after-fork hazard)
  — but explicitly setting `cores` without also setting `chainCores=1`
  silently oversubscribes both axes at once, which is a measurable (not
  catastrophic) slowdown relative to the fork-only default.
- **Cross-package C linkage**: nlmixr2bayes doesn't reimplement any
  likelihood/gradient math — it calls nlmixr2est's FOCEi/NLM C
  function-pointer tables, installed at `.onLoad()` (`R/zzz.R`:
  `.iniFoceiPtrs()`/`.iniNlmPtrs()`) from `nlmixr2est::.nlmixr2estFoceiPtrs()`/
  `.nlmixr2estNlmPtrs()`, with a runtime API-version check that refuses to
  run (rather than silently misbehave) against a mismatched nlmixr2est
  build. `src/nlmixr2bayesPtr.c` is the C shim: it re-exports the linked
  calls as the plain-C symbols the Stan-generated `external` functions
  resolve via `R_GetCCallable` (deliberately no Rcpp/Armadillo/Eigen there —
  the Stan translation unit and this package's own C++ must never share a
  C++ type across the DSO boundary).

### The rxstan bridge (hand-coded Stan)

Separate codegen and naming scheme from the above — `codegen.R`
(`rxsStanFromUi()`, `.rxsName()` mangling, `blk<N>` block ids, *not*
`.stanBlockId()`), `stan.R`/`check.R`/`registry.R` (handle registry,
`rxsIncludePath()`/`rxsStanIncludes()` for the `includes=` namespace
surgery `rstan::stan_model()` needs), `stanBuild.R` (load-time stanc probes
that turn a stanc3/rstan internals shift into an actionable refusal instead
of silent breakage), and `src/bridge.cpp`/`src/fast.cpp` (the C++ solver
policies `rx_solve()` dispatches to). Don't reach for `.stanParName()` or
`.stanBlockId()` here, or for `.rxsName()` on the likelihood-linking side —
they're not interchangeable.

## Testing conventions

`tests/testthat/helper-est.R`/`helper-link.R` define the shared small
fixtures (`.estMod`/`.linkData()`, `.linkMod`) most likelihood-linking tests
reuse — a 4-subject analytic (non-ODE) model, chosen so link tests compile
fast and the conditional density is hand-computable. Prefer extending an
existing test file over adding a new one when a fixture already fits;
`skipUnlessStan()` (`helper-rxstan.R`) is the standard skip gate for
rxstan-bridge tests that need an actual Stan compile.
