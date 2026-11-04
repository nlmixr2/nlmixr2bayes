# Compile a Stan program against the linked likelihood

Wraps \[rstan::stan_model()\] with \`allow_undefined=TRUE\` and the
\`nlmixr2bayes_lp.hpp\` header injected via \`includes=\`. Compiled
models are cached by a hash of the program, the header, and the
rstan/nlmixr2bayes versions.

## Usage

``` r
stanCompile(
  code = .stanPhase0Code(),
  cache = TRUE,
  cacheDir = NULL,
  verbose = FALSE
)
```

## Arguments

- code:

  Stan program text (default: the phase-0 eta-only program)

- cache:

  use (and populate) the on-disk cache

- cacheDir:

  cache directory (default \`tools::R_user_dir("nlmixr2bayes",
  "cache")\`)

- verbose:

  show the compiler output

## Value

an \`rstan::stanmodel\`

## Author

Matthew L Fidler
