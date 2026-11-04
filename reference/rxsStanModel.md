# Compile a Stan program against the rxstan bridge

Thin wrapper over \[rstan::stanc()\] + \[rstan::stan_model()\] that
arranges the namespace surgery \`includes=\` needs and keeps the model
name unobfuscated so the namespace is predictable.

## Usage

``` r
rxsStanModel(
  file,
  modelName = "rxstan_model",
  tag = NULL,
  pre = NULL,
  post = NULL,
  ...
)
```

## Arguments

- file:

  path to a \`.stan\` file declaring \`vector rx_solve(int, vector);\`

- modelName:

  C++-safe model name

- tag, pre, post:

  passed to \[rxsStanIncludes()\]

- ...:

  passed to \[rstan::stan_model()\]

## Value

a \`stanmodel\`

## Author

Lukas A. Widmer
