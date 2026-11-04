# Bayesian estimation via Stan linked to the rxode2/nlmixr2est likelihood

Dispatched by \`nlmixr2(model, data, est="stan", stanControl(...))\`.
The first nlmixr2 estimation method to declare \`nlmixr2Priors =
"all"\`: prior distributions come from the model's \`ini()\` block and
are required (a model without any errors, printing the lines to add).

## Usage

``` r
# S3 method for class 'stan'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  nlmixr2 estimation environment

- ...:

  passed through

## Value

an nlmixr2 fit (or an \`nlmixr2bayesCode\` when
\`stanControl(run=FALSE)\`)

## Author

Matthew L Fidler
