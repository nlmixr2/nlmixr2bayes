# NUTS estimation through Stan (sugar for est="stan" + algorithm)

\`est="nuts"\` is \`est="stan"\` with \`stanControl(algorithm="NUTS")\`
forced – naming the SAMPLING ALGORITHM (Stan's No-U-Turn full-Bayes HMC
sampler) the way \`est="advi"\` and \`est="pathfinder"\` name theirs,
instead of naming the software tool (\`est="stan"\`). \`"NUTS"\` is
already \`stanControl()\`'s default, so this is mostly a naming
convenience; asking for a variational algorithm under \`est="nuts"\` is
an error rather than a silent override. \[nutsControl()\] is its
control, and \`nlmixr2(model, data, nutsControl())\` infers
\`est="nuts"\` from it.

## Usage

``` r
# S3 method for class 'nuts'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  nlmixr2 estimation environment

- ...:

  passed through

## Value

nlmixr2 fit

## Author

Matthew L Fidler
