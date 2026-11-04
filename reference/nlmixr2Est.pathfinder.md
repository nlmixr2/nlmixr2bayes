# Pathfinder estimation through Stan (sugar for est="stan" + algorithm)

\`est="pathfinder"\` is \`est="stan"\` with
\`stanControl(algorithm="pathfinder")\` forced. \[pathfinderControl()\]
is its control, and \`nlmixr2(model, data, pathfinderControl())\` infers
\`est="pathfinder"\` from it.

## Usage

``` r
# S3 method for class 'pathfinder'
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
