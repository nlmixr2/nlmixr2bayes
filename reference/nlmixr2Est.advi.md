# ADVI estimation through Stan (sugar for est="stan" + algorithm)

\`est="advi"\` is \`est="stan"\` with
\`stanControl(algorithm="meanfield")\` unless the control already names
an ADVI variant (\`"fullrank"\` is honored; asking for
\`"NUTS"\`/\`"pathfinder"\` under \`est="advi"\` is an error rather than
a silent override). \[adviControl()\] is its control, and
\`nlmixr2(model, data, adviControl())\` infers \`est="advi"\` from it.

## Usage

``` r
# S3 method for class 'advi'
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
