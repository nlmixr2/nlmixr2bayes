# Control options for \`est="pathfinder"\` (Stan's multi-path Pathfinder)

\`pathfinderControl()\` is \[stanControl()\] with \`algorithm =
"pathfinder"\` forced – the same options, named after the INFERENCE
ALGORITHM instead of the software tool. Because
\`class(pathfinderControl())\` is \`"pathfinderControl"\`,
\`nlmixr2(model, data, pathfinderControl())\` infers
\`est="pathfinder"\` with no \`est=\` argument.

## Usage

``` r
pathfinderControl(..., algorithm = c("pathfinder", "NUTS"))
```

## Arguments

- ...:

  \[stanControl()\] options

- algorithm:

  must be \`"pathfinder"\` (\`"NUTS"\` is accepted only because it is
  \[stanControl()\]'s default; an ADVI variant errors)

## Value

a \`pathfinderControl\` list

## Details

Pathfinder is implemented in this package (rstan does not expose the
Pathfinder service), so it needs nothing beyond rstan; see
\[stanControl()\] for \`pathfinderPaths\` and the \`khat\` diagnostic.

## See also

\[stanControl()\], \[nutsControl()\], \[adviControl()\]

## Author

Matthew L Fidler

## Examples

``` r
ct <- pathfinderControl()
ct$algorithm
#> [1] "pathfinder"
```
