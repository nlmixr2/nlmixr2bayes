# Control options for \`est="advi"\` (Stan's ADVI variational Bayes)

\`adviControl()\` is \[stanControl()\] with an ADVI \`algorithm\` forced
– the same options, named after the INFERENCE ALGORITHM instead of the
software tool. Because \`class(adviControl())\` is \`"adviControl"\`,
\`nlmixr2(model, data, adviControl())\` infers \`est="advi"\` with no
\`est=\` argument.

## Usage

``` r
adviControl(..., algorithm = c("meanfield", "fullrank"))
```

## Arguments

- ...:

  \[stanControl()\] options

- algorithm:

  the ADVI variant: \`"meanfield"\` (default; independent Gaussians on
  the unconstrained scale) or \`"fullrank"\` (one full-covariance
  Gaussian). \`"NUTS"\` is accepted only because it is
  \[stanControl()\]'s default and selects \`"meanfield"\`;
  \`"pathfinder"\` errors

## Value

an \`adviControl\` list

## Details

ADVI is typically 10-100x faster than NUTS and is a fair first look, but
it is an APPROXIMATION – it understates tails and correlations; the
Pareto-k diagnostic (\`khat\`) replaces Rhat/divergences and \`khat \>
0.7\` means NUTS should be used instead.

## See also

\[stanControl()\], \[nutsControl()\], \[pathfinderControl()\]

## Author

Matthew L Fidler

## Examples

``` r
# meanfield is the default ADVI variant
ct <- adviControl()
ct$algorithm
#> [1] "meanfield"

ct <- adviControl(algorithm = "fullrank")
ct$algorithm
#> [1] "fullrank"
```
