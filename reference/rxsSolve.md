# Solve a registered handle from R

Goes through the same C entry point Stan uses, so the test suite covers
the real ABI without compiling a Stan model.

## Usage

``` r
rxsSolve(handle, p, slow = FALSE)
```

## Arguments

- handle:

  handle returned by \[rxsRegister()\]

- p:

  numeric vector of parameter values, in \`sens\` order

- slow:

  force the slow \[rxode2::rxSolve()\] route, bypassing Path A

## Value

an \`ny x (np + 1)\` matrix: value, then \`dy/dp\` per parameter

## Author

Lukas A. Widmer
