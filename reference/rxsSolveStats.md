# How often a handle has solved, and how often that failed

Counted regardless of \`quiet\`, so a silenced run still reports whether
the solver was struggling.

## Usage

``` r
rxsSolveStats(handle)
```

## Arguments

- handle:

  handle returned by \[rxsRegister()\]

## Value

named numeric vector with \`solves\` and \`failures\`

## Author

Lukas A. Widmer
