# Is Path A currently armed for this handle?

Is Path A currently armed for this handle?

## Usage

``` r
rxsFastAvailable(handle)
```

## Arguments

- handle:

  handle returned by \[rxsRegister()\]

## Value

\`TRUE\` when the cached "fast" solve is armed for \`handle\`, \`FALSE\`
when the next solve will go through the slow path

## Author

Lukas A. Widmer
