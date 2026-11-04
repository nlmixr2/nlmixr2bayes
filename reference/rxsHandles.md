# List registered handles

List registered handles

## Usage

``` r
rxsHandles()
```

## Value

an integer vector of the handles currently registered, empty when none
are

## Author

Lukas A. Widmer

## Examples

``` r
# integer(0) until rxsRegister() has been called
rxsHandles()
#> integer(0)
```
