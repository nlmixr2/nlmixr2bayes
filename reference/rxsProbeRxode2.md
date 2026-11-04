# Report whether rxode2's C function-pointer table was installed

The bridge cannot work without it, so this is the first thing to check
if solves fail in a way that looks like rxode2 is not there.

## Usage

``` r
rxsProbeRxode2()
```

## Value

named logical vector, one entry per probed rxode2 entry point

## Author

Lukas A. Widmer

## Examples

``` r
# which rxode2 entry points this build of the bridge found
rxsProbeRxode2()
```
