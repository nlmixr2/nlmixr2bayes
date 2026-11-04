# Force the next solve back through the slow path

Call this if something else in the session has solved an rxode2 model
and you want the bridge to rebuild its cached solve immediately rather
than on the next detection.

## Usage

``` r
rxsInvalidateFast()
```

## Value

\`NULL\`, invisibly; called for its effect on the cached solve

## Author

Lukas A. Widmer
