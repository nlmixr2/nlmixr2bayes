# Set the subject-parallel thread count for the linked likelihood

Stan itself runs single-threaded in nlmixr2bayes; the parallelism lives
in the rxode2/nlmixr2est subject loop inside each likelihood evaluation.
This sets how many threads that loop uses (clamped to what the loaded
solve pool supports, and to 1 when the ODE method is not thread safe).

## Usage

``` r
stanSetCores(cores = rxode2::getRxThreads())
```

## Arguments

- cores:

  number of threads (default \[rxode2::getRxThreads()\])

## Value

invisibly, the thread count set

## Author

Matthew L Fidler
