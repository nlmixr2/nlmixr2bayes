# Load the tier-0 (population-only) likelihood for a model/data pair

Wraps \[nlmixr2est::nlmObjectiveSetup()\] with \`gradient=TRUE\` and
\`scale="natural"\` (the nlm C API, nlmixr2/nlmixr2est#953), so repeated
evaluations return minus the log-likelihood and its analytic theta
gradient on the model's own scale. Shares the package-level one-linked-
problem-at-a-time lock with \[stanLinkSetup()\].

## Usage

``` r
stanPopLinkSetup(
  ui,
  data,
  rxControl = rxode2::rxControl(),
  cores = rxode2::getRxThreads(),
  print = 0L
)
```

## Arguments

- ui:

  rxode2 ui (or model function) WITHOUT etas

- data:

  estimation data

- rxControl:

  solving options

- cores:

  subject-parallel thread count for the population solves (applied
  through \`rxControl\$cores\` when that is left at its 0/"auto"
  default; an explicit \`rxControl\$cores\` wins)

- print:

  iteration-print cadence: every \`print\`-th evaluation of the linked
  objective prints an nlmixr2est-style row (0 = off); the history is
  collected by \[nlmixr2est::nlmGetParHist()\]

## Value

invisibly, a handle (ntheta, nobs, flags, initPar, setupHash)

## Author

Matthew L Fidler
