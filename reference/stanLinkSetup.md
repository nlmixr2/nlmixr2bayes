# Load the linked likelihood for a model/data pair

Wraps \[nlmixr2est::foceiLikLoad()\] with \`scale="natural"\` (so the
theta vector is directly comparable with \`ui\$iniDf\$est\`), records
the handle and a setup hash, and takes the package-level re-entrancy
lock.

## Usage

``` r
stanLinkSetup(
  ui,
  data,
  likelihood = c("focei", "foce"),
  rxControl = rxode2::rxControl(),
  thetaSens = FALSE,
  literalFix = TRUE,
  cores = rxode2::getRxThreads(),
  maxOdeRecalc = 3L,
  fallbackFD = TRUE
)
```

## Arguments

- ui:

  rxode2 ui (or model function)

- data:

  estimation data

- likelihood:

  individual likelihood type; only \`"focei"\` and \`"foce"\` are
  self-consistent value/gradient pairs (\`"focep"\` is refused by the
  capability flags)

- rxControl:

  solving options

- thetaSens:

  also wire the theta-sensitivity model (\`d/d(theta)\` of the
  conditional at fixed eta)

- literalFix:

  inline fixed population parameters as model literals (the nlmixr2est
  literal-fix hook); keep this matched with the est method's
  \`stanControl(literalFix=)\` so the linked problem and the generated
  Stan program agree on the parameter vector

- cores:

  subject-parallel thread count for each likelihood evaluation

- maxOdeRecalc:

  failed-solve tolerance relaxations before falling back (default 3, the
  tolerance-report calibration; 0 disables)

- fallbackFD:

  after the relaxations, retry with the prediction model + shi21
  central-difference eta gradients before rejecting with -Inf (default
  TRUE; FALSE rejects immediately)

## Value

invisibly, the link handle (the \`foceiLikLoad\` handle plus \`flags\`,
\`setupHash\`)

## Author

Matthew L Fidler
