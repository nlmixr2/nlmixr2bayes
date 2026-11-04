# Translate a model's \`ini()\` priors into Stan

Consumes \[rxode2::rxUiPriors()\] and \[lotri::lotriPriorDists()\] and
returns the Stan-side pieces: per-parameter bound declarations
(support-promoted) and sampling statements for population parameters,
and the classified omega-block priors for the model generator.

## Usage

``` r
stanPriors(ui)
```

## Arguments

- ui:

  an rxode2 ui (or model function)

## Value

a list:

- pop:

  data.frame over population-parameter priors: \`name\` (nlmixr2),
  \`par\` (Stan identifier), \`prior\`, \`stanName\`, \`kind\`,
  \`lower\`, \`upper\` (effective, support-promoted), \`constraint\`,
  \`truncated\`, \`members\` (comma-joined for a multivariate group,
  \`NA\` otherwise), \`statement\`

- omega:

  data.frame over omega-block priors: \`name\`, \`neta1\`, \`neta2\`,
  \`prior\`, \`stanName\`, \`kind\`

## Author

Matthew L Fidler

## Examples

``` r
# priors are read from the model's ini({}) block, and map one-to-one
# onto Stan's distributions
mod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- 0.5
    prior(tcl) ~ dnorm(0, 10)
    prior(add.sd) ~ dcauchy(0, 5)
    eta.cl ~ 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}
p <- stanPriors(rxode2::rxode2(mod))
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
p$pop$statement
#> [1] "tcl ~ normal(0, 10);"   "add_sd ~ cauchy(0, 5);"
```
