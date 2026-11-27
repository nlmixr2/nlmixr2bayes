# Build a Stan \`init\` function from a model's \`ini()\` values

Population PK models are routinely multimodal – a one-compartment oral
model can swap absorption and elimination and fit almost as well, the
"flip-flop" – so Stan's default inits, drawn uniformly on \`(-2, 2)\`
unconstrained, regularly send two chains to two different modes and
report \`Rhat\` in the twenties. Starting from the \`ini()\` estimates
avoids that.

## Usage

``` r
rxsInit(gen, jitter = 0.2)
```

## Arguments

- gen:

  the list returned by \[rxsStanFromUi()\]

- jitter:

  SD of the normal noise added to each fixed effect per chain. \`0\`
  gives every chain the same starting point, which makes \`Rhat\`
  optimistic – prefer a small positive value.

## Value

a function of \`chain_id\`, suitable for \`rstan::sampling(init = )\`

## Details

Jitter is applied to the fixed effects only. The between-subject SDs and
the residual error are left at their \`ini()\` values because they are
positively constrained and a jitter large enough to be useful is also
large enough to push them past zero; the \`z\` matrix starts at 0, which
is the prior mean of the non-centered parameterization.

## Author

Lukas A. Widmer

## Examples

``` r
# \donttest{
mod <- function() {
  ini({
    tka <- 0.5
    tcl <- 1
    tv <- 3.4
    eta.cl ~ 0.1
    add.sd <- 0.2
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}
ev <- rxode2::et(amt = 100, cmt = "depot")
ev <- as.data.frame(rxode2::et(ev, c(0.5, 1, 2, 4, 8, 12, 24)))
d <- rbind(cbind(ev, id = 1L), cbind(ev, id = 2L))
d$dv <- 1
gen <- rxsStanFromUi(mod, d)
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#>  
#>  
#>  
#>  
#> → calculate sensitivities
# one init list per chain; the fixed effects are jittered
init <- rxsInit(gen)
str(init(chain_id = 1))
#> List of 6
#>  $ tka         : num 0.694
#>  $ tcl         : num 0.799
#>  $ tv          : num 3.38
#>  $ omega_eta_cl: num 0.316
#>  $ z           : num [1, 1:2] 0 0
#>  $ add_sd      : num 0.2
# jitter = 0 starts every chain at the ini() values
str(rxsInit(gen, jitter = 0)(chain_id = 2))
#> List of 6
#>  $ tka         : num 0.5
#>  $ tcl         : num 1
#>  $ tv          : num 3.4
#>  $ omega_eta_cl: num 0.316
#>  $ z           : num [1, 1:2] 0 0
#>  $ add_sd      : num 0.2
# }
```
