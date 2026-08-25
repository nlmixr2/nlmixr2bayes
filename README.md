
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nlmixr2bayes <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
<!-- badges: end -->

A Stan interface for nlmixr2 combining two complementary approaches:

1.  **Likelihood-level linking** (`nlmixr2(model, data, est = "stan")`):
    Link Stan to the likelihoods that already exist rather than
    re-express the model in Stan’s language. rxode2 compiles the
    ODE/solving code and nlmixr2est holds the likelihood machinery, so
    Stan calls into those directly. That keeps one source of truth for
    the likelihood and avoids reimplementing solving, censoring and the
    residual-error models.

2.  **ODE-level hand-coded Stan models** (the “rxstan” bridge, part of
    this package): write your own Stan programs that declare rxode2 as
    the solver backend via `rxsRegister()`. Keep dosing events, DDEs,
    and analytic parameter sensitivities in rxode2 while writing priors,
    transformations, mixtures, and derived quantities directly in Stan.

## What linking buys

Both approaches share the same fundamental advantage: **things Stan’s
ODE interface cannot express**

- **dosing events** (bolus, infusion, `addl`, steady state, multiple
  compartments) with exact event timing
- **derivatives with respect to modeled dose handling** – an estimated
  lag time, bioavailability, duration or rate
- **delay differential equations**, via rxode2’s `delay()` / `past()` –
  through *both* approaches: a `delay()` model is an ordinary nlmixr2
  model to `est = "nuts"`, because the linked likelihood is the same one
  FOCEi solves
- **analytic parameter sensitivities** via forward sensitivity equations
  (no autodiff through the solver)
- **the residual-error and censoring machinery** (M2/M3/M4) already
  validated in nlmixr2est (likelihood-level only)

The model is compiled and linked **without engaging Stan’s threading** –
a second threading runtime layered on top of rxode2/nlmixr2est’s would
contend with it. Parallelism has two axes instead: by default (unix) the
**chains run in forked processes** (`stanControl(chainCores=)`,
inheriting the `rxode2::getRxThreads()` budget capped at the chain
count; each fork gets a copy-on-write duplicate of the linked state, and
draws are bit-identical to a sequential run), and with `chainCores=1`
the chains run sequentially with **subject-parallel OpenMP inside each
likelihood evaluation** (`stanControl(cores=)`). The rate-limiting step
in most pharmacometric models is the ODE solving; these are the two
places it parallelizes.

Two gradient tiers are provided, both supplying analytic derivatives to
Stan via `precomputed_gradients` (nothing is re-derived by autodiff):

1.  per-individual conditional gradients `d/d(eta) log p(y_i|eta_i)`
    through nlmixr2est;
2.  plus `d/d(theta)` of the conditional at fixed eta (forward
    sensitivities for non-mu structural and residual thetas; the
    mu-reference identity `d/dtheta_p = d/deta_k` for the rest, extended
    to mu-referenced covariate coefficients as
    `d/dtheta_p = cov_i * d/deta_k` when the sensitivity model does not
    already carry them), so Stan samples theta + Omega + etas jointly.
    Stan owns the full, normalized eta prior (non-centered,
    `eta = z L'`), so its autodiff composes the supplied gradients
    through the Omega parameterization for free.

The bridge is **first-order only**: `precomputed_gradients()` builds a
reverse-mode `var`, and Stan Math has no forward-mode (`fvar`)
counterpart for externally supplied partials, so anything needing
higher-order autodiff of the target through the external term –
Riemannian HMC’s metric, the embedded-Laplace (`laplace_marginal`)
machinery – is out of scope. That costs nothing in practice: NUTS, ADVI,
optimization and Pathfinder are all first-order, and a
Laplace-*marginalized* target is exactly what nlmixr2est’s own
(NONMEM-validated) FOCEi/Laplace machinery computes natively.

Beyond NUTS, `stanControl(algorithm=)` runs Stan’s ADVI variational
approximations (`"meanfield"`, `"fullrank"`) and multi-path Pathfinder
(implemented in-package against the model’s exact log density and
analytic gradient, since rstan does not expose the Pathfinder service)
on the same generated program and linked likelihood, returning the same
complete nlmixr2 fit 3-10x faster; the Pareto-khat diagnostic replaces
Rhat/ESS and warns when the approximation is not trustworthy (use NUTS
for the final answer). `est="advi"` and `est="pathfinder"` are sugar for
these, the way `"foce"`/`"focei"` are members of the focei family, each
with its own control – `nutsControl()`, `adviControl()`,
`pathfinderControl()` – so `nlmixr2(model, data, adviControl())` picks
the algorithm with no `est=` argument at all.

Prior distributions come from the `ini({})` block (lotri/rxode2), whose
distribution catalog is deliberately one-to-one with Stan’s – see
`lotri::lotriPriorDists()` and `rxode2::rxUiPriors()`. A model without
priors is refused, printing the exact `prior()` lines to add.

## Installation

You can install the development version of nlmixr2bayes from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("nlmixr2/nlmixr2bayes")
```

## Two Approaches: Automatic or Hand-Coded Stan

### 1. Likelihood-level: Automatic Stan generation via nlmixr2 dispatch

Use `nlmixr2(model, data, est = "nuts")` to automatically generate and
fit a Stan program (`"nuts"` selects the algorithm – Stan’s No-U-Turn
HMC sampler – rather than the software tool; it is sugar for
`est = "stan"` with NUTS, its default). No Stan code to write – Stan
calls rxode2’s solver backend directly and receives analytic gradients.

``` r
library(nlmixr2)
library(nlmixr2bayes)


one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    add.sd <- c(0, 0.7)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    prior(tka) ~ dnorm(0.45, 1)
    prior(tcl) ~ dnorm(1, 1)
    prior(tv) ~ dnorm(3.45, 1)
    prior(add.sd) ~ dcauchy(0, 2.5)
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

fit := nlmixr2est::nlmixr2(one.cmt, theo_sd, est = "nuts",
                           control = nutsControl(seed = 42, chains = 2,
                                                 iter = 1000))
```

(`:=` is nlmixr2save’s caching assignment – it runs the fit once, caches
it under `inst/cache`, and reloads it on later knits so the numbers on
this page are the ones this code produced. With `<-` it is an ordinary
fit. The short 2-chain run keeps the README quick to rebuild; the
defaults are 4 chains of 2000.)

It returns an ordinary nlmixr2 fit object – the sampler diagnostics are
in the header, and the posterior is summarized in the usual parameter
table:

``` r
print(fit)
#> ── nlmixr² Stan (HMC) (2 chains x 500 draws; max Rhat 1.009; min ESS 377; Wald
#> CI from posterior cov) ──
#> 
#>                                  OBJF      AIC      BIC Log-likelihood
#> FOCEi                       127.57085 384.1706 404.3502     -185.08531
#> WAIC (subject-level)         68.36909       NA       NA      -34.18455
#> LOO (leave-one-subject-out)  78.70536       NA       NA      -39.35268
#>                             Condition#(Cov) Condition#(Cor)
#> FOCEi                              18.87145         1.46836
#> WAIC (subject-level)                     NA              NA
#> LOO (leave-one-subject-out)              NA              NA
#> 
#> ── Time (sec $time): ──
#> 
#>             setup   optimize covariance preprocess postprocess table compress
#> elapsed 0.0672646 3.5103e-05 2.0814e-05      0.029       0.015 0.031    0.021
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>         Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka    0.469  0.222 47.4       1.60 (1.03, 2.47)     92.9       17.1 
#> tcl     1.01  0.103 10.2       2.75 (2.24, 3.36)     32.8       20.4 
#> tv      3.45 0.0524 1.52       31.5 (28.5, 35.0)     17.2       26.1 
#> add.sd 0.707 0.0544 7.69    0.707 (0.601, 0.814)                     
#>  
#>   Covariance Type ($covMethod): stan.posterior
#>   Fixed parameter correlations in $cor
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 18
#>   ID     TIME    DV  PRED    RES IPRED   IRES  IWRES eta.ka eta.cl   eta.v depot
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74   1.05  0.0884 -0.481 -0.0811  320.
#> 2 1      0.25  2.84  3.30 -0.465  3.86 -1.02  -1.45  0.0884 -0.481 -0.0811  207.
#> 3 1      0.57  6.57  5.90  0.673  6.80 -0.234 -0.331 0.0884 -0.481 -0.0811  118.
#> # ℹ 129 more rows
#> # ℹ 6 more variables: central <dbl>, ka <dbl>, cl <dbl>, v <dbl>, tad <dbl>,
#> #   dosenum <int>
```

On top of the usual fit contents it carries `$stanfit`,
`$posteriorSummary`, `$loo`, `$waic` and `$stanDiagnostics`.

`$stanfit` is the `rstan::stanfit` itself, so everything in the
rstan/bayesplot/posterior ecosystem works from it:

``` r
fit$stanfit
#> Inference for Stan model: nlmixr2bayes.
#> 2 chains, each with iter=1000; warmup=500; thin=1; 
#> post-warmup draws per chain=500, total post-warmup draws=1000.
#> 
#>                  mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff
#> tka              0.47    0.01 0.22   0.05   0.33   0.47   0.60   0.93   377
#> tcl              1.01    0.00 0.10   0.82   0.95   1.01   1.07   1.22   443
#> tv               3.45    0.00 0.05   3.34   3.42   3.45   3.49   3.56   576
#> add_sd           0.71    0.00 0.05   0.61   0.67   0.70   0.74   0.82   667
#> sd_eta_ka        0.76    0.01 0.20   0.48   0.62   0.73   0.85   1.30   470
#> z_eta_ka[1,1]    0.13    0.02 0.35  -0.57  -0.10   0.12   0.36   0.79   407
#> z_eta_ka[2,1]    0.27    0.02 0.37  -0.45   0.04   0.26   0.51   1.02   452
#> z_eta_ka[3,1]    0.53    0.02 0.40  -0.20   0.24   0.53   0.78   1.38   513
#> z_eta_ka[4,1]   -0.40    0.02 0.37  -1.17  -0.63  -0.40  -0.15   0.30   412
#> z_eta_ka[5,1]   -0.08    0.02 0.35  -0.75  -0.29  -0.06   0.15   0.60   419
#> z_eta_ka[6,1]   -0.55    0.02 0.41  -1.43  -0.80  -0.53  -0.28   0.18   520
#> z_eta_ka[7,1]   -1.11    0.02 0.45  -2.03  -1.40  -1.08  -0.80  -0.30   509
#> z_eta_ka[8,1]   -0.23    0.02 0.38  -0.96  -0.46  -0.22   0.03   0.50   510
#> z_eta_ka[9,1]    2.01    0.02 0.59   1.01   1.59   1.97   2.37   3.26   567
#> z_eta_ka[10,1]  -1.02    0.02 0.42  -1.89  -1.27  -1.00  -0.74  -0.27   412
#> z_eta_ka[11,1]   1.08    0.02 0.46   0.25   0.77   1.05   1.37   2.08   475
#> z_eta_ka[12,1]  -0.74    0.02 0.38  -1.58  -0.97  -0.73  -0.48  -0.04   420
#> sd_eta_cl        0.31    0.00 0.09   0.17   0.25   0.29   0.35   0.52   400
#> z_eta_cl[1,1]   -1.65    0.03 0.57  -2.82  -2.04  -1.62  -1.24  -0.59   499
#> z_eta_cl[2,1]    0.49    0.02 0.47  -0.35   0.15   0.47   0.78   1.46   624
#> z_eta_cl[3,1]    0.11    0.02 0.52  -0.84  -0.23   0.09   0.45   1.11   622
#> z_eta_cl[4,1]   -0.06    0.02 0.50  -1.10  -0.38  -0.06   0.30   0.93   688
#> z_eta_cl[5,1]   -0.50    0.02 0.45  -1.41  -0.77  -0.48  -0.19   0.32   513
#> z_eta_cl[6,1]    1.30    0.02 0.59   0.23   0.89   1.27   1.68   2.53   618
#> z_eta_cl[7,1]    0.53    0.02 0.50  -0.41   0.21   0.52   0.86   1.52   671
#> z_eta_cl[8,1]    0.55    0.02 0.51  -0.46   0.21   0.55   0.86   1.57   800
#> z_eta_cl[9,1]    0.13    0.02 0.49  -0.88  -0.18   0.11   0.45   1.13   745
#> z_eta_cl[10,1]  -1.31    0.02 0.52  -2.41  -1.63  -1.29  -0.94  -0.41   588
#> z_eta_cl[11,1]   0.99    0.02 0.52   0.08   0.64   0.96   1.33   2.11   604
#> z_eta_cl[12,1]  -0.42    0.02 0.46  -1.38  -0.71  -0.40  -0.11   0.40   700
#> sd_eta_v         0.16    0.00 0.05   0.08   0.13   0.16   0.19   0.28   416
#> z_eta_v[1,1]    -0.53    0.02 0.47  -1.51  -0.82  -0.51  -0.21   0.35   824
#> z_eta_v[2,1]     0.11    0.02 0.50  -0.88  -0.23   0.12   0.43   1.18   888
#> z_eta_v[3,1]     0.44    0.02 0.55  -0.61   0.09   0.40   0.75   1.60   732
#> z_eta_v[4,1]    -0.03    0.02 0.55  -1.11  -0.38  -0.05   0.32   1.13   926
#> z_eta_v[5,1]    -0.89    0.02 0.48  -1.82  -1.22  -0.88  -0.56   0.00   731
#> z_eta_v[6,1]     1.27    0.02 0.70  -0.09   0.81   1.25   1.72   2.62  1019
#> z_eta_v[7,1]     0.33    0.03 0.71  -1.12  -0.08   0.35   0.76   1.72   796
#> z_eta_v[8,1]     0.69    0.02 0.62  -0.47   0.28   0.67   1.08   1.95   882
#> z_eta_v[9,1]     0.10    0.02 0.47  -0.81  -0.19   0.10   0.39   1.04   842
#> z_eta_v[10,1]   -1.07    0.02 0.56  -2.19  -1.43  -1.04  -0.70   0.01   995
#> z_eta_v[11,1]    1.00    0.02 0.58  -0.02   0.59   0.94   1.37   2.24   805
#> z_eta_v[12,1]   -1.29    0.02 0.58  -2.49  -1.68  -1.25  -0.88  -0.25  1003
#> eta[1,1]         0.09    0.01 0.26  -0.43  -0.07   0.09   0.26   0.58   434
#> eta[1,2]        -0.48    0.01 0.15  -0.78  -0.57  -0.48  -0.38  -0.20   619
#> eta[1,3]        -0.08    0.00 0.07  -0.22  -0.13  -0.08  -0.03   0.06   746
#> eta[2,1]         0.19    0.01 0.27  -0.34   0.03   0.19   0.37   0.72   469
#> eta[2,2]         0.14    0.01 0.14  -0.12   0.05   0.14   0.23   0.41   647
#> eta[2,3]         0.02    0.00 0.08  -0.14  -0.04   0.02   0.06   0.17   824
#> eta[3,1]         0.38    0.01 0.28  -0.16   0.18   0.38   0.57   0.91   518
#> eta[3,2]         0.03    0.01 0.15  -0.25  -0.07   0.03   0.13   0.33   626
#> eta[3,3]         0.06    0.00 0.08  -0.08   0.01   0.07   0.11   0.21   751
#> eta[4,1]        -0.29    0.01 0.27  -0.83  -0.47  -0.29  -0.11   0.23   444
#> eta[4,2]        -0.02    0.01 0.14  -0.31  -0.11  -0.02   0.08   0.25   697
#> eta[4,3]        -0.01    0.00 0.08  -0.16  -0.06  -0.01   0.05   0.16   818
#> eta[5,1]        -0.06    0.01 0.26  -0.58  -0.22  -0.05   0.11   0.42   440
#> eta[5,2]        -0.15    0.01 0.13  -0.42  -0.23  -0.14  -0.05   0.09   532
#> eta[5,3]        -0.14    0.00 0.07  -0.29  -0.19  -0.14  -0.09   0.00   609
#> eta[6,1]        -0.40    0.01 0.29  -0.99  -0.59  -0.39  -0.21   0.17   542
#> eta[6,2]         0.38    0.01 0.16   0.08   0.27   0.38   0.49   0.72   649
#> eta[6,3]         0.20    0.00 0.11  -0.01   0.13   0.20   0.28   0.43   646
#> eta[7,1]        -0.80    0.01 0.30  -1.42  -1.00  -0.79  -0.61  -0.24   524
#> eta[7,2]         0.15    0.01 0.14  -0.13   0.06   0.15   0.25   0.44   697
#> eta[7,3]         0.05    0.00 0.11  -0.18  -0.01   0.05   0.12   0.25   758
#> eta[8,1]        -0.17    0.01 0.27  -0.72  -0.36  -0.17   0.02   0.38   520
#> eta[8,2]         0.16    0.01 0.15  -0.13   0.06   0.16   0.26   0.44   745
#> eta[8,3]         0.11    0.00 0.09  -0.07   0.04   0.11   0.17   0.28   847
#> eta[9,1]         1.46    0.01 0.38   0.79   1.22   1.43   1.66   2.34   671
#> eta[9,2]         0.04    0.01 0.14  -0.24  -0.06   0.03   0.13   0.31   696
#> eta[9,3]         0.01    0.00 0.07  -0.13  -0.03   0.01   0.06   0.17   748
#> eta[10,1]       -0.74    0.01 0.27  -1.26  -0.93  -0.74  -0.56  -0.22   447
#> eta[10,2]       -0.38    0.01 0.14  -0.69  -0.47  -0.38  -0.29  -0.11   597
#> eta[10,3]       -0.17    0.00 0.09  -0.36  -0.22  -0.16  -0.11   0.00   864
#> eta[11,1]        0.78    0.01 0.30   0.23   0.57   0.77   0.97   1.39   543
#> eta[11,2]        0.29    0.01 0.15   0.03   0.18   0.28   0.39   0.58   709
#> eta[11,3]        0.15    0.00 0.08   0.00   0.09   0.15   0.21   0.32   934
#> eta[12,1]       -0.54    0.01 0.27  -1.08  -0.71  -0.53  -0.37  -0.03   437
#> eta[12,2]       -0.12    0.01 0.13  -0.38  -0.21  -0.12  -0.03   0.13   663
#> eta[12,3]       -0.20    0.00 0.09  -0.39  -0.25  -0.19  -0.14  -0.03   763
#> theta[1]         0.47    0.01 0.22   0.05   0.33   0.47   0.60   0.93   377
#> theta[2]         1.01    0.00 0.10   0.82   0.95   1.01   1.07   1.22   443
#> theta[3]         3.45    0.00 0.05   3.34   3.42   3.45   3.49   3.56   576
#> theta[4]         0.71    0.00 0.05   0.61   0.67   0.70   0.74   0.82   667
#> L_eta_ka[1,1]    0.76    0.01 0.20   0.48   0.62   0.73   0.85   1.30   470
#> L_eta_cl[1,1]    0.31    0.00 0.09   0.17   0.25   0.29   0.35   0.52   400
#> L_eta_v[1,1]     0.16    0.00 0.05   0.08   0.13   0.16   0.19   0.28   416
#> omegaOut[1,1]    0.62    0.02 0.37   0.23   0.39   0.53   0.73   1.69   492
#> omegaOut[1,2]    0.00     NaN 0.00   0.00   0.00   0.00   0.00   0.00   NaN
#> omegaOut[1,3]    0.00     NaN 0.00   0.00   0.00   0.00   0.00   0.00   NaN
#> omegaOut[2,1]    0.00     NaN 0.00   0.00   0.00   0.00   0.00   0.00   NaN
#> omegaOut[2,2]    0.10    0.00 0.07   0.03   0.06   0.09   0.12   0.27   455
#> omegaOut[2,3]    0.00     NaN 0.00   0.00   0.00   0.00   0.00   0.00   NaN
#> omegaOut[3,1]    0.00     NaN 0.00   0.00   0.00   0.00   0.00   0.00   NaN
#> omegaOut[3,2]    0.00     NaN 0.00   0.00   0.00   0.00   0.00   0.00   NaN
#> omegaOut[3,3]    0.03    0.00 0.02   0.01   0.02   0.02   0.04   0.08   472
#> logLikSubj[1]   -2.16    0.05 1.31  -5.33  -2.80  -1.87  -1.20  -0.56   567
#> logLikSubj[2]   -6.57    0.05 1.16  -9.44  -7.17  -6.34  -5.72  -4.91   582
#> logLikSubj[3]    2.03    0.06 1.37  -1.29   1.26   2.22   3.01   4.11   446
#> logLikSubj[4]   -3.50    0.06 1.26  -6.74  -4.02  -3.18  -2.61  -2.08   466
#> logLikSubj[5]  -11.36    0.06 1.70 -15.46 -12.31 -11.18 -10.16  -8.58   953
#> logLikSubj[6]   -0.20    0.05 1.30  -3.47  -0.84   0.02   0.75   1.59   590
#> logLikSubj[7]    1.58    0.05 1.22  -1.34   0.93   1.78   2.44   3.49   505
#> logLikSubj[8]   -1.27    0.05 1.07  -3.94  -1.86  -1.00  -0.49   0.15   418
#> logLikSubj[9]   -0.24    0.05 1.36  -3.46  -0.87   0.05   0.72   1.56   645
#> logLikSubj[10]   0.87    0.06 1.39  -2.54   0.07   1.15   1.89   2.91   546
#> logLikSubj[11]   1.80    0.07 1.49  -1.54   0.92   2.03   2.85   4.19   397
#> logLikSubj[12]  -0.81    0.06 1.45  -4.46  -1.55  -0.55   0.29   1.12   525
#> lp__           -39.47    0.49 6.98 -53.94 -44.26 -39.16 -34.33 -27.29   205
#>                Rhat
#> tka            1.01
#> tcl            1.00
#> tv             1.00
#> add_sd         1.00
#> sd_eta_ka      1.00
#> z_eta_ka[1,1]  1.01
#> z_eta_ka[2,1]  1.00
#> z_eta_ka[3,1]  1.00
#> z_eta_ka[4,1]  1.01
#> z_eta_ka[5,1]  1.01
#> z_eta_ka[6,1]  1.00
#> z_eta_ka[7,1]  1.00
#> z_eta_ka[8,1]  1.00
#> z_eta_ka[9,1]  1.01
#> z_eta_ka[10,1] 1.00
#> z_eta_ka[11,1] 1.01
#> z_eta_ka[12,1] 1.01
#> sd_eta_cl      1.00
#> z_eta_cl[1,1]  1.00
#> z_eta_cl[2,1]  1.00
#> z_eta_cl[3,1]  1.00
#> z_eta_cl[4,1]  1.00
#> z_eta_cl[5,1]  1.01
#> z_eta_cl[6,1]  1.00
#> z_eta_cl[7,1]  1.00
#> z_eta_cl[8,1]  1.00
#> z_eta_cl[9,1]  1.00
#> z_eta_cl[10,1] 1.00
#> z_eta_cl[11,1] 1.00
#> z_eta_cl[12,1] 1.00
#> sd_eta_v       1.01
#> z_eta_v[1,1]   1.00
#> z_eta_v[2,1]   1.00
#> z_eta_v[3,1]   1.00
#> z_eta_v[4,1]   1.00
#> z_eta_v[5,1]   1.00
#> z_eta_v[6,1]   1.00
#> z_eta_v[7,1]   1.00
#> z_eta_v[8,1]   1.00
#> z_eta_v[9,1]   1.00
#> z_eta_v[10,1]  1.00
#> z_eta_v[11,1]  1.00
#> z_eta_v[12,1]  1.00
#> eta[1,1]       1.00
#> eta[1,2]       1.00
#> eta[1,3]       1.00
#> eta[2,1]       1.00
#> eta[2,2]       1.00
#> eta[2,3]       1.00
#> eta[3,1]       1.00
#> eta[3,2]       1.00
#> eta[3,3]       1.00
#> eta[4,1]       1.01
#> eta[4,2]       1.00
#> eta[4,3]       1.00
#> eta[5,1]       1.00
#> eta[5,2]       1.00
#> eta[5,3]       1.00
#> eta[6,1]       1.00
#> eta[6,2]       1.00
#> eta[6,3]       1.00
#> eta[7,1]       1.00
#> eta[7,2]       1.00
#> eta[7,3]       1.00
#> eta[8,1]       1.00
#> eta[8,2]       1.00
#> eta[8,3]       1.00
#> eta[9,1]       1.01
#> eta[9,2]       1.00
#> eta[9,3]       1.00
#> eta[10,1]      1.00
#> eta[10,2]      1.00
#> eta[10,3]      1.00
#> eta[11,1]      1.01
#> eta[11,2]      1.00
#> eta[11,3]      1.00
#> eta[12,1]      1.01
#> eta[12,2]      1.00
#> eta[12,3]      1.00
#> theta[1]       1.01
#> theta[2]       1.00
#> theta[3]       1.00
#> theta[4]       1.00
#> L_eta_ka[1,1]  1.00
#> L_eta_cl[1,1]  1.00
#> L_eta_v[1,1]   1.01
#> omegaOut[1,1]  1.00
#> omegaOut[1,2]   NaN
#> omegaOut[1,3]   NaN
#> omegaOut[2,1]   NaN
#> omegaOut[2,2]  1.00
#> omegaOut[2,3]   NaN
#> omegaOut[3,1]   NaN
#> omegaOut[3,2]   NaN
#> omegaOut[3,3]  1.01
#> logLikSubj[1]  1.00
#> logLikSubj[2]  1.00
#> logLikSubj[3]  1.01
#> logLikSubj[4]  1.01
#> logLikSubj[5]  1.00
#> logLikSubj[6]  1.00
#> logLikSubj[7]  1.00
#> logLikSubj[8]  1.00
#> logLikSubj[9]  1.00
#> logLikSubj[10] 1.00
#> logLikSubj[11] 1.01
#> logLikSubj[12] 1.00
#> lp__           1.01
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Aug 25 14:45:22 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

`$posteriorSummary` gives exact posterior quantiles for every monitored
parameter, with the omega blocks relabelled to the model’s own
`om.<eta>`/`cov.<eta1>.<eta2>` names. Prefer these to the Wald intervals
in `$parFixed` for anything skewed, such as a variance:

``` r
fit$posteriorSummary
#>                   mean      se_mean         sd          2.5%          25%
#> tka         0.46948906 0.0114620670 0.22247017   0.050676984   0.32747237
#> tcl         1.00991883 0.0049094682 0.10328449   0.818972569   0.94591366
#> tv          3.45136850 0.0021846576 0.05242960   3.344720241   3.41822574
#> add.sd      0.70717694 0.0021069454 0.05439601   0.609330281   0.66905165
#> om.eta.ka   0.62267844 0.0166364357 0.36918230   0.226031454   0.38906228
#> om.eta.cl   0.10232599 0.0031799534 0.06785634   0.030146072   0.06073637
#> om.eta.v    0.02905885 0.0009316727 0.02023344   0.006779386   0.01584824
#> lp__      -39.46627010 0.4875328436 6.98245697 -53.943391519 -44.25512451
#>                    50%          75%        97.5%    n_eff      Rhat
#> tka         0.46978534   0.60431376   0.92601169 376.7192 1.0066839
#> tcl         1.01198249   1.07478027   1.21642995 442.5897 1.0007022
#> tv          3.45000088   3.48572722   3.56043843 575.9521 0.9985170
#> add.sd      0.70465557   0.74078662   0.82130641 666.5421 1.0049569
#> om.eta.ka   0.53027864   0.72610224   1.69039676 492.4489 1.0009359
#> om.eta.cl   0.08655029   0.12300478   0.27034405 455.3437 0.9994274
#> om.eta.v    0.02436028   0.03707306   0.07791673 471.6421 1.0087301
#> lp__      -39.15608436 -34.33374592 -27.28621891 205.1204 1.0117524
```

`$loo` and `$waic` are the `loo` package objects, computed from the
per-subject conditional log-likelihood draws – so they are
leave-one-*subject*-out, the subject being the exchangeable unit of a
hierarchical model:

``` r
fit$loo
#> 
#> Computed from 1000 by 12 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo    -39.4 14.4
#> p_loo        27.0  2.0
#> looic        78.7 28.7
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.5, 0.9]).
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.67]   (good)     1      8.3%   53      
#>    (0.67, 1]   (bad)      6     50.0%   <NA>    
#>     (1, Inf)   (very bad) 5     41.7%   <NA>    
#> See help('pareto-k-diagnostic') for details.

fit$waic
#> 
#> Computed from 1000 by 12 log-likelihood matrix.
#> 
#>           Estimate   SE
#> elpd_waic    -34.2 14.0
#> p_waic        21.8  1.6
#> waic          68.4 27.9
#> 
#> 12 (100.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.
```

`$stanDiagnostics` is what the fit header summarizes: divergences,
tree-depth saturation, worst Rhat and smallest ESS, plus any actionable
message. An empty `$messages` means every check passed:

``` r
fit$stanDiagnostics
#> $nDivergent
#> [1] 0
#> 
#> $nMaxTreedepth
#> [1] 0
#> 
#> $maxRhat
#> [1] 1.008782
#> 
#> $minEss
#> [1] 376.7192
#> 
#> $khat
#> [1] NA
#> 
#> $messages
#> character(0)
```

### 2. ODE-level: Hand-written Stan models with the rxstan backend

Write your own Stan program and link it to rxode2’s solver via
`rxsRegister()`, part of this package (`nlmixr2bayes`). This approach
gives you full control over your Stan program while keeping dosing,
DDEs, and solver sensitivities in rxode2.

The Stan program declares rxode2 as the solver backend and is otherwise
ordinary Stan. This one ships with the package as
`system.file("stan", "pk_1cmt_oral.stan", package = "nlmixr2bayes")`,
and is what the code below compiles:

``` stan
functions {
  // declared, never defined here: the linked C++ bridge supplies it
  vector rx_solve(int handle, vector p);
}
data {
  int<lower=1> nObs;
  int<lower=1> handle;
  vector<lower=0>[nObs] cpObs;
}
parameters {
  real lka;
  real lcl;
  real lv;
  real<lower=0> sigma;
}
model {
  vector[3] p = [lka, lcl, lv]';
  vector[nObs] center = rx_solve(handle, p);
  vector[nObs] cp = center / exp(lv);

  lka ~ normal(0, 1);
  lcl ~ normal(1.4, 1);
  lv ~ normal(3.4, 1);
  sigma ~ normal(0, 1);

  log(cpObs) ~ normal(log(cp), sigma);
  target += -sum(log(cpObs));   // Jacobian of the log transform
}
```

`rx_solve()` returns ODE **states** and rxode2’s analytic derivatives of
them; the concentration `cp` is formed in Stan, which differentiates
that step itself. The rxode2 side is a plain model on the log scale,
registered with the parameters Stan will pass, in the order it packs
them:

``` r
library(nlmixr2bayes)

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
"

ev <- rxode2::et(amt = 100, cmt = "depot")
ev <- rxode2::et(ev, c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24))

h <- rxsRegister(pkModel, events = ev,
                 sens = c("lka", "lcl", "lv"), output = "center",
                 atol = 1e-10, rtol = 1e-10)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
h
#> <rxstan handle 1>
#>   parameters (3): lka, lcl, lv
#>   states     (1): center
#>   rows       : 12 obs -> 12 outputs
```

Some simulated data to fit, drawn from the same solve the bridge will
hand Stan (column 1 of `rxsSolve()` is the state, the rest its
Jacobian):

``` r
set.seed(42)
truth <- c(lka = log(1.1), lcl = log(4), lv = log(30))
sol <- rxsSolve(h, unname(truth))
cp <- sol[, 1] / exp(truth[["lv"]])
cpObs <- cp * exp(stats::rnorm(length(cp), 0, 0.15))

standata <- list(nObs = length(cpObs), handle = h, cpObs = cpObs)
```

Compile and sample. `rxsStanModel()` splices the bridge header in at
global scope so `rx_solve()` resolves at link time; from there it is
ordinary rstan:

``` r
sm <- rxsStanModel(system.file("stan", "pk_1cmt_oral.stan",
                               package = "nlmixr2bayes"),
                   modelName = "pk_1cmt_oral")

fit2 <- rstan::sampling(sm, data = standata, chains = 2, iter = 1000,
                        seed = 42, refresh = 0,
                        init = function(chain_id = 1) {
                          list(lka = 0.1, lcl = 1.4, lv = 3.4, sigma = 0.2)
                        })
print(fit2, pars = c("lka", "lcl", "lv", "sigma"))
#> Inference for Stan model: pk_1cmt_oral.
#> 2 chains, each with iter=1000; warmup=500; thin=1; 
#> post-warmup draws per chain=500, total post-warmup draws=1000.
#> 
#>       mean se_mean   sd  2.5%  25%  50%  75% 97.5% n_eff Rhat
#> lka   0.14    0.01 0.16 -0.18 0.04 0.14 0.25  0.46   401    1
#> lcl   1.25    0.00 0.05  1.15 1.22 1.25 1.29  1.35   569    1
#> lv    3.36    0.00 0.09  3.17 3.31 3.37 3.42  3.53   336    1
#> sigma 0.15    0.00 0.04  0.10 0.12 0.14 0.17  0.26   278    1
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Aug 25 17:25:51 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

The posterior sits on the simulated truth, and `rxsCheckGradient()`
confirms the analytic sensitivities match finite differences of the
model’s own log density – the check to re-run whenever the solver
tolerances are loosened:

``` r
rbind(truth = truth,
      posterior = colMeans(as.matrix(fit2, pars = c("lka", "lcl", "lv"))))
#>                  lka      lcl       lv
#> truth     0.09531018 1.386294 3.401197
#> posterior 0.14477675 1.253490 3.364142

rxsCheckGradient(fit2, c(lka = 0.1, lcl = 1.4, lv = 3.4, sigma = 0.2))
#>   par     analytic      numeric      relDiff
#> 1   1  -0.07123147  -0.07123147 2.180990e-10
#> 2   2  -1.40451837  -1.40451837 7.261012e-10
#> 3   3   0.41635454   0.41635455 3.331731e-09
#> 4   4 -13.21096698 -13.21096698 1.177218e-11
```

See `vignette("rxstan-handcoded")` and `vignette("rxstan-dde")` for the
fuller design sketch of hand-coded Stan programs.

## Is it fast?

### Likelihood-level (nlmixr2 dispatch)

|                                             | native Stan                         | nlmixr2bayes (linked)                                |
|---------------------------------------------|-------------------------------------|------------------------------------------------------|
| wall, chains sequential, two-solve path     | –                                   | 408.8 s                                              |
| wall, forked chains, two-solve path         | 142.7 s                             | 302.6 s                                              |
| wall, forked + fused single-solve entry     | 142.7 s                             | 203.8 s                                              |
| wall, forked + fused + C-side omega rebuild | 142.7 s                             | **80.8 s**                                           |
| worst bulk ESS                              | 236                                 | **286-378**                                          |
| worst ESS / s                               | 1.66                                | 0.70 -\> 1.25 -\> 1.60 -\> **4.04**                  |
| gradient evaluation                         | 1.98 ms                             | 2.80 ms via `grad_log_prob` (1.36 ms at the C level) |
| posterior means                             | agree to \< 0.01 on every parameter |                                                      |

The posteriors are the same; the linked run extracts *more* effective
samples per draw; forked chains (now the default) close most of native’s
wall-clock edge, and disabling the failure cascade
(`maxOdeRecalc=0, fallbackFD=FALSE`) changes nothing measurable – the
cascade never fires on a healthy trajectory (the forked retry-on run
reproduced the sequential retry-off draws bit-for-bit). The former
two-integrations-per-gradient gap is closed: with nlmixr2est’s combined
eta+theta sensitivity build the whole tier-2 gradient comes from ONE
solve per subject through a fused batch entry – 0.62 ms per
full-population gradient at the C level (vs 1.26 ms two-model, and
native’s 1.98 ms coupled solve), negotiated automatically when the
loaded nlmixr2est provides it. The last per-evaluation drag – fetching
Omega^-1 by evaluating rxode2’s precomputed R closure every call (~1.2
ms, more than the ODE solve) – is replaced by a C-side chol-factor
rebuild (probe-mapped at setup, verified against the R path, 0.015 ms).
Net: **2.4x native Stan’s ESS/s on the pure-ODE model** (4.04 vs 1.66,
with more effective samples per draw and bit-identical reproducibility),
and a larger win wherever the analytic `linCmt()` tier applies – a tier
a native ODE implementation cannot express. Per gradient the linked C
path is cheaper than the native solve – and for models with closed-form
solutions the `linCmt()` tier evaluates ~5-10x cheaper still, an option
a native ODE implementation does not have. ADVI completes the same fit
in ~30-40 s and Pathfinder in ~90 s (with the Pareto-k-hat diagnostic
guarding both). Reproduce with `Rscript inst/bench/native-vs-linked.R`.

## Is it right?

Because the bridge operates at the **likelihood level** (it exposes
`log p(y_i|eta_i)` and its gradients, not raw ODE solutions), the
oracles are likelihood-level. Independent oracles first – checks against
something the bridge does not share code with:

| check                                                                                                                | agreement                                                                                                                                                             |
|----------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Stan’s assembled gradient (`grad_log_prob`) vs. Richardson finite differences of its own log density                 | ~3e-8                                                                                                                                                                 |
| `d/d(eta)` and `d/d(theta)` vs. central differences of the linked conditional value                                  | ~3e-8                                                                                                                                                                 |
| quadrature marginal over the linked conditional vs. nlmixr2’s own `est="agq"` (nAGQ=101) at matched theta            | 3e-7 (independent *marginalizations* – trapezoid vs adaptive Gauss-Hermite; the integrand engine is shared, and is anchored absolutely by the hand-density row below) |
| censored (M2/M3/M4) and `ll()` conditionals vs. hand-computed textbook densities                                     | exact up to a pinned, parameter-free constant; all gradients FD-verified, including the residual-SD dependence of the censored CDF terms                              |
| prior-only sampling (likelihood stripped) vs. every declared prior                                                   | one-sample KS per parameter, including the default LKJ + half-Cauchy omega path                                                                                       |
| full posterior vs. a hand-written **native-Stan** implementation of the same model and priors (no external function) | within 3x combined MCSE on every parameter                                                                                                                            |
| the linked likelihood engine itself                                                                                  | NONMEM-validated upstream in nlmixr2est (Wang 2007 objective, per-subject ETA/CWRES)                                                                                  |

And internal-consistency checks – exact identities the implementation
must satisfy:

| check                                                                                         | agreement                                                                                              |
|-----------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------|
| linked value vs. `nlmixr2est::foceiLikRun(type="cond")`                                       | exact (same engine, one code path)                                                                     |
| mu-referenced theta-gradient column vs. the eta gradient                                      | exact                                                                                                  |
| assembled gradient under Omega^-1 swaps spanning 4 orders of magnitude                        | exact (the conditional is Omega-free by construction)                                                  |
| `log_prob` decomposes into conditional + normalized eta prior                                 | ~1e-6 (difference form, constants cancel)                                                              |
| bitwise determinism: 500+ interleaved evaluations of both batch entries at alternating points | exact (bitwise identical – the target is a pure function of state, which NUTS’s reversibility assumes) |

Plus simulate-and-recover fits, including a correlated omega block whose
generative correlation lands inside the 90% credible interval, and
FOCEi-agreement runs on `nlmixr2data::theo_sd` under weak priors.

Wrong-by-construction variants are asserted to **fail**: a sign-flipped
gradient assembly (off by `2 Omega^-1 eta`) fails the FD test, and
`likelihood="focep"` – whose value and gradient are gradients of
different functions – is refused by capability flags, with the
underlying FD mismatch itself locked in as a test upstream.

## Current scope

Supported: mixed AND population-only models (a no-eta model – e.g. a
single-subject Bayesian fit – runs as “tier 0” through nlmixr2est’s nlm
path: one external scalar `nlmixr2_pop_ll(theta)` carrying the complete
data log-likelihood and its analytic gradient, value tied to a
hand-written density and FD-verified), ODE, delay-differential
(`delay()` / `past()`) and `linCmt()` models, normal residual models
(add/prop/combined and transforms with fixed lambda), censored data
(M2/M3/M4 via CENS/LIMIT) and user-written `ll()` endpoints (both
verified: values tie to textbook densities up to a parameter-free
constant, gradients FD-verified), mu-referenced and non-mu etas,
mu-referenced covariates (subject-constant and time-varying), all 39
real-valued prior distributions in the lotri catalog (univariate,
multivariate normal families, LKJ/Wishart on omega blocks).

Mu-referenced covariates work for both subject-constant covariates (the
`cov_i * d/deta` scatter identity) and time-varying covariates (the
coefficient rides the forward-sensitivity model like any other
structural theta, the same way the other nlmixr2est methods treat a
time-varying regressor); both routes are FD-verified.

Estimated transform-both-sides lambda (Box-Cox / Yeo-Johnson) is
supported: the linked conditional supplies the transformed-scale density
with an exact d/dlambda column (FD-verified at ~1e-10), and the
generator adds the DV-transform Jacobian Stan-side as
`target += (lambda - 1) * sumLogJac` – a pure data statistic, so
lambda’s full gradient is exact end to end (the assembled target is
verified against the untransformed-scale density in difference form).

Finite mixtures (`mix()`, 2 components) are supported: the linked batch
evaluates the component-conditional likelihoods in nlmixr2est’s
component-major layout, the generator marginalizes with the Stan Users
Guide `log_sum_exp` pattern (component-specific etas whose priors factor
out; the mixing probability’s gradient is pure Stan autodiff – the
conditional is provably p-free, FD-verified), and the per-subject
membership posteriors land in `fit$env$mixProb`. The whole assembled
mixture target is FD-verified through `grad_log_prob`.

IOV is supported: `iov.x ~ v | OCC` expands (via nlmixr2est’s
preprocessing hook) into per-occasion fixed unit-variance etas scaled by
an estimated magnitude theta.

Refused with an explanatory error (not silently wrong): mixtures with
more than 2 components (for now), and the 8 discrete distributions
(every `ini({})` parameter is real-valued).

WAIC/LOO (subject-level, from the per-subject conditional log-likelihood
draws) are implemented – see the “objective-function rows” section of
the `nlmixr2bayes` vignette. See the issue tracker for what’s still
open.
