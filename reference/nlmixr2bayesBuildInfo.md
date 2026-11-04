# Probe the Stan build environment

Runs \`rstan::stanc(allow_undefined=TRUE)\` on a fixture and records
three facts about the generated C++ (cached per session):

## Usage

``` r
nlmixr2bayesBuildInfo(force = FALSE)
```

## Arguments

- force:

  re-run the probe even if cached

## Value

a list with \`ok\`, \`callShape\`, \`injectionPoint\`, \`stancVersion\`,
\`rstanVersion\`; \`ok=NA\` when rstan is not installed

## Details

1\. \`callShape\`: stanc3 still emits the external call as
\`nlmixr2_cond_all(\<expr\>, pstream\_\_)\` with no template-argument
list – the shape \`inst/include/nlmixr2bayes_lp.hpp\`'s concrete
overloads bind to. 2. \`injectionPoint\`: the first \`class \<ident\>\`
in the generated code sits inside \`namespace \<model\>\_namespace\`,
which is where \`rstan::stan_model(includes=)\` splices the header. 3.
\`stancVersion\`/\`rstanVersion\` for the record.

## Author

Matthew L Fidler

## Examples

``` r
# what this build of nlmixr2bayes can compile against
info <- nlmixr2bayesBuildInfo()
info$ok
#> [1] TRUE
info$rstanVersion
#> [1] "2.32.7"
```
