# Build the \`includes=\` string for \[rstan::stan_model()\]

rstan splices \`includes\` immediately before the model class, which is
\*inside\* \`model\_\<name\>\_namespace\`. The bridge needs global scope
(it defines a namespace, template specializations, and static
function-pointer caches), so the returned text closes that namespace,
emits the bridge at global scope, and reopens it with a \`using\`
declaration that makes \`rx_solve\` visible to the generated unqualified
call.

## Usage

``` r
rxsStanIncludes(modelName, tag = NULL, pre = NULL, post = NULL)
```

## Arguments

- modelName:

  model name passed to \[rstan::stanc()\]; must have been compiled with
  \`obfuscate_model_name = FALSE\` so the namespace is predictable.

- tag:

  optional C++ tag type selecting a \`rxstan::solve_policy\`
  specialization. Defaults to rxstan's analytic policy; pass
  \`"::rxstan::finite_diff_tag"\` to compile the same program against
  central differences instead.

- pre:

  C++ emitted at global scope \*before\* the bridge header, e.g. the
  declaration of a custom tag type.

- post:

  C++ emitted at global scope \*after\* the bridge header, e.g. a
  generated \`rxstan::solve_policy\` specialization.

## Value

a single string

## Details

Ordering matters and is why \`pre\` and \`post\` are separate: a custom
tag type must be declared \*before\* the header (the primary template is
named in \`rx_solve\`'s body, which is bound at definition time), while
a \`solve_policy\` specialization for it can only be written \*after\*
the primary template exists.

## Author

Lukas A. Widmer

## Examples

``` r
# the includes= string rstan::stan_model() needs for a bridge model
inc <- rxsStanIncludes("myModel")
nchar(inc)
#> [1] 212
```
