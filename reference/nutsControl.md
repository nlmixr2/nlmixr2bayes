# Control options for \`est="nuts"\` (Stan's NUTS sampler)

\`nutsControl()\` is \[stanControl()\] with \`algorithm = "NUTS"\`
forced – the same options, named after the SAMPLING ALGORITHM instead of
the software tool, the way nlmixr2est's \`foceControl()\` is
\`foceiControl()\` with \`interaction = FALSE\` forced. Because
\`class(nutsControl())\` is \`"nutsControl"\`, \`nlmixr2(model, data,
nutsControl())\` infers \`est="nuts"\` with no \`est=\` argument.

## Usage

``` r
nutsControl(..., algorithm = "NUTS")
```

## Arguments

- ...:

  \[stanControl()\] options

- algorithm:

  must be \`"NUTS"\` (present so a \[stanControl()\] can be converted to
  a \`nutsControl\`; a variational algorithm errors)

## Value

a \`nutsControl\` list

## Details

\`"NUTS"\` is already \[stanControl()\]'s default, so this is mostly a
naming convenience; asking for a variational algorithm here is an error
rather than a silent override.

## See also

\[stanControl()\], \[adviControl()\], \[pathfinderControl()\]

## Author

Matthew L Fidler

## Examples

``` r
# nutsControl() is stanControl() with algorithm pinned to NUTS, so
# nlmixr2(model, data, nutsControl()) needs no est= argument
ct <- nutsControl(iter = 500L)
ct$algorithm
#> [1] "NUTS"
ct$iter
#> [1] 500
```
