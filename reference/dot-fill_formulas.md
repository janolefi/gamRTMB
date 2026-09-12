# Complete and order a set of per-parameter formulas

Every parameter the family declares gets a one-sided formula in the
family's own order, `~1` where the user gave none, all sharing the outer
formula's environment so that a variable resolves from where the model
was written rather than from where each piece happened to be built.

## Usage

``` r
.fill_formulas(out, parnames, env)
```
