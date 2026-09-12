# The drawable first-order parametric terms of one parameter

Everything a `termplot` would draw: a term of order one in a single
variable, so that `poly(x, 2)` and `log(x)` count but `x:z` does not.
The intercept is not a term and is skipped; it is not a contribution to
the predictor that varies with anything.

## Usage

``` r
.para_terms(part)
```

## Value

A list of `list(k, var, label)`, `k` indexing the term within
`terms(...)`, with `var = NULL` for a term that cannot be drawn.
