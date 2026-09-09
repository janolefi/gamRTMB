# The penalized Hessian, in the design's own coefficient ordering

The engine contract: return the Hessian of the joint penalized negative
log-likelihood in the coefficients at the fitted mode, plus index
vectors giving the row of `H` for each entry of `beta` and of `b`.
Anything needing curvature —
[`edf()`](https://janolefi.github.io/gamRTMB/reference/edf.md),
smooth-term covariances — uses only this.

## Usage

``` r
.penalized_hessian(fit)
```

## Arguments

- fit:

  A `gamRTMB` fit.

## Value

`list(H, i_beta, i_b)`, or `NULL` if the engine cannot supply it.
