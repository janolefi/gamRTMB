# Covariance of the fixed-effect coefficients

The `(beta, beta)` block of the inverse joint precision, so it includes
the uncertainty in the smoothing parameters rather than conditioning on
them (mgcv's `unconditional = TRUE`). Needs a fit made with
`joint_precision = TRUE`, which is the default.

## Usage

``` r
# S3 method for class 'gamRTMB'
vcov(object, ...)
```

## Arguments

- object:

  A `gamRTMB` fit.

- ...:

  Ignored.

## Value

A labelled covariance matrix over the parametric and smooth null-space
coefficients.
