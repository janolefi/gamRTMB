# One smooth's design in joint-coefficient space

A smooth's fitted values and their standard errors are both linear forms
in the joint coefficient vector once its model matrix has been mapped
through the reparameterisation: with `Z = X Tmap`, the contribution is
`Z c` and its covariance `Z V Z'`. Mapping the design once is simpler
than mapping the coefficients and their covariance separately, and it
removes the need to assemble a block-diagonal transform for the whole
linear predictor.

## Usage

``` r
.smooth_part(fit, p, j, X)
```

## Arguments

- fit:

  A `gamRTMB` fit.

- p, j:

  Distributional parameter and smooth index.

- X:

  The smooth's model matrix, from the fit or from
  [`mgcv::PredictMat()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html).

## Value

`list(Z, coef, b, f)`: the mapped design, the coefficients it
multiplies, and the `b` and `beta` indices for locating them in a
covariance matrix.
