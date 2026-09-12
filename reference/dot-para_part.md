# One parametric term's contribution, as a linear form in the coefficients

The same object a smooth panel draws, built the same way: a design for
the term alone, evaluated along its own covariate with the others at a
typical value, times the coefficients it multiplies. Its columns are
centred at their means over the fitting data, as
[`stats::predict.lm()`](https://rdrr.io/r/stats/predict.lm.html) does
for `type = "terms"`, which fixes the free constant the intercept would
otherwise absorb and leaves a form whose standard errors are those of a
contrast rather than of an arbitrary level.

## Usage

``` r
.para_part(fit, p, k, v, values)
```

## Arguments

- fit:

  A `gamRTMB` fit.

- p, k:

  Distributional parameter and parametric term index.

- values:

  The covariate values to evaluate at.

## Value

`list(Z, coef, ib)`, with `ib` indexing `beta`.
