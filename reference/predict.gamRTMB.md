# Predictions from a gamRTMB fit

New data goes through
[`mgcv::PredictMat()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html) on
the stored `smoothCon` objects, so the basis, knots and constraint
matrices are exactly those of the fit. Smooths are never re-fitted on
new data, which would place different knots.

## Usage

``` r
# S3 method for class 'gamRTMB'
predict(
  object,
  newdata = NULL,
  type = c("link", "response", "terms", "quantile"),
  prob = c(0.05, 0.25, 0.5, 0.75, 0.95),
  se.fit = FALSE,
  ...
)
```

## Arguments

- object:

  A `gamRTMB` fit.

- newdata:

  Optional data frame. If omitted, the fitting data's in-sample linear
  predictors are returned.

- type:

  `"link"` (default) for the linear predictors, `"response"` for the
  parameters on their natural scale, `"terms"` for the contribution of
  each smooth separately, or `"quantile"` for quantiles of the fitted
  distribution. As in mgcv, `"terms"` excludes the intercept and any
  offset, which belong to the predictor rather than to a smooth.

- prob:

  Probabilities for `type = "quantile"`.

- se.fit:

  Also return standard errors; needs a fit made with
  `joint_precision = TRUE`.

- ...:

  Ignored.

## Value

For `type = "link"`/`"response"`, a named list with one vector per
distributional parameter. For `type = "terms"`, a named list with a
matrix of per-term contributions for each parameter. For
`type = "quantile"`, a matrix with one column per probability. With
`se.fit = TRUE`, a list of `fit` and `se.fit` in the same shape.

## Quantiles

`type = "quantile"` evaluates the fitted distribution's quantile
function at each observation, giving covariate-dependent quantiles —
which is much of the point of letting every parameter vary. It needs the
family to have a quantile function;
[`families()`](https://janolefi.github.io/gamRTMB/reference/families.md)
reports which do.

Their standard errors come from the delta method: the quantile is
differentiated numerically with respect to each linear predictor, and
those derivatives are combined with the predictors' joint covariance,
including the covariance *between* distributional parameters. Offered
only for a continuous response, because for a lattice or mixed one the
quantile function is a step and its derivative is not meaningful.
