# Delta-method standard errors for fitted quantiles

A quantile is a smooth function of every linear predictor at once, so
its variance needs the predictors' *joint* covariance, cross-parameter
terms included: \\Var(q) = \sum\_{k,l} (\partial q / \partial \eta_k)
(\partial q / \partial \eta_l) Cov(\eta_k, \eta_l)\\.

## Usage

``` r
.quantile_se(object, prob, forms, Vj, fx = object$fixed, h = 1e-04)
```

## Arguments

- object:

  A `gamRTMB` fit.

- prob:

  Probabilities.

- forms:

  Per-parameter linear forms from
  [`predict.gamRTMB()`](https://janolefi.github.io/gamRTMB/reference/predict.gamRTMB.md).

- Vj:

  Joint coefficient covariance from
  [`.joint_cov()`](https://janolefi.github.io/gamRTMB/reference/dot-joint_cov.md).

- fx:

  Resolved fixed arguments for the rows being predicted.

- h:

  Step for the central difference, on the link scale.

## Value

A matrix of standard errors, one column per probability.

## Details

The derivatives are taken by central difference on the family's quantile
function, which avoids needing an analytic derivative for each of the
families that has one. The covariances come from the same linear forms
the per-parameter standard errors use, so `Cov(eta_k, eta_l)` is one
row-wise product per pair.
