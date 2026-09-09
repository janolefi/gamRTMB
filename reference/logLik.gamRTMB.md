# Log-likelihood, and hence AIC and BIC

Returns the **log-likelihood** of the data at the fitted coefficients,
with `df` set to the total effective degrees of freedom. That is
deliberately not the fitting criterion: under REML the criterion is a
restricted one, is not a likelihood, and is not comparable across
different mean structures, so using it for AIC would be wrong. It
remains available as `fit$objective`.

## Usage

``` r
# S3 method for class 'gamRTMB'
logLik(object, ...)
```

## Arguments

- object:

  A `gamRTMB` fit.

- ...:

  Ignored.

## Value

An object of class `logLik`.

## Details

With `df` = total EDF this is the convention used by
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) and by GAMLSS's
GAIC, so [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) are comparable across models
fitted to the same response. Under `method = "ML"`, where effective
degrees of freedom are unavailable, `df` falls back to counting the
fixed coefficients and free smoothing parameters.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(200)); d$y <- rnorm(200, sin(2 * pi * d$x), 0.3)
fit <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d)
logLik(fit)
#> 'log Lik.' -34.98567 (df=8.392431)
AIC(fit)
#> [1] 86.7562
```
