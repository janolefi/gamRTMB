# Summarise a gamRTMB fit

Laid out like
[`mgcv::summary.gam()`](https://rdrr.io/pkg/mgcv/man/summary.gam.html),
but with one row per distributional parameter throughout, since every
parameter has its own formula, link, coefficients and smooths.

## Usage

``` r
# S3 method for class 'gamRTMB'
summary(object, ...)

# S3 method for class 'summary.gamRTMB'
print(x, ...)
```

## Arguments

- object:

  A `gamRTMB` fit.

- ...:

  Ignored.

- x:

  A `summary.gamRTMB`.

## Value

An object of class `summary.gamRTMB`.

## Details

Parametric coefficients get standard errors from
[`vcov.gamRTMB()`](https://janolefi.github.io/gamRTMB/reference/vcov.gamRTMB.md)
and approximate Gaussian (z) tests. Smooth terms are reported with their
effective degrees of freedom and smoothing parameters but **no
p-values**: mgcv's "approximate significance" rests on Wood's (2013)
test for a term's whole coefficient block, which is not implemented
here, and a naive Wald test in its place would be misleading.

## Examples

``` r
set.seed(1)
d <- data.frame(x1 = runif(300), x2 = runif(300))
d$y <- rnorm(300, sin(2 * pi * d$x1), exp(-1 + d$x2))
summary(gamRTMB(y ~ list(mean = ~ s(x1), sd = ~ s(x2)), data = d))
#> 
#> Family: norm   [dnorm from RTMB]
#> Links:  mean = identity,  sd = log
#> 
#> Formula:
#>   mean ~ s(x1)
#>     sd ~ s(x2)
#> 
#> Parametric coefficients:
#>                  Estimate Std. Error z value Pr(>|z|)    
#> mean:(Intercept)  0.05266    0.03374   1.561    0.119    
#> sd:(Intercept)   -0.45491    0.04092 -11.118   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Smooth terms:
#>        term   edf k        sp
#>  mean:s(x1) 6.377 9    0.1225
#>  sd:s(x2)   1.000 9 2.828e+08
#> 
#> Total EDF = 9.38   n = 300
#> -REML = 309.343   logLik = -289.210   AIC = 597.17
```
