
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gamRTMB

<!-- badges: start -->

[![R-CMD-check](https://github.com/janolefi/gamRTMB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janolefi/gamRTMB/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Distributional (GAMLSS-style) regression where smooth terms can enter
**every** parameter of a distribution, not just the mean.

The package is mostly glue, by design:

- **mgcv** builds the bases and penalties, and `mgcv::smooth2random()`
  reparameterises the penalized coefficients as iid Gaussian random
  effects — except where the penalty is already a sparse precision
  matrix, which keeps it and uses `RTMB::dgmrf()` instead;
- **RTMB** supplies automatic differentiation and the Laplace
  approximation;
- **RTMBdist** supplies the log-densities, each in its own **native**
  parameterisation — a skew normal is modelled as `xi`, `omega`,
  `alpha`, not as generic location/scale/shape.

Restricted maximum likelihood is the default, effective degrees of
freedom per smooth agree with mgcv, and smoothing parameters can be
shared between terms *and* between distributional parameters.

## Installation

``` r
# install.packages("pak")
pak::pak("janolefi/gamRTMB")
```

## Example

`MASS::mcycle` measures head acceleration in a simulated motorcycle
crash. The mean is famously hard to model — and the variance changes
just as much.

``` r
library(gamRTMB)
data(mcycle, package = "MASS")

fit <- gamRTMB(accel ~ list(mean = ~ s(times, bs="ad", k=30, m=4),
                            sd   = ~ s(times, k = 10)),
               data = mcycle, family = fam("norm"))
summary(fit)
#> 
#> Family: norm   [dnorm from RTMB]
#> Links:  mean = identity,  sd = log
#> 
#> Formula:
#>   mean ~ s(times, bs = "ad", k = 30, m = 4)
#>     sd ~ s(times, k = 10)
#> 
#> Parametric coefficients:
#>                  Estimate Std. Error z value Pr(>|z|)    
#> mean:(Intercept)  2.92738    9.32615   0.314    0.754    
#> sd:(Intercept)    2.57879    0.06448  39.994   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Smooth terms:
#>           term    edf  k                                sp
#>  mean:s(times) 13.992 30 0.05135,1.131e-06,1.131e-06,0.309
#>  sd:s(times)    7.228  9                          0.008796
#> 
#> Total EDF = 23.22   n = 133
#> -REML = 582.263   logLik = -531.697   AIC = 1109.83
```

Once every parameter varies, the useful output is covariate-dependent
quantiles rather than a mean and one standard error:

``` r
plot(fit, type = "quantile")
```

<img src="man/figures/README-quantiles-1.png" alt="Fitted 5th, 25th, 50th, 75th and 95th percentile curves of head acceleration against time, over the mcycle data. The curves are narrow before impact, fan out widely through it, then contract again." width="100%" />

Term plots, a worm plot of randomised quantile residuals, and everything
else are in `vignette("gamRTMB")`.

## What it does

|  |  |
|----|----|
| Families | 85, from RTMBdist plus the standard R densities; `families()` lists them, `fam()` builds one |
| Formula | `y ~ list(mean = ~ s(x), sd = ~ s(z))`, one one-sided formula per parameter, missing ones get `~1` |
| Smooths | `s()`, `te()`, `ti()`, `t2()`, `by=`, `bs="fs"`, `bs="re"`, `bs="ad"`, shrinkage bases, and `id=` to share smoothing parameters |
| Spatial | `bs="mrf"` Markov random fields over an adjacency graph, or any precision matrix via `xt=list(penalty=)`; `bs="spde"` Matern fields on an [fmesher](https://cran.r-project.org/package=fmesher) mesh. Both kept sparse, so a few thousand regions or mesh nodes is routine |
| Criterion | `REML` by default (coefficients integrated out by the same Laplace approximation), `ML`, or `qREML` — the REML criterion optimised by extended Fellner–Schall, which fits four-parameter families the other two cannot start |
| Inference | `summary()`, `vcov()`, `edf()`, `AIC()`/`BIC()`, `predict()` with standard errors |
| Diagnostics | `residuals()` gives randomised quantile residuals; `plot(type = "worm")` |
| Also | `weights`, `offset()` inside a parameter’s formula, `na.action` |

## Status

Early but checked. EDF and fitted values agree with
`mgcv::gam(family = gaulss())` to three decimals, with and without
shared smoothing parameters; `predict()` on the fitting data reproduces
the in-sample linear predictors to machine precision; and quantile
standard errors match both an analytic special case and simulation from
the joint posterior. The sparse GMRF route is a reparameterisation
rather than an approximation, and is checked against the dense one:
log-likelihood, EDF, AIC, fitted values, predictions and their standard
errors all agree, on Markov random fields and on ordinary smooths alike.
The SPDE precision matches the closed form
$\tau^2(\kappa^4 C + 2\kappa^2 G_1 + G_2)$ to 5e-16, and the fitted
range and marginal standard deviation approach their true values as the
sample grows.

Smooths that penalize one coefficient vector several times over —
`te()`, `ti()`, `bs = "ad"` — take the sparse route too, and for the
same reason a Markov random field does: `smooth2random()` needs one
variance per penalized block and has none to give. What makes that
legitimate is that the null space of $\sum_i \lambda_i S_i$ is the
intersection of the individual null spaces and so does not move with
$\lambda$, which means one corner constraint serves every penalty at
once and leaves a positive definite precision — so the generalised
determinant $|S_\lambda|_+$, and its notoriously delicate stable
evaluation, never has to be computed. EDF agree with `mgcv::gam` to 1e-3
and fitted values to 3e-3 on a `te()` over interacting and over additive
data, a `ti()`, and adaptive smooths on a jump and on `MASS::mcycle`.

Not supported: `fx = TRUE`, and smooths whose unpenalized null spaces
overlap (use a shrinkage basis, `bs = "ts"`). Not yet implemented:
smooth-term p-values and 2-D smooth plots.
