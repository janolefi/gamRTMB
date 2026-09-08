
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gamRTMB

<!-- badges: start -->

<!-- badges: end -->

Distributional (GAMLSS-style) regression where smooth terms can enter
**every** parameter of a distribution, not just the mean.

The package is mostly glue, by design:

- **mgcv** builds the bases and penalties, and `mgcv::smooth2random()`
  reparameterises the penalized coefficients as iid Gaussian random
  effects;
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

A Gaussian location-scale model: the mean varies smoothly with `x1`, and
so does the standard deviation.

``` r
library(gamRTMB)

set.seed(1)
n <- 500
d <- data.frame(x1 = runif(n), x2 = runif(n))
d$y <- rnorm(n,
             mean = sin(2 * pi * d$x1) + d$x2^2,
             sd   = exp(-1 + 0.8 * cos(2 * pi * d$x1)))

fit <- gamRTMB(y ~ list(mu    = ~ s(x1) + s(x2),
                        sigma = ~ s(x1)),
               family = gaussian_ls(), data = d)
fit
#> gamRTMB fit
#>   family:    gaussian_ls (mu/identity, sigma/log)
#>   criterion: REML   engine: laplace
#>   converged: TRUE   -logLik: 258.4182   max|grad|: 9.58e-06
#>   coefficients: 5 fixed (incl. null spaces), 24 penalized; 3 smoothing parameters
```

Effective degrees of freedom per smooth, as in a `gam` summary:

``` r
edf(fit)
#>   parameter  term      edf k     sp id
#> 1        mu s(x1) 7.482610 9 0.1192   
#> 2        mu s(x2) 3.688019 9  9.331   
#> 3     sigma s(x1) 5.710145 9 0.1814
```

Predictions come back as one vector per distributional parameter. New
data is pushed through `mgcv::PredictMat()` on the stored smooths, so
the basis and knots are exactly those of the fit:

``` r
grid <- data.frame(x1 = seq(0, 1, length.out = 5), x2 = 0.5)
predict(fit, newdata = grid, type = "response")
#> $mu
#> [1] -0.07480724  1.23104173  0.25907607 -0.78599712  0.25564047
#> 
#> $sigma
#> [1] 0.9104810 0.4027244 0.1655772 0.3676115 0.8616064
```

### Any RTMBdist family

`rtmbdist_family()` derives the parameter names, links and starting
values from the density itself, so most of RTMBdist works without a
hand-written family. Here a skew normal with a smooth on all three
parameters at once:

``` r
fam <- rtmbdist_family("skewnorm2")
fam
#> gamRTMB family: skewnorm2  [dskewnorm2]
#>   modelled: mean (identity), sd (log), alpha (identity)

set.seed(2)
n <- 1500
d2 <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
d2$y <- RTMBdist::rskewnorm2(n,
  mean  = sin(2 * pi * d2$x1),
  sd    = exp(-0.2 + 0.5 * cos(2 * pi * d2$x2)),
  alpha = 2.5 * sin(2 * pi * d2$x3))

fit2 <- gamRTMB(y ~ list(mean = ~ s(x1), sd = ~ s(x2), alpha = ~ s(x3)),
                family = fam, data = d2)
edf(fit2)
#>   parameter  term      edf k       sp id
#> 1      mean s(x1) 7.429360 9   0.1494   
#> 2        sd s(x2) 5.522424 9   0.6772   
#> 3     alpha s(x3) 4.317055 9 0.003831
```

Parameters that are known data rather than quantities to model — the
number of binomial trials, truncation bounds — are declared separately
and are not reachable from the formula:

``` r
rtmbdist_family("betabinom", fixed = list(size = "trials"))
#> gamRTMB family: betabinom  [dbetabinom]
#>   modelled: shape1 (log), shape2 (log)
#>   fixed:    size
```

### Shared smoothing parameters

`s(..., id = )` ties smooths to one smoothing parameter, exactly as in
mgcv (including the linked bases, so the shared variance means the same
thing for each term). Unlike mgcv, the group may span distributional
parameters:

``` r
fit3 <- gamRTMB(y ~ list(mu    = ~ s(x1, id = "sh"),
                         sigma = ~ s(x1, id = "sh")),
                family = gaussian_ls(), data = d)
fit3
#> gamRTMB fit
#>   family:    gaussian_ls (mu/identity, sigma/log)
#>   criterion: REML   engine: laplace
#>   converged: TRUE   -logLik: 423.9595   max|grad|: 8.13e-08
#>   coefficients: 4 fixed (incl. null spaces), 16 penalized; 2 smoothing parameters (1 free, 1 tied by id)
edf(fit3)
#>   parameter  term      edf k     sp id
#> 1        mu s(x1) 6.585786 9 0.1523 sh
#> 2     sigma s(x1) 5.964221 9 0.1523 sh
```

## Confidence bands

Ask for the joint precision matrix at fit time and `predict()` will
return standard errors, per term or for the whole linear predictor:

``` r
fitb <- gamRTMB(y ~ list(mu = ~ s(x1) + s(x2), sigma = ~ s(x1)),
                family = gaussian_ls(), data = d, joint_precision = TRUE)

g <- data.frame(x1 = seq(0, 1, length.out = 200), x2 = 0.5)
tm <- predict(fitb, newdata = g, type = "terms", se.fit = TRUE)

par(mar = c(4, 4, 1, 1))
plot(g$x1, tm$mu$fit[, 1], type = "l", lwd = 2, ylim = c(-2, 2),
     xlab = "x1", ylab = "s(x1)")
polygon(c(g$x1, rev(g$x1)),
        c(tm$mu$fit[, 1] + 2 * tm$mu$se[, 1],
          rev(tm$mu$fit[, 1] - 2 * tm$mu$se[, 1])),
        col = adjustcolor("steelblue", 0.25), border = NA)
lines(g$x1, sin(2 * pi * g$x1) - mean(sin(2 * pi * d$x1)),
      col = 2, lty = 2, lwd = 2)
legend("topright", c("fitted", "truth"), col = c(1, 2), lty = c(1, 2), bty = "n")
```

<img src="man/figures/README-bands-1.png" width="100%" />

## Status and scope

Early but checked. EDF and fitted values agree with
`mgcv::gam(family = gaulss())` to three decimals, with and without
shared smoothing parameters, and `predict()` on the fitting data
reproduces the in-sample linear predictors to machine precision.

Supported: `s()`, `t2()`, `by =` variables, `bs = "fs"`, `bs = "re"`,
and `id =`. Not supported: `te()` (mgcv itself declines
`smooth2random(type = 2)` for it — use `t2()`) and `fx = TRUE` (no
penalized part to make random).

An extended Fellner–Schall fitting engine is planned as an alternative
to the Laplace one; `dev/NOTES-fellner-schall.md` records the derivation
and the seams it needs.
