# QQ and worm plots

A worm plot is the detrended QQ plot, so both are the same picture with
the theoretical quantile subtracted or not. The reference is the
identity line rather than a fitted one, because these residuals should
be standard normal and not merely normal.

## Usage

``` r
.plot_res(x, nsim, detrend, se, band.col, ...)
```
