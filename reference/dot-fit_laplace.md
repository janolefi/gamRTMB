# Laplace engine

Declares the coefficients random so that RTMB supplies the marginal
criterion, then optimises the variance components with `nlminb`.

## Usage

``` r
.fit_laplace(
  nll,
  pars,
  design,
  method,
  joint_precision,
  silent,
  control,
  family
)
```

## Details

Smoothing parameters tied by an `id` are collapsed through
[`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)'s
`map`: several entries of `log_sigma` become one estimated value. An
`NA` level there would instead fix an entry at its starting value, which
is how a user-specified smoothing parameter would be implemented.
