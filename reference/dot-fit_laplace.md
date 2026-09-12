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
  family,
  inner_control = list(),
  repars = NULL
)
```

## Arguments

- inner_control:

  Passed to
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)'s
  `inner.control`.

- repars:

  Function of a `sigma_frac` returning a fresh starting parameter list,
  used to walk the ladder. `NULL` disables the retries.

## Details

Smoothing parameters tied by an `id` are collapsed through
[`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)'s
`map`: several entries of `log_sigma` become one estimated value. An
`NA` level there would instead fix an entry at its starting value, which
is how a user-specified smoothing parameter would be implemented.

## Getting started at all

Before any of that, the marginal objective has to be finite at the
starting values, and on the harder families it often is not. That is
handled in two steps, both of which are about failing cheaply rather
than about finding a better start:
[`.probe_finite()`](https://janolefi.github.io/gamRTMB/reference/dot-probe_finite.md)
asks the question with a short inner iteration cap, and a `FALSE` sends
the caller to the next rung of `.sigma_frac_ladder()`. If every rung
fails, the fit proceeds from the original starting values under the full
cap, so the ladder can only add fits, never remove one.
