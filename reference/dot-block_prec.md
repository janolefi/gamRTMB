# A penalized block's precision matrix, at given parameters

The one place that knows how a block's parameters become a precision, so
that the objective and
[`edf()`](https://janolefi.github.io/gamRTMB/reference/edf.md) cannot
drift apart. Written to work both on an RTMB tape and on plain numbers.

## Usage

``` r
.block_prec(bl, theta)
```

## Arguments

- bl:

  A block from
  [`.build_design()`](https://janolefi.github.io/gamRTMB/reference/dot-build_design.md).

- theta:

  That block's parameters, `bl$ntheta` of them.

## Details

`"iid"` is the
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
basis, where the penalty is the identity. `"gmrf"` keeps a fixed sparse
precision and scales it. `"multi"` combines several penalty matrices
through mgcv's `L` convention, \\\lambda = \exp(L\theta)\\ and \\Q =
\sum_i \lambda_i S_i\\, which is how a Matern SPDE writes
\\\tau^2(\kappa^4 C + 2\kappa^2 G_1 + G_2)\\ with \\\theta = (\log\tau,
\log\kappa)\\.
