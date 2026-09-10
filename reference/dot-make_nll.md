# Build the joint negative log-likelihood

The penalized coefficients carry a mean-zero Gaussian prior, one
variance per penalized block. Which Gaussian depends on the route the
block took through
[`.build_design()`](https://janolefi.github.io/gamRTMB/reference/dot-build_design.md):
after
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
the penalty is the identity and the prior is iid \\N(0, \sigma_k^2)\\,
while a block that kept its own sparse penalty gets \\N(0, \sigma_k^2
Q_k^{-1})\\ through
[`RTMB::dgmrf()`](https://rdrr.io/pkg/RTMB/man/MVgauss.html). The two
are the same model written two ways; keeping \\Q_k\\ sparse is what lets
a Markov random field over many regions stay affordable, since
`smooth2random`'s rotation would fill it in.

## Usage

``` r
.make_nll(design, family, y, fx = list(), w = NULL)
```

## Arguments

- design:

  From
  [`.build_design()`](https://janolefi.github.io/gamRTMB/reference/dot-build_design.md).

- family:

  A `gamRTMB_family`.

- y:

  Response.

- fx:

  Resolved fixed arguments.

- w:

  Prior weights, or `NULL` for unweighted.

## Value

A function of a parameter list `list(beta, b, log_sigma)`.

## Details

Null-space coefficients live in `beta` and are never given a prior.
Nothing is shared across smooths or across distributional parameters
unless an `id` says so.

Prior weights multiply each observation's log-density contribution, as
in [`stats::glm()`](https://rdrr.io/r/stats/glm.html). Offsets are added
to the relevant parameter's linear predictor.
