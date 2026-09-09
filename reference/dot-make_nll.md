# Build the joint negative log-likelihood

The penalized coefficients carry an iid \\N(0, \sigma_k^2)\\ prior, one
variance per penalized block. Null-space coefficients live in `beta` and
are never given a prior. Nothing is shared across smooths or across
distributional parameters unless an `id` says so.

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

Prior weights multiply each observation's log-density contribution, as
in [`stats::glm()`](https://rdrr.io/r/stats/glm.html). Offsets are added
to the relevant parameter's linear predictor.
