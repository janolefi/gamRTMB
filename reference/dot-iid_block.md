# Build an iid block from a smooth, the `smooth2random` way

The other half of
[`.gmrf_block()`](https://janolefi.github.io/gamRTMB/reference/dot-gmrf_block.md),
with the same shape of answer, so that
[`.build_design()`](https://janolefi.github.io/gamRTMB/reference/dot-build_design.md)
can pick a route and then stop caring which it picked. Here the
reparameterisation makes the penalty the identity, so the block needs no
precision matrix and exactly one variance.

## Usage

``` r
.iid_block(sm)
```

## Arguments

- sm:

  A `smoothCon` object built with `absorb.cons = TRUE`.
