# Map reparameterised coefficients back to a smooth's own basis

[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
returns penalized design blocks (`rand`), an unpenalized null-space
block (`Xf`), and the transformation back to the smooth's original
constrained basis: \$\$\beta = U (D \cdot c(b\_{rand},
b\_{fixed})\[idx\])\$\$ with `idx = c(rind, offset + seq_len(n_fixed))`.

## Usage

``` r
.reconstruct_map(re)
```

## Details

For ordinary smooths `rind` is the identity, but for `bs = "fs"` it is a
genuine permutation: the random blocks are penalty-major while the
coefficients are level-major. Assuming `cbind(rand, Xf)` therefore gives
silently wrong coefficients for factor-smooth interactions, so the map
is built once here and reused. Returning it as an explicit matrix means
the same object serves the point reconstruction and the delta-method
covariance of a smooth.

Verified exact (to 6e-16) for `s()`, `s(bs = "cr")`, `t2()`, `by =`
factors, `bs = "fs"` and `bs = "re"`.
