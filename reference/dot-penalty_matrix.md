# The penalty, assembled in the Hessian's own ordering

\\S\\ is block diagonal with one block per penalized block of
coefficients: \\\sigma_k^{-2} I\\ for a block that went through
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html),
and \\\sigma_k^{-2} Q_k\\ for one that kept a sparse penalty. Built
sparse so that `solve(H, S)` stays a sparse solve rather than a full
inverse.

## Usage

``` r
.penalty_matrix(design, ls, ph)
```
