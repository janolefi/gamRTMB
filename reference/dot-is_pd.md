# Is a sparse symmetric matrix positive definite?

[`Matrix::chol()`](https://rdrr.io/pkg/Matrix/man/chol-methods.html)
rather than
[`Matrix::Cholesky()`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html):
the latter computes an LDL' factorisation, which exists perfectly well
for an indefinite matrix and returns with nothing worse than a CHOLMOD
warning, so it answers a different question from the one being asked
here. [`chol()`](https://rdrr.io/r/base/chol.html) fails on an
indefinite matrix, which is the answer wanted.

## Usage

``` r
.is_pd(h)
```

## Arguments

- h:

  A symmetric sparse matrix, or `NULL`.

## Value

`TRUE`, `FALSE`, or `NA` if `h` is missing or has non-finite entries, so
the question cannot be put.
