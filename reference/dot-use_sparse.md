# Is this penalty worth keeping sparse?

Two conditions. The smooth must have exactly one penalty, since several
penalties on one coefficient block is what the whole `smooth2random`
machinery exists to handle. And the penalty must be big and sparse
enough that the rotation would actually cost something: a 10-coefficient
P-spline penalty is technically banded, but nothing is gained by
treating it specially, and the well-travelled route is the safer one.

## Usage

``` r
.use_sparse(S, mode)
```

## Arguments

- S:

  The penalty matrix.

- mode:

  `"auto"`, `"always"` or `"never"`.
