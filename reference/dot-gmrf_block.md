# Build a sparse GMRF block from a smooth

Returns the same three things the
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
route returns – a design matrix on the penalized coefficients, an
unpenalized null-space matrix for `beta`, and the map back to the
smooth's own basis – so that everything downstream is unchanged, plus
the pieces the prior needs.

## Usage

``` r
.gmrf_block(sm)
```

## Arguments

- sm:

  A `smoothCon` object built with `absorb.cons = FALSE`.

## Details

Null-space columns whose contribution the parameter's intercept already
covers are left out, which is what the sum-to-zero constraint achieves
on the other route. A pivoted QR of `cbind(1, X N)` finds them: for a
connected Markov random field the whole null space is the constant, so
no free column survives and the intercept carries the level; for a field
with an island, one contrast between the two components survives, which
is right, because their levels really are separately free.

A smooth with an `L` matrix keeps all of its penalty matrices and gets
`ncol(L)` parameters instead of one. Its penalty is assumed proper – an
SPDE precision is, for any positive range – so no constraint arises.
