# The null space of a penalty, and which coefficients to drop for it

An intrinsic GMRF – an ICAR field, a random walk – has a rank-deficient
precision: its density depends only on differences between coefficients,
so `k` directions are unpenalized. Both routes through the package have
to account for them, and they differ in how.

## Usage

``` r
.null_space(S, k)
```

## Arguments

- S:

  The penalty matrix.

- k:

  Dimension of its null space, from `smoothCon`'s `null.space.dim`.

## Value

`list(N, drop)`: a `q` by `k` null-space basis, and the `k` indices
whose coefficients are fixed at zero.

## Details

`mgcv` absorbs a sum-to-zero constraint into the basis with a dense
orthogonal rotation. That is what turns an ICAR penalty from 1% nonzero
into 100% nonzero, and it is why a Markov random field over a few
hundred regions is slow. Here the penalized part is corner-constrained
instead: `k` coefficients are fixed at zero and dropped, which leaves a
sparse submatrix, and the null space is handed to `beta` as unpenalized
columns exactly as
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
would.

The two are the same model. Write \\b = b_0 + N\gamma\\ with \\\gamma =
(N\_{drop})^{-1} b\_{drop}\\, so that \\b_0\\ vanishes on the dropped
indices. The decomposition is a bijection whenever \\N\_{drop}\\ is
invertible – which is what the pivoted QR below picks the dropped
indices to guarantee, and is the familiar fact that deleting one node of
a connected graph Laplacian leaves a positive definite matrix. Then \\Xb
= X\_{keep} b\_{0,keep} + (XN)\gamma\\ and, because \\SN = 0\\, \\b'Sb =
b\_{0,keep}' Q b\_{0,keep}\\. Same span, same penalty, sparse \\Q\\.

A proper penalty has no null space and gets no constraint: the prior
itself identifies every direction. That is the case a Matern/SPDE
precision lands in.

## Finding the null space

`k` comes from `mgcv`, which knows it by construction. Only a basis has
to be found here, and there is one case worth doing cheaply, because it
is every `bs = "mrf"` and every first-order random walk: when the
penalty's rows sum to zero the constant lies in the null space, and if
the graph then has exactly `k` connected components their indicators
*are* the null space. Reading those off the sparsity pattern costs a
traversal where an eigendecomposition would cost `q^3`. When the counts
disagree – a second difference penalty has a linear null direction as
well as a constant one – the general route is taken instead. A Cholesky
cannot stand in for this test: CHOLMOD and
[`chol()`](https://rdrr.io/r/base/chol.html) both factor a penalty whose
smallest eigenvalue is -8e-18 without complaint.
