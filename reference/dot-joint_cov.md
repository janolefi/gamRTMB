# Joint covariance of the coefficients

The `(beta, b)` block of the inverse joint precision, which includes the
uncertainty in the smoothing parameters (mgcv's `unconditional = TRUE`).

## Usage

``` r
.joint_cov(fit)
```

## Arguments

- fit:

  A `gamRTMB` fit.

## Value

`list(V, ib, ir)`: the coefficient covariance, and the positions of the
`beta` and `b` entries within it.

## Details

Computed as a Schur complement rather than by inverting the whole
matrix. Writing the joint precision over coefficients `c` and log
smoothing parameters `s` as `[[Qcc, Qcs], [Qsc, Qss]]`, the block needed
is \$\$\[Q^{-1}\]\_{cc} = (Q\_{cc} - Q\_{cs} Q\_{ss}^{-1}
Q\_{sc})^{-1},\$\$ which is algebraically identical to inverting the
whole thing but isolates the awkward part.

Boundary smoothing parameters make this delicate in two separate ways,
and both are handled here because a plain
[`solve()`](https://rdrr.io/r/base/solve.html) of the whole matrix fails
on either, taking the standard errors, bands and summary with it.

**A flat smoothing parameter.** A term shrunk onto its null space leaves
the criterion flat in its own `log_sigma`, so the joint precision is
genuinely singular — but only in `Qss`: with two such terms the
offending eigenvalues were 3e-12 and 2e-07, loading on `log_sigma` with
weight 1.00, while the coefficient block stayed invertible. A
pseudo-inverse of `Qss` drops exactly those directions, which amounts to
treating a boundary smoothing parameter as known rather than estimated —
the conditional treatment, for that parameter only. Every other one
still contributes.

**Scaling.** The same boundary puts precision entries of order 1e19 next
to entries of order 1, and at that dynamic range double precision loses
positive definiteness outright: a factor-smooth model measured an
eigenvalue of -7.5e3 in a matrix whose largest was 2e19. So the matrix
is first scaled to a unit diagonal, which leaves only the correlation
structure to invert, and the result is unscaled afterwards.
