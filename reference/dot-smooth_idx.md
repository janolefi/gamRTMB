# Where one smooth's coefficients live in the joint vectors

Four places need the same thing: which entries of `b` and which entries
of `beta` belong to smooth `j` of parameter `p` (its penalized blocks,
then its null-space columns). Gathering it once keeps
[`edf()`](https://janolefi.github.io/gamRTMB/reference/edf.md), the
coefficient reconstruction and both covariance helpers in step.

## Usage

``` r
.smooth_idx(design, p, j)
```

## Value

`list(b, f)`, indices into the `b` and `beta` vectors.
