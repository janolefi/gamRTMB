# Build the design for every distributional parameter

Build the design for every distributional parameter

## Usage

``` r
.build_design(par_formulas, data, parnames, knots = NULL, sparse = "auto")
```

## Arguments

- par_formulas:

  List of one-sided formulas, one per parameter.

- data:

  Model frame.

- parnames:

  The family's modelled parameters.

- knots:

  Passed to
  [`mgcv::smoothCon()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html).

- sparse:

  Whether a smooth with a single sparse penalty may skip
  [`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
  and keep that penalty; see
  [`.gmrf_block()`](https://janolefi.github.io/gamRTMB/reference/dot-gmrf_block.md).

## Value

A design object: per-parameter parts (parametric matrix, terms, xlevels,
smooths), the stacked fixed-effect matrices, index bookkeeping into
`beta` and `b`, the penalized blocks, and the variance-component
grouping.

## Identifiability

A smooth must not collide with its parameter's intercept. On the
`smooth2random` route `smoothCon(absorb.cons = TRUE)` handles that with
the sum-to-zero constraint; on the sparse route the constraint would
destroy the sparsity, so the term is corner-constrained instead and the
intercept takes the level – see
[`.null_space()`](https://janolefi.github.io/gamRTMB/reference/dot-null_space.md).
Either way the fixed-effect design is checked for rank at the end, since
two smooths can still share a null space between them.

## Null spaces

Only the penalized blocks become random effects. The unpenalized
null-space columns are appended to that parameter's fixed-effect matrix
and are never pooled into the random-effect sum.

## Offsets

An [`offset()`](https://rdrr.io/r/stats/offset.html) term inside a
parameter's formula adds a fixed, known contribution to that parameter's
linear predictor. Offsets are per-parameter because that is the only
meaningful reading in a distributional model: `sd = ~ offset(log(s))`
says something quite different from the same term on `mean`. The term is
recomputed from `newdata` when predicting.

## Shared smoothing parameters

`s(..., id = )` does two things in mgcv and both are reproduced. *Linked
bases*: the basis is built from the pooled covariate values of the whole
id group and then evaluated on each term's own data, via
`smoothCon(spec, pooled, n = nrow(data), dataX = data)`. This matters
beyond knot placement, because `scale.penalty` rescales each penalty by
a norm of its own model matrix; without pooling, two id-linked smooths
get different `S.scale` and one shared variance would mean different
amounts of smoothing for each. *Shared smoothing parameter*: blocks in a
group get the same `sig_key`, which the engine turns into one variance
component.

id groups are resolved globally across distributional parameters, so
`s(x, id = 1)` on `mu` and on `sigma` share one smoothness. A gam
formula cannot express this, having only one linear predictor.
