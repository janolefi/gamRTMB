# Is the marginal objective finite, cheaply?

The first marginal evaluation is where a badly posed inner problem shows
up, and it is an expensive place to discover it: the inner Newton runs
all the way to its iteration cap before handing back a non-finite value.
On the four-parameter Box-Cox families that is over two minutes to learn
that the fit will not start. Probing with a short cap answers the same
question in seconds, and costs nothing on a model that was going to
work.

## Usage

``` r
.probe_finite(obj)
```

## Details

A probe is only ever an accelerator: a `FALSE` sends the caller to the
next rung of
[.sigma_frac_ladder](https://janolefi.github.io/gamRTMB/reference/dot-sigma_frac_ladder.md),
and if every rung fails the fit proceeds from the original starting
values under the full cap, exactly as it would have. So a probe that is
wrong about a slow-but-sound inner solve costs a few seconds, never a
fit.
