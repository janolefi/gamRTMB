# Starting parameter values

Intercepts start on the link scale from the family's
[`start()`](https://rdrr.io/r/stats/start.html); coefficients start at
zero.

## Usage

``` r
.init_pars(design, family, y, frac = 0.05, start = NULL)
```

## Details

For the variance components, cold-starting every log-sigma at zero is
slow and can wander on flat marginal surfaces. Instead pick \\\sigma_k\\
so that the term's implied prior standard deviation, \\\sigma_k
\sqrt{mean(rowSums(X_r^2))}\\, is `frac` of the rough scale of that
parameter's linear predictor, which the family supplies. Blocks tied by
an `id` get a common value, since only one of them survives the mapping.

`frac` errs on the smooth side deliberately. Too flexible a start lets
the inner Newton solve push a parameter out of the family's support,
which shows up as a non-finite marginal objective; too clamped a start
leaves the smooth pinned to its null space, where the REML gradient in
log-sigma is nearly zero and the outer optimiser stalls.

On well-behaved families the two failure modes bracket a wide, flat
optimum, and `frac` may as well not exist: over `norm`, `gamma2` and
`beta2` every value from 0.2 down to 0.001 converges on all fifteen fits
and agrees to seven significant figures. On four-parameter families it
matters, but **not monotonically, and no value dominates** – 0.01 does
worst, with 0.2 and 0.001 on either side of it doing better, and the
spread between families is far larger than the spread across `frac`. So
0.05 is a default rather than an optimum, kept because the evidence for
moving it is four fits out of twenty-five. The small values additionally
converge to a worse optimum more often, which is the
pinned-to-null-space mode above.

Since no single value serves, a start that leaves the objective
non-finite is retried over
[.sigma_frac_ladder](https://janolefi.github.io/gamRTMB/reference/dot-sigma_frac_ladder.md)
rather than left to the user to guess. See `dev/NOTES-sigma-frac.md` for
the measurements, and `dev/bench-sigma-frac.R` to reproduce them.
