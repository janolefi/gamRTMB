# Fallback variance-component starts

Tried in this order when the first marginal evaluation is not finite.
They are not a refinement of one another, and the ladder is not a search
for a better default: over eight families and five seeds the number of
fits reaching a converged solution is not monotone in `frac`, no value
dominates, and the spread between families is far larger than the spread
across `frac`. What the measurements support is that no single value
serves, so the point of the ladder is coverage.

## Usage

``` r
.sigma_frac_ladder
```

## Format

A numeric vector of `sigma_frac` values.

## Details

Ordered by that table: the two that did best on count, then 0.001 last,
because it converges most often and misses the optimum most often – a
start that clamped leaves the smooth pinned near its null space.

See
[`.init_pars()`](https://janolefi.github.io/gamRTMB/reference/dot-init_pars.md)
for what `frac` means, `dev/NOTES-sigma-frac.md` for the measurements
and `dev/bench-sigma-frac.R` to reproduce them.
