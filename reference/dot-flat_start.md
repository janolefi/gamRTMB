# Detect a distributional parameter that starts at a useless stationary point

A parameter whose intercept has an identically zero score cannot move,
and a smooth on it then drifts on a flat surface instead of failing
loudly. `dskewnorm2`'s `alpha` does this at `alpha = 0`; see
[`fam()`](https://janolefi.github.io/gamRTMB/reference/fam.md).

## Usage

``` r
.flat_start(grad, eval_at, pfull, is_beta, design)
```

## Arguments

- grad:

  Joint gradient at the starting values.

- eval_at:

  Function of a full parameter vector returning the objective.

- pfull:

  The full starting parameter vector.

- is_beta:

  Logical index of the `beta` entries within `pfull`.

- design:

  The design object.

## Value

A message describing the affected parameters, or `NULL`.

## Details

A zero score is not on its own a problem: an intercept started at its
own marginal optimum has one too, and that is exactly where it should be
(`dnbinom2`'s `mu` started at `log(mean(y))` has a zero score and fits
perfectly). The two are told apart by probing rather than by curvature —
perturb the intercept either way and see whether the objective actually
falls. A stationary point that can be improved on by stepping away from
it is the bad kind.
