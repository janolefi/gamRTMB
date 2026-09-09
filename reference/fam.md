# Build a family object

Reads a density's [`formals()`](https://rdrr.io/r/base/formals.html) and
derives everything needed to model it: the native parameter names in the
density's own order, a link per parameter, which arguments are data
rather than parameters, and starting values. Most of RTMBdist works with
no hand-written family; use
[`families()`](https://janolefi.github.io/gamRTMB/reference/families.md)
to see what is available.

## Usage

``` r
fam(
  dist,
  links = NULL,
  fixed = NULL,
  support = NULL,
  start = NULL,
  eta_scale = NULL
)

# S3 method for class 'gamRTMB_family'
print(x, ...)
```

## Arguments

- dist:

  Density name, with or without the leading `d` (`"skewnorm2"` or
  `"dskewnorm2"`).

- links:

  Named character vector overriding the derived links.

- fixed:

  Named list of values for fixed arguments; each is a constant or the
  name of a column of the data.

- support:

  Response support: `"continuous"`, `"lattice"` (integer) or `"mixed"`
  (continuous with atoms). Derived from
  [.lattice_families](https://janolefi.github.io/gamRTMB/reference/dot-lattice_families.md)
  and the presence of an inflation parameter; override it for a family
  those rules get wrong. Only
  [`residuals.gamRTMB()`](https://janolefi.github.io/gamRTMB/reference/residuals.gamRTMB.md)
  uses it.

- start, eta_scale:

  Optional replacements for the starting-value and
  linear-predictor-scale heuristics.

- x:

  A `gamRTMB_family`.

- ...:

  Ignored.

## Value

An object of class `gamRTMB_family`.

## Where densities come from

RTMBdist is searched first, then a curated list of the standard
densities that RTMB makes AD-aware (`norm`, `pois`, `binom`, `gamma`,
`exp`, `lnorm`, `weibull`, `cauchy`, `logis`, `t`, `chisq`). Parameter
names are always the density's own: a skew normal is `xi`, `omega`,
`alpha`, and a Gaussian is `mean`, `sd`.

## Modelled versus fixed arguments

Modelled parameters get a formula, a link and smooths. Fixed arguments
are known data or constants — the number of binomial trials, truncation
bounds, a numerical `eps` — and are passed straight to the density,
unreachable from the formula interface. Supply one with
`fixed = list(size = "trials")` naming a column of the data, or a
constant.

## Starting values

Location parameters get data-driven starts and the rest keep the
density's own default, with one exception. A shape parameter entering an
already mean/sd standardised density can have an **identically zero
score** at the symmetric point: in `dskewnorm2` the direct effect of
`alpha` on the log density cancels exactly against the shift in the
internal location needed to hold `mean` and `sd` fixed, so the
derivative is zero for every observation at `alpha = 0` (measured at
~1e-15, i.e. exactly zero). The default `alpha = 0` is therefore a
perfectly neutral value and a useless starting point: the optimiser has
no descent direction, and because that coefficient block has no
curvature either, the inner Newton solve is singular and the Laplace
approximation is undefined. Such parameters start off the symmetric
point instead, with the sign taken from the sample skewness — starting
at `-0.5` on right-skewed data gets stuck just as badly as starting at
`0`.

## See also

[`families()`](https://janolefi.github.io/gamRTMB/reference/families.md)
for what is available,
[`gamRTMB()`](https://janolefi.github.io/gamRTMB/reference/gamRTMB.md)
to fit.

## Examples

``` r
fam("norm")
#> gamRTMB family: norm  [dnorm, RTMB]
#>   modelled: mean (identity), sd (log)
#>   support:  continuous
fam("gamma2")
#> gamRTMB family: gamma2  [dgamma2, RTMBdist]
#>   modelled: mean (log), sd (log)
#>   support:  continuous
fam("skewnorm2")
#> gamRTMB family: skewnorm2  [dskewnorm2, RTMBdist]
#>   modelled: mean (identity), sd (log), alpha (identity)
#>   support:  continuous
fam("betabinom", fixed = list(size = "trials"))
#> gamRTMB family: betabinom  [dbetabinom, RTMBdist]
#>   modelled: shape1 (log), shape2 (log)
#>   fixed:    size
#>   support:  lattice
#>   no residuals or quantiles: the density has no CDF or quantile function
```
