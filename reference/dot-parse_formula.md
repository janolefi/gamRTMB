# Parse a distributional model formula

The response is taken from the left-hand side of the outer formula and
each distributional parameter gets a one-sided formula from the
[`list()`](https://rdrr.io/r/base/list.html) on the right, in the style
of `brms::bf()`. Any parameter the family declares but the user does not
mention is given `~1`, which is why the family spec has to exist before
formula processing.

## Usage

``` r
.parse_formula(formula, parnames)
```

## Arguments

- formula:

  e.g. `y ~ list(mean = ~ s(x1) + s(x2), sd = ~ s(x1))`.

- parnames:

  Character vector of the family's modelled parameters.

## Value

A list with the response expression and one formula per parameter,
ordered as `parnames`.
