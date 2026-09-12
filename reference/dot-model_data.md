# Assemble the model data

Collects every variable the model touches — via
[`mgcv::interpret.gam()`](https://rdrr.io/pkg/mgcv/man/interpret.gam.html)'s
`fake.formula`, which reports the variables inside `s()` terms as well
as the parametric and offset ones — applies `na.action` to that set, and
returns the data with incomplete rows dropped, together with the
response and prior weights aligned to it.

## Usage

``` r
.model_data(
  response,
  par_formulas,
  data,
  weights,
  na.action,
  env = parent.frame()
)
```

## Arguments

- response:

  The response expression (LHS of the outer formula).

- par_formulas:

  One-sided formulas, one per distributional parameter.

- data:

  A data frame, or `NULL` to look the variables up in `env`.

- weights:

  Evaluated prior weights, or `NULL`.

- na.action:

  Missing-data action, e.g.
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html).

- env:

  Environment for the variables when `data` is `NULL`.

## Value

`list(data, y, weights, dropped)`.

## Details

With a `data` argument, model variables must be columns of it. R would
otherwise let some of them come from the calling environment, where
dropping rows for missing values could silently misalign them against
the response. Without one, everything comes from `env` and is rectangled
here, before any row is dropped, so the same guarantee holds.
