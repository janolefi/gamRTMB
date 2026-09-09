# Assemble the model data

Collects every variable the model touches — via
[`mgcv::interpret.gam()`](https://rdrr.io/pkg/mgcv/man/interpret.gam.html)'s
`fake.formula`, which reports the variables inside `s()` terms as well
as the parametric and offset ones — applies `na.action` to that set, and
returns the data with incomplete rows dropped, together with the
response and prior weights aligned to it.

## Usage

``` r
.model_data(response, par_formulas, data, weights, na.action)
```

## Arguments

- response:

  The response expression (LHS of the outer formula).

- par_formulas:

  One-sided formulas, one per distributional parameter.

- data:

  A data frame.

- weights:

  Evaluated prior weights, or `NULL`.

- na.action:

  Missing-data action, e.g.
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html).

## Value

`list(data, y, weights, dropped)`.

## Details

Model variables must be columns of `data`. R would otherwise let them
come from the calling environment, where dropping rows for missing
values could silently misalign them against the response.
