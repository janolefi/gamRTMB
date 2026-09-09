# What a term plot needs, recorded at design time

A one-dimensional smooth of a numeric covariate can be drawn on its own,
so record which covariate it uses and, for a `by=` smooth, the value of
the `by` variable that
[`mgcv::PredictMat()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html) will
want: the smooth's own factor level, or 1 for a numeric `by`. The
covariate values themselves come from the fit's stored model frame.

## Usage

``` r
.plot_spec(sm, data)
```

## Details

Returns `NULL` for anything not drawable as a single curve — tensor
products, random effects, factor-smooth interactions — which the plot
method reports rather than drawing wrongly.
