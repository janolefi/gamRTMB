# Rotated conditional densities at a few covariate values

The fitted density of the response, turned on its side and drawn at
chosen positions along a covariate — the shape the quantile fan only
summarises. It uses the family's own log-density, so it works for every
family, not only those with a quantile function.

## Usage

``` r
.plot_density(x, xvar, at, ngrid, band.col, ...)
```

## Details

Each density opens to the left of its position line and is scaled to a
common width rather than a common height. A common height would be more
faithful — a concentrated distribution really does have a taller density
— but on data where the spread changes by a factor of 40 it makes the
wide ones invisible, and the shape is the point.
