# Locate the CDF or quantile function matching a density

`p<dist>` or `q<dist>` in the density's own namespace, else in stats.
Neither has to be AD-compatible: residuals and quantiles are computed
after fitting, on plain numerics.

## Usage

``` r
.find_pfun(prefix, dist, source, modelled)
```

## Arguments

- prefix:

  `"p"` for the CDF, `"q"` for the quantile function.

- dist, source:

  Density name and the package it came from.

- modelled:

  The modelled parameter names.

## Value

A function of `(x, theta, fx)`, or `NULL`.

## Details

A few RTMBdist versions take extra arguments the density does not
(`ncp`, `method`, `from`, `tol`), which is harmless. The case that
matters is the reverse: if a modelled parameter is missing from the
arguments the function cannot be evaluated faithfully, so nothing is
offered rather than something with a parameter silently dropped.
