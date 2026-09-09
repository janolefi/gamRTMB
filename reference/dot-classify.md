# Derive a family's structure from its density

Everything that can be read off
[`formals()`](https://rdrr.io/r/base/formals.html) plus the correction
tables, and nothing that needs the data. Shared by
[`fam()`](https://janolefi.github.io/gamRTMB/reference/fam.md) and
[`families()`](https://janolefi.github.io/gamRTMB/reference/families.md),
so that listing families does not mean catching errors from constructing
them.

## Usage

``` r
.classify(dist, fixed = NULL)
```

## Value

`NULL` if no such density; otherwise a list with the modelled
parameters, their links, dropped (derived) arguments, the fixed
arguments and which of them the user must still supply, and default
starts.
