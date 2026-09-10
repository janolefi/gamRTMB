# Connected components of a penalty's sparsity pattern

A breadth-first sweep over the off-diagonal nonzeros, returning a
component label per coefficient. Isolated coefficients – an island
region with no neighbours – come back as components of their own, which
is the right answer: their level is free too.

## Usage

``` r
.components(S, tol)
```
