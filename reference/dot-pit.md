# Put PIT values on the requested scale, keeping them finite

[`qnorm()`](https://rdrr.io/r/stats/Normal.html) at exactly 0 or 1 is
infinite, which would drop the very observations a diagnostic most wants
to show; clamping keeps them finite and visibly extreme instead.

## Usage

``` r
.pit(u, type)
```
