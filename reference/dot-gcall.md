# Call a graphics function with our defaults under the caller's `...`

Passing both a default and the caller's value for the same argument to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) is an error
("matched by multiple actual arguments"), so they are merged first and
the caller wins.

## Usage

``` r
.gcall(fun, defaults, dots)
```
