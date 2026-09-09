# Colours for a fan of quantile curves

One hue throughout, with the median in black: a fan of quantiles is one
ordered family, not a set of unrelated series, so different colours and
line types would imply distinctions that are not there. Opacity carries
the ordering instead, fading outwards from the median, which keeps a
dense fan legible and reads as the density it approximates.

## Usage

``` r
.fan_cols(prob, hue = "#0B5D9E")
```

## Arguments

- prob:

  Probabilities, in the order they will be drawn.

## Value

A character vector of colours.
