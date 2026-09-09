# Starting parameter values

Intercepts start on the link scale from the family's
[`start()`](https://rdrr.io/r/stats/start.html); coefficients start at
zero.

## Usage

``` r
.init_pars(design, family, y, frac = 0.2, start = NULL)
```

## Details

For the variance components, cold-starting every log-sigma at zero is
slow and can wander on flat marginal surfaces. Instead pick \\\sigma_k\\
so that the term's implied prior standard deviation, \\\sigma_k
\sqrt{mean(rowSums(X_r^2))}\\, is `frac` of the rough scale of that
parameter's linear predictor, which the family supplies. Blocks tied by
an `id` get a common value, since only one of them survives the mapping.
