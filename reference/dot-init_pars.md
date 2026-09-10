# Starting parameter values

Intercepts start on the link scale from the family's
[`start()`](https://rdrr.io/r/stats/start.html); coefficients start at
zero.

## Usage

``` r
.init_pars(design, family, y, frac = 0.05, start = NULL)
```

## Details

For the variance components, cold-starting every log-sigma at zero is
slow and can wander on flat marginal surfaces. Instead pick \\\sigma_k\\
so that the term's implied prior standard deviation, \\\sigma_k
\sqrt{mean(rowSums(X_r^2))}\\, is `frac` of the rough scale of that
parameter's linear predictor, which the family supplies. Blocks tied by
an `id` get a common value, since only one of them survives the mapping.

`frac` errs on the smooth side deliberately. Too flexible a start lets
the inner Newton solve push a parameter out of the family's support,
which shows up as a non-finite marginal objective; too clamped a start
leaves the smooth pinned to its null space, where the REML gradient in
log-sigma is nearly zero and the outer optimiser stalls. Over six
simulated designs the two failure modes bracket a wide, flat optimum: on
well-behaved families (Gaussian, gamma, beta, t) any `frac` between 0.5
and 0.005 reaches the same optimum to two decimals, while on
four-parameter families it matters, with 0.05 converging on 172/180 fits
against 168/180 at 0.2 and collapsing below 0.02. Hence the default. The
cost of the smaller start is 20-40%% more outer iterations on the models
that never had trouble.
