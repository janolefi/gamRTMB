# Effective degrees of freedom per smooth

TMB does not hand these over the way mgcv's PIRLS does, so they are
derived explicitly. At the fitted smoothing parameters the penalized
Hessian of the joint negative log-likelihood in the coefficients is
\$\$H = H\_{data} + S, \quad S = diag(0 \text{ for fixed}, 1/\sigma_k^2
\text{ for block } k),\$\$ with \\S\\ block diagonal – \\\sigma_k^{-2}
I\\ in the
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
basis, \\\sigma_k^{-2} Q_k\\ for a block that kept its own sparse
penalty. Wood's effective degrees of freedom, \\tr((X'WX + S)^{-1}
X'WX)\\, therefore generalise to \$\$F = H^{-1} H\_{data} = I - H^{-1}
S, \quad edf_j = 1 - \[H^{-1} S\]\_{jj},\$\$ so only the diagonal of
\\H^{-1} S\\ is needed, and that comes from a sparse solve rather than a
full inverse. Null-space coefficients contribute exactly 1 and penalized
ones between 0 and 1.

## Usage

``` r
edf(object, ...)

# S3 method for class 'gamRTMB'
edf(object, ...)
```

## Arguments

- object:

  A `gamRTMB` fit, made with `method = "REML"`.

- ...:

  Ignored.

## Value

A data frame with one row per smooth: parameter, term label, EDF, basis
dimension, smoothing parameter(s) and `id`. The total EDF over all
coefficients is attached as attribute `"edf.total"`.

## Details

Checked against
[`mgcv::gaulss()`](https://rdrr.io/pkg/mgcv/man/gaulss.html), which
agrees to three decimals on both untied and `id`-tied models.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(200)); d$y <- rnorm(200, sin(2 * pi * d$x), 0.3)
edf(gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d))
#>   parameter term      edf k      sp id
#> 1      mean s(x) 6.392431 7 0.09817   
```
