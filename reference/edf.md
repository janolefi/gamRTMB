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

## When the EDF do not exist

All of that assumes \\H\_{data}\\ is positive semi-definite, which is
what makes \\F = H^{-1} H\_{data}\\ a projection and puts every
\\edf_j\\ in \\\[0, 1\]\\. At a point the optimiser never converged to
it need not be, and the solve still returns numbers: a four-parameter
Box-Cox fit on `film90` gives EDF near \\-6000\\ for a rank-9 basis.

Note that it is \\H\_{data}\\ and not \\H\\ that has to be checked. The
penalty can and does rescue the sum: on that same fit \\H\\ is positive
definite while \\H\_{data} = H - S\\ has two negative eigenvalues, so a
test on \\H\\ passes and the EDF are still nonsense. Rather than
factorise a second matrix, the \\edf_j\\ are checked against the \\\[0,
1\]\\ they are guaranteed to lie in – the same statement, and already
computed.

A negative EDF is not a small inaccuracy to report with a caveat; it
means the quantity does not exist at this point. So the column is `NA`
instead, with a warning pointing at `max_grad`. See
[`.inner_indefinite()`](https://janolefi.github.io/gamRTMB/reference/dot-inner_indefinite.md)
for how a fit gets into that state.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(200)); d$y <- rnorm(200, sin(2 * pi * d$x), 0.3)
edf(gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d))
#>   parameter term      edf k      sp id
#> 1      mean s(x) 6.392431 7 0.09817   
```
