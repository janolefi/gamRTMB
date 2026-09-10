# Matern SPDE smooth

A smooth term for a spatially continuous Gaussian random field, usable
as `s(x, y, bs = "spde", xt = list(mesh = mesh))` in any `gamRTMB`
formula and therefore on any parameter of any family.

## Usage

``` r
# S3 method for class 'spde.smooth.spec'
smooth.construct(object, data, knots)

# S3 method for class 'spde.smooth'
Predict.matrix(object, data)
```

## Arguments

- object, data, knots:

  As
  [`mgcv::smooth.construct()`](https://rdrr.io/pkg/mgcv/man/smooth.construct.html).

- object, data:

  As
  [`mgcv::Predict.matrix()`](https://rdrr.io/pkg/mgcv/man/Predict.matrix.html).

## Value

A `smoothCon` object of class `spde.smooth`, with sparse `X` and sparse
penalties.

## Building a mesh

The mesh is yours to build and to pass in, because a good one depends on
the domain, the data and the range you expect, and no default can know
those.
[`fmesher::fm_mesh_2d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)
takes locations, a maximum triangle edge length inside the domain and in
the outer extension, and a cutoff below which nearby points are merged:

    mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.1, 0.3),
                                cutoff = 0.05, offset = c(0.1, 0.3))
    f <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh))),
                 data = d)

The outer extension matters: without it the field's variance is inflated
at the boundary. A rule of thumb is `max.edge` no larger than a third of
the range you expect, and an offset about equal to the range.

For one dimension pass a
[`fmesher::fm_mesh_1d()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_1d.html),
or leave `xt` out and a mesh with `k` evenly spaced knots is built over
the range of the covariate.

## Parameterisation

The precision is \$\$Q(\tau, \kappa) = \tau^2(\kappa^4 C + 2\kappa^2
G_1 + G_2)\$\$ with \\C\\, \\G_1\\, \\G_2\\ the finite element matrices
from
[`fmesher::fm_fem()`](https://inlabru-org.github.io/fmesher/reference/fm_fem.html).
For `alpha = 2` in two dimensions the Matern smoothness is \\\nu = 1\\,
the range is \\\rho = \sqrt{8\nu}/\kappa = 2\sqrt2/\kappa\\ and the
marginal standard deviation is \\\sigma = 1/(\tau\kappa\sqrt{4\pi})\\.

The two estimated parameters are \\\log\sigma\\ and \\\log\rho\\, not
\\\log\tau\\ and \\\log\kappa\\, and they are reported as `sd` and
`range` by
[`edf()`](https://janolefi.github.io/gamRTMB/reference/edf.md) and named
in [`coef()`](https://rdrr.io/r/stats/coef.html). The change is exact –
it is a linear map of the log parameters, absorbed into the constants in
front of the three matrices – and it is worth making for two reasons.
The parameters mean something on their own, and the SPDE likelihood has
a long ridge along constant marginal variance which in \\(\tau,
\kappa)\\ runs diagonally and in \\(\sigma, \rho)\\ runs along an axis,
which the outer optimiser finds much easier.

The three matrices are carried through mgcv's `L` convention – `L` holds
the powers of each parameter in each penalty, so that \\\lambda =
\exp(L\theta)\\ – which is what tells `gamRTMB` to give the block two
parameters and a
[`RTMB::dgmrf()`](https://rdrr.io/pkg/RTMB/man/MVgauss.html) prior
instead of one variance.
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
cannot represent such a term at all, so an SPDE smooth always takes the
sparse route regardless of the `sparse` argument.

Both parameters are estimated freely, with no prior on them. A single
realisation of a spatial field identifies them only loosely, and a field
whose range approaches the size of the domain is close to improper: the
range then runs off and the marginal standard deviation follows it, with
the fitted surface barely changing. If that happens, the fit is still
usable but the two parameters are not, and a shorter mesh `max.edge`
will not help – it is the data that are uninformative. INLA's PC priors
exist for this.

The field is proper for any positive \\\kappa\\, so no identifiability
constraint is imposed, as in INLA. It is not centred, and its overall
level is only weakly separated from the intercept when the range is
large.

## References

Lindgren, F., Rue, H. and Lindstrom, J. (2011). An explicit link between
Gaussian fields and Gaussian Markov random fields. *JRSS-B* 73, 423-498.

## Examples

``` r
# \donttest{
if (requireNamespace("fmesher", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(x = runif(300), y = runif(300))
  d$z <- rnorm(300, sin(6 * d$x) * cos(6 * d$y), 0.3)
  mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.15, 0.4),
                              cutoff = 0.06, offset = c(0.1, 0.3))
  fit <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde",
                                     xt = list(mesh = mesh)), sd = ~ 1),
                 data = d)
  edf(fit)          # the sp column reports the field's sd and range
}
#>   parameter   term      edf   k                    sp id
#> 1      mean s(x,y) 50.61775 239 sd=0.8243,range=1.226   
# }
```
