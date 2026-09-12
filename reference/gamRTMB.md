# Fit a distributional GAM

Smooth terms on every parameter of a distribution, fitted by combining
mgcv's basis and penalty construction with RTMB's automatic
differentiation and Laplace approximation, over the log-densities in
RTMBdist.

## Usage

``` r
gamRTMB(
  formula,
  family = fam("norm"),
  data = NULL,
  weights = NULL,
  na.action = stats::na.omit,
  knots = NULL,
  method = c("REML", "ML"),
  engine = c("laplace", "efs"),
  sigma_frac = 0.05,
  sparse = c("auto", "never", "always"),
  joint_precision = TRUE,
  start = NULL,
  silent = TRUE,
  control = list(),
  inner_control = list()
)

# S3 method for class 'gamRTMB'
nobs(object, ...)

# S3 method for class 'gamRTMB'
coef(object, ...)

# S3 method for class 'gamRTMB'
fitted(object, ...)
```

## Arguments

- formula:

  A two-sided formula whose right-hand side is a
  [`list()`](https://rdrr.io/r/base/list.html) of per-parameter
  formulas, or a plain right-hand side for the family's first parameter.

- family:

  A `gamRTMB_family`, from
  [`fam()`](https://janolefi.github.io/gamRTMB/reference/fam.md). See
  [`families()`](https://janolefi.github.io/gamRTMB/reference/families.md).

- data:

  A data frame holding every model variable, or `NULL` (the default) to
  take them from the environment of `formula`.

- weights:

  Optional prior weights, evaluated in `data`. As in
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html), each observation's
  log-density contribution is multiplied by its weight.

- na.action:

  How to treat missing values in any model variable;
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html) by default,
  which drops those rows and reports how many in the fit's summary.

- knots:

  Passed to
  [`mgcv::smoothCon()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html).

- method:

  `"REML"` (default) or `"ML"`.

- engine:

  Fitting engine; only `"laplace"` is implemented.

- sigma_frac:

  Tuning constant for the variance-component starting values: each
  smooth starts contributing this fraction of its parameter's
  linear-predictor scale. See
  [`.init_pars()`](https://janolefi.github.io/gamRTMB/reference/dot-init_pars.md).
  Raise it if a fit converges to an over-smooth solution. If the
  objective is not finite here, a few other values are tried
  automatically before giving up (see `.sigma_frac_ladder()`) and the
  one used is reported; passing this argument explicitly does not switch
  that off, but passing `start` does.

- sparse:

  How to treat a smooth whose single penalty is already sparse – a
  Markov random field, a random walk, a supplied GMRF precision.
  `"auto"` (default) keeps the penalty and gives the block a
  [`RTMB::dgmrf()`](https://rdrr.io/pkg/RTMB/man/MVgauss.html) prior
  when it has at least 50 coefficients and is at most 20% nonzero, and
  sends everything else through
  [`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
  as usual. `"never"` is the old behaviour; `"always"` takes the sparse
  route for every single-penalty smooth, which is mainly useful for
  checking that the two agree. See
  [`.gmrf_block()`](https://janolefi.github.io/gamRTMB/reference/dot-gmrf_block.md).

- joint_precision:

  Ask
  [`RTMB::sdreport()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)
  for the joint precision matrix, which
  [`predict()`](https://rdrr.io/r/stats/predict.html) needs for standard
  errors. On by default; turn it off to save time and memory on large
  models where bands are not wanted.

- start:

  Optional named list overriding entries of the starting parameter list
  (`beta`, `b`, `log_sigma`).

- silent:

  Passed to
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

- control:

  Passed to [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- inner_control:

  Passed to
  [`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)'s
  `inner.control`, which governs the inner Newton solve over the
  coefficients rather than the outer optimisation of the smoothing
  parameters. `list(maxit = ...)` is the entry worth reaching for; see
  `dev/NOTES-inner-method.md` for why `inner.method` is not exposed.

- object:

  A `gamRTMB` fit.

- ...:

  Ignored.

## Value

An object of class `gamRTMB`.

## Methods (by generic)

- `nobs(gamRTMB)`: Number of observations actually used.

- `coef(gamRTMB)`: Coefficients, as a list of the fixed (`beta`) and
  penalized (`b`) vectors plus the log variance components. `beta` and
  `log_sigma` are named `parameter:term`.

- `fitted(gamRTMB)`: Fitted values of every distributional parameter, on
  the response scale.

## Formula

The response is the left-hand side of the outer formula and each
distributional parameter gets a one-sided formula:
`y ~ list(mean = ~ s(x1) + s(x2), sd = ~ s(x1))`. Parameters the family
declares but the formula omits are given `~1`. Parameter names are the
density's own (`xi`, `omega`, `alpha` for a skew normal), not generic
location/scale/shape labels.

A plain right-hand side is shorthand for modelling the family's *first*
parameter and leaving the rest constant, so `y ~ s(x)` means
`y ~ list(mean = ~ s(x))` for `fam("norm")`. That is the ordinary gam
formula, and it is what makes the one-parameter case read like one.

## Data

`data` is optional. Without it the model variables are looked up where
the formula was written, as in
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) or
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html); they are
collected into a data frame first, so everything downstream —
`na.action`, [`predict()`](https://rdrr.io/r/stats/predict.html), the
plots — sees the same rectangle either way. Variables found this way
must all have the same length.

## Sparse penalties

A smooth whose penalty is a sparse precision matrix – `bs = "mrf"` over
an adjacency graph, a random walk, or any precision supplied through
`xt = list(penalty = )` – can skip
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html).
That rotation makes the coefficients iid, which is convenient but fills
the penalty in completely; keeping it instead and giving the block an
\\N(0, \sigma^2 Q^{-1})\\ prior through
[`RTMB::dgmrf()`](https://rdrr.io/pkg/RTMB/man/MVgauss.html) is the same
model at a fraction of the cost. The `sparse` argument controls this.

An intrinsic field is corner-constrained rather than sum-to-zero
constrained, since the latter is what destroys the sparsity. The two are
equivalent up to a constant absorbed by the intercept, so such a term
needs its parameter to have one. See
[`.null_space()`](https://janolefi.github.io/gamRTMB/reference/dot-null_space.md).

## REML

With `method = "REML"` the mean-structure coefficients join the random
vector alongside the spline coefficients, so the same Laplace
approximation integrates out both. This is the bias and stability
correction of Wood (2011); it is not a sparsity argument, since mgcv's
bases have global support and are dense either way. ML keeps them as
fixed effects, which also means
[`edf()`](https://janolefi.github.io/gamRTMB/reference/edf.md) is
unavailable.

## Engines

`engine = "laplace"` hands the smoothing parameters to `nlminb` and lets
RTMB supply the REML criterion and its gradient. `engine = "efs"` is
reserved for an extended Fellner-Schall fit, which would avoid the
third-derivative term in that gradient at the cost of owning its own
inner optimisation; it is not implemented, and the seams it needs are
documented in `dev/NOTES-fellner-schall.md`.

## Supported smooths

`s()`, `t2()`, `by =` variables, `bs = "fs"` and `bs = "re"` all
reconstruct exactly through
[`.reconstruct_map()`](https://janolefi.github.io/gamRTMB/reference/dot-reconstruct_map.md).
`te()` is not supported, because mgcv itself declines
`smooth2random(type = 2)` for it and directs you to `t2()`. `fx = TRUE`
is rejected, having no penalized part. `s(..., id = )` shares one
smoothing parameter across a group of smooths, including across
distributional parameters.

## Examples

``` r
set.seed(1)
d <- data.frame(x1 = runif(200), x2 = runif(200))
d$y <- rnorm(200, sin(2 * pi * d$x1), exp(-1 + d$x2))
fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)),
               data = d)
fit
#> gamRTMB fit
#>   family:    norm (mean/identity, sd/log)
#>   criterion: REML   engine: laplace
#>   converged: TRUE   -REML: 195.2955   max|grad|: 4.2e-05
#>   observations: 200
#>   coefficients: 4 fixed (incl. null spaces), 12 penalized; 2 smoothing parameters
edf(fit)
#>   parameter  term      edf k        sp id
#> 1      mean s(x1) 5.302745 7     0.121   
#> 2        sd s(x2) 1.000000 7 3.658e+08   

## a plain right-hand side models the family's first parameter
gamRTMB(y ~ s(x1, k = 8), data = d)
#> gamRTMB fit
#>   family:    norm (mean/identity, sd/log)
#>   criterion: REML   engine: laplace
#>   converged: TRUE   -REML: 205.3397   max|grad|: 1.89e-08
#>   observations: 200
#>   coefficients: 3 fixed (incl. null spaces), 6 penalized; 1 smoothing parameter
```
