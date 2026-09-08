#' Fit a distributional GAM
#'
#' Smooth terms on every parameter of a distribution, fitted by combining
#' \pkg{mgcv}'s basis and penalty construction with \pkg{RTMB}'s automatic
#' differentiation and Laplace approximation, over the log-densities in
#' \pkg{RTMBdist}.
#'
#' @section Formula:
#' The response is the left-hand side of the outer formula and each
#' distributional parameter gets a one-sided formula:
#' `y ~ list(mu = ~ s(x1) + s(x2), sigma = ~ s(x1))`. Parameters the family
#' declares but the formula omits are given `~1`. Parameter names are the
#' density's own (`xi`, `omega`, `alpha` for a skew normal), not generic
#' location/scale/shape labels.
#'
#' @section REML:
#' With `method = "REML"` the mean-structure coefficients join the random
#' vector alongside the spline coefficients, so the same Laplace
#' approximation integrates out both. This is the bias and stability
#' correction of Wood (2011); it is not a sparsity argument, since mgcv's
#' bases have global support and are dense either way. ML keeps them as fixed
#' effects, which also means [edf()] is unavailable.
#'
#' @section Engines:
#' `engine = "laplace"` hands the smoothing parameters to `nlminb` and lets
#' RTMB supply the REML criterion and its gradient. `engine = "efs"` is
#' reserved for an extended Fellner-Schall fit, which would avoid the
#' third-derivative term in that gradient at the cost of owning its own inner
#' optimisation; it is not implemented, and the seams it needs are documented
#' in `dev/NOTES-fellner-schall.md`.
#'
#' @section Supported smooths:
#' `s()`, `t2()`, `by =` variables, `bs = "fs"` and `bs = "re"` all
#' reconstruct exactly through [.reconstruct_map()]. `te()` is not supported,
#' because mgcv itself declines `smooth2random(type = 2)` for it and directs
#' you to `t2()`. `fx = TRUE` is rejected, having no penalized part.
#' `s(..., id = )` shares one smoothing parameter across a group of smooths,
#' including across distributional parameters.
#'
#' @param formula A two-sided formula whose right-hand side is a `list()` of
#'   per-parameter formulas.
#' @param family A `gamRTMB_family`, from [rtmbdist_family()] or
#'   [gaussian_ls()].
#' @param data A data frame.
#' @param knots Passed to [mgcv::smoothCon()].
#' @param method `"REML"` (default) or `"ML"`.
#' @param engine Fitting engine; only `"laplace"` is implemented.
#' @param sigma_frac Tuning constant for the variance-component starting
#'   values; see [.init_log_sigma()].
#' @param joint_precision Ask [RTMB::sdreport()] for the joint precision
#'   matrix, needed for standard errors on smooth terms. Costs more on larger
#'   models.
#' @param start Optional named list overriding entries of the starting
#'   parameter list (`beta`, `b`, `log_sigma`).
#' @param silent Passed to [RTMB::MakeADFun()].
#' @param control Passed to [stats::nlminb()].
#' @return An object of class `gamRTMB`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = runif(200), x2 = runif(200))
#' d$y <- rnorm(200, sin(2 * pi * d$x1), exp(-1 + d$x2))
#' fit <- gamRTMB(y ~ list(mu = ~ s(x1, k = 8), sigma = ~ s(x2, k = 8)),
#'                family = gaussian_ls(), data = d)
#' fit
#' edf(fit)
#' @export
gamRTMB <- function(formula, family = gaussian_ls(), data,
                    knots = NULL, method = c("REML", "ML"),
                    engine = c("laplace", "efs"), sigma_frac = 0.2,
                    joint_precision = FALSE, start = NULL, silent = TRUE,
                    control = list()) {
  method <- match.arg(method)
  engine <- match.arg(engine)
  if (!inherits(family, "gamRTMB_family"))
    stop("`family` must be a gamRTMB_family, e.g. gaussian_ls() or ",
         "rtmbdist_family(\"gamma2\")")
  if (missing(data) || !is.data.frame(data))
    stop("`data` must be a data frame")

  pf <- .parse_formula(formula, family$parnames)
  y <- eval(pf$response, data, environment(formula))
  design <- .build_design(pf$par_formulas, data, family$parnames, knots = knots)
  fx <- .resolve_fixed(family, data, length(y))
  pars <- .init_pars(design, family, y, sigma_frac, start)
  nll <- .make_nll(design, family, y, fx)

  fit <- switch(engine,
    laplace = .fit_laplace(nll, pars, design, method, joint_precision, silent,
                           control, family),
    efs     = .fit_efs(nll, pars, design, method, family))

  fit$family <- family
  fit$method <- method
  fit$engine <- engine
  fit$design <- design
  fit$formula <- formula
  fit$par_formulas <- pf$par_formulas
  fit$y <- y
  fit$fixed <- fx
  structure(fit, class = "gamRTMB")
}

#' Laplace engine
#'
#' Declares the coefficients random so that RTMB supplies the marginal
#' criterion, then optimises the variance components with `nlminb`.
#'
#' Smoothing parameters tied by an `id` are collapsed through
#' [RTMB::MakeADFun()]'s `map`: several entries of `log_sigma` become one
#' estimated value. An `NA` level there would instead fix an entry at its
#' starting value, which is how a user-specified smoothing parameter would be
#' implemented.
#'
#' @keywords internal
.fit_laplace <- function(nll, pars, design, method, joint_precision, silent,
                          control, family) {
  random <- if (method == "REML") c("beta", "b") else "b"
  map <- if (design$nsigma_free < design$nsigma)
    list(log_sigma = design$sig_group) else list()

  obj <- RTMB::MakeADFun(nll, pars, random = random, map = map, silent = silent)

  ## start diagnostics, on the joint objective rather than the marginal one
  pfull <- obj$env$par
  g0 <- tryCatch(obj$env$f(pfull, order = 1), error = function(e) NULL)
  msg <- if (!is.null(g0))
    .flat_start(g0, function(p) obj$env$f(p, order = 0), pfull,
                names(pfull) == "beta", design, design$parnames) else NULL

  ## A flat direction can be worse than leaving a parameter stuck: with no
  ## curvature in that coefficient block either, the inner Newton solve is
  ## singular and the Laplace approximation is undefined, so the marginal
  ## objective comes back non-finite. Name the parameter rather than blaming
  ## the response's support.
  if (!is.finite(obj$fn(obj$par)))
    stop("the objective is not finite at the starting values. ",
         if (!is.null(msg)) msg else
           paste0("Check that the response is in the support of family '",
                  family$family, "', and consider passing start = list(beta = ...)."))
  if (!is.null(msg)) warning(msg, call. = FALSE)

  ctl <- utils::modifyList(list(eval.max = 2000, iter.max = 1000), control)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = ctl)
  sdr <- RTMB::sdreport(obj, getJointPrecision = joint_precision)

  pl <- obj$env$parList(par = obj$env$last.par.best)
  list(obj = obj, opt = opt, sdr = sdr,
       coefficients = list(beta = pl$beta, b = pl$b),
       log_sigma = pl$log_sigma,               # full length, not the mapped one
       objective = opt$objective,
       convergence = opt$convergence == 0,
       max_grad = max(abs(obj$gr(opt$par))))
}

#' Extended Fellner-Schall engine (not implemented)
#'
#' @keywords internal
.fit_efs <- function(nll, pars, design, method, family) {
  stop("engine = \"efs\" is not implemented yet.\n",
       "It is not an optimiser setting but a separate fit: nothing may be ",
       "declared random, so the Laplace machinery is never built, and the ",
       "engine owns both an inner Newton loop over the coefficients and the ",
       "multiplicative smoothing-parameter update. In this parameterisation ",
       "that update reduces to lambda_k <- edf_k / sum(b_k^2). ",
       "See dev/NOTES-fellner-schall.md.", call. = FALSE)
}
