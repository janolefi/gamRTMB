## The objective. Deliberately a plain closure over the design: it says
## nothing about how it will be optimised, and in particular nothing about
## which coefficients are treated as random effects. That is the engine's
## choice, which is what keeps a non-Laplace engine (see
## dev/NOTES-fellner-schall.md) a matter of adding a fitting function rather
## than rewriting the model.

#' Build the joint negative log-likelihood
#'
#' The penalized coefficients carry an iid \eqn{N(0, \sigma_k^2)} prior, one
#' variance per penalized block. Null-space coefficients live in `beta` and
#' are never given a prior. Nothing is shared across smooths or across
#' distributional parameters unless an `id` says so.
#'
#' Prior weights multiply each observation's log-density contribution, as in
#' [stats::glm()]. Offsets are added to the relevant parameter's linear
#' predictor.
#'
#' @param design From [.build_design()].
#' @param family A `gamRTMB_family`.
#' @param y Response.
#' @param fx Resolved fixed arguments.
#' @param w Prior weights, or `NULL` for unweighted.
#' @return A function of a parameter list `list(beta, b, log_sigma)`.
#' @keywords internal
.make_nll <- function(design, family, y, fx = list(), w = NULL) {
  parnames  <- design$parnames
  Xfix      <- design$Xfix
  beta_idx  <- design$beta_idx
  blocks    <- design$blocks
  par_blocks <- design$par_blocks
  offsets   <- lapply(design$parts, `[[`, "offset")
  linkinv   <- lapply(family$links, function(l) .links[[l]]$linkinv)
  names(linkinv) <- parnames
  logdens   <- family$logdens

  function(pv) {
    RTMB::getAll(pv)
    yo <- RTMB::OBS(y)
    jnll <- 0

    for (k in seq_along(blocks))
      jnll <- jnll - sum(dnorm(b[blocks[[k]]$idx], 0, exp(log_sigma[k]),
                               log = TRUE))

    theta <- list()
    for (p in parnames) {
      eta <- as.vector(Xfix[[p]] %*% beta[beta_idx[[p]]]) + offsets[[p]]
      for (k in par_blocks[[p]])
        eta <- eta + as.vector(blocks[[k]]$X %*% b[blocks[[k]]$idx])
      theta[[p]] <- linkinv[[p]](eta)
    }
    ld <- logdens(yo, theta, fx)
    jnll - if (is.null(w)) sum(ld) else sum(w * ld)
  }
}

#' Starting parameter values
#'
#' Intercepts start on the link scale from the family's `start()`; coefficients
#' start at zero.
#'
#' For the variance components, cold-starting every log-sigma at zero is slow
#' and can wander on flat marginal surfaces. Instead pick \eqn{\sigma_k} so
#' that the term's implied prior standard deviation,
#' \eqn{\sigma_k \sqrt{mean(rowSums(X_r^2))}}, is `frac` of the rough scale
#' of that parameter's linear predictor, which the family supplies. Blocks
#' tied by an `id` get a common value, since only one of them survives the
#' mapping.
#'
#' @keywords internal
.init_pars <- function(design, family, y, frac = 0.2, start = NULL) {
  beta0 <- numeric(design$nbeta)
  s0 <- family$start(y)
  for (p in design$parnames) {
    ii <- design$beta_idx[[p]]
    k <- match("(Intercept)", attr(ii, "labels"))
    if (!is.na(k)) beta0[ii[k]] <- s0[[p]]
  }
  sc <- family$eta_scale(y)
  ls0 <- vapply(design$blocks, function(bl) {
    rs <- sqrt(mean(rowSums(bl$X^2)))
    log(max(frac * sc[[bl$par]] / max(rs, 1e-8), 1e-4))
  }, numeric(1))
  pars <- list(beta = beta0, b = numeric(design$nb),
               log_sigma = stats::ave(ls0, design$sig_group))
  if (!is.null(start)) pars[names(start)] <- start
  pars
}

#' Detect a distributional parameter that starts at a useless stationary point
#'
#' A parameter whose intercept has an identically zero score cannot move, and
#' a smooth on it then drifts on a flat surface instead of failing loudly.
#' `dskewnorm2`'s `alpha` does this at `alpha = 0`; see [fam()].
#'
#' A zero score is not on its own a problem: an intercept started at its own
#' marginal optimum has one too, and that is exactly where it should be
#' (`dnbinom2`'s `mu` started at `log(mean(y))` has a zero score and fits
#' perfectly). The two are told apart by probing rather than by curvature —
#' perturb the intercept either way and see whether the objective actually
#' falls. A stationary point that can be improved on by stepping away from it
#' is the bad kind.
#'
#' @param grad Joint gradient at the starting values.
#' @param eval_at Function of a full parameter vector returning the objective.
#' @param pfull The full starting parameter vector.
#' @param is_beta Logical index of the `beta` entries within `pfull`.
#' @param design The design object.
#' @return A message describing the affected parameters, or `NULL`.
#' @keywords internal
.flat_start <- function(grad, eval_at, pfull, is_beta, design) {
  if (!any(is.finite(grad))) return(NULL)
  tol <- 1e-8 * max(abs(grad[is.finite(grad)]), 1)
  f0 <- eval_at(pfull)
  gb <- grad[is_beta]
  flat <- character(0)
  for (p in design$parnames) {
    ii <- design$beta_idx[[p]]
    k <- match("(Intercept)", attr(ii, "labels"))
    if (is.na(k) || !is.finite(gb[ii[k]]) || abs(gb[ii[k]]) >= tol) next
    j <- which(is_beta)[ii[k]]
    drop <- vapply(c(-0.25, 0.25), function(h) {
      pp <- pfull; pp[j] <- pp[j] + h
      f0 - tryCatch(eval_at(pp), error = function(e) Inf)
    }, numeric(1))
    if (any(is.finite(drop) & drop > 1e-6 * max(abs(f0), 1))) flat <- c(flat, p)
  }
  if (!length(flat)) return(NULL)
  paste0("the score for '", paste(flat, collapse = "', '"), "' is numerically ",
         "zero at the starting values, at a point that is not optimal, so ",
         if (length(flat) > 1) "these parameters cannot" else "this parameter cannot",
         " move. Pass a different intercept via start = list(beta = ...), or a ",
         "start() in the family spec.")
}

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
#' `y ~ list(mean = ~ s(x1) + s(x2), sd = ~ s(x1))`. Parameters the family
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
#' @param family A `gamRTMB_family`, from [fam()]. See [families()].
#' @param data A data frame. Every model variable must be a column of it.
#' @param weights Optional prior weights, evaluated in `data`. As in
#'   [stats::glm()], each observation's log-density contribution is multiplied
#'   by its weight.
#' @param na.action How to treat missing values in any model variable;
#'   [stats::na.omit()] by default, which drops those rows and reports how
#'   many in the fit's summary.
#' @param knots Passed to [mgcv::smoothCon()].
#' @param method `"REML"` (default) or `"ML"`.
#' @param engine Fitting engine; only `"laplace"` is implemented.
#' @param sigma_frac Tuning constant for the variance-component starting
#'   values; see [.init_pars()].
#' @param joint_precision Ask [RTMB::sdreport()] for the joint precision
#'   matrix, which [predict()] needs for standard errors. On by default; turn
#'   it off to save time and memory on large models where bands are not
#'   wanted.
#' @param start Optional named list overriding entries of the starting
#'   parameter list (`beta`, `b`, `log_sigma`).
#' @param silent Passed to [RTMB::MakeADFun()].
#' @param control Passed to [stats::nlminb()].
#' @return An object of class `gamRTMB`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = runif(200), x2 = runif(200))
#' d$y <- rnorm(200, sin(2 * pi * d$x1), exp(-1 + d$x2))
#' fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)),
#'                data = d)
#' fit
#' edf(fit)
#' @export
gamRTMB <- function(formula, family = fam("norm"), data, weights = NULL,
                    na.action = stats::na.omit, knots = NULL,
                    method = c("REML", "ML"),
                    engine = c("laplace", "efs"), sigma_frac = 0.2,
                    joint_precision = TRUE, start = NULL, silent = TRUE,
                    control = list()) {
  method <- match.arg(method)
  engine <- match.arg(engine)
  if (!inherits(family, "gamRTMB_family"))
    stop("`family` must be a gamRTMB_family, e.g. fam(\"norm\") or ",
         "fam(\"gamma2\"); see families()")
  if (missing(data) || !is.data.frame(data))
    stop("`data` must be a data frame")

  w <- eval(substitute(weights), data, parent.frame())
  pf <- .parse_formula(formula, family$parnames)
  md <- .model_data(pf$response, pf$par_formulas, data, w, na.action)
  data <- md$data; y <- md$y; w <- md$weights
  design <- .build_design(pf$par_formulas, data, family$parnames, knots = knots)
  fx <- .resolve_fixed(family, data, length(y))
  pars <- .init_pars(design, family, y, sigma_frac, start)
  nll <- .make_nll(design, family, y, fx, w)

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
  fit$data <- data
  fit$weights <- w
  fit$dropped <- md$dropped
  fit$na.action <- na.action
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
                names(pfull) == "beta", design) else NULL

  ## A finite objective with a non-finite gradient is not a data problem: it
  ## means the density's derivative is broken at these parameter values, which
  ## no starting value or optimiser setting can rescue. Say that, rather than
  ## letting nlminb fail obscurely. (RTMBdist's dtruncnorm does this for an
  ## infinite bound: the value is right and d/d(sd) is NaN.)
  if (!is.null(g0) && is.finite(obj$env$f(pfull, order = 0)) && anyNA(g0))
    stop("the objective is finite at the starting values but its gradient is ",
         "not, so family '", family$family, "' cannot be differentiated here. ",
         "This is a property of the density rather than of the data: check ",
         "for infinite fixed arguments (bounds), and pass large finite values ",
         "instead if so.", call. = FALSE)

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

  ## With no smooths at all there is nothing for the outer optimiser to do:
  ## every coefficient is already handled by the inner Laplace problem.
  ctl <- utils::modifyList(list(eval.max = 2000, iter.max = 1000), control)
  ## nlminb warns whenever its line search probes a point where the objective
  ## is not finite, which is routine and self-correcting: it backtracks and
  ## carries on. Whether the fit actually worked is reported by `convergence`
  ## and `max_grad`, so that intermediate complaint is muted rather than left
  ## looking like a failure.
  ##
  ## A non-finite *gradient* is different: nlminb raises an error and stops.
  ## That happens when a smoothing parameter reaches a region where the
  ## family's own derivatives break down, and it can happen after real
  ## progress has been made. TMB has kept the best point it saw, so the fit is
  ## returned from there and flagged as unconverged, rather than thrown away.
  opt <- if (!length(obj$par))
    list(par = obj$par, objective = obj$fn(obj$par), convergence = 0L,
         message = "no smoothing parameters to estimate")
  else tryCatch(
    withCallingHandlers(
      stats::nlminb(obj$par, obj$fn, obj$gr, control = ctl),
      warning = function(w) {
        if (grepl("NA/NaN function evaluation", conditionMessage(w)))
          invokeRestart("muffleWarning")
      }),
    error = function(e) {
      warning("the outer optimiser stopped early: ", conditionMessage(e),
              ". The best point reached is returned, but the fit has not ",
              "converged -- check `fit$convergence`, and treat the smoothing ",
              "parameters and any standard errors with suspicion. This ",
              "usually means a smoothing parameter ran into a region where ",
              "family '", family$family, "' cannot be differentiated.",
              call. = FALSE)
      best <- obj$env$last.par.best
      pf <- tryCatch(best[obj$env$lfixed()], error = function(e2) obj$par)
      val <- obj$env$value.best
      list(par = pf, convergence = 1L, message = conditionMessage(e),
           objective = if (length(val) == 1L && is.finite(val)) val else NA_real_)
    })
  sdr <- RTMB::sdreport(obj, getJointPrecision = joint_precision)

  pl <- obj$env$parList(par = obj$env$last.par.best)
  list(obj = obj, opt = opt, sdr = sdr,
       coefficients = list(beta = pl$beta, b = pl$b),
       log_sigma = pl$log_sigma,               # full length, not the mapped one
       objective = opt$objective,
       convergence = opt$convergence == 0,
       max_grad = if (!length(obj$par)) 0 else
         tryCatch(max(abs(obj$gr(opt$par))), error = function(e) NA_real_))
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
