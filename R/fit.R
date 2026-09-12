## The objective. Deliberately a plain closure over the design: it says
## nothing about how it will be optimised, and in particular nothing about
## which coefficients are treated as random effects. That is the engine's
## choice, which is what keeps a non-Laplace engine (see
## dev/NOTES-fellner-schall.md) a matter of adding a fitting function rather
## than rewriting the model.

#' Build the joint negative log-likelihood
#'
#' The penalized coefficients carry a mean-zero Gaussian prior, one variance
#' per penalized block. Which Gaussian depends on the route the block took
#' through [.build_design()]: after [mgcv::smooth2random()] the penalty is the
#' identity and the prior is iid \eqn{N(0, \sigma_k^2)}, while a block that
#' kept its own sparse penalty gets \eqn{N(0, \sigma_k^2 Q_k^{-1})} through
#' [RTMB::dgmrf()]. The two are the same model written two ways; keeping
#' \eqn{Q_k} sparse is what lets a Markov random field over many regions stay
#' affordable, since `smooth2random`'s rotation would fill it in.
#'
#' Null-space coefficients live in `beta` and are never given a prior.
#' Nothing is shared across smooths or across distributional parameters
#' unless an `id` says so.
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

    for (k in seq_along(blocks)) {
      bl <- blocks[[k]]
      th <- log_sigma[bl$theta_idx]
      jnll <- jnll - if (identical(bl$kind, "iid"))
        sum(dnorm(b[bl$idx], 0, exp(th[1L]), log = TRUE))
      else
        RTMB::dgmrf(b[bl$idx], 0, .block_prec(bl, th), log = TRUE)
    }

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
#' `frac` errs on the smooth side deliberately. Too flexible a start lets the
#' inner Newton solve push a parameter out of the family's support, which
#' shows up as a non-finite marginal objective; too clamped a start leaves the
#' smooth pinned to its null space, where the REML gradient in log-sigma is
#' nearly zero and the outer optimiser stalls.
#'
#' On well-behaved families the two failure modes bracket a wide, flat
#' optimum, and `frac` may as well not exist: over `norm`, `gamma2` and
#' `beta2` every value from 0.2 down to 0.001 converges on all fifteen fits
#' and agrees to seven significant figures. On four-parameter families it
#' matters, but **not monotonically, and no value dominates** -- 0.01 does
#' worst, with 0.2 and 0.001 on either side of it doing better, and the
#' spread between families is far larger than the spread across `frac`. So
#' 0.05 is a default rather than an optimum, kept because the evidence for
#' moving it is four fits out of twenty-five. The small values additionally
#' converge to a worse optimum more often, which is the pinned-to-null-space
#' mode above.
#'
#' Since no single value serves, a start that leaves the objective non-finite
#' is retried over [.sigma_frac_ladder()] rather than left to the user to
#' guess. See `dev/NOTES-sigma-frac.md` for the measurements, and
#' `dev/bench-sigma-frac.R` to reproduce them.
#'
#' @keywords internal
.init_pars <- function(design, family, y, frac = 0.05, start = NULL) {
  beta0 <- numeric(design$nbeta)
  s0 <- family$start(y)
  for (p in design$parnames) {
    ii <- design$beta_idx[[p]]
    k <- match("(Intercept)", attr(ii, "labels"))
    if (!is.na(k)) beta0[ii[k]] <- s0[[p]]
  }
  sc <- family$eta_scale(y)
  ls0 <- as.numeric(unlist(lapply(design$blocks, function(bl) {
    rs <- sqrt(mean(Matrix::rowSums(bl$X^2)))
    s1 <- log(max(frac * sc[[bl$par]] * bl$qscale / max(rs, 1e-8), 1e-4))
    ## A block with its own parameterisation is not a variance, so the rule
    ## above does not apply to it. Its constructor supplies the start.
    if (bl$ntheta == 1L) s1 else bl$theta_start
  })))
  pars <- list(beta = beta0, b = numeric(design$nb),
               log_sigma = stats::ave(ls0, design$sig_group))
  if (!is.null(start)) pars[names(start)] <- start
  pars
}

## Alternative variance-component starts, tried in this order when the first
## marginal evaluation is not finite. They are not a refinement of one
## another: over eight families and five seeds the number of fits reaching a
## converged solution is not monotone in `frac` and no value dominates, so the
## point of the ladder is coverage rather than a better default. Ordered by
## that table, with 0.001 last because it converges most often and misses the
## optimum most often. See .init_pars() for what `frac` means, and
## dev/NOTES-sigma-frac.md for the measurements.
.sigma_frac_ladder <- c(0.005, 0.2, 0.001)

## Inner iteration cap for the probe below. A well-posed inner solve reaches
## its mode in a few dozen Newton steps from these starting values, so this is
## generous; a problem that needs more than this is one the cap is looking for.
.probe_maxit <- 100L

#' Is the marginal objective finite, cheaply?
#'
#' The first marginal evaluation is where a badly posed inner problem shows
#' up, and it is an expensive place to discover it: the inner Newton runs all
#' the way to its iteration cap before handing back a non-finite value. On the
#' four-parameter Box-Cox families that is over two minutes to learn that the
#' fit will not start. Probing with a short cap answers the same question in
#' seconds, and costs nothing on a model that was going to work.
#'
#' A probe is only ever an accelerator: a `FALSE` sends the caller to the next
#' rung of [.sigma_frac_ladder], and if every rung fails the fit proceeds from
#' the original starting values under the full cap, exactly as it would have.
#' So a probe that is wrong about a slow-but-sound inner solve costs a few
#' seconds, never a fit.
#'
#' @keywords internal
.probe_finite <- function(obj) {
  keep <- obj$env$inner.control$maxit
  on.exit(obj$env$inner.control$maxit <- keep, add = TRUE)
  obj$env$inner.control$maxit <- min(.probe_maxit, keep)
  v <- tryCatch(suppressWarnings(obj$fn(obj$par)),
                error = function(e) NA_real_)
  length(v) == 1L && is.finite(v)
}

#' Label every coefficient the inner problem solves over
#'
#' In the order [RTMB::MakeADFun()] lays them out, which is the order of the
#' parameter list: all of `beta`, then all of `b`. Under `"ML"` only `b` is
#' random, so only those are named.
#'
#' @keywords internal
.random_labels <- function(design, method) {
  bl <- unlist(lapply(design$blocks, function(z)
    paste0(z$par, ":", z$label, ".", seq_len(z$q))))
  if (method == "REML") c(.beta_labels(design), bl) else bl
}

#' Diagnose a non-finite marginal objective as non-concavity
#'
#' A non-finite marginal objective is usually read as the response leaving the
#' family's support, and the error message used to say so. That is the wrong
#' diagnosis for a whole class of models, and a misleading one: the Laplace
#' approximation needs the log determinant of the inner Hessian, so an
#' *indefinite* Hessian produces exactly the same symptom with the data
#' entirely inside the support.
#'
#' The Box-Cox power exponential is the case in hand. Its log-likelihood is
#' concave in the four intercepts alone -- an intercept-only fit converges in
#' half a second -- but adding a single covariate column to `mu` puts three
#' negative eigenvalues into the inner Hessian, every one of them a direction
#' mixing that column with the `sigma`, `nu` and `tau` intercepts. No starting
#' value repairs it: sweeping `tau` from 2 to 9 and `nu` from 1 to 2.5 never
#' gets below two negative directions. It is a property of the family's
#' parameterisation, not of the start.
#'
#' Under `"REML"` the negative curvature lands in `beta`, which is declared
#' random and so passes through the inner solve carrying no prior to convexify
#' it. The penalized blocks are not the problem -- their Gaussian prior leaves
#' them comfortably positive definite. Hence the suggestion of a basis with no
#' null space, which is what puts those columns under a penalty.
#'
#' @param obj The `MakeADFun` object, evaluated at its starting values.
#' @param design The design object.
#' @param method `"REML"` or `"ML"`.
#' @param famname The family's name, for the message.
#' @return A sentence describing the negative curvature, or `NULL` if the
#'   Hessian is unavailable or positive definite.
#' @keywords internal
.inner_indefinite <- function(obj, design, method, famname) {
  h <- tryCatch(obj$env$spHess(obj$env$par, random = TRUE),
                error = function(e) NULL)
  ## A Cholesky is the cheap question ("is this positive definite?"); the
  ## eigen decomposition is only needed to name the directions, and is dense,
  ## so it is reserved for problems small enough to afford it.
  if (!isFALSE(.is_pd(h))) return(NULL)
  base <- paste0("the inner Hessian is not positive definite at the starting ",
                 "values, so the Laplace approximation's log determinant is ",
                 "undefined there")
  if (ncol(h) > 1000L) return(paste0(base, "."))
  e <- tryCatch(eigen(as.matrix(h), symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) return(paste0(base, "."))
  neg <- which(e$values < 0)
  if (!length(neg)) return(NULL)
  lab <- .random_labels(design, method)
  ## The flattest negative direction is the informative one: the deepest is
  ## dominated by whichever coefficient happens to carry the most curvature,
  ## which is usually an intercept and says nothing. `eigen()` sorts
  ## decreasing, so the flattest negative is the first of them.
  v <- abs(e$vectors[, neg[1L]])
  top <- utils::head(order(-v), 3L)
  paste0(base, " (", length(neg), " negative ",
         if (length(neg) > 1L) "directions" else "direction",
         ", the flattest dominated by ",
         paste0("`", lab[top], "`", collapse = ", "),
         "). This is non-concavity of family '", famname,
         "' in its own parameterisation rather than a bad starting value, and ",
         "under REML it is the unpenalized coefficients that carry it. A ",
         "basis with no null space puts them under a penalty: try bs = \"cs\" ",
         "or bs = \"ts\", or method = \"ML\".")
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
#' A plain right-hand side is shorthand for modelling the family's *first*
#' parameter and leaving the rest constant, so `y ~ s(x)` means
#' `y ~ list(mean = ~ s(x))` for `fam("norm")`. That is the ordinary gam
#' formula, and it is what makes the one-parameter case read like one.
#'
#' @section Data:
#' `data` is optional. Without it the model variables are looked up where the
#' formula was written, as in [stats::glm()] or [mgcv::gam()]; they are
#' collected into a data frame first, so everything downstream — `na.action`,
#' [predict()], the plots — sees the same rectangle either way. Variables
#' found this way must all have the same length.
#'
#' @section Sparse penalties:
#' A smooth whose penalty is a sparse precision matrix -- `bs = "mrf"` over an
#' adjacency graph, a random walk, or any precision supplied through
#' `xt = list(penalty = )` -- can skip [mgcv::smooth2random()]. That rotation
#' makes the coefficients iid, which is convenient but fills the penalty in
#' completely; keeping it instead and giving the block an
#' \eqn{N(0, \sigma^2 Q^{-1})} prior through [RTMB::dgmrf()] is the same
#' model at a fraction of the cost. The `sparse` argument controls this.
#'
#' An intrinsic field is corner-constrained rather than sum-to-zero
#' constrained, since the latter is what destroys the sparsity. The two are
#' equivalent up to a constant absorbed by the intercept, so such a term needs
#' its parameter to have one. See [.null_space()].
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
#'   per-parameter formulas, or a plain right-hand side for the family's first
#'   parameter.
#' @param family A `gamRTMB_family`, from [fam()]. See [families()].
#' @param data A data frame holding every model variable, or `NULL` (the
#'   default) to take them from the environment of `formula`.
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
#'   values: each smooth starts contributing this fraction of its parameter's
#'   linear-predictor scale. See [.init_pars()]. Raise it if a fit converges
#'   to an over-smooth solution. If the objective is not finite here, a few
#'   other values are tried automatically before giving up (see
#'   [.sigma_frac_ladder()]) and the one used is reported; passing this
#'   argument explicitly does not switch that off, but passing `start` does.
#' @param sparse How to treat a smooth whose single penalty is already sparse
#'   -- a Markov random field, a random walk, a supplied GMRF precision.
#'   `"auto"` (default) keeps the penalty and gives the block a
#'   [RTMB::dgmrf()] prior when it has at least 50 coefficients and is at
#'   most 20% nonzero, and sends everything else through
#'   [mgcv::smooth2random()] as usual. `"never"` is the old behaviour;
#'   `"always"` takes the sparse route for every single-penalty smooth, which
#'   is mainly useful for checking that the two agree. See [.gmrf_block()].
#' @param joint_precision Ask [RTMB::sdreport()] for the joint precision
#'   matrix, which [predict()] needs for standard errors. On by default; turn
#'   it off to save time and memory on large models where bands are not
#'   wanted.
#' @param start Optional named list overriding entries of the starting
#'   parameter list (`beta`, `b`, `log_sigma`).
#' @param silent Passed to [RTMB::MakeADFun()].
#' @param control Passed to [stats::nlminb()].
#' @param inner_control Passed to [RTMB::MakeADFun()]'s `inner.control`, which
#'   governs the inner Newton solve over the coefficients rather than the
#'   outer optimisation of the smoothing parameters. `list(maxit = ...)` is
#'   the entry worth reaching for; see `dev/NOTES-inner-method.md` for why
#'   `inner.method` is not exposed.
#' @return An object of class `gamRTMB`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = runif(200), x2 = runif(200))
#' d$y <- rnorm(200, sin(2 * pi * d$x1), exp(-1 + d$x2))
#' fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)),
#'                data = d)
#' fit
#' edf(fit)
#'
#' ## a plain right-hand side models the family's first parameter
#' gamRTMB(y ~ s(x1, k = 8), data = d)
#' @export
gamRTMB <- function(formula, family = fam("norm"), data = NULL, weights = NULL,
                    na.action = stats::na.omit, knots = NULL,
                    method = c("REML", "ML"),
                    engine = c("laplace", "efs"), sigma_frac = 0.05,
                    sparse = c("auto", "never", "always"),
                    joint_precision = TRUE, start = NULL, silent = TRUE,
                    control = list(), inner_control = list()) {
  method <- match.arg(method)
  engine <- match.arg(engine)
  sparse <- match.arg(sparse)
  if (!inherits(family, "gamRTMB_family"))
    stop("`family` must be a gamRTMB_family, e.g. fam(\"norm\") or ",
         "fam(\"gamma2\"); see families()")
  if (!is.null(data) && !is.data.frame(data))
    stop("`data` must be a data frame, or NULL to use the formula's environment")

  ## `parent.frame()` is taken here rather than passed along lazily: inside a
  ## promise it would resolve against whatever frame forced it.
  caller <- parent.frame()
  fenv <- if (is.null(environment(formula))) caller else environment(formula)

  w <- if (is.null(data)) eval(substitute(weights), caller)
       else eval(substitute(weights), data, caller)
  pf <- .parse_formula(formula, family$parnames)
  md <- .model_data(pf$response, pf$par_formulas, data, w, na.action, fenv)
  data <- md$data; y <- md$y; w <- md$weights
  design <- .build_design(pf$par_formulas, data, family$parnames,
                          knots = knots, sparse = sparse)
  fx <- .resolve_fixed(family, data, length(y))
  pars <- .init_pars(design, family, y, sigma_frac, start)
  nll <- .make_nll(design, family, y, fx, w)

  ## The ladder needs to be able to rebuild the starting values at another
  ## `sigma_frac`; an explicit `start` is the user's and is never overwritten,
  ## so the retries are switched off in that case.
  repars <- if (is.null(start))
    function(frac) .init_pars(design, family, y, frac, NULL) else NULL

  fit <- switch(engine,
    laplace = .fit_laplace(nll, pars, design, method, joint_precision, silent,
                           control, family, inner_control, repars),
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
#' @section Getting started at all:
#' Before any of that, the marginal objective has to be finite at the starting
#' values, and on the harder families it often is not. That is handled in two
#' steps, both of which are about failing cheaply rather than about finding a
#' better start: [.probe_finite()] asks the question with a short inner
#' iteration cap, and a `FALSE` sends the caller to the next rung of
#' [.sigma_frac_ladder()]. If every rung fails, the fit proceeds from the
#' original starting values under the full cap, so the ladder can only add
#' fits, never remove one.
#'
#' @param repars Function of a `sigma_frac` returning a fresh starting
#'   parameter list, used to walk the ladder. `NULL` disables the retries.
#' @param inner_control Passed to [RTMB::MakeADFun()]'s `inner.control`.
#' @keywords internal
.fit_laplace <- function(nll, pars, design, method, joint_precision, silent,
                          control, family, inner_control = list(),
                          repars = NULL) {
  random <- if (method == "REML") c("beta", "b") else "b"
  map <- if (design$nsigma_free < design$nsigma)
    list(log_sigma = design$sig_group) else list()

  ## Merged onto MakeADFun's own default rather than replacing it: handing it
  ## a bare list() would drop `maxit` and leave the inner Newton running on
  ## newton()'s much smaller formal default, which is a silent change of
  ## behaviour for every fit that passes nothing. Pinned here rather than read
  ## back out of RTMB so that an upstream change cannot move it either.
  ic <- utils::modifyList(list(maxit = 1000L), inner_control)
  build <- function(p) RTMB::MakeADFun(nll, p, random = random, map = map,
                                       silent = silent, inner.control = ic)
  obj <- build(pars)

  ## Walk the ladder only if the first start does not work. The probe is
  ## skipped when there is nothing to walk to, so a fit with `repars = NULL`
  ## or an explicit `start` behaves exactly as it did before.
  if (!is.null(repars) && !.probe_finite(obj)) {
    for (fr in .sigma_frac_ladder) {
      cand <- tryCatch(build(repars(fr)), error = function(e) NULL)
      if (is.null(cand) || !.probe_finite(cand)) next
      message("the marginal objective was not finite at the default ",
              "variance-component start; refitting with sigma_frac = ", fr,
              ". Pass sigma_frac explicitly to pin this down.")
      obj <- cand
      attr(obj, "sigma_frac") <- fr
      break
    }
  }

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
  ##
  ## Indefinite curvature is the other way to get here, and used to be
  ## reported as a support problem, which is both wrong and a hard thing to
  ## recover from as a user. Ask the Hessian before guessing; see
  ## [.inner_indefinite()].
  if (!is.finite(obj$fn(obj$par)))
    stop("the objective is not finite at the starting values. ",
         ## `.or_else` is lazy in its second argument, so the Hessian is only
         ## factorised when there is no cheaper explanation to give.
         .or_else(msg, .or_else(
           .inner_indefinite(obj, design, method, family$family),
           paste0("Check that the response is in the support of family '",
                  family$family, "', and consider passing start = ",
                  "list(beta = ...), or a different sigma_frac."))),
         call. = FALSE)
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
       sigma_frac_used = attr(obj, "sigma_frac"),
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
