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
#' @param design From [.build_design()].
#' @param family A `gamRTMB_family`.
#' @param y Response.
#' @param fx Resolved fixed arguments.
#' @return A function of a parameter list `list(beta, b, log_sigma)`.
#' @keywords internal
.make_nll <- function(design, family, y, fx = list()) {
  parnames  <- design$parnames
  Xfix      <- design$Xfix
  beta_idx  <- design$beta_idx
  blocks    <- design$blocks
  par_blocks <- design$par_blocks
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
      eta <- as.vector(Xfix[[p]] %*% beta[beta_idx[[p]]])
      for (k in par_blocks[[p]])
        eta <- eta + as.vector(blocks[[k]]$X %*% b[blocks[[k]]$idx])
      theta[[p]] <- linkinv[[p]](eta)
    }
    jnll - sum(logdens(yo, theta, fx))
  }
}

#' Starting values for the variance components
#'
#' Cold-starting every log-sigma at zero is slow and can wander on flat
#' marginal surfaces. Instead pick \eqn{\sigma_k} so that the term's implied
#' prior standard deviation, \eqn{\sigma_k \sqrt{mean(rowSums(X_r^2))}}, is
#' `frac` of the rough scale of that parameter's linear predictor, which the
#' family supplies.
#'
#' @keywords internal
.init_log_sigma <- function(design, family, y, frac = 0.2) {
  sc <- family$eta_scale(y)
  vapply(design$blocks, function(bl) {
    rs <- sqrt(mean(rowSums(bl$X^2)))
    log(max(frac * sc[[bl$par]] / max(rs, 1e-8), 1e-4))
  }, numeric(1))
}

#' Assemble the full starting parameter list
#'
#' Intercepts start on the link scale from the family's `start()`; everything
#' else starts at zero. Blocks tied by an `id` are given a common starting
#' value, since only one of them survives the mapping.
#'
#' @keywords internal
.init_pars <- function(design, family, y, sigma_frac = 0.2, start = NULL) {
  beta0 <- numeric(design$nbeta)
  s0 <- family$start(y)
  for (p in design$parnames) {
    ii <- design$beta_idx[[p]]
    k <- match("(Intercept)", attr(ii, "labels"))
    if (!is.na(k)) beta0[ii[k]] <- s0[[p]]
  }
  ls0 <- stats::ave(.init_log_sigma(design, family, y, sigma_frac),
                    design$sig_group)
  pars <- list(beta = beta0, b = numeric(design$nb), log_sigma = ls0)
  if (!is.null(start)) pars[names(start)] <- start
  pars
}

#' Detect a distributional parameter that starts at a useless stationary point
#'
#' A parameter whose intercept has an identically zero score cannot move, and
#' a smooth on it then drifts on a flat surface instead of failing loudly.
#' `dskewnorm2`'s `alpha` does this at `alpha = 0`; see [rtmbdist_family()].
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
#' @param design,parnames Design and parameter names.
#' @return A message describing the affected parameters, or `NULL`.
#' @keywords internal
.flat_start <- function(grad, eval_at, pfull, is_beta, design, parnames) {
  tol <- 1e-8 * max(abs(grad), 1)
  f0 <- eval_at(pfull)
  gb <- grad[is_beta]
  flat <- character(0)
  for (p in parnames) {
    ii <- design$beta_idx[[p]]
    k <- match("(Intercept)", attr(ii, "labels"))
    if (is.na(k) || abs(gb[ii[k]]) >= tol) next
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
