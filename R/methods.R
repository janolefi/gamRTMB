#' @export
print.gamRTMB <- function(x, ...) {
  D <- x$design
  cat("gamRTMB fit\n")
  cat("  family:    ", x$family$family, " (",
      paste(sprintf("%s/%s", x$family$parnames, x$family$links), collapse = ", "),
      ")\n", sep = "")
  cat("  criterion: ", x$method, "   engine: ", x$engine, "\n", sep = "")
  cat("  converged: ", x$convergence, "   -logLik: ",
      sprintf("%.4f", x$objective), "   max|grad|: ",
      sprintf("%.3g", x$max_grad), "\n", sep = "")
  cat("  observations: ", D$n,
      if (isTRUE(x$dropped > 0)) paste0(" (", x$dropped, " dropped: missing)") else "",
      if (!is.null(x$weights)) ", prior weights" else "", "\n", sep = "")
  cat("  coefficients: ", D$nbeta, " fixed (incl. null spaces), ", D$nb,
      " penalized; ", D$nsigma, " smoothing parameter",
      if (D$nsigma != 1L) "s" else "",
      if (D$nsigma_free < D$nsigma)
        paste0(" (", D$nsigma_free, " free, ", D$nsigma - D$nsigma_free,
               " tied by id)") else "",
      "\n", sep = "")
  invisible(x)
}

#' Labels for the fixed-effect coefficients
#'
#' `"<parameter>:<coefficient>"`, in the order of the `beta` vector.
#'
#' @keywords internal
.beta_labels <- function(design) {
  unlist(lapply(design$parnames, function(p)
    paste0(p, ":", attr(design$beta_idx[[p]], "labels"))))
}

#' Covariance of the fixed-effect coefficients
#'
#' The `(beta, beta)` block of the inverse joint precision, so it includes the
#' uncertainty in the smoothing parameters rather than conditioning on them
#' (mgcv's `unconditional = TRUE`). Needs a fit made with
#' `joint_precision = TRUE`, which is the default.
#'
#' @param object A `gamRTMB` fit.
#' @param ... Ignored.
#' @return A labelled covariance matrix over the parametric and smooth
#'   null-space coefficients.
#' @export
vcov.gamRTMB <- function(object, ...) {
  Vj <- .joint_cov(object)
  V <- Vj$V[Vj$ib, Vj$ib, drop = FALSE]
  lab <- .beta_labels(object$design)
  dimnames(V) <- list(lab, lab)
  V
}

#' Log-likelihood, and hence AIC and BIC
#'
#' Returns the \strong{log-likelihood} of the data at the fitted
#' coefficients, with `df` set to the total effective degrees of freedom.
#' That is deliberately not the fitting criterion: under REML the criterion is
#' a restricted one, is not a likelihood, and is not comparable across
#' different mean structures, so using it for AIC would be wrong. It remains
#' available as `fit$objective`.
#'
#' With `df` = total EDF this is the convention used by [mgcv::gam()] and by
#' GAMLSS's GAIC, so `AIC()` and `BIC()` are comparable across models fitted
#' to the same response. Under `method = "ML"`, where effective degrees of
#' freedom are unavailable, `df` falls back to counting the fixed
#' coefficients and free smoothing parameters.
#'
#' @param object A `gamRTMB` fit.
#' @param ... Ignored.
#' @return An object of class `logLik`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(200)); d$y <- rnorm(200, sin(2 * pi * d$x), 0.3)
#' fit <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d)
#' logLik(fit)
#' AIC(fit)
#' @export
logLik.gamRTMB <- function(object, ...) {
  theta <- stats::predict(object, type = "response")
  ld <- object$family$logdens(object$y, theta, object$fixed)
  ll <- if (is.null(object$weights)) sum(ld) else sum(object$weights * ld)
  df <- tryCatch(attr(edf(object), "edf.total"), error = function(e)
    object$design$nbeta + object$design$nsigma_free)
  structure(ll, df = df, nobs = object$design$n, class = "logLik")
}

#' @param object A `gamRTMB` fit.
#' @param ... Ignored.
#' @describeIn gamRTMB Number of observations actually used.
#' @export
nobs.gamRTMB <- function(object, ...) object$design$n

#' Summarise a gamRTMB fit
#'
#' Laid out like [mgcv::summary.gam()], but with one row per distributional
#' parameter throughout, since every parameter has its own formula, link,
#' coefficients and smooths.
#'
#' Parametric coefficients get standard errors from [vcov.gamRTMB()] and
#' approximate Gaussian (z) tests. Smooth terms are reported with their
#' effective degrees of freedom and smoothing parameters but \strong{no
#' p-values}: mgcv's "approximate significance" rests on Wood's (2013) test
#' for a term's whole coefficient block, which is not implemented here, and a
#' naive Wald test in its place would be misleading.
#'
#' @param object A `gamRTMB` fit.
#' @param ... Ignored.
#' @return An object of class `summary.gamRTMB`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = runif(300), x2 = runif(300))
#' d$y <- rnorm(300, sin(2 * pi * d$x1), exp(-1 + d$x2))
#' summary(gamRTMB(y ~ list(mean = ~ s(x1), sd = ~ s(x2)), data = d))
#' @export
summary.gamRTMB <- function(object, ...) {
  D <- object$design
  est <- object$coefficients$beta
  lab <- .beta_labels(D)
  keep <- !grepl("\\.null[0-9]+$", lab)          # null-space bases are not
                                                # parametric terms
  se <- tryCatch(sqrt(diag(vcov(object)))[keep], error = function(e) NULL)
  cf <- if (is.null(se)) {
    matrix(est[keep], ncol = 1L, dimnames = list(lab[keep], "Estimate"))
  } else {
    z <- est[keep] / se
    cbind(Estimate = est[keep], `Std. Error` = se, `z value` = z,
          `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  }

  e <- tryCatch(edf(object), error = function(err) NULL)
  ll <- tryCatch(stats::logLik(object), error = function(err) NULL)
  ## right-pad the parameter names so the tildes line up
  wp <- max(nchar(D$parnames))
  fml <- vapply(D$parnames, function(p)
    paste(formatC(p, width = wp), deparse(object$par_formulas[[p]][[2L]]),
          sep = " ~ "), "")

  structure(list(family = object$family, method = object$method,
                 engine = object$engine, formulas = fml,
                 coefficients = cf, smooth = e,
                 n = D$n, dropped = object$dropped,
                 weighted = !is.null(object$weights),
                 objective = object$objective, convergence = object$convergence,
                 logLik = ll, edf.total = if (!is.null(e)) attr(e, "edf.total"),
                 aic = if (!is.null(ll)) stats::AIC(ll)),
            class = "summary.gamRTMB")
}

#' @param x A `summary.gamRTMB`.
#' @rdname summary.gamRTMB
#' @export
print.summary.gamRTMB <- function(x, ...) {
  cat("\nFamily: ", x$family$family, "   [", x$family$dist, " from ",
      x$family$source, "]\n", sep = "")
  cat("Links:  ", paste(sprintf("%s = %s", x$family$parnames, x$family$links),
                        collapse = ",  "), "\n", sep = "")
  cat("\nFormula:\n")
  cat(paste0("  ", x$formulas, collapse = "\n"), "\n", sep = "")

  cat("\nParametric coefficients:\n")
  if (ncol(x$coefficients) == 1L) {
    print(x$coefficients, digits = 4)
    cat("(standard errors need joint_precision = TRUE)\n")
  } else {
    stats::printCoefmat(x$coefficients, digits = 4, signif.stars = TRUE,
                        has.Pvalue = TRUE, P.values = TRUE)
  }

  if (!is.null(x$smooth) && nrow(x$smooth)) {
    cat("\nSmooth terms:\n")
    tm <- paste0(x$smooth$parameter, ":", x$smooth$term)
    sm <- data.frame(term = formatC(tm, width = max(nchar(tm)), flag = "-"),
                     edf = round(x$smooth$edf, 3), k = x$smooth$k,
                     sp = x$smooth$sp)
    if (any(nzchar(x$smooth$id))) sm$id <- x$smooth$id   # only when used
    print(sm, row.names = FALSE, right = TRUE)
  }

  cat("\n")
  if (!is.null(x$edf.total))
    cat("Total EDF = ", sprintf("%.2f", x$edf.total), "   ", sep = "")
  cat("n = ", x$n, sep = "")
  if (isTRUE(x$dropped > 0))
    cat("  (", x$dropped, " observation", if (x$dropped > 1) "s", " deleted: missingness)",
        sep = "")
  if (x$weighted) cat("  [prior weights]")
  cat("\n")
  cat("-", x$method, " = ", sprintf("%.3f", x$objective), sep = "")
  if (!is.null(x$logLik))
    cat("   logLik = ", sprintf("%.3f", as.numeric(x$logLik)),
        "   AIC = ", sprintf("%.2f", x$aic), sep = "")
  cat("\n")
  if (!isTRUE(x$convergence)) cat("** the optimiser did not converge **\n")
  invisible(x)
}

#' Predictions from a gamRTMB fit
#'
#' New data goes through [mgcv::PredictMat()] on the stored `smoothCon`
#' objects, so the basis, knots and constraint matrices are exactly those of
#' the fit. Smooths are never re-fitted on new data, which would place
#' different knots.
#'
#' @param object A `gamRTMB` fit.
#' @param newdata Optional data frame. If omitted, the fitting data's
#'   in-sample linear predictors are returned.
#' @param type `"link"` (default) for the linear predictors, `"response"` for
#'   the parameters on their natural scale, `"terms"` for the contribution of
#'   each smooth separately. As in mgcv, `"terms"` excludes the intercept and
#'   any offset, which belong to the predictor rather than to a smooth.
#' @param se.fit Also return standard errors; needs a fit made with
#'   `joint_precision = TRUE`.
#' @param ... Ignored.
#' @return For `type = "link"`/`"response"`, a named list with one vector per
#'   distributional parameter (or, with `se.fit = TRUE`, a list of `fit` and
#'   `se.fit`). For `type = "terms"`, a named list with a matrix of per-term
#'   contributions for each parameter.
#' @export
predict.gamRTMB <- function(object, newdata = NULL,
                            type = c("link", "response", "terms"),
                            se.fit = FALSE, ...) {
  type <- match.arg(type)
  D <- object$design
  Vj <- if (se.fit) .joint_cov(object) else NULL
  out <- list(); se <- list(); trm <- list()

  for (p in D$parnames) {
    P <- D$parts[[p]]
    ## an offset() term is rebuilt from newdata, never carried over
    if (is.null(newdata)) {
      Xpara <- P$Xpara; off <- P$offset
    } else {
      tt <- stats::delete.response(P$terms)
      mf <- stats::model.frame(tt, newdata, xlev = P$xlev)
      Xpara <- stats::model.matrix(tt, mf)
      off <- stats::model.offset(mf)
      if (is.null(off)) off <- 0
    }
    npara <- ncol(P$Xpara)
    bfull <- object$coefficients$beta[D$beta_idx[[p]]]
    eta <- as.vector(Xpara %*% bfull[seq_len(npara)]) + off

    nm <- vapply(P$smooths, `[[`, "", "label")
    tm <- matrix(0, nrow(Xpara), length(P$smooths), dimnames = list(NULL, nm))
    tse <- tm
    Xs_all <- vector("list", length(P$smooths))
    for (j in seq_along(P$smooths)) {
      s <- P$smooths[[j]]
      Xs <- if (is.null(newdata)) s$sm$X else mgcv::PredictMat(s$sm, newdata)
      Xs_all[[j]] <- Xs
      tm[, j] <- as.vector(Xs %*% .smooth_beta(object, p, j))
      eta <- eta + tm[, j]
      if (se.fit) {
        Vb <- .smooth_vcov(object, p, j, Vj)
        tse[, j] <- sqrt(pmax(rowSums((Xs %*% Vb) * Xs), 0))
      }
    }
    out[[p]] <- if (type == "response")
      .links[[object$family$links[[p]]]]$linkinv(eta) else eta
    trm[[p]] <- list(fit = tm, se = if (se.fit) tse else NULL)
    if (se.fit) {
      Xall <- do.call(cbind, c(list(Xpara), Xs_all))
      V <- .eta_vcov(object, p, Vj)
      se[[p]] <- sqrt(pmax(rowSums((Xall %*% V) * Xall), 0))
    }
  }
  if (type == "terms") return(trm)
  if (se.fit) list(fit = out, se.fit = se) else out
}

#' @param object A `gamRTMB` fit.
#' @param ... Ignored.
#' @describeIn gamRTMB Coefficients, as a list of the fixed (`beta`) and
#'   penalized (`b`) vectors plus the log variance components.
#' @export
coef.gamRTMB <- function(object, ...) {
  c(object$coefficients, list(log_sigma = object$log_sigma))
}

#' @describeIn gamRTMB Fitted values of every distributional parameter, on the
#'   response scale.
#' @export
fitted.gamRTMB <- function(object, ...) {
  stats::predict(object, type = "response")
}

