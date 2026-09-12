#' @export
print.gamRTMB <- function(x, ...) {
  D <- x$design
  cat("gamRTMB fit\n")
  cat("  family:    ", x$family$family, " (",
      paste(sprintf("%s/%s", x$family$parnames, x$family$links), collapse = ", "),
      ")\n", sep = "")
  cat("  criterion: ", x$method,
      if (identical(x$method, "aREML")) "  (extended Fellner-Schall)" else "",
      "\n", sep = "")
  ## The Fellner-Schall gradient drops a third-derivative term and so does not
  ## reach zero at the optimum; naming it differently keeps it from being read
  ## as a stationarity measure. See [.fit_efs()].
  cat("  converged: ", x$convergence, "   -", .criterion_label(x$method), ": ",
      sprintf("%.4f", x$objective),
      if (identical(x$method, "aREML")) "   max|FS grad|: "
      else "   max|grad|: ",
      sprintf("%.3g", x$max_grad), "\n", sep = "")
  if (!isTRUE(x$convergence) && nzchar(.or_else(x$opt$message, "")))
    cat("  optimiser: ", x$opt$message, "\n", sep = "")
  ## A repaired data Hessian means the criterion reported above is not the one
  ## the model posed, and neither are the EDF; see [.psd_repair()].
  if (isTRUE(x$psd_repairs > 0))
    cat("  note:      the data Hessian was repaired at ", x$psd_repairs,
        " iteration", if (x$psd_repairs != 1L) "s",
        if (isTRUE(x$psd_repaired_at_mode)) ", including the last" else
          " (but not the last)", "\n", sep = "")
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

#' What to call the reported criterion value
#'
#' `"aREML"` optimises the REML criterion -- the approximation is in the
#' gradient, not in the quantity -- so its value is a REML value and is
#' directly comparable with one from `method = "REML"`. The header line says
#' which route produced it, so labelling the number for what it is costs no
#' ambiguity and gains a comparison.
#'
#' @keywords internal
.criterion_label <- function(method)
  if (identical(method, "aREML")) "REML" else method

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
#' Under `method = "REML"` or `"ML"`, the `(beta, beta)` block of the inverse
#' joint precision, so it includes the uncertainty in the smoothing parameters
#' rather than conditioning on them (mgcv's `unconditional = TRUE`). Needs a
#' fit made with `joint_precision = TRUE`, which is the default.
#'
#' \strong{Under `method = "aREML"` it is conditional instead} -- mgcv's
#' `unconditional = FALSE` -- and so a little too narrow. There is no joint
#' precision to take a block of: its smoothing-parameter part comes from
#' differentiating the marginal criterion twice in the smoothing parameters,
#' which is exactly the term the extended Fellner-Schall method exists to
#' avoid computing. What that criterion has is the penalized Hessian, and its
#' inverse is the covariance given the smoothing parameters it settled on.
#' Everything downstream inherits this: the standard errors in
#' [summary.gamRTMB()], the bands from [predict.gamRTMB()] and the term plots.
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
                 formulas = fml,
                 coefficients = cf, smooth = e,
                 n = D$n, dropped = object$dropped,
                 weighted = !is.null(object$weights),
                 objective = object$objective, convergence = object$convergence,
                 logLik = ll,
                 ## edf() has already decided this by returning NA, and has
                 ## warned; reading it back keeps the judgement in one place
                 edf_defined = is.null(e) || !all(is.na(e$edf)),
                 edf.total = if (!is.null(e)) attr(e, "edf.total"),
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
  if (!is.null(x$edf.total) && !is.na(x$edf.total))
    cat("Total EDF = ", sprintf("%.2f", x$edf.total), "   ", sep = "")
  cat("n = ", x$n, sep = "")
  if (isTRUE(x$dropped > 0))
    cat("  (", x$dropped, " observation", if (x$dropped > 1) "s", " deleted: missingness)",
        sep = "")
  if (x$weighted) cat("  [prior weights]")
  cat("\n")
  cat("-", .criterion_label(x$method), " = ", sprintf("%.3f", x$objective),
      sep = "")
  if (!is.null(x$logLik)) {
    cat("   logLik = ", sprintf("%.3f", as.numeric(x$logLik)), sep = "")
    ## AIC counts the EDF as its parameter count, so it goes when they do
    if (!is.null(x$aic) && !is.na(x$aic))
      cat("   AIC = ", sprintf("%.2f", x$aic), sep = "")
  }
  cat("\n")
  if (!isTRUE(x$convergence)) cat("** the optimiser did not converge **\n")
  if (isFALSE(x$edf_defined))
    cat("** the effective degrees of freedom are undefined at these values, ",
        "not merely\n   imprecise: the data Hessian is not positive ",
        "semi-definite here **\n", sep = "")
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
#' @section Quantiles:
#' `type = "quantile"` evaluates the fitted distribution's quantile function
#' at each observation, giving covariate-dependent quantiles — which is much
#' of the point of letting every parameter vary. It needs the family to have a
#' quantile function; [families()] reports which do.
#'
#' Their standard errors come from the delta method. The quantile is
#' differentiated numerically with respect to each linear predictor —
#' \emph{on the link scale}, so the link's own derivative is absorbed into
#' the difference and never has to be supplied — and those derivatives are
#' combined with the predictors' joint covariance, cross-parameter terms
#' included. Checked against the analytic normal case to a relative 1e-9.
#'
#' Offered only for a continuous response, because for a lattice or mixed one
#' the quantile function is a step and its derivative is not meaningful.
#'
#' The result is a standard error on the \strong{response} scale, since a
#' quantile is a value of the response and has no link of its own. A symmetric
#' interval built from it can therefore cross a boundary of the support — for
#' a positive response, `q - 2 * se` can be negative in the lower tail where
#' the quantile is small and its uncertainty is not. Read such an interval as
#' a local measure of precision rather than a range of plausible values, or
#' take percentiles of quantiles simulated from the joint posterior instead.
#'
#' @param type `"link"` (default) for the linear predictors, `"response"` for
#'   the parameters on their natural scale, `"terms"` for the contribution of
#'   each smooth separately, or `"quantile"` for quantiles of the fitted
#'   distribution. As in mgcv, `"terms"` excludes the intercept and any
#'   offset, which belong to the predictor rather than to a smooth.
#' @param prob Probabilities for `type = "quantile"`.
#' @param se.fit Also return standard errors; needs a fit made with
#'   `joint_precision = TRUE`.
#' @param ... Ignored.
#' @return For `type = "link"`/`"response"`, a named list with one vector per
#'   distributional parameter. For `type = "terms"`, a named list with a matrix
#'   of per-term contributions for each parameter. For `type = "quantile"`, a
#'   matrix with one column per probability. With `se.fit = TRUE`, a list of
#'   `fit` and `se.fit` in the same shape.
#' @export
predict.gamRTMB <- function(object, newdata = NULL,
                            type = c("link", "response", "terms", "quantile"),
                            prob = c(0.05, 0.25, 0.5, 0.75, 0.95),
                            se.fit = FALSE, ...) {
  type <- match.arg(type)
  D <- object$design
  Vj <- if (se.fit) .joint_cov(object) else NULL
  out <- list(); se <- list(); trm <- list(); forms <- list()

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
    Z_all <- vector("list", length(P$smooths))
    cols <- integer(0)
    for (j in seq_along(P$smooths)) {
      s <- P$smooths[[j]]
      Xs <- if (is.null(newdata)) s$sm$X else mgcv::PredictMat(s$sm, newdata)
      sp <- .smooth_part(object, p, j, Xs)
      Z_all[[j]] <- sp$Z
      tm[, j] <- as.vector(sp$Z %*% sp$coef)
      eta <- eta + tm[, j]
      if (se.fit) {
        ii <- c(Vj$ir[sp$b], Vj$ib[sp$f])
        cols <- c(cols, ii)
        tse[, j] <- .qform_se(sp$Z, Vj$V[ii, ii, drop = FALSE])
      }
    }
    out[[p]] <- if (type == "link") eta
      else .links[[object$family$links[[p]]]]$linkinv(eta)
    trm[[p]] <- list(fit = tm, se = if (se.fit) tse else NULL)
    ## the whole predictor is the same linear form, stacked
    L <- do.call(cbind, c(list(Xpara), Z_all))
    forms[[p]] <- list(eta = eta, L = L,
                       ii = if (se.fit)
                         c(Vj$ib[D$beta_idx[[p]][seq_len(npara)]], cols))
    if (se.fit)
      se[[p]] <- .qform_se(L, Vj$V[forms[[p]]$ii, forms[[p]]$ii, drop = FALSE])
  }

  if (type == "terms") return(trm)
  if (type == "quantile") {
    fam <- object$family
    if (is.null(fam$qf))
      stop("family '", fam$family, "' has no quantile function in ",
           fam$source, ", so quantile predictions are not available. ",
           "families() reports which families support them.")
    nn <- length(out[[1L]])
    ## a fixed argument taken from the data (binomial trials, say) belongs to
    ## the rows being predicted, not to the rows that were fitted
    fx <- if (is.null(newdata)) object$fixed
          else .resolve_fixed(object$family, newdata, nn)
    q <- vapply(prob, function(pr) fam$qf(rep(pr, nn), out, fx), numeric(nn))
    dim(q) <- c(nn, length(prob))          # vapply drops to a vector when n = 1
    dimnames(q) <- list(NULL, paste0("q", prob))
    if (!se.fit) return(q)
    if (fam$support != "continuous")
      stop("standard errors for quantiles need a continuous response; for a ",
           substr(fam$support, 1, 20), " one the quantile function is a step")
    sq <- .quantile_se(object, prob, forms, Vj, fx = fx)
    dim(sq) <- dim(q); dimnames(sq) <- dimnames(q)
    return(list(fit = q, se.fit = sq))
  }
  if (se.fit) list(fit = out, se.fit = se) else out
}

#' Delta-method standard errors for fitted quantiles
#'
#' A quantile is a smooth function of every linear predictor at once, so its
#' variance needs the predictors' \emph{joint} covariance, cross-parameter
#' terms included: \eqn{Var(q) = \sum_{k,l} (\partial q / \partial \eta_k)
#' (\partial q / \partial \eta_l) Cov(\eta_k, \eta_l)}.
#'
#' The derivatives are taken by central difference on the family's quantile
#' function, which avoids needing an analytic derivative for each of the
#' families that has one. The covariances come from the same linear forms the
#' per-parameter standard errors use, so `Cov(eta_k, eta_l)` is one row-wise
#' product per pair.
#'
#' @param object A `gamRTMB` fit.
#' @param prob Probabilities.
#' @param forms Per-parameter linear forms from [predict.gamRTMB()].
#' @param Vj Joint coefficient covariance from [.joint_cov()].
#' @param fx Resolved fixed arguments for the rows being predicted.
#' @param h Step for the central difference, on the link scale.
#' @return A matrix of standard errors, one column per probability.
#' @keywords internal
.quantile_se <- function(object, prob, forms, Vj, fx = object$fixed,
                         h = 1e-4) {
  fam <- object$family
  pn <- object$design$parnames
  K <- length(pn); n <- length(forms[[1L]]$eta)
  linkinv <- lapply(fam$links, function(l) .links[[l]]$linkinv)
  etas <- lapply(forms, `[[`, "eta")

  ## Cov(eta_k, eta_l), per observation
  cv <- lapply(seq_len(K), function(k) lapply(seq_len(K), function(l)
    rowSums((forms[[k]]$L %*% Vj$V[forms[[k]]$ii, forms[[l]]$ii, drop = FALSE]) *
              forms[[l]]$L)))

  qat <- function(e, pr) fam$qf(rep(pr, n),
    stats::setNames(lapply(seq_len(K), function(k) linkinv[[k]](e[[k]])), pn), fx)

  vapply(prob, function(pr) {
    g <- lapply(seq_len(K), function(k) {
      up <- dn <- etas
      up[[k]] <- up[[k]] + h; dn[[k]] <- dn[[k]] - h
      (qat(up, pr) - qat(dn, pr)) / (2 * h)
    })
    v <- numeric(n)
    for (k in seq_len(K)) for (l in seq_len(K))
      v <- v + g[[k]] * g[[l]] * cv[[k]][[l]]
    sqrt(pmax(v, 0))
  }, numeric(n))
}

#' @param object A `gamRTMB` fit.
#' @param ... Ignored.
#' @describeIn gamRTMB Coefficients, as a list of the fixed (`beta`) and
#'   penalized (`b`) vectors plus the log variance components. `beta` and
#'   `log_sigma` are named `parameter:term`.
#' @export
coef.gamRTMB <- function(object, ...) {
  b <- object$coefficients$beta
  names(b) <- .beta_labels(object$design)
  ls <- object$log_sigma
  names(ls) <- unlist(lapply(object$design$blocks, function(z)
    paste0(z$par, ":", z$label,
           if (z$ntheta > 1L) paste0(":", z$theta_names) else "")))
  list(beta = b, b = object$coefficients$b, log_sigma = ls)
}

#' @describeIn gamRTMB Fitted values of every distributional parameter, on the
#'   response scale.
#' @export
fitted.gamRTMB <- function(object, ...) {
  stats::predict(object, type = "response")
}

