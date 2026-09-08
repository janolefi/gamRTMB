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
  cat("  coefficients: ", D$nbeta, " fixed (incl. null spaces), ", D$nb,
      " penalized; ", D$nsigma, " smoothing parameter",
      if (D$nsigma != 1L) "s" else "",
      if (D$nsigma_free < D$nsigma)
        paste0(" (", D$nsigma_free, " free, ", D$nsigma - D$nsigma_free,
               " tied by id)") else "",
      "\n", sep = "")
  invisible(x)
}

#' Summarise a gamRTMB fit
#'
#' @param object A `gamRTMB` fit.
#' @param ... Ignored.
#' @return An object of class `summary.gamRTMB`.
#' @export
summary.gamRTMB <- function(object, ...) {
  e <- tryCatch(edf(object), error = function(err) NULL)
  para <- lapply(object$design$parnames, function(p) {
    ii <- object$design$beta_idx[[p]]
    lab <- attr(ii, "labels")
    keep <- !grepl("\\.null[0-9]+$", lab)
    data.frame(parameter = p, coefficient = lab[keep],
               estimate = object$coefficients$beta[ii][keep], row.names = NULL)
  })
  structure(list(family = object$family, method = object$method,
                 engine = object$engine, objective = object$objective,
                 convergence = object$convergence, n = object$design$n,
                 parametric = do.call(rbind, para), edf = e),
            class = "summary.gamRTMB")
}

#' @param x A `summary.gamRTMB`.
#' @rdname summary.gamRTMB
#' @export
print.summary.gamRTMB <- function(x, ...) {
  cat("gamRTMB: ", x$family$family, " | ", x$method, " | n = ", x$n, "\n", sep = "")
  cat("parameters: ",
      paste(sprintf("%s (%s)", x$family$parnames, x$family$links),
            collapse = ", "), "\n\n", sep = "")
  cat("Parametric coefficients (link scale):\n")
  print(x$parametric, row.names = FALSE, digits = 4)
  if (!is.null(x$edf)) {
    cat("\nSmooth terms:\n")
    print(x$edf, row.names = FALSE, digits = 4)
    cat("\ntotal EDF: ", sprintf("%.2f", attr(x$edf, "edf.total")),
        "   -logLik: ", sprintf("%.4f", x$objective), "\n", sep = "")
  } else cat("\n-logLik: ", sprintf("%.4f", x$objective), "\n", sep = "")
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
#'   each smooth separately.
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
    Xpara <- if (is.null(newdata)) P$Xpara else
      stats::model.matrix(stats::delete.response(P$terms), newdata, xlev = P$xlev)
    npara <- ncol(P$Xpara)
    bfull <- object$coefficients$beta[D$beta_idx[[p]]]
    eta <- as.vector(Xpara %*% bfull[seq_len(npara)])

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

#' @describeIn gamRTMB The value of the fitting criterion. Note that under
#'   `method = "REML"` this is the restricted criterion, not the likelihood,
#'   so it is not comparable across different fixed-effect structures.
#' @export
logLik.gamRTMB <- function(object, ...) {
  structure(-object$objective,
            df = object$design$nsigma_free,
            nobs = object$design$n,
            criterion = object$method,
            class = "logLik")
}
