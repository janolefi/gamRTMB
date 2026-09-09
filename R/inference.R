## Everything downstream of the fit goes through .penalized_hessian(), so that
## a second engine only has to supply that one accessor rather than have
## edf() and vcov() reach into its internals.

#' The penalized Hessian, in the design's own coefficient ordering
#'
#' The engine contract: return the Hessian of the joint penalized negative
#' log-likelihood in the coefficients at the fitted mode, plus index vectors
#' giving the row of `H` for each entry of `beta` and of `b`. Anything needing
#' curvature — [edf()], smooth-term covariances — uses only this.
#'
#' @param fit A `gamRTMB` fit.
#' @return `list(H, i_beta, i_b)`, or `NULL` if the engine cannot supply it.
#' @keywords internal
.penalized_hessian <- function(fit) {
  if (fit$engine != "laplace" || is.null(fit$obj)) return(NULL)
  if (fit$method != "REML") return(NULL)   # beta is not in the random vector
  obj <- fit$obj
  H <- obj$env$spHess(obj$env$last.par.best, random = TRUE)
  nm <- names(obj$env$par[obj$env$random])
  list(H = H, i_beta = which(nm == "beta"), i_b = which(nm == "b"))
}

#' Effective degrees of freedom per smooth
#'
#' TMB does not hand these over the way mgcv's PIRLS does, so they are derived
#' explicitly. At the fitted smoothing parameters the penalized Hessian of the
#' joint negative log-likelihood in the coefficients is
#' \deqn{H = H_{data} + S, \quad S = diag(0 \text{ for fixed}, 1/\sigma_k^2
#' \text{ for block } k),}
#' the penalty being exactly diagonal in the [mgcv::smooth2random()] basis.
#' Wood's effective degrees of freedom, \eqn{tr((X'WX + S)^{-1} X'WX)},
#' therefore generalise to
#' \deqn{F = H^{-1} H_{data} = I - H^{-1} S, \quad
#'       edf_j = 1 - s_j [H^{-1}]_{jj},}
#' so only the diagonal of \eqn{H^{-1}} is needed: null-space coefficients
#' contribute exactly 1 and penalized ones between 0 and 1.
#'
#' Checked against [mgcv::gaulss()], which agrees to three decimals on both
#' untied and `id`-tied models.
#'
#' @param object A `gamRTMB` fit, made with `method = "REML"`.
#' @param ... Ignored.
#' @return A data frame with one row per smooth: parameter, term label, EDF,
#'   basis dimension, smoothing parameter(s) and `id`. The total EDF over all
#'   coefficients is attached as attribute `"edf.total"`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(200)); d$y <- rnorm(200, sin(2 * pi * d$x), 0.3)
#' edf(gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d))
#' @export
edf <- function(object, ...) UseMethod("edf")

#' @rdname edf
#' @export
edf.gamRTMB <- function(object, ...) {
  ph <- .penalized_hessian(object)
  if (is.null(ph))
    stop("effective degrees of freedom need the coefficients in the random ",
         "vector, i.e. method = \"REML\" with engine = \"laplace\"")
  D <- object$design
  dH <- diag(solve(as.matrix(ph$H)))
  ls <- object$log_sigma

  sj <- numeric(length(dH))                       # the diagonal penalty
  for (k in seq_along(D$blocks))
    sj[ph$i_b[D$blocks[[k]]$idx]] <- exp(-2 * ls[k])
  edf_all <- 1 - sj * dH

  rows <- list()
  for (p in D$parnames) for (j in seq_along(D$parts[[p]]$smooths)) {
    s <- D$parts[[p]]$smooths[[j]]
    idx <- .smooth_idx(D, p, j)
    ii <- c(ph$i_b[idx$b], ph$i_beta[idx$f])
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = p, term = s$label, edf = sum(edf_all[ii]),
      k = ncol(s$sm$X),
      sp = paste(sprintf("%.4g", exp(-2 * ls[s$block_ids])), collapse = ","),
      id = if (is.na(s$id)) "" else s$id, row.names = NULL)
  }
  res <- if (length(rows)) do.call(rbind, rows) else
    data.frame(parameter = character(), term = character(), edf = numeric())
  attr(res, "edf.total") <- sum(edf_all)
  res
}

#' Where one smooth's coefficients live in the joint vectors
#'
#' Four places need the same thing: which entries of `b` and which entries of
#' `beta` belong to smooth `j` of parameter `p` (its penalized blocks, then
#' its null-space columns). Gathering it once keeps [edf()], the coefficient
#' reconstruction and both covariance helpers in step.
#'
#' @return `list(b, f)`, indices into the `b` and `beta` vectors.
#' @keywords internal
.smooth_idx <- function(design, p, j) {
  s <- design$parts[[p]]$smooths[[j]]
  list(b = unlist(lapply(s$block_ids, function(k) design$blocks[[k]]$idx)),
       f = if (length(s$f_local)) design$beta_idx[[p]][s$f_local] else integer(0))
}

#' Joint covariance of all coefficients
#'
#' From [RTMB::sdreport()]'s joint precision, which includes the uncertainty
#' in the smoothing parameters. Requires `joint_precision = TRUE` at fit time.
#'
#' @keywords internal
.joint_cov <- function(fit) {
  jp <- fit$sdr$jointPrecision
  if (is.null(jp))
    stop("standard errors on smooth terms need a fit made with ",
         "joint_precision = TRUE")
  nm <- colnames(jp)
  list(V = solve(as.matrix(jp)), ib = which(nm == "beta"), ir = which(nm == "b"))
}

#' Covariance of one smooth's coefficients, in its own basis
#'
#' The same map that reconstructs the coefficients propagates their
#' covariance: `Tmap %*% V %*% t(Tmap)`.
#'
#' @keywords internal
.smooth_vcov <- function(fit, p, j, Vj) {
  Tm <- fit$design$parts[[p]]$smooths[[j]]$Tmap
  idx <- .smooth_idx(fit$design, p, j)
  ii <- c(Vj$ir[idx$b], Vj$ib[idx$f])
  Tm %*% Vj$V[ii, ii, drop = FALSE] %*% t(Tm)
}

#' Covariance of a whole linear predictor's coefficients
#'
#' Ordered as `(parametric, smooth 1, smooth 2, ...)` in the original bases,
#' matching how [predict.gamRTMB()] stacks the model matrices.
#'
#' @keywords internal
.eta_vcov <- function(fit, p, Vj) {
  D <- fit$design; P <- D$parts[[p]]
  npara <- ncol(P$Xpara)
  cols <- Vj$ib[D$beta_idx[[p]][seq_len(npara)]]
  Tlist <- list(diag(1, npara))
  for (j in seq_along(P$smooths)) {
    idx <- .smooth_idx(D, p, j)
    cols <- c(cols, Vj$ir[idx$b], Vj$ib[idx$f])
    Tlist[[length(Tlist) + 1L]] <- P$smooths[[j]]$Tmap
  }
  Tb <- as.matrix(Matrix::bdiag(Tlist))
  Tb %*% Vj$V[cols, cols, drop = FALSE] %*% t(Tb)
}

#' Coefficients of one smooth, in its own (constrained) basis
#'
#' @keywords internal
.smooth_beta <- function(fit, p, j) {
  idx <- .smooth_idx(fit$design, p, j)
  cf <- c(fit$coefficients$b[idx$b], fit$coefficients$beta[idx$f])
  as.vector(fit$design$parts[[p]]$smooths[[j]]$Tmap %*% cf)
}
