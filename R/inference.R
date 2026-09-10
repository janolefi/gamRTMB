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
#' with \eqn{S} block diagonal -- \eqn{\sigma_k^{-2} I} in the
#' [mgcv::smooth2random()] basis, \eqn{\sigma_k^{-2} Q_k} for a block that
#' kept its own sparse penalty. Wood's effective degrees of freedom,
#' \eqn{tr((X'WX + S)^{-1} X'WX)}, therefore generalise to
#' \deqn{F = H^{-1} H_{data} = I - H^{-1} S, \quad
#'       edf_j = 1 - [H^{-1} S]_{jj},}
#' so only the diagonal of \eqn{H^{-1} S} is needed, and that comes from a
#' sparse solve rather than a full inverse. Null-space coefficients contribute
#' exactly 1 and penalized ones between 0 and 1.
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
  ls <- object$log_sigma
  edf_all <- 1 - Matrix::diag(Matrix::solve(ph$H, .penalty_matrix(D, ls, ph)))

  rows <- list()
  for (p in D$parnames) for (j in seq_along(D$parts[[p]]$smooths)) {
    s <- D$parts[[p]]$smooths[[j]]
    idx <- .smooth_idx(D, p, j)
    ii <- c(ph$i_b[idx$b], ph$i_beta[idx$f])
    rows[[length(rows) + 1L]] <- data.frame(
      parameter = p, term = s$label, edf = sum(edf_all[ii]),
      k = ncol(s$sm$X),
      sp = paste(unlist(lapply(D$blocks[s$block_ids], function(bl)
        if (identical(bl$theta_names, "sd"))
          sprintf("%.4g", exp(-2 * ls[bl$theta_idx]))
        else sprintf("%s=%.4g", bl$theta_names, exp(ls[bl$theta_idx])))),
        collapse = ","),
      id = if (is.na(s$id)) "" else s$id, row.names = NULL)
  }
  res <- if (length(rows)) do.call(rbind, rows) else
    data.frame(parameter = character(), term = character(), edf = numeric())
  attr(res, "edf.total") <- sum(edf_all)
  res
}

#' The penalty, assembled in the Hessian's own ordering
#'
#' \eqn{S} is block diagonal with one block per penalized block of
#' coefficients: \eqn{\sigma_k^{-2} I} for a block that went through
#' [mgcv::smooth2random()], and \eqn{\sigma_k^{-2} Q_k} for one that kept a
#' sparse penalty. Built sparse so that `solve(H, S)` stays a sparse solve
#' rather than a full inverse.
#'
#' @keywords internal
.penalty_matrix <- function(design, ls, ph) {
  np <- nrow(ph$H)
  ii <- jj <- integer(0); xx <- numeric(0)
  for (bl in design$blocks) {
    r <- ph$i_b[bl$idx]
    if (identical(bl$kind, "iid")) {
      ii <- c(ii, r); jj <- c(jj, r)
      xx <- c(xx, rep(exp(-2 * ls[bl$theta_idx]), length(r)))
    } else {
      ## A precision is stored symmetric, and a symmetric sparse matrix keeps
      ## only one triangle; going through "generalMatrix" first gets both.
      Q <- .block_prec(bl, ls[bl$theta_idx])
      tq <- Matrix::summary(as(as(Q, "generalMatrix"), "TsparseMatrix"))
      ii <- c(ii, r[tq$i]); jj <- c(jj, r[tq$j]); xx <- c(xx, tq$x)
    }
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(np, np))
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

#' Joint covariance of the coefficients
#'
#' The `(beta, b)` block of the inverse joint precision, which includes the
#' uncertainty in the smoothing parameters (mgcv's `unconditional = TRUE`).
#'
#' Computed as a Schur complement rather than by inverting the whole matrix.
#' Writing the joint precision over coefficients `c` and log smoothing
#' parameters `s` as `[[Qcc, Qcs], [Qsc, Qss]]`, the block needed is
#' \deqn{[Q^{-1}]_{cc} = (Q_{cc} - Q_{cs} Q_{ss}^{-1} Q_{sc})^{-1},}
#' which is algebraically identical to inverting the whole thing but isolates
#' the awkward part.
#'
#' Boundary smoothing parameters make this delicate in two separate ways, and
#' both are handled here because a plain `solve()` of the whole matrix fails
#' on either, taking the standard errors, bands and summary with it.
#'
#' \strong{A flat smoothing parameter.} A term shrunk onto its null space
#' leaves the criterion flat in its own `log_sigma`, so the joint precision is
#' genuinely singular — but only in `Qss`: with two such terms the offending
#' eigenvalues were 3e-12 and 2e-07, loading on `log_sigma` with weight 1.00,
#' while the coefficient block stayed invertible. A pseudo-inverse of `Qss`
#' drops exactly those directions, which amounts to treating a boundary
#' smoothing parameter as known rather than estimated — the conditional
#' treatment, for that parameter only. Every other one still contributes.
#'
#' \strong{Scaling.} The same boundary puts precision entries of order 1e19
#' next to entries of order 1, and at that dynamic range double precision
#' loses positive definiteness outright: a factor-smooth model measured an
#' eigenvalue of -7.5e3 in a matrix whose largest was 2e19. So the matrix is
#' first scaled to a unit diagonal, which leaves only the correlation
#' structure to invert, and the result is unscaled afterwards.
#'
#' @param fit A `gamRTMB` fit.
#' @return `list(V, ib, ir)`: the coefficient covariance, and the positions of
#'   the `beta` and `b` entries within it.
#' @keywords internal
.joint_cov <- function(fit) {
  jp <- fit$sdr$jointPrecision
  if (is.null(jp))
    stop("standard errors on smooth terms need a fit made with ",
         "joint_precision = TRUE")
  Q <- as.matrix(jp)
  nm <- colnames(Q)
  ## guard before the square root: at this dynamic range the diagonal itself
  ## can come back negative
  dg <- diag(Q); dg[!is.finite(dg) | dg <= 0] <- 1
  sc <- sqrt(dg)
  Q <- Q / outer(sc, sc)                          # unit diagonal
  ic <- which(nm %in% c("beta", "b"))
  is <- which(nm == "log_sigma")
  Sc <- Q[ic, ic, drop = FALSE]
  if (length(is)) {
    Qcs <- Q[ic, is, drop = FALSE]
    e <- eigen(Q[is, is, drop = FALSE], symmetric = TRUE)
    keep <- e$values > max(e$values, 0) * 1e-10
    if (any(keep)) {
      U <- e$vectors[, keep, drop = FALSE]
      Sc <- Sc - Qcs %*% (U %*% (t(U) / e$values[keep])) %*% t(Qcs)
    }
  }
  sub <- nm[ic]
  list(V = solve(Sc) / outer(sc[ic], sc[ic]),      # and unscaled again
       ib = which(sub == "beta"), ir = which(sub == "b"))
}

#' One smooth's design in joint-coefficient space
#'
#' A smooth's fitted values and their standard errors are both linear forms in
#' the joint coefficient vector once its model matrix has been mapped through
#' the reparameterisation: with `Z = X Tmap`, the contribution is `Z c` and its
#' covariance `Z V Z'`. Mapping the design once is simpler than mapping the
#' coefficients and their covariance separately, and it removes the need to
#' assemble a block-diagonal transform for the whole linear predictor.
#'
#' @param fit A `gamRTMB` fit.
#' @param p,j Distributional parameter and smooth index.
#' @param X The smooth's model matrix, from the fit or from
#'   [mgcv::PredictMat()].
#' @return `list(Z, coef, b, f)`: the mapped design, the coefficients it
#'   multiplies, and the `b` and `beta` indices for locating them in a
#'   covariance matrix.
#' @keywords internal
.smooth_part <- function(fit, p, j, X) {
  idx <- .smooth_idx(fit$design, p, j)
  list(Z = X %*% fit$design$parts[[p]]$smooths[[j]]$Tmap,
       coef = c(fit$coefficients$b[idx$b], fit$coefficients$beta[idx$f]),
       b = idx$b, f = idx$f)
}

#' Standard errors of a linear form in the coefficients
#' @keywords internal
.qform_se <- function(Z, V) sqrt(pmax(rowSums((Z %*% V) * Z), 0))

#' Randomised quantile (pseudo) residuals
#'
#' Residuals by the probability integral transform: if the fitted
#' distribution is right, \eqn{u_i = F(y_i; \hat\theta_i)} is uniform and
#' \eqn{r_i = \Phi^{-1}(u_i)} is standard normal, so a QQ plot of `r` checks
#' the whole distributional assumption rather than just the mean.
#'
#' @section Discrete and mixed responses:
#' Where the response has atoms, \eqn{F} is a step function and \eqn{u} cannot
#' be uniform, so the residual is randomised within the step
#' (Dunn & Smyth 1996):
#' \deqn{u_i = F(y_i^-) + v_i (F(y_i) - F(y_i^-)), \quad v_i \sim U(0,1).}
#' The left limit \eqn{F(y^-)} comes from the family's declared support, since
#' it cannot be obtained reliably any other way:
#' \describe{
#'   \item{continuous}{\eqn{F(y^-) = F(y)} and no randomisation happens.}
#'   \item{lattice}{\eqn{F(y^-) = F(y-1)}, evaluated at integer arguments
#'     only, and taken as 0 at \eqn{y = 0} rather than evaluating the CDF
#'     below its support.}
#'   \item{mixed}{\eqn{F(y^-) = F(y) - p(y)} at an atom, where the density
#'     returns the atom's mass, and \eqn{F(y)} elsewhere.}
#' }
#' Nudging the argument instead (\eqn{F(y-\delta)}) is not safe: CDF
#' implementations disagree about whether they floor a non-integer argument,
#' and some reject one outright.
#'
#' @section Randomisation:
#' With `randomise = TRUE` (the default, and what you want for a QQ plot) the
#' residuals are not a deterministic function of the fit; set a seed for
#' reproducibility, or inspect several draws. With `randomise = FALSE` a
#' discrete or mixed response returns the bounding interval per observation
#' instead of a point.
#'
#' Prior weights are ignored: a residual belongs to a row of the data, not to
#' the replicate count a weight stands for.
#'
#' @param object A `gamRTMB` fit.
#' @param type `"quantile"` for normal-scale residuals (the default), or
#'   `"uniform"` for the PIT values themselves.
#' @param randomise Randomise within the step for a discrete or mixed
#'   response. Ignored for a continuous one.
#' @param ... Ignored.
#' @return A numeric vector, or with `randomise = FALSE` and a non-continuous
#'   response a data frame of `lower`, `upper` and `mid`. Values are clamped
#'   away from 0 and 1 so that extreme observations stay finite and visible.
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996) Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics} 5, 236--244.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(300))
#' d$y <- rpois(300, exp(1 + sin(2 * pi * d$x)))
#' fit <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8)), family = fam("pois"),
#'                data = d)
#' r <- residuals(fit)
#' qqnorm(r); qqline(r)
#' @export
residuals.gamRTMB <- function(object, type = c("quantile", "uniform"),
                              randomise = TRUE, ...) {
  type <- match.arg(type)
  fam <- object$family
  if (is.null(fam$cdf))
    stop("family '", fam$family, "' has no CDF in ", fam$source,
         ", so quantile residuals are not available. families() reports which ",
         "families support them.")
  y <- object$y
  theta <- stats::predict(object, type = "response")
  fx <- object$fixed
  hi <- fam$cdf(y, theta, fx)

  lo <- switch(fam$support,
    continuous = hi,
    ## integer arguments only, and 0 at the bottom of the support, so the CDF
    ## is never asked for a value below it
    lattice = {
      l <- numeric(length(y)); pos <- y > 0
      if (any(pos))
        l[pos] <- fam$cdf(y[pos] - 1, .at(theta, pos), .at(fx, pos))
      l
    },
    ## at an atom the density returns its mass, not a density
    mixed = {
      l <- hi; at <- y %in% fam$atoms
      if (any(at))
        l[at] <- hi[at] - exp(fam$logdens(y[at], .at(theta, at), .at(fx, at)))
      l
    })
  hi <- pmin(pmax(hi, 0), 1)
  lo <- pmin(pmax(lo, 0), hi)

  if (fam$support == "continuous") return(.pit(hi, type))
  if (randomise) return(.pit(stats::runif(length(y), lo, hi), type))
  data.frame(lower = .pit(lo, type), upper = .pit(hi, type),
             mid = .pit((lo + hi) / 2, type))
}

#' Subset the per-observation entries of a parameter list
#' @keywords internal
.at <- function(l, i) lapply(l, function(z) if (length(z) == 1L) z else z[i])

#' Put PIT values on the requested scale, keeping them finite
#'
#' `qnorm()` at exactly 0 or 1 is infinite, which would drop the very
#' observations a diagnostic most wants to show; clamping keeps them finite
#' and visibly extreme instead.
#'
#' @keywords internal
.pit <- function(u, type) {
  eps <- .Machine$double.eps
  u <- pmin(pmax(u, eps), 1 - eps)
  if (type == "uniform") u else stats::qnorm(u)
}
