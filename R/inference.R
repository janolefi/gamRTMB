## Everything downstream of the fit goes through .penalized_hessian(), so that
## a second fitting routine only has to supply that one accessor rather than
## have edf() and vcov() reach into its internals.

#' The penalized Hessian, in the design's own coefficient ordering
#'
#' The contract each fitting routine fills: return the Hessian of the joint
#' penalized negative log-likelihood in the coefficients at the fitted mode,
#' plus index vectors giving the row of `H` for each entry of `beta` and of
#' `b`. Anything needing curvature — [edf()], smooth-term covariances — uses
#' only this.
#'
#' @param fit A `gamRTMB` fit.
#' @return `list(H, i_beta, i_b)`, or `NULL` when the fit did not form it.
#' @keywords internal
.penalized_hessian <- function(fit) {
  ## "aREML" forms this matrix itself, in the parameter list's own ordering,
  ## `beta` then `b` -- which is the ordering this function promises.
  if (identical(fit$method, "aREML"))
    return(if (is.null(fit$H)) NULL else
      list(H = fit$H, i_beta = seq_len(fit$design$nbeta),
           i_b = fit$design$nbeta + seq_len(fit$design$nb)))
  if (is.null(fit$obj)) return(NULL)
  if (fit$method != "REML") return(NULL)   # beta is not in the random vector
  obj <- fit$obj
  H <- obj$env$spHess(obj$env$last.par.best, random = TRUE)
  nm <- names(obj$env$par[obj$env$random])
  list(H = H, i_beta = which(nm == "beta"), i_b = which(nm == "b"))
}

## Round-off allowance on an edf_j, which is mathematically confined to
## [0, 1]. Generous by the standards of the sparse solve that produces it and
## nowhere near the scale of a real violation: the fits that break this are
## out by thousands, not by parts in a million.
.edf_tol <- 1e-6

#' Is a sparse symmetric matrix positive definite?
#'
#' `Matrix::chol()` rather than `Matrix::Cholesky()`: the latter computes an
#' LDL' factorisation, which exists perfectly well for an indefinite matrix
#' and returns with nothing worse than a CHOLMOD warning, so it answers a
#' different question from the one being asked here. `chol()` fails on an
#' indefinite matrix, which is the answer wanted.
#'
#' @param h A symmetric sparse matrix, or `NULL`.
#' @return `TRUE`, `FALSE`, or `NA` if `h` is missing or has non-finite
#'   entries, so the question cannot be put.
#' @keywords internal
.is_pd <- function(h) {
  if (is.null(h)) return(NA)
  ## The stored values, not `as.numeric(h)`: the latter expands a sparse
  ## matrix to a dense vector, and this is asked once per outer iteration by
  ## method = "aREML".
  x <- if (methods::is(h, "sparseMatrix")) h@x else as.numeric(h)
  if (anyNA(x)) return(NA)
  !inherits(tryCatch(Matrix::chol(h), error = function(e) e,
                     warning = function(w) w), "condition")
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
#' @section When the EDF do not exist:
#' All of that assumes \eqn{H_{data}} is positive semi-definite, which is what
#' makes \eqn{F = H^{-1} H_{data}} a projection and puts every \eqn{edf_j} in
#' \eqn{[0, 1]}. At a point the optimiser never converged to it need not be,
#' and the solve still returns numbers: a four-parameter Box-Cox fit on
#' `film90` gives EDF near \eqn{-6000} for a rank-9 basis.
#'
#' Note that it is \eqn{H_{data}} and not \eqn{H} that has to be checked.
#' The penalty can and does rescue the sum: on that same fit \eqn{H} is
#' positive definite while \eqn{H_{data} = H - S} has two negative
#' eigenvalues, so a test on \eqn{H} passes and the EDF are still nonsense.
#' Rather than factorise a second matrix, the \eqn{edf_j} are checked against
#' the \eqn{[0, 1]} they are guaranteed to lie in -- the same statement, and
#' already computed.
#'
#' A negative EDF is not a small inaccuracy to report with a caveat; it means
#' the quantity does not exist at this point. So the column is `NA` instead,
#' with a warning pointing at `max_grad`. See [.inner_indefinite()] for how a
#' fit gets into that state.
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
    stop("effective degrees of freedom need the penalized Hessian over every ",
         "coefficient, and method = \"ML\" does not form it: the unpenalized ",
         "coefficients are fixed effects there rather than part of the random ",
         "vector. Use method = \"REML\" or \"aREML\".")
  D <- object$design
  ls <- object$log_sigma
  edf_all <- 1 - Matrix::diag(Matrix::solve(ph$H, .penalty_matrix(D, ls, ph)))
  ## Every edf_j lies in [0, 1] when the EDF exist at all; the tolerance is
  ## for the sparse solve's round-off, not for genuinely out-of-range values,
  ## which run to thousands rather than to fractions.
  if (any(edf_all < -.edf_tol | edf_all > 1 + .edf_tol)) {
    warning("the effective degrees of freedom are not defined at these ",
            "values: ", sum(edf_all < -.edf_tol | edf_all > 1 + .edf_tol),
            " of ", length(edf_all), " coefficients fall outside [0, 1], so ",
            "the data Hessian is not positive semi-definite here and ",
            "H^-1 H_data is not a projection. They are reported as NA. This ",
            "fit has not converged -- check `max_grad`; see ?gamRTMB for the ",
            "starting-value and basis options.", call. = FALSE)
    edf_all[] <- NA_real_
  }

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
#' uncertainty in the smoothing parameters (mgcv's `unconditional = TRUE`) --
#' under `method = "REML"` or `"ML"`. `"aREML"` has no joint precision and
#' returns the inverse penalized Hessian, which conditions on the smoothing
#' parameters instead; see the branch at the top of the function and
#' [vcov.gamRTMB()].
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
  ## "aREML" never builds an `sdreport`, and could not fill this in from one
  ## if it did: the joint precision's smoothing-parameter block comes from
  ## differentiating the marginal criterion twice, which is the term
  ## Fellner-Schall exists to avoid. What it has is the penalized Hessian,
  ## whose inverse is the covariance *conditional* on the fitted smoothing
  ## parameters -- mgcv's `unconditional = FALSE`. Intervals from it are a
  ## little too narrow for the same reason mgcv's conditional ones are.
  if (identical(fit$method, "aREML")) {
    ph <- .penalized_hessian(fit)
    if (is.null(ph))
      stop("this fit did not keep its penalized Hessian, so the coefficient ",
           "covariance is not available")
    return(list(V = as.matrix(Matrix::solve(ph$H)),
                ib = ph$i_beta, ir = ph$i_b))
  }
  jp <- fit$sdr$jointPrecision
  if (is.null(jp))
    stop("standard errors on smooth terms need a fit made with ",
         "joint_precision = TRUE")
  Q <- as.matrix(jp)
  ## A model with no smooths has nothing for the outer optimiser to estimate,
  ## and RTMB leaves the joint precision unnamed when the fixed-effect vector
  ## is empty. The block order is the parameter list's, so the names can be
  ## recovered from it; without this a purely parametric fit has no standard
  ## errors at all.
  nm <- colnames(Q)
  if (is.null(nm)) nm <- names(fit$obj$env$par)
  if (length(nm) != ncol(Q))
    stop("the joint precision does not match the parameter vector")
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
