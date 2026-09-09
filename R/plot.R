## Plotting. One entry point, three workers: term plots, a QQ plot and a worm
## plot. Base graphics throughout, no colour unless something is being
## overlaid, and `bty = "n"` as a default the caller can still override.

.or_else <- function(a, b) if (is.null(a)) b else a

#' Panel grid for n plots
#' @keywords internal
.panel_grid <- function(n) { nc <- ceiling(sqrt(n)); c(ceiling(n / nc), nc) }

#' Colour-blind safe palette
#'
#' Okabe-Ito, in base R since 4.0, so no dependency. Used only where series
#' are overlaid; a single-series plot needs no colour at all. The pale yellow
#' is dropped: it is fine for fills but too low-contrast for a line on white.
#'
#' @keywords internal
.pal <- function(n) {
  p <- grDevices::palette.colors(8, "Okabe-Ito")[c(1, 2, 3, 4, 6, 7, 8)]
  rep_len(p, max(n, 1L))[seq_len(max(n, 1L))]
}

#' Call a graphics function with our defaults under the caller's `...`
#'
#' Passing both a default and the caller's value for the same argument to
#' `plot()` is an error ("matched by multiple actual arguments"), so they are
#' merged first and the caller wins.
#'
#' @keywords internal
.gcall <- function(fun, defaults, dots) do.call(fun, utils::modifyList(defaults, dots))

#' A light polygon for a pointwise interval
#' @keywords internal
.band <- function(x, lo, hi, col) {
  ok <- is.finite(x) & is.finite(lo) & is.finite(hi)
  if (any(ok))
    graphics::polygon(c(x[ok], rev(x[ok])), c(hi[ok], rev(lo[ok])),
                      col = col, border = NA)
}

#' Plot a gamRTMB fit
#'
#' @section Term plots (`type = "terms"`):
#' One panel per smooth, showing its contribution to that distributional
#' parameter's linear predictor — the scale on which terms are additive — with
#' a **pointwise** \eqn{\pm 2} standard error band. The band comes from
#' [vcov.gamRTMB()]'s joint covariance, so it includes the uncertainty in the
#' smoothing parameters (mgcv's `unconditional = TRUE`); it needs a fit made
#' with `joint_precision = TRUE`, which is the default. Panels are titled
#' `parameter: term`, since terms belong to different parameters.
#'
#' Only one-dimensional smooths of a numeric covariate are drawn. Tensor
#' products, random effects (`bs = "re"`) and factor-smooth interactions
#' (`bs = "fs"`) are reported and skipped rather than drawn misleadingly.
#'
#' @section Quantile plots (`type = "quantile"`):
#' Fitted quantiles of the response against one covariate, over the observed
#' data. This is the picture that letting every parameter vary is *for*: the
#' curves fan out and contract as the fitted spread and shape change, which no
#' mean-only model can show. Other covariates are held at a typical value —
#' the median for a numeric one, the modal level for a factor — so with more
#' than one covariate the points and the curves do not condition on quite the
#' same thing.
#'
#' A band is drawn for the central curve only, and only for a continuous
#' response. Bands on every quantile would be unreadable, and it is the middle
#' of the distribution whose estimation uncertainty one usually wants next to
#' the spread the other curves already show. Note that this band is the
#' central curve's own uncertainty and is not a prediction interval — the
#' outer quantiles are that.
#'
#' @section Diagnostics (`type = "qq"`, `type = "worm"`):
#' Both use the randomised quantile residuals of [residuals.gamRTMB()]. The QQ
#' plot references the identity line, because these residuals should be
#' standard normal rather than merely normal. The worm plot is its detrended
#' version — deviation from the theoretical quantile against that quantile —
#' which makes it much easier to see *where* a distribution is wrong, with the
#' usual pointwise 95% band.
#'
#' For a discrete or mixed response the residuals are randomised, so a single
#' panel is one realisation; `nsim > 1` overlays draws so that the
#' randomisation is visible instead of hidden. For a continuous response the
#' draws are identical and `nsim` has no effect.
#'
#' @section Graphical arguments:
#' `...` is passed to the underlying [plot()] call, so `main`, `xlim`, `ylim`,
#' `cex` and friends work as usual; `col`, `lwd` and `pch` are applied to the
#' lines and points that are drawn on top. `bty = "n"` is the default and can
#' be overridden like any other argument.
#'
#' @param x A `gamRTMB` fit.
#' @param type `"terms"` (default), `"qq"` or `"worm"`.
#' @param select Which smooths to draw: an integer index, or a pattern matched
#'   against the `parameter: term` labels (e.g. `"sd"` or `"s(x1)"`). Defaults
#'   to all of them.
#' @param se Draw the interval band.
#' @param rug Add a rug of the observed covariate values.
#' @param band.col Fill for the interval band. A solid light grey by default,
#'   which renders the same on every device.
#' @param prob Probabilities for `type = "quantile"`.
#' @param xvar Covariate for the x-axis of a quantile plot. Defaults to the
#'   covariate of the first smooth.
#' @param nsim Randomisation draws to overlay in a QQ or worm plot.
#' @param ask Draw one panel per page, waiting between them, instead of
#'   fitting everything onto one page.
#' @param n Grid resolution for term curves.
#' @param ... Passed to [plot()].
#' @return Invisibly, a list of the plotted data: one data frame per term for
#'   `type = "terms"`, or the residual quantiles for the diagnostics, so any
#'   panel can be rebuilt by hand.
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = runif(300), x2 = runif(300))
#' d$y <- rnorm(300, sin(2 * pi * d$x1), exp(-1 + d$x2))
#' fit <- gamRTMB(y ~ list(mean = ~ s(x1), sd = ~ s(x2)), data = d)
#' plot(fit)
#' plot(fit, select = "sd")
#' plot(fit, type = "worm")
#' @export
plot.gamRTMB <- function(x, type = c("terms", "qq", "worm", "quantile"),
                         select = NULL, prob = c(0.05, 0.25, 0.5, 0.75, 0.95),
                         xvar = NULL, se = TRUE, rug = TRUE,
                         band.col = "grey85", nsim = 1, ask = FALSE,
                         n = 200, ...) {
  type <- match.arg(type)
  switch(type,
    terms    = .plot_terms(x, select, se, rug, band.col, ask, n, ...),
    qq       = .plot_res(x, nsim, FALSE, se, band.col, ...),
    worm     = .plot_res(x, nsim, TRUE,  se, band.col, ...),
    quantile = .plot_quantile(x, xvar, prob, se, band.col, n, ...))
}

#' The emptiest corner, for a legend
#'
#' Counts what is drawn in each corner of the plotting region and returns the
#' least crowded, so a five-curve legend does not land on top of the curves.
#' `xs` and `ys` must be matched coordinates of everything drawn.
#'
#' @keywords internal
.empty_corner <- function(xs, ys, yl) {
  fx <- (xs - min(xs)) / max(diff(range(xs)), 1e-12)
  fy <- (ys - yl[1L]) / max(diff(yl), 1e-12)
  n <- c(topleft = sum(fx < .35 & fy > .65), topright = sum(fx > .65 & fy > .65),
         bottomleft = sum(fx < .35 & fy < .35),
         bottomright = sum(fx > .65 & fy < .35))
  names(n)[which.min(n)]
}

#' Fitted quantiles against one covariate
#' @keywords internal
.plot_quantile <- function(x, xvar, prob, se, band.col, ngrid, ...) {
  D <- x$design
  if (is.null(xvar)) {
    v <- unlist(lapply(D$parnames, function(p)
      lapply(D$parts[[p]]$smooths, function(s) s$plot1d$var)))
    if (!length(v))
      stop("no smooth covariate to plot against; name one with xvar =")
    xvar <- v[[1L]]
  }
  if (!xvar %in% names(x$data))
    stop("`xvar` must name a column of the model data: ",
         paste(names(x$data), collapse = ", "))
  xv <- x$data[[xvar]]
  if (!is.numeric(xv)) stop("`xvar` must be a numeric covariate")

  ## other covariates at a typical value, so one curve per probability
  nd <- x$data[rep(1L, ngrid), , drop = FALSE]
  for (v in setdiff(names(nd), xvar)) {
    cl <- x$data[[v]]
    nd[[v]] <- if (is.factor(cl))
      factor(rep(names(which.max(table(cl))), ngrid), levels = levels(cl))
      else rep(stats::median(cl), ngrid)
  }
  nd[[xvar]] <- seq(min(xv), max(xv), length.out = ngrid)

  band <- se && x$family$support == "continuous"
  qq <- stats::predict(x, newdata = nd, type = "quantile", prob = prob,
                       se.fit = band)
  Q <- if (band) qq$fit else qq
  S <- if (band) qq$se.fit else NULL
  mid <- which.min(abs(prob - 0.5))

  dots <- list(...)
  pch <- .or_else(dots$pch, 20); cx <- .or_else(dots$cex, 0.5)
  dots$pch <- dots$cex <- NULL
  ## Quantiles come in symmetric pairs, so code them by distance from the
  ## median: same colour and line type for a pair, black and thickest for the
  ## middle. That stays readable in greyscale as well as in colour.
  ## rounded, or 0.05 and 0.95 fail to pair up in floating point
  dist <- round(abs(prob - 0.5), 8)
  lvl <- match(dist, sort(unique(dist)))
  cols <- .or_else(dots$col, .pal(max(lvl)))[lvl]
  dots$col <- NULL
  ## step curves for a lattice response, where quantiles really are steps
  ltype <- if (x$family$support == "continuous") "l" else "s"
  yl <- range(c(x$y, Q, if (band) c(Q[, mid] - 2 * S[, mid],
                                    Q[, mid] + 2 * S[, mid])), finite = TRUE)

  .gcall(graphics::plot,
         list(x = xv, y = x$y, type = "n", bty = "n", ylim = yl,
              xlab = xvar, ylab = deparse(x$formula[[2L]]),
              main = paste0("fitted quantiles: ", x$family$family)), dots)
  if (band) .band(nd[[xvar]], Q[, mid] - 2 * S[, mid], Q[, mid] + 2 * S[, mid],
                  band.col)
  graphics::points(xv, x$y, pch = pch, cex = cx,
                   col = grDevices::adjustcolor("black", 0.4))
  for (k in seq_along(prob))
    graphics::lines(nd[[xvar]], Q[, k], type = ltype, col = cols[k],
                    lty = lvl[k], lwd = if (k == mid) 3 else 2)
  graphics::legend(.empty_corner(c(xv, rep(nd[[xvar]], ncol(Q))),
                                 c(x$y, as.vector(Q)), yl),
                   legend = paste0(100 * prob, "%"), bty = "n", col = cols,
                   lty = lvl, lwd = 2, cex = 0.85)
  out <- data.frame(nd[[xvar]], Q)
  names(out) <- c(xvar, colnames(Q))
  if (band) out$se.mid <- S[, mid]
  invisible(out)
}

#' @keywords internal
.plot_terms <- function(x, select, se, rug, band.col, ask, ngrid, ...) {
  D <- x$design
  tl <- list()
  for (p in D$parnames) for (j in seq_along(D$parts[[p]]$smooths)) {
    s <- D$parts[[p]]$smooths[[j]]
    tl[[length(tl) + 1L]] <- list(par = p, j = j, s = s,
                                  label = paste0(p, ": ", s$label))
  }
  if (!length(tl)) stop("the model has no smooth terms to plot")

  if (!is.null(select)) {
    keep <- if (is.numeric(select)) {
      if (any(select < 1 | select > length(tl)))
        stop("`select` must index the ", length(tl), " smooth terms")
      as.integer(select)
    } else grep(select, vapply(tl, `[[`, "", "label"), fixed = FALSE)
    if (!length(keep))
      stop("`select` matched none of: ",
           paste(vapply(tl, `[[`, "", "label"), collapse = ", "))
    tl <- tl[keep]
  }

  drawable <- !vapply(tl, function(t) is.null(t$s$plot1d), TRUE)
  if (any(!drawable))
    message("not drawable as a single curve, skipped: ",
            paste(vapply(tl[!drawable], `[[`, "", "label"), collapse = ", "))
  tl <- tl[drawable]
  if (!length(tl)) stop("no one-dimensional smooth terms left to plot")

  Vj <- if (se) tryCatch(.joint_cov(x), error = function(e) NULL) else NULL
  if (se && is.null(Vj))
    message("no interval band: the joint covariance is unavailable ",
            "(refit with joint_precision = TRUE)")

  ## build every curve first, so a failure cannot leave a half-drawn page
  curves <- lapply(tl, function(t) {
    ps <- t$s$plot1d
    xv <- x$data[[ps$var]]
    xg <- seq(min(xv), max(xv), length.out = ngrid)
    nd <- stats::setNames(list(xg), ps$var)
    if (!is.null(ps$by)) nd[[ps$by]] <- rep(ps$by_val, ngrid)
    Xs <- mgcv::PredictMat(t$s$sm, as.data.frame(nd, stringsAsFactors = FALSE))
    sp <- .smooth_part(x, t$par, t$j, Xs)
    sef <- if (is.null(Vj)) rep(NA_real_, ngrid) else {
      ii <- c(Vj$ir[sp$b], Vj$ib[sp$f])
      .qform_se(sp$Z, Vj$V[ii, ii, drop = FALSE])
    }
    data.frame(x = xg, fit = as.vector(sp$Z %*% sp$coef), se = sef)
  })
  names(curves) <- vapply(tl, `[[`, "", "label")

  ## one page by default; `ask` steps through instead. Either way the
  ## caller's device settings are put back.
  if (ask) {
    old <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(old), add = TRUE)
  } else {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = .panel_grid(length(tl)), mar = c(4, 4, 2, 1))
  }

  dots <- list(...)
  col <- .or_else(dots$col, "black"); lwd <- .or_else(dots$lwd, 2)
  dots$col <- dots$lwd <- NULL
  fill <- if (identical(col, "black")) band.col
          else grDevices::adjustcolor(col, alpha.f = 0.2)

  for (i in seq_along(tl)) {
    cv <- curves[[i]]; ps <- tl[[i]]$s$plot1d
    hi <- cv$fit + 2 * cv$se; lo <- cv$fit - 2 * cv$se
    yl <- if (all(is.na(cv$se))) range(cv$fit) else range(c(lo, hi), finite = TRUE)
    .gcall(graphics::plot,
           list(x = cv$x, y = cv$fit, type = "n", bty = "n", ylim = yl,
                xlab = ps$var, ylab = tl[[i]]$s$label, main = names(curves)[i]),
           dots)
    if (!all(is.na(cv$se))) .band(cv$x, lo, hi, fill)
    graphics::lines(cv$x, cv$fit, col = col, lwd = lwd)
    if (rug) graphics::rug(x$data[[ps$var]], col = grDevices::adjustcolor(col, 0.4))
  }
  invisible(curves)
}

#' QQ and worm plots
#'
#' A worm plot is the detrended QQ plot, so both are the same picture with the
#' theoretical quantile subtracted or not. The reference is the identity line
#' rather than a fitted one, because these residuals should be standard
#' normal and not merely normal.
#'
#' @keywords internal
.plot_res <- function(x, nsim, detrend, se, band.col, ...) {
  ## one column per randomisation draw; identical columns if continuous
  R <- vapply(seq_len(max(1L, as.integer(nsim))),
              function(i) sort(as.numeric(stats::residuals(x))),
              numeric(x$design$n))
  nn <- nrow(R)
  p <- stats::ppoints(nn); th <- stats::qnorm(p)
  Y <- if (detrend) R - th else R
  ## pointwise 95% band for the deviation of a normal order statistic
  bw <- 1.96 * sqrt(p * (1 - p) / nn) / stats::dnorm(th)

  dots <- list(...)
  pch <- .or_else(dots$pch, 20); cx <- .or_else(dots$cex, 0.5)
  cols <- if (ncol(R) > 1L) .pal(ncol(R)) else .or_else(dots$col, "grey30")
  dots$pch <- dots$cex <- dots$col <- NULL
  yl <- range(c(Y, if (detrend && se) c(-bw, bw)), finite = TRUE)
  .gcall(graphics::plot,
         list(x = th, y = Y[, 1], type = "n", bty = "n", ylim = yl,
              xlab = "theoretical quantile",
              ylab = if (detrend) "deviation" else "quantile residual",
              main = if (detrend) "Worm plot" else "Normal QQ"), dots)
  if (detrend && se) .band(th, -bw, bw, band.col)
  if (detrend) graphics::abline(h = 0, col = "grey60", lwd = 2)
  else graphics::abline(0, 1, col = "grey60", lwd = 2)
  for (k in seq_len(ncol(R)))
    graphics::points(th, Y[, k], pch = pch, cex = cx,
                     col = if (length(cols) > 1L) cols[k] else cols)
  out <- data.frame(theoretical = th, Y)
  nmY <- if (detrend) "deviation" else "residual"
  names(out)[-1L] <- if (ncol(Y) == 1L) nmY else paste0(nmY, seq_len(ncol(Y)))
  invisible(out)
}
