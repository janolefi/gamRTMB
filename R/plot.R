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
#' One panel per term, showing its contribution to that distributional
#' parameter's linear predictor — the scale on which terms are additive — with
#' a **pointwise** \eqn{\pm 2} standard error band. The band comes from
#' [vcov.gamRTMB()]'s joint covariance. Under `method = "REML"` or `"ML"` that
#' includes the uncertainty in the smoothing parameters (mgcv's
#' `unconditional = TRUE`) and needs a fit made with `joint_precision = TRUE`,
#' which is the default; under `method = "aREML"` it conditions on them, so the
#' band is a little too narrow. Panels are titled
#' `parameter: term`, since terms belong to different parameters.
#'
#' Parametric terms get panels too, as in `gamlss::term.plot` and unlike
#' [mgcv::plot.gam()], whose `all.terms` defaults to `FALSE`. A term is a
#' contribution to a parameter's predictor whether or not it is penalized, and
#' what a distributional model is for is usually how those contributions
#' differ between parameters; a parameter carrying no smooth at all should
#' still show what does act on it. Set `all.terms = FALSE` for the mgcv
#' behaviour. A parametric term's columns are centred at their means over the
#' fitting data, as in [stats::predict.lm()] with `type = "terms"`, so a panel
#' shows the term's variation rather than the level the intercept already
#' carries. A categorical term is drawn as an estimate and interval per level
#' rather than as a curve.
#'
#' Only one-dimensional terms in a single variable are drawn: on the smooth
#' side that excludes tensor products, random effects (`bs = "re"`) and
#' factor-smooth interactions (`bs = "fs"`), and on the parametric side
#' interactions. They are reported and skipped rather than drawn misleadingly.
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
#' @section Conditional densities (`type = "density"`):
#' The fitted density of the response turned on its side and drawn at a few
#' covariate values, over the data — the shape that the quantile fan only
#' summarises, and the clearest way to see a changing skewness or a changing
#' spread. It uses the family's log-density, so unlike the quantile plot it
#' works for every family.
#'
#' Densities are scaled to a common width, not a common height, so their
#' shapes are comparable; a lattice response gets a spike per integer rather
#' than a filled outline.
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
#' @param select Which terms to draw: an integer index, or a pattern matched
#'   against the `parameter: term` labels (e.g. `"sd"` or `"s(x1)"`). Defaults
#'   to all of them.
#' @param all.terms Give parametric terms their own panels in a term plot,
#'   alongside the smooths. `TRUE` by default; see the section above.
#' @param se Draw the interval band.
#' @param rug Add a rug of the observed covariate values.
#' @param band.col Fill for the interval band. A solid light grey by default,
#'   which renders the same on every device.
#' @param prob Probabilities for `type = "quantile"`. The default is a fine
#'   fan; a short vector such as `c(0.1, 0.5, 0.9)` also gets a legend.
#' @param xvar Covariate for the x-axis of a quantile or density plot.
#'   Defaults to the covariate of the first smooth.
#' @param at Covariate values at which to draw conditional densities.
#'   Defaults to five, evenly spaced and inset from the ends.
#' @param nsim Randomisation draws to overlay in a QQ or worm plot.
#' @param ask Draw one panel per page, waiting between them, instead of
#'   fitting everything onto one page.
#' @param n Grid resolution for term curves.
#' @param ... Passed to [plot()].
#' @return Invisibly, the plotted data, so any panel can be rebuilt by hand:
#'   one data frame per term for `type = "terms"`, the residual quantiles for
#'   the diagnostics, the fitted quantiles for `type = "quantile"`, and a long
#'   data frame of `at`, `y` and `density` for `type = "density"`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = runif(300), x2 = runif(300))
#' d$y <- rnorm(300, sin(2 * pi * d$x1), exp(-1 + d$x2))
#' fit <- gamRTMB(y ~ list(mean = ~ s(x1), sd = ~ s(x2)), data = d)
#' plot(fit)
#' plot(fit, select = "sd")
#' plot(fit, type = "worm")
#' @export
plot.gamRTMB <- function(x, type = c("terms", "qq", "worm", "quantile",
                                     "density"),
                         select = NULL, prob = seq(0.05, 0.95, by = 0.05),
                         xvar = NULL, at = NULL, se = TRUE, rug = TRUE,
                         band.col = "grey85", nsim = 1, ask = FALSE,
                         n = 200, all.terms = TRUE, ...) {
  type <- match.arg(type)
  switch(type,
    terms    = .plot_terms(x, select, se, rug, band.col, ask, n, all.terms,
                           ...),
    qq       = .plot_res(x, nsim, FALSE, se, band.col, ...),
    worm     = .plot_res(x, nsim, TRUE,  se, band.col, ...),
    quantile = .plot_quantile(x, xvar, prob, se, band.col, n, ...),
    density  = .plot_density(x, xvar, at, n, band.col, ...))
}

#' Colours for a fan of quantile curves
#'
#' One hue throughout, with the median in black: a fan of quantiles is one
#' ordered family, not a set of unrelated series, so different colours and
#' line types would imply distinctions that are not there. Opacity carries the
#' ordering instead, fading outwards from the median, which keeps a dense fan
#' legible and reads as the density it approximates.
#'
#' @param prob Probabilities, in the order they will be drawn.
#' @return A character vector of colours.
#' @keywords internal
.fan_cols <- function(prob, hue = "#0B5D9E") {
  d <- abs(prob - 0.5) / 0.5                       # 0 at the median, 1 at the ends
  out <- grDevices::adjustcolor(rep(hue, length(prob)),
                                alpha.f = 1)
  vapply(seq_along(prob), function(k)
    if (isTRUE(all.equal(d[k], 0))) "black"
    else grDevices::adjustcolor(hue, alpha.f = max(0.9 - 0.6 * d[k], 0.2)),
    "")
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

#' The covariate a response plot varies along
#'
#' Defaults to the covariate of the first drawable smooth.
#'
#' @keywords internal
.resolve_xvar <- function(fit, xvar) {
  if (is.null(xvar)) {
    v <- unlist(lapply(fit$design$parnames, function(p)
      lapply(fit$design$parts[[p]]$smooths, function(s) s$plot1d$var)))
    if (!length(v))
      stop("no smooth covariate to plot against; name one with xvar =")
    xvar <- v[[1L]]
  }
  if (!xvar %in% names(fit$data))
    stop("`xvar` must name a column of the model data: ",
         paste(names(fit$data), collapse = ", "))
  if (!is.numeric(fit$data[[xvar]]))
    stop("`xvar` must be a numeric covariate")
  xvar
}

#' New data varying one covariate, the others at a typical value
#'
#' The median for a numeric covariate and the modal level for a factor, so
#' that a response plot shows one curve per probability rather than a bundle.
#' With more than one covariate the fitted curves and the plotted points
#' therefore do not condition on quite the same thing.
#'
#' @keywords internal
.newdata_along <- function(fit, xvar, values) {
  nd <- fit$data[rep(1L, length(values)), , drop = FALSE]
  for (v in setdiff(names(nd), xvar))
    nd[[v]] <- rep(.typical(fit$data[[v]]), length(values))
  nd[[xvar]] <- values
  nd
}

#' Rotated conditional densities at a few covariate values
#'
#' The fitted density of the response, turned on its side and drawn at chosen
#' positions along a covariate — the shape the quantile fan only summarises.
#' It uses the family's own log-density, so it works for every family, not
#' only those with a quantile function.
#'
#' Each density opens to the left of its position line and is scaled to a
#' common width rather than a common height. A common height would be more
#' faithful — a concentrated distribution really does have a taller density —
#' but on data where the spread changes by a factor of 40 it makes the wide
#' ones invisible, and the shape is the point.
#'
#' @keywords internal
.plot_density <- function(x, xvar, at, ngrid, band.col, ...) {
  xvar <- .resolve_xvar(x, xvar)
  xv <- x$data[[xvar]]
  if (is.null(at)) {
    ## inset from the ends, or half of each density hangs off the panel
    r <- range(xv); pad <- 0.08 * diff(r)
    at <- seq(r[1L] + pad, r[2L] - pad, length.out = 5L)
  }
  nd <- .newdata_along(x, xvar, at)
  theta <- stats::predict(x, newdata = nd, type = "response")
  fx <- .resolve_fixed(x$family, nd, length(at))

  lattice <- x$family$support != "continuous"
  yl <- range(x$y, finite = TRUE)
  ## widest a density may be, in covariate units
  w <- if (length(at) > 1L) 0.3 * min(diff(sort(at))) else 0.12 * diff(range(xv))

  dfun <- function(yv, i) {
    d <- exp(x$family$logdens(yv, .at(theta, i), .at(fx, i)))
    d[!is.finite(d)] <- 0
    d
  }
  ## Two passes for a continuous response: locate each conditional's support
  ## on a coarse grid spanning the data, then re-evaluate finely over just
  ## that stretch. A single global grid renders a narrow conditional -- an sd
  ## of 1.5 where the response spans 250 -- as three or four points.
  curves <- lapply(seq_along(at), function(i) {
    if (lattice) {
      yv <- seq(floor(yl[1L]), ceiling(yl[2L]))
      d <- dfun(yv, i)
      k <- which(d > 0.005 * max(d, 0))
      if (length(k)) { k <- seq.int(min(k), max(k)); yv <- yv[k]; d <- d[k] }
    } else {
      y0 <- seq(yl[1L], yl[2L], length.out = ngrid)
      d0 <- dfun(y0, i)
      k <- which(d0 > 0.005 * max(d0, 0))
      r <- if (length(k)) range(y0[k]) else yl
      yv <- seq(r[1L], r[2L], length.out = ngrid)
      d <- dfun(yv, i)
    }
    if (max(d) > 0) d <- d / max(d)
    list(y = yv, d = d)
  })

  dots <- list(...)
  pch <- .or_else(dots$pch, 20); cx <- .or_else(dots$cex, 0.5)
  col <- .or_else(dots$col, "#0B5D9E")
  dots$pch <- dots$cex <- dots$col <- NULL
  ## include the drawn densities, or an outline that reaches past the data
  ## runs off the panel edge and looks like a rendering fault
  ylim <- range(c(yl, unlist(lapply(curves, `[[`, "y"))), finite = TRUE)
  .gcall(graphics::plot,
         list(x = xv, y = x$y, type = "n", bty = "n",
              xlim = range(c(xv, at - w)), ylim = ylim, xlab = xvar,
              ylab = deparse(x$formula[[2L]]),
              main = paste0("conditional densities: ", x$family$family)), dots)
  ## a full-height rule marks each position; the density outline is trimmed to
  ## where it has mass, so the two do not compete
  graphics::abline(v = at, col = grDevices::adjustcolor("black", 0.45), lwd = 1)
  graphics::points(xv, x$y, pch = pch, cex = cx,
                   col = grDevices::adjustcolor("black", 0.35))
  ## the fitted median, where the family has a quantile function, for a sense
  ## of where the centre runs between the positions
  if (!is.null(x$family$qf)) {
    xg <- seq(min(xv), max(xv), length.out = ngrid)
    med <- stats::predict(x, newdata = .newdata_along(x, xvar, xg),
                          type = "quantile", prob = 0.5)
    graphics::lines(xg, med[, 1], col = grDevices::adjustcolor("black", 0.8),
                    lwd = 1, type = if (lattice) "s" else "l")
  }
  ## Outline only: a fill competes with the data points it sits over.
  for (i in seq_along(at)) {
    yy <- curves[[i]]$y; xx <- at[i] - w * curves[[i]]$d
    if (lattice) graphics::segments(at[i], yy, xx, yy, col = col, lwd = 1.5)
    else graphics::lines(xx, yy, col = col, lwd = 1.4)
  }
  invisible(do.call(rbind, lapply(seq_along(at), function(i)
    data.frame(at = at[i], y = curves[[i]]$y, density = curves[[i]]$d))))
}

#' Fitted quantiles against one covariate
#' @keywords internal
.plot_quantile <- function(x, xvar, prob, se, band.col, ngrid, ...) {
  xvar <- .resolve_xvar(x, xvar)
  xv <- x$data[[xvar]]
  nd <- .newdata_along(x, xvar, seq(min(xv), max(xv), length.out = ngrid))

  band <- se && x$family$support == "continuous"
  qq <- stats::predict(x, newdata = nd, type = "quantile", prob = prob,
                       se.fit = band)
  Q <- if (band) qq$fit else qq
  S <- if (band) qq$se.fit else NULL
  mid <- which.min(abs(prob - 0.5))

  dots <- list(...)
  pch <- .or_else(dots$pch, 20); cx <- .or_else(dots$cex, 0.5)
  dots$pch <- dots$cex <- NULL
  cols <- .or_else(dots$col, .fan_cols(prob))
  dots$col <- NULL
  ## step curves for a lattice response, where quantiles really are steps
  ltype <- if (x$family$support == "continuous") "l" else "s"
  yl <- range(c(x$y, Q, if (band) c(Q[, mid] - 2 * S[, mid],
                                    Q[, mid] + 2 * S[, mid])), finite = TRUE)

  .gcall(graphics::plot,
         list(x = xv, y = x$y, type = "n", bty = "n", ylim = yl,
              xlab = xvar, ylab = deparse(x$formula[[2L]]),
              main = paste0("fitted quantiles: ", x$family$family)), dots)
  graphics::points(xv, x$y, pch = pch, cex = cx,
                   col = grDevices::adjustcolor("black", 0.4))
  for (k in seq_along(prob))
    if (k != mid)
      graphics::lines(nd[[xvar]], Q[, k], type = ltype, col = cols[k], lwd = 1)
  ## the band goes over the fan and under the median line: drawn underneath, a
  ## fine fan buries it completely
  if (band) .band(nd[[xvar]], Q[, mid] - 2 * S[, mid], Q[, mid] + 2 * S[, mid],
                  grDevices::adjustcolor(band.col, 0.75))
  graphics::lines(nd[[xvar]], Q[, mid], type = ltype, col = cols[mid], lwd = 2.5)
  ## a legend earns its place only for a handful of curves; a dense fan reads
  ## off the fading itself
  if (length(prob) <= 6L)
    graphics::legend(.empty_corner(c(xv, rep(nd[[xvar]], ncol(Q))),
                                   c(x$y, as.vector(Q)), yl),
                     legend = paste0(100 * prob, "%"), bty = "n", col = cols,
                     lwd = c(1, 2.5)[1 + (seq_along(prob) == mid)], cex = 0.85)
  out <- data.frame(nd[[xvar]], Q)
  names(out) <- c(xvar, colnames(Q))
  if (band) out$se.mid <- S[, mid]
  invisible(out)
}

#' A typical value of a covariate
#'
#' The median of a numeric one and the modal level of a categorical one, for
#' holding the rest of a term plot's covariates still.
#'
#' @keywords internal
.typical <- function(z) {
  if (is.numeric(z)) return(stats::median(z, na.rm = TRUE))
  u <- names(which.max(table(z)))
  if (is.factor(z)) factor(u, levels = levels(z))
  else if (is.logical(z)) as.logical(u) else u
}

#' The values a term plot's x-axis runs over
#'
#' A grid across the observed range of a numeric covariate, and every observed
#' level of a categorical one.
#'
#' @keywords internal
.term_grid <- function(z, ngrid) {
  if (is.numeric(z)) seq(min(z), max(z), length.out = ngrid)
  else if (is.factor(z)) factor(levels(z), levels = levels(z))
  else sort(unique(z))
}

#' The drawable first-order parametric terms of one parameter
#'
#' Everything a `termplot` would draw: a term of order one in a single
#' variable, so that `poly(x, 2)` and `log(x)` count but `x:z` does not. The
#' intercept is not a term and is skipped; it is not a contribution to the
#' predictor that varies with anything.
#'
#' @return A list of `list(k, var, label)`, `k` indexing the term within
#'   `terms(...)`, with `var = NULL` for a term that cannot be drawn.
#' @keywords internal
.para_terms <- function(part) {
  tt <- part$terms
  labs <- attr(tt, "term.labels")
  ord <- attr(tt, "order")
  asg <- attr(part$Xpara, "assign")
  lapply(seq_along(labs), function(k) {
    v <- all.vars(stats::reformulate(labs[k])[[2L]])
    list(k = k, label = labs[k],
         var = if (ord[k] == 1L && length(v) == 1L && any(asg == k)) v)
  })
}

#' One parametric term's contribution, as a linear form in the coefficients
#'
#' The same object a smooth panel draws, built the same way: a design for the
#' term alone, evaluated along its own covariate with the others at a typical
#' value, times the coefficients it multiplies. Its columns are centred at
#' their means over the fitting data, as [stats::predict.lm()] does for
#' `type = "terms"`, which fixes the free constant the intercept would
#' otherwise absorb and leaves a form whose standard errors are those of a
#' contrast rather than of an arbitrary level.
#'
#' @param fit A `gamRTMB` fit.
#' @param p,k Distributional parameter and parametric term index.
#' @param values The covariate values to evaluate at.
#' @return `list(Z, coef, ib)`, with `ib` indexing `beta`.
#' @keywords internal
.para_part <- function(fit, p, k, v, values) {
  P <- fit$design$parts[[p]]
  tt <- stats::delete.response(P$terms)
  nd <- lapply(all.vars(tt), function(u)
    if (identical(u, v)) values else rep(.typical(fit$data[[u]]), length(values)))
  names(nd) <- all.vars(tt)
  mf <- stats::model.frame(tt, as.data.frame(nd, stringsAsFactors = FALSE),
                           xlev = P$xlev)
  X <- stats::model.matrix(tt, mf, contrasts.arg = attr(P$Xpara, "contrasts"))
  cols <- which(attr(P$Xpara, "assign") == k)
  Z <- X[, cols, drop = FALSE]
  Z <- sweep(Z, 2L, colMeans(P$Xpara[, cols, drop = FALSE]))
  ib <- fit$design$beta_idx[[p]][cols]
  list(Z = Z, coef = fit$coefficients$beta[ib], ib = ib)
}

#' @keywords internal
.plot_terms <- function(x, select, se, rug, band.col, ask, ngrid, all.terms,
                        ...) {
  D <- x$design
  tl <- list()
  add <- function(...) tl[[length(tl) + 1L]] <<- list(...)
  for (p in D$parnames) {
    if (all.terms) for (t in .para_terms(D$parts[[p]]))
      add(par = p, kind = "para", k = t$k, var = t$var,
          label = paste0(p, ": ", t$label), ylab = t$label)
    for (j in seq_along(D$parts[[p]]$smooths)) {
      s <- D$parts[[p]]$smooths[[j]]
      add(par = p, kind = "smooth", j = j, s = s, var = s$plot1d$var,
          label = paste0(p, ": ", s$label), ylab = s$label)
    }
  }
  if (!length(tl))
    stop(if (all.terms)
      "the model has no terms to plot: every distributional parameter is an intercept."
    else
      "the model has no smooth terms to plot; all.terms = TRUE would draw its parametric ones.")
  labs <- function(z) vapply(z, `[[`, "", "label")

  if (!is.null(select)) {
    keep <- if (is.numeric(select)) {
      if (any(select < 1 | select > length(tl)))
        stop("`select` must index the ", length(tl), " terms")
      as.integer(select)
    } else grep(select, labs(tl), fixed = FALSE)
    if (!length(keep))
      stop("`select` matched none of: ", paste(labs(tl), collapse = ", "))
    tl <- tl[keep]
  }

  drawable <- !vapply(tl, function(t) is.null(t$var), TRUE)
  if (any(!drawable))
    message("not drawable as a single panel, skipped: ",
            paste(labs(tl[!drawable]), collapse = ", "))
  tl <- tl[drawable]
  if (!length(tl)) stop("no one-dimensional terms left to plot")

  Vj <- if (se) tryCatch(.joint_cov(x), error = function(e) NULL) else NULL
  if (se && is.null(Vj))
    message("no interval band: the joint covariance is unavailable ",
            "(refit with joint_precision = TRUE)")

  ## build every curve first, so a failure cannot leave a half-drawn page
  curves <- lapply(tl, function(t) {
    xv <- x$data[[t$var]]
    xg <- if (identical(t$kind, "smooth")) .term_grid(xv, ngrid)
          else .term_grid(xv, min(ngrid, 100L))
    sp <- if (identical(t$kind, "para")) {
      q <- .para_part(x, t$par, t$k, t$var, xg)
      list(Z = q$Z, coef = q$coef, ii = Vj$ib[q$ib])
    } else {
      nd <- stats::setNames(list(xg), t$var)
      ps <- t$s$plot1d
      if (!is.null(ps$by)) nd[[ps$by]] <- rep(ps$by_val, length(xg))
      Xs <- mgcv::PredictMat(t$s$sm,
                             as.data.frame(nd, stringsAsFactors = FALSE))
      q <- .smooth_part(x, t$par, t$j, Xs)
      list(Z = q$Z, coef = q$coef, ii = c(Vj$ir[q$b], Vj$ib[q$f]))
    }
    sef <- if (is.null(Vj)) rep(NA_real_, length(xg))
           else .qform_se(sp$Z, Vj$V[sp$ii, sp$ii, drop = FALSE])
    data.frame(x = xg, fit = as.vector(sp$Z %*% sp$coef), se = sef)
  })
  names(curves) <- labs(tl)

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
    cv <- curves[[i]]; t <- tl[[i]]
    hi <- cv$fit + 2 * cv$se; lo <- cv$fit - 2 * cv$se
    yl <- if (all(is.na(cv$se))) range(cv$fit) else range(c(lo, hi), finite = TRUE)
    ## a categorical term is a handful of levels, not a curve: draw the
    ## estimate and its interval at each level rather than joining them, which
    ## would imply an ordering the levels need not have
    if (is.numeric(cv$x)) {
      .gcall(graphics::plot,
             list(x = cv$x, y = cv$fit, type = "n", bty = "n", ylim = yl,
                  xlab = t$var, ylab = t$ylab, main = names(curves)[i]), dots)
      if (!all(is.na(cv$se))) .band(cv$x, lo, hi, fill)
      graphics::lines(cv$x, cv$fit, col = col, lwd = lwd)
      if (rug) graphics::rug(x$data[[t$var]],
                             col = grDevices::adjustcolor(col, 0.4))
    } else {
      at <- seq_along(cv$x)
      ## a point sitting on the frame reads as clipped, and a tight interval
      ## puts every one of them there; a curve running to the edge does not
      yl <- yl + c(-1, 1) * 0.04 * max(diff(yl), 1e-8)
      .gcall(graphics::plot,
             list(x = at, y = cv$fit, type = "n", bty = "n", ylim = yl,
                  xlim = c(0.5, length(at) + 0.5), xaxt = "n",
                  xlab = t$var, ylab = t$ylab, main = names(curves)[i]), dots)
      graphics::axis(1, at = at, labels = as.character(cv$x), lwd = 0,
                     lwd.ticks = 1)
      if (!all(is.na(cv$se)))
        graphics::segments(at, lo, at, hi, col = col, lwd = lwd)
      graphics::points(at, cv$fit, pch = .or_else(dots$pch, 19), col = col)
    }
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
