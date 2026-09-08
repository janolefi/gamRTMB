## The mgcv-facing half of the package: build bases and penalties with
## smoothCon, reparameterise them with smooth2random, and record enough to
## reconstruct and predict later. No RTMB here.

#' Copy a basis specification onto another term
#'
#' Local stand-in for the unexported `mgcv:::clone.smooth.spec`: the first
#' member of an `id` group defines the basis, and each later term gets that
#' definition back with its own identity (term names, label, `by`, `xt`).
#'
#' @keywords internal
.clone_spec <- function(base, spec) {
  out <- base
  out$term <- spec$term; out$label <- spec$label; out$by <- spec$by
  if (inherits(base, c("tensor.smooth.spec", "t2.smooth.spec"))) {
    if (!is.null(spec$margin))
      for (i in seq_along(out$margin)) {
        out$margin[[i]]$term  <- spec$margin[[i]]$term
        out$margin[[i]]$label <- spec$margin[[i]]$label
        out$margin[[i]]$xt    <- spec$margin[[i]]$xt
      }
  } else out$xt <- spec$xt
  out
}

#' Map reparameterised coefficients back to a smooth's own basis
#'
#' [mgcv::smooth2random()] returns penalized design blocks (`rand`), an
#' unpenalized null-space block (`Xf`), and the transformation back to the
#' smooth's original constrained basis:
#' \deqn{\beta = U (D \cdot c(b_{rand}, b_{fixed})[idx])}
#' with `idx = c(rind, offset + seq_len(n_fixed))`.
#'
#' For ordinary smooths `rind` is the identity, but for `bs = "fs"` it is a
#' genuine permutation: the random blocks are penalty-major while the
#' coefficients are level-major. Assuming `cbind(rand, Xf)` therefore gives
#' silently wrong coefficients for factor-smooth interactions, so the map is
#' built once here and reused. Returning it as an explicit matrix means the
#' same object serves the point reconstruction and the delta-method
#' covariance of a smooth.
#'
#' Verified exact (to 6e-16) for `s()`, `s(bs = "cr")`, `t2()`, `by =`
#' factors, `bs = "fs"` and `bs = "re"`.
#'
#' @keywords internal
.reconstruct_map <- function(re) {
  nb <- sum(vapply(re$rand, ncol, 1L))
  nf <- ncol(re$Xf)
  idx <- c(re$rind, if (nf) nb + seq_len(nf))
  Tm <- matrix(0, length(idx), nb + nf)
  Tm[cbind(seq_along(idx), idx)] <- 1
  Tm <- re$trans.D * Tm
  if (!is.null(re$trans.U)) Tm <- re$trans.U %*% Tm
  Tm
}

#' Build the design for every distributional parameter
#'
#' @section Identifiability:
#' `smoothCon(absorb.cons = TRUE)` applies the sum-to-zero constraint, so
#' smooths cannot collide with the intercept.
#'
#' @section Null spaces:
#' Only the penalized blocks become random effects. The unpenalized
#' null-space columns are appended to that parameter's fixed-effect matrix and
#' are never pooled into the random-effect sum.
#'
#' @section Shared smoothing parameters:
#' `s(..., id = )` does two things in mgcv and both are reproduced.
#' \emph{Linked bases}: the basis is built from the pooled covariate values of
#' the whole id group and then evaluated on each term's own data, via
#' `smoothCon(spec, pooled, n = nrow(data), dataX = data)`. This matters
#' beyond knot placement, because `scale.penalty` rescales each penalty by a
#' norm of its own model matrix; without pooling, two id-linked smooths get
#' different `S.scale` and one shared variance would mean different amounts of
#' smoothing for each. \emph{Shared smoothing parameter}: blocks in a group
#' get the same `sig_key`, which the engine turns into one variance component.
#'
#' id groups are resolved globally across distributional parameters, so
#' `s(x, id = 1)` on `mu` and on `sigma` share one smoothness. A gam formula
#' cannot express this, having only one linear predictor.
#'
#' @param par_formulas List of one-sided formulas, one per parameter.
#' @param data Model frame.
#' @param parnames The family's modelled parameters.
#' @param knots Passed to [mgcv::smoothCon()].
#' @return A design object: per-parameter parts (parametric matrix, terms,
#'   xlevels, smooths), the stacked fixed-effect matrices, index bookkeeping
#'   into `beta` and `b`, the penalized blocks, and the variance-component
#'   grouping.
#' @keywords internal
.build_design <- function(par_formulas, data, parnames, knots = NULL) {
  n <- nrow(data)

  ## pass 1 -- interpret every formula and flatten the smooth specifications,
  ## so that id groups can be found across distributional parameters
  gps <- list(); specs <- list()
  for (p in parnames) {
    gps[[p]] <- mgcv::interpret.gam(par_formulas[[p]])
    for (j in seq_along(gps[[p]]$smooth.spec))
      specs[[length(specs) + 1L]] <- list(par = p, j = j,
                                          spec = gps[[p]]$smooth.spec[[j]])
  }
  ids <- vapply(specs, function(z)
    if (is.null(z$spec$id)) NA_character_ else as.character(z$spec$id), "")

  ## pass 2 -- within each id group, clone the basis specification and pool
  ## the covariate values
  pooled <- list()
  for (idv in unique(ids[!is.na(ids)])) {
    grp <- which(ids == idv)
    base <- specs[[grp[1L]]]$spec
    if (length(grp) > 1L)
      for (g in grp[-1L]) specs[[g]]$spec <- .clone_spec(base, specs[[g]]$spec)
    pooled[[idv]] <- lapply(seq_along(base$term), function(k)
      do.call(cbind, lapply(grp, function(g)
        mgcv::get.var(specs[[g]]$spec$term[k], data, vecMat = FALSE))))
  }
  spec_of <- function(p, j) {
    k <- which(vapply(specs, function(z) z$par == p && z$j == j, TRUE))
    list(spec = specs[[k]]$spec, id = ids[k])
  }

  ## pass 3 -- construct the smooths
  parts <- list()
  for (p in parnames) {
    gp <- gps[[p]]
    tt <- stats::terms(gp$pf, data = data)
    Xp <- stats::model.matrix(tt, data)
    xlev <- stats::.getXlevels(tt, stats::model.frame(tt, data))

    smooths <- list()
    for (j in seq_along(gp$smooth.spec)) {
      sj <- spec_of(p, j)
      scl <- if (is.na(sj$id))
        mgcv::smoothCon(sj$spec, data = data, knots = knots,
                        absorb.cons = TRUE, null.space.penalty = FALSE)
      else {
        pd <- pooled[[sj$id]]; names(pd) <- sj$spec$term
        mgcv::smoothCon(sj$spec, data = pd, knots = knots, absorb.cons = TRUE,
                        n = n, dataX = data, null.space.penalty = FALSE)
      }
      for (sm in scl) {
        if (isTRUE(sm$fixed))
          stop("fx = TRUE smooths are not supported: there is no penalized ",
               "part to make random")
        re <- mgcv::smooth2random(sm, "", type = 2)
        smooths[[length(smooths) + 1L]] <- list(
          sm = sm, re = re, Tmap = .reconstruct_map(re),
          label = sm$label, id = sj$id,
          Xr = lapply(re$rand, as.matrix), Xf = re$Xf)
      }
    }
    parts[[p]] <- list(Xpara = Xp, terms = tt, xlev = xlev, smooths = smooths)
  }

  ## index bookkeeping over the two coefficient vectors ---------------------
  beta_idx <- list(); blocks <- list(); par_blocks <- list(); Xfix <- list()
  nbeta <- 0L; nb <- 0L
  for (p in parnames) {
    P <- parts[[p]]
    Xf_all <- P$Xpara
    labels <- colnames(P$Xpara)
    for (s in P$smooths)
      if (ncol(s$Xf)) {
        Xf_all <- cbind(Xf_all, s$Xf)
        labels <- c(labels, paste0(s$label, ".null", seq_len(ncol(s$Xf))))
      }
    Xfix[[p]] <- Xf_all
    beta_idx[[p]] <- nbeta + seq_len(ncol(Xf_all))
    attr(beta_idx[[p]], "labels") <- labels
    nbeta <- nbeta + ncol(Xf_all)

    off <- ncol(P$Xpara); ids_p <- integer(0)
    for (j in seq_along(P$smooths)) {
      s <- P$smooths[[j]]
      parts[[p]]$smooths[[j]]$f_local <-
        if (ncol(s$Xf)) off + seq_len(ncol(s$Xf)) else integer(0)
      off <- off + ncol(s$Xf)
      loc <- integer(0)
      for (k in seq_along(s$Xr)) {
        q <- ncol(s$Xr[[k]])
        blocks[[length(blocks) + 1L]] <- list(
          par = p, smooth = j, penalty = k, q = q,
          label = if (length(s$Xr) > 1L) paste0(s$label, ".s", k) else s$label,
          ## one variance component per penalty per smooth per parameter,
          ## unless an id ties this penalty across a group
          sig_key = if (is.na(s$id)) paste(p, j, k, sep = ":")
                    else paste0("id", s$id, ":", k),
          idx = nb + seq_len(q), X = s$Xr[[k]])
        loc <- c(loc, length(blocks)); nb <- nb + q
      }
      parts[[p]]$smooths[[j]]$block_ids <- loc
      ids_p <- c(ids_p, loc)
    }
    par_blocks[[p]] <- ids_p
  }

  keys <- vapply(blocks, `[[`, "", "sig_key")
  sig_group <- factor(keys, levels = unique(keys))

  structure(list(
    parts = parts, parnames = parnames, n = n, Xfix = Xfix,
    beta_idx = beta_idx, blocks = blocks, par_blocks = par_blocks,
    nbeta = nbeta, nb = nb, nsigma = length(blocks),
    sig_group = sig_group, nsigma_free = nlevels(sig_group)
  ), class = "gamRTMB_design")
}
