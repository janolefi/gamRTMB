## Turning a formula and a data frame into bases, penalties and index
## bookkeeping. This is the whole mgcv-facing half of the package; no RTMB.

#' Parse a distributional model formula
#'
#' The response is taken from the left-hand side of the outer formula and each
#' distributional parameter gets a one-sided formula from the `list()` on the
#' right, in the style of `brms::bf()`. Any parameter the family declares but
#' the user does not mention is given `~1`, which is why the family spec has
#' to exist before formula processing.
#'
#' @param formula e.g. `y ~ list(mean = ~ s(x1) + s(x2), sd = ~ s(x1))`.
#' @param parnames Character vector of the family's modelled parameters.
#' @return A list with the response expression and one formula per parameter,
#'   ordered as `parnames`.
#' @keywords internal
.parse_formula <- function(formula, parnames) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be two-sided, e.g. y ~ list(mean = ~ s(x))")
  rhs <- formula[[3L]]
  if (!(is.call(rhs) && identical(rhs[[1L]], as.name("list"))))
    stop("the right-hand side must be a list() of per-parameter formulas, ",
         "e.g. y ~ list(mean = ~ s(x), sd = ~ 1)")
  args <- as.list(rhs)[-1L]
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args)))))
    stop("every element of the list() must be named after a distributional parameter")
  unknown <- setdiff(names(args), parnames)
  if (length(unknown))
    stop("unknown distributional parameter(s) ", paste(unknown, collapse = ", "),
         "; this family models ", paste(parnames, collapse = ", "))
  env <- environment(formula)
  out <- lapply(args, function(a) {
    f <- eval(a, env)
    if (!inherits(f, "formula")) stop("each list() element must be a formula")
    if (length(f) == 3L) stop("per-parameter formulas must be one-sided")
    environment(f) <- env
    f
  })
  for (p in setdiff(parnames, names(out))) out[[p]] <- stats::as.formula("~1", env)
  list(response = formula[[2L]], par_formulas = out[parnames])
}

#' Assemble the model data
#'
#' Collects every variable the model touches — via [mgcv::interpret.gam()]'s
#' `fake.formula`, which reports the variables inside `s()` terms as well as
#' the parametric and offset ones — applies `na.action` to that set, and
#' returns the data with incomplete rows dropped, together with the response
#' and prior weights aligned to it.
#'
#' Model variables must be columns of `data`. R would otherwise let them come
#' from the calling environment, where dropping rows for missing values could
#' silently misalign them against the response.
#'
#' @param response The response expression (LHS of the outer formula).
#' @param par_formulas One-sided formulas, one per distributional parameter.
#' @param data A data frame.
#' @param weights Evaluated prior weights, or `NULL`.
#' @param na.action Missing-data action, e.g. [stats::na.omit()].
#' @return `list(data, y, weights, dropped)`.
#' @keywords internal
.model_data <- function(response, par_formulas, data, weights, na.action) {
  vars <- unique(c(all.vars(response),
                   unlist(lapply(par_formulas, function(f)
                     all.vars(mgcv::interpret.gam(f)$fake.formula)))))
  extra <- setdiff(vars, names(data))
  if (length(extra))
    stop("every model variable must be a column of `data`; not found: ",
         paste(extra, collapse = ", "))

  n0 <- nrow(data)
  if (!is.null(weights)) {
    if (!is.numeric(weights))
      stop("`weights` must be numeric")
    if (length(weights) == 1L) weights <- rep(weights, n0)
    if (length(weights) != n0)
      stop("`weights` has length ", length(weights), ", need 1 or ", n0)
    if (any(weights < 0, na.rm = TRUE)) stop("`weights` must be non-negative")
  }

  mf <- stats::model.frame(stats::reformulate(vars), data = data,
                           na.action = na.action)
  drop <- attr(mf, "na.action")
  if (length(drop)) {
    data <- data[-drop, , drop = FALSE]
    if (!is.null(weights)) weights <- weights[-drop]
  }
  if (!is.null(weights) && anyNA(weights))
    stop("`weights` contains missing values")
  if (!nrow(data)) stop("no complete observations left after na.action")

  ## every model variable is a column of `data` (checked above), so the
  ## response needs no enclosing environment beyond base
  list(data = data, y = eval(response, data, baseenv()),
       weights = weights, dropped = length(drop))
}

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

#' What a term plot needs, recorded at design time
#'
#' A one-dimensional smooth of a numeric covariate can be drawn on its own, so
#' record which covariate it uses and, for a `by=` smooth, the value of the
#' `by` variable that [mgcv::PredictMat()] will want: the smooth's own factor
#' level, or 1 for a numeric `by`. The covariate values themselves come from
#' the fit's stored model frame.
#'
#' Returns `NULL` for anything not drawable as a single curve — tensor
#' products, random effects, factor-smooth interactions — which the plot
#' method reports rather than drawing wrongly.
#'
#' @keywords internal
.plot_spec <- function(sm, data) {
  if (sm$dim != 1L || length(sm$term) != 1L) return(NULL)
  if (inherits(sm, c("random.effect", "fs.interaction"))) return(NULL)
  xv <- data[[sm$term]]
  if (!is.numeric(xv)) return(NULL)
  has_by <- !is.null(sm$by) && !identical(as.character(sm$by), "NA")
  by_val <- NULL
  if (has_by) {
    bv <- data[[sm$by]]
    by_val <- if (is.factor(bv)) factor(sm$by.level, levels = levels(bv)) else 1
  }
  list(var = sm$term, by = if (has_by) sm$by else NULL, by_val = by_val)
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
#' @section Offsets:
#' An `offset()` term inside a parameter's formula adds a fixed, known
#' contribution to that parameter's linear predictor. Offsets are
#' per-parameter because that is the only meaningful reading in a
#' distributional model: `sd = ~ offset(log(s))` says something quite
#' different from the same term on `mean`. The term is recomputed from
#' `newdata` when predicting.
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
#' @param sparse Whether a smooth with a single sparse penalty may skip
#'   [mgcv::smooth2random()] and keep that penalty; see [.gmrf_block()].
#' @return A design object: per-parameter parts (parametric matrix, terms,
#'   xlevels, smooths), the stacked fixed-effect matrices, index bookkeeping
#'   into `beta` and `b`, the penalized blocks, and the variance-component
#'   grouping.
#' @keywords internal
.build_design <- function(par_formulas, data, parnames, knots = NULL,
                          sparse = "auto") {
  n <- nrow(data)

  ## pass 1 -- interpret every formula and flatten the smooth specifications,
  ## so that id groups can be found across distributional parameters
  gps <- lapply(par_formulas, mgcv::interpret.gam)
  at <- do.call(rbind, lapply(parnames, function(p)
    if (length(gps[[p]]$smooth.spec))
      data.frame(par = p, j = seq_along(gps[[p]]$smooth.spec)) else NULL))
  ids <- if (is.null(at)) character(0) else
    vapply(seq_len(nrow(at)), function(i) {
      idv <- gps[[at$par[i]]]$smooth.spec[[at$j[i]]]$id
      if (is.null(idv)) NA_character_ else as.character(idv)
    }, "")

  ## pass 2 -- within each id group, clone the basis specification and pool
  ## the covariate values
  ## the cloned specs are written back where pass 3 will find them, so the
  ## id can simply be read off the spec rather than tracked alongside it
  pooled <- list()
  for (idv in unique(ids[!is.na(ids)])) {
    grp <- which(ids == idv)
    base <- gps[[at$par[grp[1L]]]]$smooth.spec[[at$j[grp[1L]]]]
    for (g in grp[-1L])
      gps[[at$par[g]]]$smooth.spec[[at$j[g]]] <-
        .clone_spec(base, gps[[at$par[g]]]$smooth.spec[[at$j[g]]])
    pooled[[idv]] <- lapply(seq_along(base$term), function(k)
      do.call(cbind, lapply(grp, function(g)
        mgcv::get.var(gps[[at$par[g]]]$smooth.spec[[at$j[g]]]$term[k],
                      data, vecMat = FALSE))))
  }

  ## pass 3 -- construct the smooths
  parts <- list()
  for (p in parnames) {
    gp <- gps[[p]]
    tt <- stats::terms(gp$pf, data = data)
    pmf <- stats::model.frame(tt, data)
    Xp <- stats::model.matrix(tt, data)
    xlev <- stats::.getXlevels(tt, pmf)
    off <- stats::model.offset(pmf)
    if (is.null(off)) off <- 0 else if (length(off) == 1L) off <- rep(off, n)

    smooths <- list()
    for (j in seq_along(gp$smooth.spec)) {
      spec <- gp$smooth.spec[[j]]
      idv <- if (is.null(spec$id)) NA_character_ else as.character(spec$id)
      build <- function(absorb) if (is.na(idv))
        mgcv::smoothCon(spec, data = data, knots = knots,
                        absorb.cons = absorb, null.space.penalty = FALSE)
      else {
        pd <- pooled[[idv]]; names(pd) <- spec$term
        mgcv::smoothCon(spec, data = pd, knots = knots, absorb.cons = absorb,
                        n = n, dataX = data, null.space.penalty = FALSE)
      }
      ## The sparse route needs the penalty as the basis constructor wrote
      ## it, so the term is built unconstrained first and only rebuilt with
      ## the constraint absorbed if it turns out not to qualify. Absorbing
      ## the constraint is what would densify the penalty, so there is no
      ## way to make this decision from the constrained object.
      scl <- build(sparse == "never")
      use_sp <- sparse != "never" && length(scl[[1L]]$S) == 1L &&
        .use_sparse(scl[[1L]]$S[[1L]], sparse)
      if (!use_sp && sparse != "never") scl <- build(TRUE)
      for (sm in scl) {
        if (isTRUE(sm$fixed))
          stop("fx = TRUE smooths are not supported: there is no penalized ",
               "part to make random")
        B <- if (use_sp) .gmrf_block(sm) else {
          re <- mgcv::smooth2random(sm, "", type = 2)
          list(Xr = lapply(re$rand, as.matrix), Xf = re$Xf,
               Tmap = .reconstruct_map(re),
               Q = vector("list", length(re$rand)), intrinsic = FALSE)
        }
        smooths[[length(smooths) + 1L]] <- list(
          sm = sm, Tmap = B$Tmap, label = sm$label, id = idv,
          Xr = B$Xr, Xf = B$Xf, Q = B$Q,
          sparse = use_sp, intrinsic = B$intrinsic,
          plot1d = .plot_spec(sm, data))
      }
    }
    parts[[p]] <- list(Xpara = Xp, terms = tt, xlev = xlev, offset = off,
                       smooths = smooths)
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
        Qk <- s$Q[[k]]
        blocks[[length(blocks) + 1L]] <- list(
          par = p, smooth = j, penalty = k, q = q,
          label = if (length(s$Xr) > 1L) paste0(s$label, ".s", k) else s$label,
          ## one variance component per penalty per smooth per parameter,
          ## unless an id ties this penalty across a group
          sig_key = if (is.na(s$id)) paste(p, j, k, sep = ":")
                    else paste0("id", s$id, ":", k),
          idx = nb + seq_len(q), X = s$Xr[[k]], Q = Qk,
          ## `sigma` scales an iid coefficient directly, but a GMRF
          ## coefficient only through the penalty: its conditional standard
          ## deviation is sigma / sqrt(Q_ii). Recording a typical Q_ii lets
          ## one starting-value rule serve both.
          qscale = if (is.null(Qk)) 1 else
            sqrt(mean(Matrix::diag(Qk))))
        loc <- c(loc, length(blocks)); nb <- nb + q
      }
      parts[[p]]$smooths[[j]]$block_ids <- loc
      ids_p <- c(ids_p, loc)
    }
    par_blocks[[p]] <- ids_p
  }

  ## A corner-constrained intrinsic field pins one coefficient per connected
  ## component at zero instead of centring the whole term, so the level it
  ## gives up has to land somewhere. With an intercept that is exactly what
  ## happens and the fit is identical to the sum-to-zero version; without one,
  ## the constraint silently forces the dropped region's effect to zero, which
  ## is a different model and not the one anybody meant.
  for (p in parnames) {
    if (!any(vapply(parts[[p]]$smooths, function(s) isTRUE(s$intrinsic), NA)))
      next
    if (!"(Intercept)" %in% colnames(parts[[p]]$Xpara))
      stop("'", p, "' has an intrinsic GMRF smooth but no intercept. Such a ",
           "term is identified only up to a constant per connected component, ",
           "which the intercept absorbs; without one, the constraint would ",
           "pin one region's effect at zero instead. Add an intercept, or fit ",
           "this term with sparse = \"never\".", call. = FALSE)
  }

  ## Overlapping unpenalized null spaces make the fixed-effect design rank
  ## deficient. A smooth's null space is a low-order polynomial in its own
  ## covariate, so `x + s(x)` fits the x main effect twice; a tensor product
  ## carries main effects for each margin, so `t2(x, z) + s(z)` does too; and
  ## a by-factor smooth carries one per level, so `s(x) + s(x, by = g)` does
  ## as well.
  ##
  ## mgcv tolerates this because its fitted values stay identifiable even
  ## when the individual coefficients do not. Here it is fatal rather than
  ## untidy: a flat direction in the coefficients makes the penalized Hessian
  ## singular, so the Laplace approximation's log-determinant is undefined and
  ## the objective comes back non-finite. Caught here, where the cause is
  ## still visible, rather than as a puzzling remark about the response.
  for (p in parnames) {
    X <- Xfix[[p]]
    if (ncol(X) > 1L && qr(X)$rank < ncol(X)) {
      lab <- attr(beta_idx[[p]], "labels")
      stop("the fixed-effect design for '", p, "' is rank deficient, so the ",
           "model is not identified: ", paste(lab, collapse = ", "), ".\n",
           "Two smooths (or a smooth and a parametric term) share an ",
           "unpenalized null space. Drop the duplicate: `x + s(x)` should be ",
           "just `s(x)`, `t2(x, z) + s(z)` just `t2(x, z)`, and ",
           "`s(x) + s(x, by = g)` is better written `s(x, g, bs = \"fs\")`, ",
           "which is fully penalized and needs no separate main effect.")
    }
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
