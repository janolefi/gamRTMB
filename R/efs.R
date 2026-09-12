## method = "aREML": the REML criterion optimised by extended Fellner-Schall
## rather than by handing the smoothing parameters to nlminb. The reason it is
## a separate fitting routine rather than an optimiser setting is in
## dev/NOTES-fellner-schall.md: the REML outer gradient differentiates log|H|
## with respect to the smoothing parameters, which needs third derivatives of
## the log-likelihood in the coefficients. Fellner-Schall drops exactly that
## term, and the saving is only real if the Laplace machinery is never built
## -- so nothing is declared `random` here, and this file owns both halves of
## the fit: an inner solve over the coefficients and an outer multiplicative
## update of the smoothing parameters.
##
## Everything below is written against `design$blocks`, never against a
## particular kind of block. A block's penalty is
##   Q_k(theta) = sum_i exp(a_ki(theta)) S_ki,   a_ki linear in theta,
## which covers all three routes the design can produce (see
## [.block_penalties()]), so a smooth with several penalty matrices on one
## coefficient vector -- a `te()`, an adaptive smooth, an SPDE -- needs
## nothing added here.

#' A block's penalty matrices and their exponents
#'
#' The engine's view of a penalized block, and the only place that knows how
#' a block's parameters turn into a precision for a *numeric* theta:
#' \deqn{Q_k(\theta) = \sum_i \exp(A_{i\cdot}\theta) S_{ki}.}
#' [.block_prec()] is the same statement written for the AD tape, where the
#' switch on `kind` has to stay so that `dnorm` and `dgmrf` see the shapes
#' they expect. The two must agree; `test-efs.R` checks that they do.
#'
#' The `A` matrix is mgcv's `L` convention generalised to every block:
#' `\eqn{\lambda = \exp(L\theta)}`. An `"iid"` or `"gmrf"` block has one
#' penalty matrix and one parameter, the log standard deviation, so
#' \eqn{\lambda = \exp(-2\theta)} and `A` is the 1 by 1 matrix `-2`.
#'
#' @param bl A block from [.build_design()].
#' @return `list(S, A)`: a list of `nS` sparse penalty matrices, and an `nS`
#'   by `bl$ntheta` matrix of exponent coefficients.
#' @keywords internal
.block_penalties <- function(bl) {
  switch(bl$kind,
    iid   = list(S = list(.sparse_I(bl$q)), A = matrix(-2, 1L, 1L)),
    gmrf  = list(S = list(bl$Q), A = matrix(-2, 1L, 1L)),
    multi = list(S = bl$Smats, A = as.matrix(bl$L)),
    stop("unknown block kind ", sQuote(bl$kind)))
}

#' Force a matrix into the symmetric sparse storage the solves want
#'
#' `forceSymmetric()` alone returns a `dsy`/`dsC` of whatever storage it was
#' given; `Matrix::Cholesky()` and the sparse solves want `dsCMatrix`.
#'
#' @keywords internal
.sym_sparse <- function(M)
  as(as(Matrix::forceSymmetric(M), "symmetricMatrix"), "CsparseMatrix")

#' A sparse identity, in the storage the trace helpers expect
#' @keywords internal
.sparse_I <- function(q)
  Matrix::sparseMatrix(i = seq_len(q), j = seq_len(q), x = 1, dims = c(q, q))

#' A sparse penalty as row/column/value triplets
#'
#' `S` is symmetric and may be stored as one triangle, so it goes through
#' `"generalMatrix"` first to get both. Done once per penalty matrix at setup,
#' because the pattern is fixed by the design and only the weights move.
#'
#' @keywords internal
.tri_S <- function(S) {
  ts <- Matrix::summary(as(as(S, "generalMatrix"), "TsparseMatrix"))
  list(ij = cbind(ts$i, ts$j), x = ts$x)
}

#' Trace of a dense inverse block against a sparse penalty
#'
#' \eqn{tr(V S) = \sum_{rs} V_{rs} S_{sr}}, summed over `S`'s nonzeros only,
#' from the triplets [.tri_S()] prepared.
#'
#' @keywords internal
.tr_VS <- function(V, tri)
  if (!length(tri$x)) 0 else sum(V[tri$ij] * tri$x)

#' Quadratic form in a sparse penalty
#' @keywords internal
.quad_S <- function(x, S) as.numeric(crossprod(x, as.numeric(S %*% x)))

## Above this many coefficients the positive-definiteness repair stops being
## able to afford a dense eigendecomposition and falls back to a ridge. Same
## threshold as [.inner_indefinite()] uses for the same reason.
.efs_dense_max <- 1000L

## Floor on the numerator and denominator of the multiplicative update. Both
## are non-negative by construction -- the numerator is a block's penalized
## effective degrees of freedom, the denominator a quadratic form in a
## positive semi-definite penalty -- and both go to zero together as a term
## is shrunk onto its null space. The floor keeps 0/0 from deciding the step;
## the step cap then decides it instead.
.efs_floor <- 1e-12

#' Control parameters for the extended Fellner-Schall engine
#'
#' @format A named list of defaults, merged with the user's `control`.
#' @keywords internal
.efs_defaults <- list(
  maxit       = 200L,   # outer iterations
  tol         = 1e-8,   # relative change in the criterion, accepted steps
  ## Looser than `tol`, and necessarily so: that compares successive accepted
  ## values of a monotone sequence, where the only floor is arithmetic, while
  ## this compares a *rejected* candidate at a step already halved to nearly
  ## nothing. Measured rather than chosen -- on a Markov random field the gap
  ## it has to forgive sits between 1e-6 and 1e-5 relative, and it is the
  ## inner solve's own accuracy rather than curvature, since halving the step
  ## twenty-five more times does not shrink it. At 1e-8 a converged fit was
  ## reported as a failure.
  stall_tol   = 1e-5,
  ## A sufficient condition for convergence, never a necessary one: the EFS
  ## gradient drops a third-derivative term and so does not reach zero at the
  ## optimum. The criterion's own flatness is the test that always applies.
  gtol        = 1e-3,   # max |outer gradient|, the quick way to converge
  max_step    = 2.5,    # cap on |delta theta| per outer step
  max_halve   = 12L,    # step halvings before the outer step is abandoned
  inner_method = "bfgs",  # or "newton" (TMB::newton); see [.efs_inner_bfgs()]
  inner_maxit = 100L,   # TMB::newton iterations
  inner_tol   = 1e-8,   # its gradient tolerance, relative to the data scale
  inner_maxit_bfgs = 1000L,
  inner_reltol = 1e-10,
  ## 0 silent, 1 one line per outer iteration, 2 the inner Newton as well.
  ## `silent = FALSE` sets 1; see [.fit_efs()].
  trace       = 0L
)

#' Is a symmetric sparse matrix positive semi-definite, to round-off?
#'
#' Asked of the *data* Hessian \eqn{H - S}, not of \eqn{H}. The penalty can
#' and does rescue \eqn{H}'s definiteness while leaving \eqn{H - S}
#' indefinite, and it is \eqn{H - S} that Wood & Fasiolo's Theorem 1 needs:
#' without it the multiplicative update is not guaranteed to increase the
#' criterion, and the effective degrees of freedom are not in \eqn{[0, 1]}
#' either (see [edf()]).
#'
#' Positive *semi*-definite is the question, since a data Hessian is
#' routinely singular -- a coefficient the data say nothing about -- and that
#' is not a problem. So the test is a Cholesky of \eqn{H - S + \delta I}:
#' it succeeds exactly when every eigenvalue exceeds \eqn{-\delta}, which is
#' "positive semi-definite up to round-off" with \eqn{\delta} setting the
#' allowance.
#'
#' @keywords internal
.psd_ok <- function(Hd, delta)
  !isFALSE(.is_pd(Hd + Matrix::Diagonal(ncol(Hd), delta)))

#' Repair an indefinite data Hessian
#'
#' Two routes, and which one is used is a question of what is affordable
#' rather than of what is right. Clamping the negative eigenvalues to zero is
#' the nearest positive semi-definite matrix in the Frobenius norm, and it is
#' what mgcv does, but it is dense and \eqn{O(q^3)}. Inflating the diagonal
#' until a Cholesky succeeds keeps the sparsity that makes a large model
#' affordable, at the price of also moving the directions that were fine.
#'
#' Either way the fit is no longer at the criterion it set out to optimise,
#' so the number of repaired iterations is counted and reported.
#'
#' @param Hd The data Hessian, symmetric sparse.
#' @param delta Round-off allowance, as in [.psd_ok()].
#' @return The repaired matrix, or `NULL` if even the ridge could not
#'   rescue it.
#' @keywords internal
.psd_repair <- function(Hd, delta) {
  if (ncol(Hd) <= .efs_dense_max) {
    e <- tryCatch(eigen(as.matrix(Hd), symmetric = TRUE),
                  error = function(e) NULL)
    if (!is.null(e)) {
      v <- pmax(e$values, 0)
      M <- e$vectors %*% (v * t(e$vectors))
      return(as(as(Matrix::Matrix((M + t(M)) / 2, sparse = TRUE),
                   "symmetricMatrix"), "CsparseMatrix"))
    }
  }
  r <- delta
  for (i in seq_len(40L)) {
    cand <- Hd + Matrix::Diagonal(ncol(Hd), r)
    if (!isFALSE(.is_pd(cand))) return(cand)
    r <- r * 4
  }
  NULL
}

#' Cholesky of a penalized Hessian, with a ridge if it is needed
#'
#' Even a positive semi-definite data Hessian can leave `H = H_data + S`
#' singular: the penalty is zero on the unpenalized coefficients, so a
#' direction the data do not identify is not rescued by it. That is a model
#' the user should hear about rather than a numerical accident, but it also
#' has to be got past to say anything at all, so a ridge is added and
#' recorded.
#'
#' @return A `CHMfactor`, or `NULL` if no ridge sufficed.
#' @keywords internal
.efs_chol <- function(H, delta) {
  r <- 0
  for (i in seq_len(40L)) {
    M <- if (r > 0) H + Matrix::Diagonal(ncol(H), r) else H
    L <- tryCatch(Matrix::Cholesky(M, LDL = FALSE, perm = TRUE),
                  error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(L)) return(L)
    r <- if (r > 0) r * 4 else delta
  }
  NULL
}

#' Inner solve: the coefficients at fixed smoothing parameters
#'
#' Maximises the penalized log-likelihood over `c = (beta, b)` with
#' `log_sigma` held where the caller left it. Both are solved for, and both
#' are then integrated out by the criterion: `"aREML"` is the REML criterion,
#' so the mean-structure coefficients are random here exactly as they are
#' under `"REML"`.
#'
#' A dispatcher over two solvers, neither of them written here:
#' [.efs_inner_bfgs()] is the default and [.efs_inner_newton()] is
#' [TMB::newton()]. Both of those topics record what was measured and why the
#' default is what it is.
#'
#' @param p Coefficient vector. The smoothing parameters are not in it: they
#'   reach the tape as data, and the caller has already put them where it
#'   wants them.
#' @return `list(p, f, g, converged, iter)`.
#' @keywords internal
.efs_inner <- function(p, fn, gr, hfun, ctl, tag = "") {
  switch(ctl$inner_method,
    bfgs   = .efs_inner_bfgs(p, fn, gr, ctl, tag),
    newton = .efs_inner_newton(p, fn, gr, hfun, ctl, tag),
    stop("`inner_method` must be \"bfgs\" or \"newton\"", call. = FALSE))
}

#' Inner solve by TMB's own Newton
#'
#' [TMB::newton()] is the inner solver the Laplace approximation itself uses.
#' This engine had a hand-rolled imitation of it for a while; that is gone,
#' and this is what `inner_method = "newton"` now means. Two things the real
#' one does that the imitation did not.
#'
#' It factorises **once symbolically** and updates the factor in place at
#' every step (`updateCholesky`), where the hand-rolled loop builds a fresh
#' `Matrix::Cholesky()` from scratch each time. On a problem whose Hessian
#' keeps its sparsity pattern -- which is every problem here, since the
#' pattern comes from the design -- that is most of the per-step cost.
#'
#' And `smartsearch` regularises adaptively: it carries a parameter in
#' \eqn{[0, 1]} interpolating between the Newton step and a gradient step,
#' raising it when the Cholesky succeeds and the objective falls and lowering
#' it when either fails. That is a trust region in all but name, and it is the
#' right answer to an indefinite Hessian. The hand-rolled version instead
#' escalated a ridge by factors of four until a Cholesky happened to succeed,
#' which gives a step of arbitrary length -- the defect that made a
#' four-parameter Box-Cox fit crawl.
#'
#' `TMB::newton()` signals failure by `stop()`, including a deliberate
#' "Newton drop out" after repeated failed attempts, so the call is wrapped:
#' a failure here is a fact about this candidate's smoothing parameters for
#' the outer step control to act on, not an error for the user.
#'
#' @keywords internal
.efs_inner_newton <- function(p, fn, gr, hfun, ctl, tag = "") {
  say <- function(...) if (isTRUE(ctl$trace >= 2L)) message(...)
  grv <- function(q) as.vector(gr(q))
  hev <- function(q) .sym_sparse(hfun(q))
  o <- tryCatch(TMB::newton(p, fn, grv, hev, trace = 0, silent = TRUE,
                            maxit = ctl$inner_maxit,
                            tol = ctl$inner_tol * max(abs(fn(p)), 1)),
                error = function(e) NULL)
  if (is.null(o)) {
    say(sprintf("   inner%s  TMB::newton did not reach a mode", tag))
    return(list(p = p, f = fn(p), g = NA_real_, converged = FALSE, iter = 0L))
  }
  gmax <- max(abs(o$gradient))
  say(sprintf("   inner%s  tmb %d iterations  f = %.6f  max|g| = %.3g",
              tag, o$iterations, o$value, gmax))
  list(p = unname(o$par), f = o$value, g = gmax, converged = TRUE,
       iter = o$iterations)
}

#' Inner solve by BFGS
#'
#' The default, and the alternative to [.efs_inner_newton()]. BFGS keeps a
#' positive definite inverse-Hessian approximation by construction, so its
#' direction always has a sensible scale, and `optim`'s line search enforces a
#' curvature condition as well as sufficient decrease. Crucially it never
#' forms the Hessian during the solve, so an inner problem whose curvature is
#' indefinite costs it nothing: that only matters at the end, where
#' \eqn{\log|H|} is needed once.
#'
#' @section Why this is the default rather than Newton:
#' Not speed -- Newton wins on speed everywhere it works. It is that the
#' four-parameter families this engine exists for have an inner Hessian that
#' is indefinite by the family's own parameterisation, not by the starting
#' values (see [.inner_indefinite()]), and a Newton method's response to
#' negative curvature is a correction that assumes the curvature is local.
#'
#' On `bcpe` over the full `film90`, four smooths: BFGS takes 42 seconds and
#' converges in 12 outer iterations, matching `LaMa::qreml()` on the same
#' model to 0.4 on the log-likelihood and 0.2 on the total EDF.
#' [TMB::newton()] returns in 138 seconds **reporting success**, at a
#' criterion 4500 nats worse, with every smooth collapsed onto its null space
#' (EDF 1, 1, 1, 1). That is not a tolerance artefact: 1e-8 and 1e-10 give the
#' same answer. Its `smartsearch` regularises toward a gradient step when the
#' Hessian is not positive definite, the mode it settles on is not the
#' penalized likelihood's, and the outer step reads that mode as calling for
#' infinite smoothing on the very first iteration.
#'
#' A wrong answer delivered confidently is worse than a slow one, so the
#' robust solver is the default and the fast one is asked for.
#'
#' @section What it costs where Newton would do:
#' A dense \eqn{p \times p} approximation against a sparse Hessian with a
#' factorisation updated in place. To the same criterion and the same EDF in
#' every row:
#'
#' ```
#'                      newton    bfgs
#'   gaussian, 2 smooths  0.51s   1.26s
#'   gamma                0.43s   0.62s
#'   mrf, 196 regions     0.13s   0.34s
#'   spde, 267 nodes      1.46s   4.06s
#' ```
#'
#' So `inner_method = "newton"` is worth reaching for when the coefficient
#' vector is large and the family is ordinary, and is a trap when the family
#' is not.
#'
#' `dev/NOTES-inner-method.md` measured something else again -- TMB's
#' `inner.method = "BFGS"` driven from inside the Laplace approximation's own
#' machinery, where the criterion needs the inner mode to high accuracy at
#' every one of many calls, and where warm starts made it drift. It says
#' nothing about a single solve per outer iteration whose result the outer
#' step control then checks.
#'
#' @keywords internal
.efs_inner_bfgs <- function(p, fn, gr, ctl, tag = "") {
  say <- function(...) if (isTRUE(ctl$trace >= 2L)) message(...)
  grv <- function(q) as.vector(gr(q))
  o <- tryCatch(stats::optim(p, fn, grv, method = "BFGS",
                             control = list(maxit = ctl$inner_maxit_bfgs,
                                            reltol = ctl$inner_reltol)),
                error = function(e) NULL)
  if (is.null(o)) return(list(p = p, f = fn(p), g = NA_real_,
                              converged = FALSE, iter = 0L))
  gmax <- max(abs(grv(o$par)))
  say(sprintf("   inner%s  bfgs %d gradient calls  f = %.6f  max|g| = %.3g",
              tag, o$counts[2L], o$value, gmax))
  list(p = o$par, f = o$value, g = gmax,
       converged = o$convergence == 0L, iter = unname(o$counts[1L]))
}

#' The exponent map from free smoothing parameters to penalty weights
#'
#' Collects every (block, penalty matrix) pair in the design into one table,
#' and records how each pair's exponent \eqn{a_{ki}} depends on the *free*
#' parameter vector -- free meaning after `id`-tying has collapsed several
#' entries of `log_sigma` into one, which [.build_design()] records in
#' `sig_group`.
#'
#' @return `list(pairs, Jac, group, nfree)`: one row of `pairs` per (block,
#'   penalty) pair, `Jac[row, l] = da_row / drho_l`, and `group[m]` the free
#'   parameter that entry `m` of `log_sigma` belongs to.
#' @keywords internal
.efs_pairs <- function(design) {
  grp <- as.integer(design$sig_group)
  nfree <- design$nsigma_free
  pairs <- list(); rows <- list()
  for (k in seq_along(design$blocks)) {
    bl <- design$blocks[[k]]
    pen <- .block_penalties(bl)
    for (i in seq_along(pen$S)) {
      r <- numeric(nfree)
      for (j in seq_len(bl$ntheta))
        r[grp[bl$theta_idx[j]]] <- r[grp[bl$theta_idx[j]]] + pen$A[i, j]
      pairs[[length(pairs) + 1L]] <- list(block = k, S = pen$S[[i]],
                                          tri = .tri_S(pen$S[[i]]),
                                          nS = length(pen$S))
      rows[[length(rows) + 1L]] <- r
    }
  }
  Jac <- if (length(rows)) do.call(rbind, rows) else matrix(0, 0L, nfree)
  ## Which free parameters have penalty weights that can move one at a time,
  ## and which are tied to others through an `L` matrix -- the two get
  ## different steps, for the reason in [.efs_step()]. A parameter is tied if
  ## any of its penalty matrices shares a block with another, or if its
  ## exponent coefficients differ; and tying spreads, since an `id` can put a
  ## tied parameter and a free one on the same penalty matrix.
  tied <- logical(nfree)
  if (nrow(Jac)) {
    ## A penalty matrix whose weight is its block's whole penalty, and which
    ## no other free parameter touches, can be moved on its own.
    alone  <- rowSums(Jac != 0) == 1L
    single <- vapply(pairs, function(z) z$nS == 1L, NA)
    for (l in seq_len(nfree)) {
      ii <- which(Jac[, l] != 0)
      if (!length(ii)) next
      tied[l] <- !(all(single[ii]) && all(alone[ii]) &&
                   all(Jac[ii, l] == Jac[ii[1L], l]))
    }
    ## and tying spreads along shared penalty matrices
    while (any(tied)) {
      rc <- rowSums(Jac[, tied, drop = FALSE] != 0) > 0
      grown <- colSums(Jac[rc, , drop = FALSE] != 0) > 0
      if (identical(unname(grown), unname(tied))) break
      tied <- grown
    }
  }
  ## Which blocks combine several penalty matrices, so that [.efs_step()] can
  ## ask without rebuilding the penalty list -- which, for an "iid" block,
  ## means materialising a q by q sparse identity it will not use.
  multi <- logical(length(design$blocks))
  for (pr in pairs) if (pr$nS > 1L) multi[pr$block] <- TRUE
  list(pairs = pairs, Jac = Jac, group = grp, nfree = nfree, tied = tied,
       multi = multi)
}

#' One EFS step: the numerators, denominators, gradient and proposal
#'
#' For penalty matrix \eqn{S_{ki}} with weight \eqn{\lambda_{ki}} in block
#' \eqn{k}, at the inner mode,
#' \deqn{n_{ki} = tr(Q_k^{-1}\lambda_{ki}S_{ki}) - tr(H^{-1}\lambda_{ki}S_{ki}),
#'       \quad d_{ki} = b_k'(\lambda_{ki}S_{ki})b_k,}
#' and \eqn{\partial V/\partial a_{ki} = (d_{ki} - n_{ki})/2}. The first trace
#' is \eqn{q_k} whenever the block has a single penalty matrix, because then
#' \eqn{Q_k = \lambda_{ki}S_{ki}}; only a block combining several needs a
#' solve. For an `"iid"` block the whole thing collapses to
#' \eqn{n_k = q_k - \lambda_k \sum_j (H^{-1})_{jj}}, which is the block's
#' penalized effective degrees of freedom, and \eqn{d_k = \lambda_k||b_k||^2}.
#'
#' @section The step:
#' Wood & Fasiolo's update is multiplicative in the penalty weights,
#' \eqn{\lambda \leftarrow \lambda\, n/d}, and that is what is used wherever
#' the weights are free to move one at a time: a free parameter all of whose
#' penalty matrices carry the same exponent coefficient. That is every
#' ordinary smooth and every `id` group of them, coefficient \eqn{-2}, with
#' the numerators and denominators pooled over the group. It is exact, and
#' always downhill -- \eqn{sign(\Delta\theta) = -sign(\partial V/\partial
#' \theta)} follows directly.
#'
#' A block that combines several penalty matrices through an `L` matrix -- an
#' SPDE, and a `te()` or an adaptive smooth when those arrive -- does not have
#' free weights: \eqn{\lambda_i = \exp(L_{i\cdot}\theta)} ties them together.
#' Applying the multiplicative update weight by weight and projecting the
#' result onto \eqn{\theta} looks reasonable and is wrong: it has a fixed
#' point where the projected log-ratios vanish, which is not where the
#' gradient vanishes, and the iteration converges neatly to a point that is
#' not stationary. On an SPDE field it stalled 1.3 nats short with
#' \eqn{|\partial V/\partial\theta| = 3.9}.
#'
#' What those parameters get instead is a step in the metric the
#' multiplicative update induces. Linearising \eqn{\log(n/d)} about
#' \eqn{n = d} gives \eqn{2(n-d)/(n+d)}, so the exact update is, to first
#' order, \eqn{-M^{-1}g} with
#' \deqn{M = L' \, diag((n_i + d_i)/4) \, L.}
#' `M` is positive definite, so the step is always downhill; it reproduces the
#' multiplicative update exactly to first order on the free-weight case; and
#' its fixed point is \eqn{g = 0}, which is the one wanted. Far from the
#' optimum the logarithmic form takes larger, better steps, which is why the
#' free-weight case keeps it rather than using the metric throughout.
#'
#' @param fac Cholesky factor of the penalized Hessian, from the criterion
#'   evaluation that produced this point -- the same factorisation, not a
#'   second one.
#' @param np Its dimension.
#' @param rows Positions of each block's coefficients within it.
#' @return `list(num, den, grad, step)`.
#' @keywords internal
.efs_step <- function(design, pmap, theta, bcoef, fac, np, rows, ctl) {
  npair <- length(pmap$pairs)
  rho <- .efs_free(theta, pmap)
  num <- den <- numeric(npair)
  if (!npair) return(list(num = num, den = den, grad = numeric(pmap$nfree),
                          step = numeric(pmap$nfree)))
  ## The block inverse of H, one block at a time: q_k right-hand sides against
  ## the factor already in hand. Solving for the whole of H^-1 at once would
  ## be one call instead of several but would hold an n by nb dense matrix,
  ## which is the thing this engine exists to avoid.
  Vk <- vector("list", length(design$blocks))
  Qk <- vector("list", length(design$blocks))
  for (k in seq_along(design$blocks)) {
    bl <- design$blocks[[k]]
    r <- rows[[k]]
    E <- Matrix::sparseMatrix(i = r, j = seq_along(r), x = 1,
                              dims = c(np, length(r)))
    Vk[[k]] <- as.matrix(Matrix::solve(fac, E, system = "A"))[r, , drop = FALSE]
    ## Only a block that combines several penalty matrices needs its own
    ## precision: with one, tr(Q^-1 lambda S) is q_k and no solve arises.
    if (pmap$multi[k]) Qk[[k]] <- .block_prec(bl, theta[bl$theta_idx])
  }

  lam <- as.vector(exp(pmap$Jac %*% rho))
  for (r in seq_len(npair)) {
    pr <- pmap$pairs[[r]]
    bl <- design$blocks[[pr$block]]
    den[r] <- lam[r] * .quad_S(bcoef[bl$idx], pr$S)
    tr_Q <- if (pr$nS == 1L) bl$q else
      lam[r] * sum(Matrix::diag(Matrix::solve(Qk[[pr$block]],
                                              as(pr$S, "generalMatrix"))))
    num[r] <- tr_Q - lam[r] * .tr_VS(Vk[[pr$block]], pr$tri)
  }

  grad <- as.vector(crossprod(pmap$Jac, den - num)) / 2
  nf <- pmax(num, .efs_floor); df <- pmax(den, .efs_floor)

  step <- numeric(pmap$nfree)
  for (l in which(!pmap$tied)) {
    ii <- which(pmap$Jac[, l] != 0)
    if (!length(ii)) next
    step[l] <- log(sum(nf[ii]) / sum(df[ii])) / pmap$Jac[ii[1L], l]
  }
  if (any(pmap$tied)) {
    J <- pmap$Jac[, pmap$tied, drop = FALSE]
    M <- crossprod(J, (nf + df) / 4 * J)
    M <- M + diag(1e-8 * max(abs(diag(M)), 1), ncol(M))
    step[pmap$tied] <- -as.vector(solve(M, grad[pmap$tied]))
  }
  big <- abs(step) > ctl$max_step
  if (any(big)) step[big] <- sign(step[big]) * ctl$max_step
  list(num = num, den = den, grad = grad, step = step)
}

#' The free smoothing parameters on the scale [edf()] reports them
#'
#' So that a progress line and the fitted summary show the same numbers. A
#' block parameterised by a log standard deviation reports
#' \eqn{\lambda = \sigma^{-2}}; one with its own parameterisation -- an SPDE's
#' range, say -- reports \eqn{\exp(\theta)}, since there is no single variance
#' to invert. Which of the two a free parameter gets is fixed by the design,
#' so it is worked out once and reused.
#'
#' @return A character vector, one entry per free parameter, naming the scale.
#' @keywords internal
.efs_sp_kind <- function(design, pmap) {
  vapply(seq_len(pmap$nfree), function(l) {
    m <- match(l, pmap$group)
    for (bl in design$blocks)
      if (m %in% bl$theta_idx)
        return(if (identical(bl$theta_names, "sd")) "sd" else "raw")
    "raw"
  }, "")
}

## `@rdname` takes the file name, which roxygen munges for a leading dot.
#' @rdname dot-efs_sp_kind
#' @keywords internal
.efs_sp <- function(rho, kind) ifelse(kind == "sd", exp(-2 * rho), exp(rho))

#' Between `log_sigma` and the free smoothing parameters
#'
#' `id`-tying makes several entries of `log_sigma` one estimated value, which
#' [.build_design()] records in `sig_group`. The engine optimises over the
#' free vector and expands it back before every evaluation; entries in a group
#' are equal by construction, so collapsing is just a lookup of the first.
#'
#' @keywords internal
.efs_free <- function(theta, pmap) theta[match(seq_len(pmap$nfree), pmap$group)]

## `@rdname` takes the file name, which roxygen munges for a leading dot.
#' @rdname dot-efs_free
#' @keywords internal
.efs_expand <- function(rho, pmap) rho[pmap$group]

#' Extended Fellner-Schall engine
#'
#' Nothing is declared `random`, so [RTMB::MakeADFun()] builds the plain joint
#' objective and none of the Laplace machinery. The engine then alternates
#'
#' 1. an inner Newton solve for the coefficients at fixed smoothing
#'    parameters ([.efs_inner()]), and
#' 2. an outer multiplicative update of the smoothing parameters
#'    ([.efs_step()]), safeguarded by halving the step in the log parameters
#'    until the criterion actually falls.
#'
#' The objective is built once. `log_sigma` stays an ordinary parameter of it
#' and is simply held fixed for the inner solve, rather than being mapped out
#' -- mapping bakes the value into the tape, and the tape would then have to
#' be rebuilt every outer iteration.
#'
#' @section The criterion:
#' The joint negative log-likelihood already carries the Gaussian priors with
#' their normalising constants, so
#' \deqn{V(\theta) = f(\hat c, \theta) + \tfrac12\log|H| -
#'       \tfrac{p}{2}\log 2\pi}
#' over all `p` of `(beta, b)` -- the same set `method = "REML"` puts in
#' `random =`, and therefore the same criterion, which is why the two report
#' comparable numbers and why this one is called approximate *REML*.
#'
#' There is deliberately no `"ML"` counterpart. It would be approximate twice
#' over -- the Fellner-Schall gradient drops a term, and profiling the
#' unpenalized coefficients at the penalized likelihood's mode rather than at
#' the maximiser of \eqn{V} would drop another -- and it measured ten times
#' further from the Laplace engine's answer than `"aREML"` does from
#' `"REML"`'s. Approximate REML is a criterion people ask for; approximate ML
#' with a second approximation inside it is not.
#'
#' @section What is dropped, and what it costs:
#' \eqn{dV/d\theta} has a term in \eqn{dH/d\hat c \cdot d\hat c/d\theta},
#' which needs third derivatives of the log-likelihood in the coefficients.
#' EFS drops it. The reported gradient is therefore the Fellner-Schall
#' gradient, exact for a Gaussian likelihood and approximate otherwise, and
#' `max_grad` should be read as a diagnostic rather than as a stationarity
#' certificate. The *criterion* is not approximated: the step control
#' evaluates \eqn{V} itself, by re-solving the inner problem and
#' re-factorising \eqn{H}.
#'
#' @section Positive definiteness:
#' Wood & Fasiolo's Theorem 1 guarantees the update increases the criterion
#' only when the data Hessian \eqn{H - S} is positive semi-definite, and it is
#' also what puts the effective degrees of freedom in \eqn{[0, 1]}. It is
#' checked and, if necessary, repaired at every outer iteration; see
#' [.psd_repair()]. A fit that needed repairs says so.
#'
#' @param nll The joint objective, from [.make_nll()].
#' @param pars Starting values, from [.init_pars()].
#' @param control Overrides for [.efs_defaults].
#' @keywords internal
.fit_efs <- function(nll, pars, design, family, silent = TRUE,
                     control = list()) {
  ctl <- utils::modifyList(.efs_defaults, control)
  ## `silent = FALSE` is what a user reaches for to see a fit work, and on
  ## this engine there is no `MakeADFun` for it to reach. So it means the same
  ## thing here as it does there -- show me what is happening -- while an
  ## explicit `control$trace` still wins, in either direction.
  ##
  ## A level rather than a flag, because on a slow family the outer iterations
  ## are not where the time goes: one inner Newton solve on a four-parameter
  ## Box-Cox fit can take longer than every outer step of a Gaussian model put
  ## together, and a progress report that says nothing until the first one
  ## returns is no use on exactly the fits that need one.
  ctl$trace <- if (is.null(control$trace)) as.integer(!isTRUE(silent))
               else if (is.logical(ctl$trace)) as.integer(isTRUE(ctl$trace))
               else as.integer(ctl$trace)
  nbeta <- design$nbeta; nb <- design$nb
  ib <- nbeta + seq_len(nb)
  ## Every coefficient is integrated out -- "aREML" is the REML criterion,
  ## which is what `method = "REML"` puts in `random =` too -- so the matrix
  ## the log determinant needs is the whole penalized Hessian, with no
  ## sub-block to take.
  np <- nbeta + nb

  ## The smoothing parameters are held outside the tape and read back into it
  ## through [RTMB::DataEval()], so that the taped function is of the
  ## coefficients alone. Two things come of that. The tape never has to
  ## differentiate the priors' log determinants, which is what makes a
  ## second-order sparse Jacobian possible at all: a `dgmrf` block's
  ## log|Q(theta)| already has an inverse subset in its first derivative, and
  ## RTMB declines to take a sparse second one. And the smoothing parameters
  ## can still be moved between outer iterations without re-taping, which
  ## mapping them out with `MakeADFun(map = )` would not allow -- so the
  ## objective really is built once.
  ##
  ## `MakeTape` rather than `MakeADFun` because a plain tape is all this engine
  ## wants: none of the parameter mapping, reporting or random-effect
  ## apparatus is used, and a tape is what `jacfun()` composes. `force.update()`
  ## is belt and braces -- an un-reordered tape re-reads its `DataEval` nodes
  ## on every pass as things stand -- but it is the documented way to say the
  ## data has moved, and each derived tape carries its own copy of the node, so
  ## each is told separately.
  ##
  ## A model with no penalized terms has nothing to move and gets no node.
  tenv <- new.env(parent = emptyenv())
  tenv$theta <- pars$log_sigma
  fc <- if (design$nsigma)
    function(x) nll(list(beta = x[seq_len(nbeta)], b = x[ib],
                         log_sigma = RTMB::DataEval(function() tenv$theta)))
  else
    function(x) nll(list(beta = x[seq_len(nbeta)], b = x[ib],
                         log_sigma = numeric(0)))
  Fv <- RTMB::MakeTape(fc, c(pars$beta, pars$b))
  Fg <- Fv$jacfun()
  Fh <- Fg$jacfun(sparse = TRUE)
  set_theta <- if (design$nsigma) function(th) {
    tenv$theta <- th
    Fv$force.update(); Fg$force.update(); Fh$force.update()
    invisible(th)
  } else function(th) invisible(th)

  p <- c(pars$beta, pars$b)
  fn <- function(q) tryCatch(suppressWarnings(Fv(q)),
                             error = function(e) NA_real_)

  if (!is.finite(fn(p))) {
    g0 <- tryCatch(as.vector(Fg(p)), error = function(e) NULL)
    msg <- if (!is.null(g0))
      .flat_start(g0, fn, p, seq_along(p) <= nbeta, design) else NULL
    stop("the penalized objective is not finite at the starting values. ",
         .or_else(msg, paste0("Check that the response is in the support of ",
                              "family '", family$family, "', and consider ",
                              "passing start = list(beta = ...), or a ",
                              "different sigma_frac.")), call. = FALSE)
  }

  pmap <- .efs_pairs(design)
  rho <- .efs_free(pars$log_sigma, pmap)
  ## .penalty_matrix() reads only the dimension and the index vectors from
  ## this, none of which move, so it is built once rather than per evaluation.
  ph <- list(H = Matrix::Diagonal(np), i_beta = seq_len(nbeta), i_b = ib)
  ## Positions of each block's coefficients within the penalized Hessian; `b`
  ## follows `beta` in the coefficient vector.
  rows <- lapply(design$blocks, function(bl) nbeta + bl$idx)

  ## One evaluation of the criterion: solve the inner problem at these
  ## smoothing parameters, repair the data Hessian if it needs it, and add the
  ## log determinant. Everything the outer step needs comes back with it, so
  ## the Hessian is factorised once per evaluation rather than once per use.
  evaluate <- function(rho, p, tag = "") {
    theta <- set_theta(.efs_expand(rho, pmap))
    inner <- .efs_inner(p, fn, Fg, Fh, ctl, tag)
    if (!is.finite(inner$f)) return(NULL)
    H <- .sym_sparse(Fh(inner$p))
    S <- .sym_sparse(.penalty_matrix(design, theta, ph))
    delta <- 1e-8 * max(abs(Matrix::diag(H)), 1)
    Hd <- H - S
    fixed <- FALSE
    if (!.psd_ok(Hd, delta)) {
      Hd <- .psd_repair(Hd, delta)
      if (is.null(Hd)) return(NULL)
      H <- Hd + S
      fixed <- TRUE
    }
    fac <- .efs_chol(H, delta)
    if (is.null(fac)) return(NULL)
    V <- inner$f + sum(log(Matrix::diag(fac))) / 2 - np / 2 * log(2 * pi)
    if (!is.finite(V)) return(NULL)
    list(p = inner$p, theta = theta, f = inner$f, V = V, H = H,
         fac = fac, repaired = fixed)
  }

  ## Progress. Announced *before* the first inner solve rather than after it:
  ## that solve is the longest single step in the fit on a hard family, and a
  ## trace that begins with the criterion has already made the user wait
  ## through the part they most wanted to watch.
  say <- function(...) if (isTRUE(ctl$trace >= 1L)) message(...)
  say(sprintf("efs      %d observations, %d + %d coefficients, %d smoothing %s",
              design$n, nbeta, nb, pmap$nfree,
              if (pmap$nfree == 1L) "parameter" else "parameters"))
  say(if (ctl$trace >= 2L) "efs      first inner solve:"
      else "efs      first inner solve (control = list(trace = 2) to watch it)")

  cur <- evaluate(rho, p)
  if (is.null(cur))
    stop("the aREML criterion is not finite at the starting values: the ",
         "inner problem did not reach a usable mode. Try a different ",
         "sigma_frac, or start = list(beta = ...).", call. = FALSE)

  nrep <- as.integer(cur$repaired)
  trace <- list(V = cur$V, theta = list(cur$theta))
  conv <- FALSE
  msg <- "maximum number of outer iterations reached"
  gmax <- NA_real_
  halved <- 0L

  ## Progress. One line per accepted outer step, plus where it started and why
  ## it stopped, so the trajectory can be read rather than just its endpoint.
  ## `sp` is what edf() will report, and is truncated rather than allowed to
  ## wrap: a model with twenty smooths is watched through its criterion, not
  ## its parameter vector.
  kind <- .efs_sp_kind(design, pmap)
  fmt_sp <- function(rho) {
    v <- .efs_sp(rho, kind)
    if (!length(v)) return("")
    paste0("  sp ", paste(signif(utils::head(v, 6L), 4L), collapse = " "),
           if (length(v) > 6L) " ..." else "")
  }
  say(sprintf("efs   0  -REML = %.6f%s   [start]%s", cur$V, fmt_sp(rho),
              if (isTRUE(cur$repaired)) " repaired" else ""))

  ## With no smoothing parameters to estimate the inner solve is the whole
  ## fit, and there is no outer loop to run.
  if (!pmap$nfree) {
    conv <- TRUE; msg <- "no smoothing parameters to estimate"; gmax <- 0
  }

  for (k in seq_len(if (pmap$nfree) ctl$maxit else 0L)) {
    st <- .efs_step(design, pmap, cur$theta, cur$p[ib], cur$fac, np, rows, ctl)
    gmax <- max(abs(st$grad))
    if (gmax < ctl$gtol) {
      conv <- TRUE; msg <- "outer gradient below tolerance"; break
    }

    ## Wood & Fasiolo's safeguard: halve the step in the log smoothing
    ## parameters until the criterion does not increase. The criterion is the
    ## real one, re-solved and re-factorised at the candidate, so a step that
    ## looks good to the dropped-third-derivative gradient and is not gets
    ## caught here.
    a <- 1
    acc <- NULL
    best <- Inf
    for (h in seq_len(ctl$max_halve)) {
      cand <- evaluate(rho + a * st$step, cur$p, sprintf(" %d.%d", k, h))
      if (!is.null(cand)) {
        if (cand$V <= cur$V) { acc <- cand; break }
        best <- min(best, cand$V)
      }
      a <- a / 2
      halved <- halved + 1L
    }
    if (is.null(acc)) {
      ## Either the criterion is flat here, which is convergence, or the step
      ## direction is being defeated by something -- an inner solve that moved
      ## when it should not have, most likely. What tells them apart is how
      ## much the rejected candidates actually cost: a halved step that makes
      ## the criterion worse by less than the convergence tolerance is a step
      ## into the flat bottom, not a failure.
      ##
      ## Not the gradient. EFS drops a third-derivative term from it, so it
      ## does not go to zero even at the optimum -- on a Markov random field
      ## it settles at 0.3 while the criterion is within 0.001 of the
      ## Laplace engine's. An absolute tolerance on it can only ever be a
      ## sufficient condition for convergence, never a necessary one.
      ##
      ## The separation this has to make is wide: a genuinely flat optimum
      ## costs around 1e-6 relative, while the one case that really does
      ## wander -- a data Hessian repaired at every iteration, so that the
      ## criterion is not a fixed function -- was losing 0.1 per step.
      conv <- is.finite(best) && best - cur$V < ctl$stall_tol * (abs(cur$V) + 1)
      msg <- if (conv) "the criterion is flat in the smoothing parameters"
             else "no step in the smoothing parameters improved the criterion"
      break
    }
    dV <- cur$V - acc$V
    rho <- rho + a * st$step
    cur <- acc
    nrep <- nrep + as.integer(cur$repaired)
    trace$V <- c(trace$V, cur$V)
    trace$theta[[length(trace$theta) + 1L]] <- cur$theta
    say(sprintf(
      "efs %3d  -REML = %.6f  dV = %.3g  max|g| = %.3g  step = %.3g%s%s",
                k, cur$V, dV, gmax, a, fmt_sp(rho),
                if (isTRUE(cur$repaired)) "  [repaired]" else ""))
    if (dV < ctl$tol * (abs(cur$V) + ctl$tol)) {
      conv <- TRUE; msg <- "relative change in the criterion below tolerance"
      break
    }
  }

  say(sprintf("efs      %s after %d iteration%s (%s)",
              if (conv) "converged" else "STOPPED", length(trace$V) - 1L,
              if (length(trace$V) == 2L) "" else "s", msg))

  ## Only a repair at the *final* point changes what is reported. Repairs
  ## along the way are routine -- the starting values are often a place where
  ## a four-parameter family has negative curvature, and the iteration walks
  ## away from it -- so those are counted and printed, not warned about.
  if (isTRUE(cur$repaired))
    warning("the data Hessian is not positive semi-definite at the fitted ",
            "values and was repaired, so the reported criterion and the ",
            "effective degrees of freedom are those of the repaired problem ",
            "rather than of the one that was posed. See ?gamRTMB for the ",
            "basis options; a basis with no null space often fixes it.",
            call. = FALSE)
  if (!conv)
    warning("the extended Fellner-Schall iteration did not converge: ", msg,
            ". Treat the smoothing parameters with suspicion.", call. = FALSE)

  list(tapes = list(value = Fv, grad = Fg, hessian = Fh, theta = tenv),
       opt = list(par = rho, objective = cur$V,
                  convergence = if (conv) 0L else 1L,
                  message = msg, iterations = length(trace$V) - 1L,
                  halvings = halved),
       coefficients = list(beta = cur$p[seq_len(nbeta)], b = cur$p[ib]),
       log_sigma = cur$theta,
       H = cur$H, penalized_loglik = -cur$f,
       objective = cur$V, convergence = conv, max_grad = gmax,
       psd_repairs = nrep, psd_repaired_at_mode = isTRUE(cur$repaired),
       efs_trace = trace)
}
