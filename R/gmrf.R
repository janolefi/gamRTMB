## The sparse route. A smooth whose penalty is sparse can skip
## `mgcv::smooth2random()` entirely: instead of rotating the coefficients into
## an iid basis -- which destroys the sparsity that made the penalty cheap in
## the first place -- keep the penalty as it is and give the block a
## `dgmrf(b, 0, Q / sigma^2)` prior. The two are the same model; only the
## parameterisation differs.

#' Is this smooth's penalty worth keeping sparse?
#'
#' Two ways to qualify. A smooth with a single penalty qualifies if that
#' penalty is big and sparse enough that `smooth2random`'s rotation would
#' actually cost something -- a 10-coefficient P-spline penalty is technically
#' banded, but nothing is gained by treating it specially and the
#' well-travelled route is the safer one.
#'
#' A smooth carrying an `L` matrix qualifies outright, whatever its size. `L`
#' is mgcv's way of saying that several penalty matrices combine into one
#' through fewer smoothing parameters than there are matrices, as
#' \eqn{\lambda = \exp(L\theta)}. `smooth2random` cannot represent that at
#' all -- it needs one variance per penalized block -- so such a smooth has
#' nowhere else to go. The SPDE smooth is the case in hand: three finite
#' element matrices, two parameters.
#'
#' @param sm A `smoothCon` object.
#' @param mode `"auto"`, `"always"` or `"never"`.
#' @keywords internal
.use_sparse <- function(sm, mode) {
  if (mode == "never") return(FALSE)
  if (!is.null(sm$L)) return(TRUE)
  if (length(sm$S) != 1L) return(FALSE)
  if (mode == "always") return(TRUE)
  ncol(sm$S[[1L]]) >= 50L && mean(sm$S[[1L]] != 0) <= 0.2
}

#' The null space of a penalty, and which coefficients to drop for it
#'
#' An intrinsic GMRF -- an ICAR field, a random walk -- has a rank-deficient
#' precision: its density depends only on differences between coefficients, so
#' `k` directions are unpenalized. Both routes through the package have to
#' account for them, and they differ in how.
#'
#' `mgcv` absorbs a sum-to-zero constraint into the basis with a dense
#' orthogonal rotation. That is what turns an ICAR penalty from 1% nonzero
#' into 100% nonzero, and it is why a Markov random field over a few hundred
#' regions is slow. Here the penalized part is corner-constrained instead:
#' `k` coefficients are fixed at zero and dropped, which leaves a sparse
#' submatrix, and the null space is handed to `beta` as unpenalized columns
#' exactly as [mgcv::smooth2random()] would.
#'
#' The two are the same model. Write \eqn{b = b_0 + N\gamma} with
#' \eqn{\gamma = (N_{drop})^{-1} b_{drop}}, so that \eqn{b_0} vanishes on the
#' dropped indices. The decomposition is a bijection whenever \eqn{N_{drop}}
#' is invertible -- which is what the pivoted QR below picks the dropped
#' indices to guarantee, and is the familiar fact that deleting one node of a
#' connected graph Laplacian leaves a positive definite matrix. Then
#' \eqn{Xb = X_{keep} b_{0,keep} + (XN)\gamma} and, because \eqn{SN = 0},
#' \eqn{b'Sb = b_{0,keep}' Q b_{0,keep}}. Same span, same penalty, sparse
#' \eqn{Q}.
#'
#' A proper penalty has no null space and gets no constraint: the prior itself
#' identifies every direction. That is the case a Matern/SPDE precision lands
#' in.
#'
#' @section Finding the null space:
#' `k` comes from `mgcv`, which knows it by construction. Only a basis has to
#' be found here, and there is one case worth doing cheaply, because it is
#' every `bs = "mrf"` and every first-order random walk: when the penalty's
#' rows sum to zero the constant lies in the null space, and if the graph then
#' has exactly `k` connected components their indicators *are* the null space.
#' Reading those off the sparsity pattern costs a traversal where an
#' eigendecomposition would cost `q^3`. When the counts disagree -- a second
#' difference penalty has a linear null direction as well as a constant one --
#' the general route is taken instead. A Cholesky cannot stand in for this
#' test: CHOLMOD and `chol()` both factor a penalty whose smallest eigenvalue
#' is -8e-18 without complaint.
#'
#' @param S The penalty matrix.
#' @param k Dimension of its null space, from `smoothCon`'s `null.space.dim`.
#' @return `list(N, drop)`: a `q` by `k` null-space basis, and the `k` indices
#'   whose coefficients are fixed at zero.
#' @keywords internal
.null_space <- function(S, k) {
  q <- ncol(S)
  if (k <= 0L) return(list(N = matrix(0, q, 0L), drop = integer(0)))
  tol <- 1e-8 * max(abs(S))
  N <- NULL
  if (all(abs(Matrix::rowSums(S)) < tol)) {
    comp <- .components(S, tol)
    if (max(comp) == k)
      N <- outer(comp, seq_len(k), "==") * 1
  }
  if (is.null(N)) {
    if (q > 2000L)
      stop("the penalty for this smooth is singular in a way that needs an ",
           "eigendecomposition to resolve, and it has ", q, " columns. ",
           "Supply a positive definite penalty (a proper CAR or an SPDE ",
           "precision), or fit with sparse = \"never\".", call. = FALSE)
    e <- eigen(as.matrix(S), symmetric = TRUE)
    kk <- sum(e$values < 1e-8 * max(e$values))
    if (!kk) return(list(N = matrix(0, q, 0L), drop = integer(0)))
    N <- e$vectors[, seq(q, q - kk + 1L), drop = FALSE]
  }
  list(N = N, drop = sort(qr(t(N))$pivot[seq_len(ncol(N))]))
}

#' Connected components of a penalty's sparsity pattern
#'
#' A breadth-first sweep over the off-diagonal nonzeros, returning a component
#' label per coefficient. Isolated coefficients -- an island region with no
#' neighbours -- come back as components of their own, which is the right
#' answer: their level is free too.
#'
#' @keywords internal
.components <- function(S, tol) {
  q <- ncol(S)
  A <- Matrix::which(abs(S) > tol, arr.ind = TRUE)
  A <- A[A[, 1L] != A[, 2L], , drop = FALSE]
  nbrs <- split(A[, 2L], factor(A[, 1L], levels = seq_len(q)))
  comp <- integer(q); nc <- 0L
  for (i in seq_len(q)) {
    if (comp[i]) next
    nc <- nc + 1L; front <- i; comp[i] <- nc
    while (length(front)) {
      nxt <- unique(unlist(nbrs[front], use.names = FALSE))
      front <- nxt[comp[nxt] == 0L]
      comp[front] <- nc
    }
  }
  comp
}

#' Build a sparse GMRF block from a smooth
#'
#' Returns the same three things the [mgcv::smooth2random()] route returns --
#' a design matrix on the penalized coefficients, an unpenalized null-space
#' matrix for `beta`, and the map back to the smooth's own basis -- so that
#' everything downstream is unchanged, plus the pieces the prior needs.
#'
#' Null-space columns whose contribution the parameter's intercept already
#' covers are left out, which is what the sum-to-zero constraint achieves on
#' the other route. A pivoted QR of `cbind(1, X N)` finds them: for a
#' connected Markov random field the whole null space is the constant, so no
#' free column survives and the intercept carries the level; for a field with
#' an island, one contrast between the two components survives, which is
#' right, because their levels really are separately free.
#'
#' A smooth with an `L` matrix keeps all of its penalty matrices and gets
#' `ncol(L)` parameters instead of one. Its penalty is assumed proper -- an
#' SPDE precision is, for any positive range -- so no constraint arises.
#'
#' @param sm A `smoothCon` object built with `absorb.cons = FALSE`.
#' @keywords internal
.gmrf_block <- function(sm) {
  spm <- function(M) as(as(Matrix::Matrix(M, sparse = TRUE), "generalMatrix"),
                        "CsparseMatrix")
  q <- ncol(sm$X)

  if (!is.null(sm$L)) {
    nm <- if (!is.null(sm$theta.names)) sm$theta.names else
      paste0("theta", seq_len(ncol(sm$L)))
    return(list(Xr = list(spm(sm$X)), Xf = matrix(0, nrow(sm$X), 0L),
                Tmap = Matrix::Diagonal(q), intrinsic = FALSE,
                spec = list(list(kind = "multi", Smats = lapply(sm$S, spm),
                                 L = as.matrix(sm$L), ntheta = ncol(sm$L),
                                 theta_names = nm,
                                 theta_start = if (!is.null(sm$theta.start))
                                   sm$theta.start else rep(0, ncol(sm$L))))))
  }

  S <- Matrix::Matrix(sm$S[[1L]], sparse = TRUE)
  ns <- .null_space(S, as.integer(sm$null.space.dim))
  keep <- setdiff(seq_len(q), ns$drop)
  Q <- as(as(S[keep, keep, drop = FALSE], "symmetricMatrix"), "CsparseMatrix")

  V <- as.matrix(sm$X %*% ns$N)
  free <- if (!ncol(V)) integer(0) else {
    qq <- qr(cbind(1, V))
    if (qq$rank > 1L) sort(qq$pivot[seq(2L, qq$rank)] - 1L) else integer(0)
  }
  Tm <- cbind(Matrix::sparseMatrix(i = keep, j = seq_along(keep), x = 1,
                                   dims = c(q, length(keep))),
              Matrix::Matrix(ns$N[, free, drop = FALSE], sparse = TRUE))
  list(Xr = list(spm(sm$X[, keep, drop = FALSE])),
       Xf = V[, free, drop = FALSE], Tmap = Tm,
       intrinsic = length(ns$drop) > 0L,
       spec = list(list(kind = "gmrf", Q = Q, ntheta = 1L,
                        theta_names = "sd")))
}

#' A penalized block's precision matrix, at given parameters
#'
#' The one place that knows how a block's parameters become a precision, so
#' that the objective and [edf()] cannot drift apart. Written to work both on
#' an RTMB tape and on plain numbers.
#'
#' `"iid"` is the [mgcv::smooth2random()] basis, where the penalty is the
#' identity. `"gmrf"` keeps a fixed sparse precision and scales it. `"multi"`
#' combines several penalty matrices through mgcv's `L` convention,
#' \eqn{\lambda = \exp(L\theta)} and \eqn{Q = \sum_i \lambda_i S_i}, which
#' is how a Matern SPDE writes
#' \eqn{\tau^2(\kappa^4 C + 2\kappa^2 G_1 + G_2)} with
#' \eqn{\theta = (\log\tau, \log\kappa)}.
#'
#' @param bl A block from [.build_design()].
#' @param theta That block's parameters, `bl$ntheta` of them.
#' @keywords internal
.block_prec <- function(bl, theta) {
  switch(bl$kind,
    iid  = Matrix::Diagonal(bl$q, exp(-2 * theta[1L])),
    gmrf = bl$Q * exp(-2 * theta[1L]),
    multi = {
      Q <- bl$Smats[[1L]] * exp(sum(bl$L[1L, ] * theta))
      for (i in seq_along(bl$Smats)[-1L])
        Q <- Q + bl$Smats[[i]] * exp(sum(bl$L[i, ] * theta))
      Q
    })
}
