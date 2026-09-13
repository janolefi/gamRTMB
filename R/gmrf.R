## The sparse route. A smooth whose penalty is sparse can skip
## `mgcv::smooth2random()` entirely: instead of rotating the coefficients into
## an iid basis -- which destroys the sparsity that made the penalty cheap in
## the first place -- keep the penalty as it is and give the block a
## `dgmrf(b, 0, Q / sigma^2)` prior. The two are the same model; only the
## parameterisation differs.
##
## The same route is the only one open to a smooth carrying several penalty
## matrices on one coefficient vector -- a `te()`, a `ti()`, an adaptive
## smooth, an SPDE -- since `smooth2random()` needs one variance per penalized
## block and has none to give. Everything below is written against a *list* of
## penalties for that reason; a list of one is the ordinary case.

#' A matrix in the sparse storage the blocks are built from
#' @keywords internal
.spm <- function(M)
  as(as(Matrix::Matrix(M, sparse = TRUE), "generalMatrix"), "CsparseMatrix")

#' Is this smooth's penalty worth keeping sparse?
#'
#' Three ways to qualify, two of them outright.
#'
#' A smooth carrying an `L` matrix qualifies whatever its size. `L` is mgcv's
#' way of saying that several penalty matrices combine into one through fewer
#' smoothing parameters than there are matrices, as
#' \eqn{\lambda = \exp(L\theta)}. `smooth2random` cannot represent that at
#' all -- it needs one variance per penalized block -- so such a smooth has
#' nowhere else to go. The SPDE smooth is the case in hand: three finite
#' element matrices, two parameters.
#'
#' A smooth with several penalty matrices whose supports *overlap* qualifies
#' for the same reason, and overlap is the whole of the test. `t2()` also has
#' several penalties, but by construction each acts on its own disjoint set of
#' coefficients, which is exactly what lets `smooth2random` split it into one
#' iid block per penalty -- so `t2()` keeps the well-travelled route. A
#' `te()`, a `ti()` and an adaptive smooth all penalise the same coefficients
#' several times over and cannot be split, which is why mgcv itself declines
#' them (*"te smooths not useable with gamm4"*, *"Can not convert this smooth
#' class to a random effect"*).
#'
#' A smooth with a single penalty is the discretionary case: it qualifies if
#' that penalty is big and sparse enough that `smooth2random`'s rotation would
#' actually cost something -- a 10-coefficient P-spline penalty is technically
#' banded, but nothing is gained by treating it specially and the
#' well-travelled route is the safer one.
#'
#' @param sm A `smoothCon` object.
#' @param mode `"auto"`, `"always"` or `"never"`.
#' @keywords internal
.use_sparse <- function(sm, mode) {
  if (mode == "never") return(FALSE)
  if (!is.null(sm$L)) return(TRUE)
  if (length(sm$S) > 1L) return(.penalties_overlap(sm$S))
  if (length(sm$S) != 1L) return(FALSE)
  if (mode == "always") return(TRUE)
  ncol(sm$S[[1L]]) >= 50L && mean(sm$S[[1L]] != 0) <= 0.2
}

#' Do several penalty matrices act on a shared coefficient?
#'
#' The one question that separates a `t2()`, which [mgcv::smooth2random()]
#' splits into independent iid blocks, from a `te()` or an adaptive smooth,
#' which it cannot split at all. Answered from the column supports rather than
#' from the smooth's class, so a basis nobody has thought of yet gets the
#' right route for the right reason.
#'
#' @param S A list of penalty matrices, all the same size.
#' @keywords internal
.penalties_overlap <- function(S) {
  seen <- logical(ncol(S[[1L]]))
  for (M in S) {
    used <- as.vector(Matrix::colSums(abs(M)) > 0)
    if (any(seen & used)) return(TRUE)
    seen <- seen | used
  }
  FALSE
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
#' @section Several penalties at once:
#' Everything above survives the step to a `te()`, a `ti()` or an adaptive
#' smooth unchanged, and one fact is why. For \eqn{\lambda_i > 0} and
#' \eqn{S_i \succeq 0},
#' \deqn{v'S_\lambda v = 0 \iff S_i v = 0 \ \forall i,}
#' so the null space of \eqn{S_\lambda = \sum_i \lambda_i S_i} is the
#' *intersection* of the individual null spaces and **does not depend on
#' \eqn{\lambda}**. It is a structural constant, and mgcv already knows it as
#' `null.space.dim`.
#'
#' That matters twice over. The columns of `N` lie in every \eqn{null(S_i)},
#' so \eqn{S_i N = 0} holds for each penalty separately and one
#' reparameterisation serves all of them at once. And \eqn{Q(\lambda) =
#' S_\lambda[keep, keep]} is positive definite for every positive
#' \eqn{\lambda}: if \eqn{\tilde v} (padded with zeros on `drop`) has
#' \eqn{\tilde v'S_\lambda\tilde v = 0} then \eqn{\tilde v = Nc}, and
#' vanishing on `drop` forces \eqn{N_{drop}c = 0} and so \eqn{c = 0}.
#'
#' So [RTMB::dgmrf()] is handed an honestly positive definite matrix and its
#' Cholesky log determinant is the quantity wanted: with \eqn{T = [P\ N]},
#' \deqn{\log\det Q(\lambda) = \log|S_\lambda|_+ + 2\log|\det T|,}
#' and the second term depends on `N` alone, not on \eqn{\lambda}. The
#' generalised determinant and its notoriously delicate stable evaluation
#' (Wood, 2011, Appendix B) never have to be computed.
#'
#' @section Finding the null space:
#' `k` comes from `mgcv`, which knows it by construction, and is taken as
#' authoritative -- a numerical rank test is exactly what this is avoiding.
#' Only a basis has to be found, and there is one case worth doing cheaply,
#' because it is every `bs = "mrf"` and every first-order random walk: when
#' the penalty's rows sum to zero the constant lies in the null space, and if
#' the graph then has exactly `k` connected components their indicators *are*
#' the null space. Reading those off the sparsity pattern costs a traversal
#' where an eigendecomposition would cost `q^3`. When the counts disagree -- a
#' second difference penalty has a linear null direction as well as a constant
#' one -- the general route is taken instead. A Cholesky cannot stand in for
#' this test: CHOLMOD and `chol()` both factor a penalty whose smallest
#' eigenvalue is -8e-18 without complaint.
#'
#' With several penalties the sum they are searched in is **balanced**,
#' \eqn{\sum_i S_i/\|S_i\|_F}, rather than taken at \eqn{\lambda = 1}. Any
#' positive weights give the same null space in exact arithmetic; balanced
#' ones give it in floating point too, since a penalty allowed to dominate
#' rounds the others' null directions away. It is also the only scaling at
#' which a fixed relative tolerance means the same thing for every term.
#'
#' The eigenvalue either side of the cut is then checked rather than assumed.
#' The gap is wide when the structure is what mgcv says it is -- 2.3e-4
#' against 4.9e-17 for a `te()`, 3.6e-5 against 3.3e-17 for an adaptive
#' smooth -- so a narrow one means the penalties are not the ones this
#' reduction was derived for, and saying so beats silently fitting a different
#' model.
#'
#' @param S A list of penalty matrices, all `q` by `q`.
#' @param k Dimension of their common null space, from `smoothCon`'s
#'   `null.space.dim`.
#' @param label The smooth's label, for error messages.
#' @return `list(N, drop)`: a `q` by `k` null-space basis, and the `k` indices
#'   whose coefficients are fixed at zero.
#' @references
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and marginal
#' likelihood estimation of semiparametric generalized linear models.
#' \emph{JRSS-B} 73, 3-36.
#' @keywords internal
.null_space <- function(S, k, label = "this smooth") {
  q <- ncol(S[[1L]])
  if (k <= 0L) return(list(N = matrix(0, q, 0L), drop = integer(0)))
  B <- .balanced_sum(S)
  tol <- 1e-8 * max(abs(B))

  N <- NULL
  if (all(abs(Matrix::rowSums(B)) < tol)) {
    comp <- .components(B, tol)
    if (max(comp) == k) N <- outer(comp, seq_len(k), "==") * 1
  }
  if (is.null(N)) {
    if (q > 2000L)
      stop("the penalty for ", sQuote(label), " is singular in a way that ",
           "needs an eigendecomposition to resolve, and it has ", q,
           " columns. Supply a positive definite penalty (a proper CAR or an ",
           "SPDE precision), or fit with sparse = \"never\".", call. = FALSE)
    e <- eigen(as.matrix(B), symmetric = TRUE)
    .check_null_gap(e$values, k, label)
    N <- e$vectors[, seq(q, q - k + 1L), drop = FALSE]
  }
  list(N = N, drop = sort(qr(t(N))$pivot[seq_len(ncol(N))]))
}

#' Several penalties on one comparable scale
#'
#' \eqn{\sum_i S_i/\|S_i\|_F}. Used wherever the penalties have to be looked
#' at together without one of them setting the scale -- see [.null_space()]
#' for the null-space search, and [.penalty_spec()] for the starting values.
#' mgcv's own `scale.penalty` already brings them within a factor of two of
#' each other, so this is insurance rather than a correction.
#'
#' @keywords internal
.balanced_sum <- function(S) {
  w <- vapply(S, function(M)
    1 / max(sqrt(sum(M^2)), .Machine$double.eps), numeric(1))
  Reduce(`+`, Map(`*`, S, w))
}

#' Is the null space where mgcv says it is?
#'
#' Both sides of the cut are checked against the same relative tolerance: the
#' smallest eigenvalue that is being *kept* must be clearly nonzero, and the
#' largest one being dropped clearly zero. See [.null_space()] for why the
#' answer is not simply read off a rank test instead.
#'
#' @param ev Eigenvalues of the balanced penalty sum, decreasing.
#' @keywords internal
.check_null_gap <- function(ev, k, label) {
  q <- length(ev); tol <- 1e-8 * max(ev)
  bad <- function(what, v)
    stop("the penalties of ", sQuote(label), " do not have the ",
         k, "-dimensional null space mgcv reports for them: ", what,
         " eigenvalue of the balanced penalty sum is ", format(v, digits = 3),
         " against a tolerance of ", format(tol, digits = 3),
         ", so the reduction to a positive definite precision is not safe ",
         "here. Fit this term with sparse = \"never\" if it has a single ",
         "penalty, and otherwise report it.", call. = FALSE)
  if (q > k && ev[q - k] <= tol) bad("the smallest retained", ev[q - k])
  if (ev[q - k + 1L] >= tol) bad("the largest discarded", ev[q - k + 1L])
  invisible(TRUE)
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
#' One path, whatever the smooth brought with it. Its penalties are reduced by
#' their common null space ([.null_space()]), which is a no-op for a proper
#' precision, a corner constraint for an intrinsic field, and exactly the same
#' corner constraint for a `te()` or an adaptive smooth -- since \eqn{S_iN = 0}
#' holds for every penalty at once, one reduction serves all of them.
#'
#' Null-space columns whose contribution the parameter's intercept already
#' covers are left out, which is what the sum-to-zero constraint achieves on
#' the other route. A pivoted QR of `cbind(1, X N)` finds them: for a
#' connected Markov random field the whole null space is the constant, so no
#' free column survives and the intercept carries the level; for a field with
#' an island, one contrast between the two components survives, which is
#' right, because their levels really are separately free; and for a `te()`
#' the main effects survive and the constant does not, which is mgcv's own
#' account of a tensor product's null space.
#'
#' @param sm A `smoothCon` object built with `absorb.cons = FALSE`.
#' @keywords internal
.gmrf_block <- function(sm) {
  q <- ncol(sm$X)
  S <- lapply(sm$S, .spm)
  ns <- .null_space(S, as.integer(sm$null.space.dim), sm$label)
  keep <- setdiff(seq_len(q), ns$drop)
  Sk <- lapply(S, function(M) .sym_sparse(M[keep, keep, drop = FALSE]))

  V <- as.matrix(sm$X %*% ns$N)
  free <- if (!ncol(V)) integer(0) else {
    qq <- qr(cbind(1, V))
    if (qq$rank > 1L) sort(qq$pivot[seq(2L, qq$rank)] - 1L) else integer(0)
  }
  Tm <- cbind(Matrix::sparseMatrix(i = keep, j = seq_along(keep), x = 1,
                                   dims = c(q, length(keep))),
              Matrix::Matrix(ns$N[, free, drop = FALSE], sparse = TRUE))
  list(Xr = list(.spm(sm$X[, keep, drop = FALSE])),
       Xf = V[, free, drop = FALSE], Tmap = Tm,
       intrinsic = length(ns$drop) > 0L,
       spec = list(.penalty_spec(sm, Sk)))
}

#' How a block's parameters become its penalty
#'
#' The prior on every penalized block is
#' \deqn{Q_k(\theta) = \sum_i \exp(A_{i\cdot}\theta)\, S_{ki},}
#' and `A` -- a linear map from the parameters actually estimated to the log
#' penalty weights -- is the only thing that differs between kinds of smooth.
#' mgcv writes it `L`, for the blocks it has a name for. This is where a
#' smooth's penalties are turned into that description; [.block_prec()] and
#' [.block_penalties()] are the two readings of it, on the tape and off.
#'
#' | penalties | \eqn{S_i} | `A` | \eqn{\theta} |
#' | --- | --- | --- | --- |
#' | one, no `L` | the penalty itself | \eqn{-2} | \eqn{\log\sigma} |
#' | several, with `L` | \eqn{C, G_1, G_2} | mgcv's `L` | \eqn{(\log\sigma, \log\rho)} |
#' | several, no `L` | \eqn{S_1 \ldots S_m} | \eqn{I_m} | \eqn{\log\lambda} |
#'
#' All three are the same `"multi"` block. A single sparse penalty scaled by a
#' variance is not a different kind of object from several combined through an
#' `L`; it is that object with \eqn{L = [-2]} and one matrix, and giving it a
#' kind of its own bought a branch in every function that switches on one.
#' Only `"iid"` is genuinely apart, and only because its penalty is the
#' identity: it needs no matrix at all, and gets `dnorm` rather than `dgmrf`.
#'
#' @section Reporting and starting values:
#' `sp_pow` records the diagonal of `A` when `A` is diagonal, so that one
#' penalty weight belongs to one parameter and can be reported as
#' \eqn{\lambda_i = \exp(sp\_pow_i\,\theta_i)}. An SPDE's `L` is not
#' diagonal -- three matrices, two parameters -- so there is no such
#' \eqn{\lambda} to report and its own named parameters are shown instead;
#' `sp_pow` is `NULL` there and that is what says so.
#'
#' The same field carries the starting values. Scaling a whole block's penalty
#' by \eqn{\sigma^{-2}} is a shift of \eqn{-2/sp\_pow} in \eqn{\theta},
#' which is \eqn{+1} for a log standard deviation and \eqn{-2} for a log
#' penalty weight, so [.init_pars()]'s one rule covers both without knowing
#' which it has. A block whose `L` is not diagonal keeps the absolute
#' `theta.start` its constructor chose.
#'
#' `qscale` is the matching piece: a typical diagonal entry of the penalty at
#' \eqn{\theta = 0}, which is what turns \eqn{\sigma} into the term's implied
#' prior standard deviation. Balanced across penalties, for the reason in
#' [.balanced_sum()].
#'
#' @param sm The `smoothCon` object.
#' @param Sk Its penalties, already reduced to the retained coefficients.
#' @keywords internal
.penalty_spec <- function(sm, Sk) {
  m <- length(Sk)
  if (is.null(sm$L)) {
    ## No `L`. One penalty is a variance to scale it by; several are mgcv's
    ## own convention of one free log penalty weight each, which is `A = I`.
    ## Keeping those free is not only simplest, it is what makes the
    ## Fellner-Schall update exact -- see [.efs_step()].
    L  <- if (m == 1L) matrix(-2, 1L, 1L) else diag(m)
    nm <- if (m == 1L) "sd" else as.character(seq_len(m))
    ## A typical diagonal entry of the penalty at theta = 0. Balanced across
    ## several, for the reason in [.balanced_sum()]; taken as it stands when
    ## there is only one, since then the scale *is* that penalty's own.
    qsc <- sqrt(mean(Matrix::diag(if (m == 1L) Sk[[1L]] else .balanced_sum(Sk))))
  } else {
    L  <- as.matrix(sm$L)
    nm <- if (!is.null(sm$theta.names)) sm$theta.names else
      paste0("theta", seq_len(ncol(L)))
    qsc <- 1
  }
  list(kind = "multi", Smats = Sk, L = L, ntheta = ncol(L), theta_names = nm,
       theta_start = if (!is.null(sm$theta.start)) sm$theta.start else
         numeric(ncol(L)),
       sp_pow = if (nrow(L) == ncol(L) && all(L[row(L) != col(L)] == 0))
         diag(L) else NULL,
       qscale = qsc)
}

#' Build an iid block from a smooth, the `smooth2random` way
#'
#' The other half of [.gmrf_block()], with the same shape of answer, so that
#' [.build_design()] can pick a route and then stop caring which it picked.
#' Here the reparameterisation makes the penalty the identity, so the block
#' needs no precision matrix and exactly one variance.
#'
#' @param sm A `smoothCon` object built with `absorb.cons = TRUE`.
#' @keywords internal
.iid_block <- function(sm) {
  re <- mgcv::smooth2random(sm, "", type = 2)
  list(Xr = lapply(re$rand, as.matrix), Xf = re$Xf,
       Tmap = .reconstruct_map(re), intrinsic = FALSE,
       spec = rep(list(list(kind = "iid", ntheta = 1L, theta_names = "sd",
                            theta_start = 0, sp_pow = -2, qscale = 1)),
                  length(re$rand)))
}

#' A block's smoothing parameters, as reported
#'
#' The one rendering of `sp_pow`; see [.penalty_spec()] for what it records
#' and [.efs_sp_kind()] for the same question asked of a free parameter
#' rather than of a block. A weight that belongs to a single parameter is
#' shown as that weight; anything else is shown as the named parameter it is.
#'
#' @param bl A block from [.build_design()].
#' @param ls The full `log_sigma` vector.
#' @keywords internal
.sp_show <- function(bl, ls) {
  th <- ls[bl$theta_idx]
  if (is.null(bl$sp_pow)) sprintf("%s=%.4g", bl$theta_names, exp(th))
  else sprintf("%.4g", exp(bl$sp_pow * th))
}

#' A penalized block's precision matrix, at given parameters
#'
#' The one place that knows how a block's parameters become a precision, so
#' that the objective and [edf()] cannot drift apart. Written to work both on
#' an RTMB tape and on plain numbers.
#'
#' `"iid"` is the [mgcv::smooth2random()] basis, where the penalty is the
#' identity and no matrix is needed. `"multi"` is every other block: it
#' combines its penalty matrices through mgcv's `L` convention,
#' \eqn{\lambda = \exp(L\theta)} and \eqn{Q = \sum_i \lambda_i S_i}, which
#' is how a Matern SPDE writes
#' \eqn{\tau^2(\kappa^4 C + 2\kappa^2 G_1 + G_2)} with
#' \eqn{\theta = (\log\tau, \log\kappa)}, and -- with `L` the identity --
#' how a `te()`, a `ti()` or an adaptive smooth carries one free
#' \eqn{\log\lambda} per penalty. See [.penalty_spec()].
#'
#' @param bl A block from [.build_design()].
#' @param theta That block's parameters, `bl$ntheta` of them.
#' @keywords internal
.block_prec <- function(bl, theta) {
  if (identical(bl$kind, "iid"))
    return(Matrix::Diagonal(bl$q, exp(-2 * theta[1L])))
  Q <- bl$Smats[[1L]] * exp(sum(bl$L[1L, ] * theta))
  for (i in seq_along(bl$Smats)[-1L])
    Q <- Q + bl$Smats[[i]] * exp(sum(bl$L[i, ] * theta))
  Q
}
