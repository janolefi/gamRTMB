## A Matern SPDE field as an mgcv smooth. The construction is Lindgren, Rue
## and Lindstrom (2011): a Gaussian field with a Matern covariance is the
## solution of a stochastic PDE, and a finite element approximation on a
## triangulation turns that into a GMRF whose precision is sparse and known in
## closed form. `fmesher` builds the mesh and the finite element matrices;
## nothing here needs INLA.

#' Matern SPDE smooth
#'
#' A smooth term for a spatially continuous Gaussian random field, usable as
#' `s(x, y, bs = "spde", xt = list(mesh = mesh))` in any `gamRTMB` formula and
#' therefore on any parameter of any family.
#'
#' @section Building a mesh:
#' The mesh is yours to build and to pass in, because a good one depends on
#' the domain, the data and the range you expect, and no default can know
#' those. [fmesher::fm_mesh_2d()] takes locations, a maximum triangle edge
#' length inside the domain and in the outer extension, and a cutoff below
#' which nearby points are merged:
#'
#' ```
#' mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.1, 0.3),
#'                             cutoff = 0.05, offset = c(0.1, 0.3))
#' f <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh))),
#'              data = d)
#' ```
#'
#' The outer extension matters: without it the field's variance is inflated at
#' the boundary. A rule of thumb is `max.edge` no larger than a third of the
#' range you expect, and an offset about equal to the range.
#'
#' For one dimension pass a [fmesher::fm_mesh_1d()], or leave `xt` out and a
#' mesh with `k` evenly spaced knots is built over the range of the covariate.
#'
#' @section Parameterisation:
#' The precision is
#' \deqn{Q(\tau, \kappa) = \tau^2(\kappa^4 C + 2\kappa^2 G_1 + G_2)}
#' with \eqn{C}, \eqn{G_1}, \eqn{G_2} the finite element matrices from
#' [fmesher::fm_fem()]. For `alpha = 2` in two dimensions the Matern
#' smoothness is \eqn{\nu = 1}, the range is
#' \eqn{\rho = \sqrt{8\nu}/\kappa = 2\sqrt2/\kappa} and the marginal
#' standard deviation is \eqn{\sigma = 1/(\tau\kappa\sqrt{4\pi})}.
#'
#' The two estimated parameters are \eqn{\log\sigma} and \eqn{\log\rho},
#' not \eqn{\log\tau} and \eqn{\log\kappa}, and they are reported as `sd`
#' and `range` by [edf()] and named in `coef()`. The change is exact -- it is
#' a linear map of the log parameters, absorbed into the constants in front of
#' the three matrices -- and it is worth making for two reasons. The
#' parameters mean something on their own, and the SPDE likelihood has a long
#' ridge along constant marginal variance which in \eqn{(\tau, \kappa)} runs
#' diagonally and in \eqn{(\sigma, \rho)} runs along an axis, which the outer
#' optimiser finds much easier.
#'
#' The three matrices are carried through mgcv's `L` convention -- `L` holds
#' the powers of each parameter in each penalty, so that
#' \eqn{\lambda = \exp(L\theta)} -- which is what tells `gamRTMB` to give the
#' block two parameters and a [RTMB::dgmrf()] prior instead of one variance.
#' `mgcv::smooth2random()` cannot represent such a term at all, so an SPDE
#' smooth always takes the sparse route regardless of the `sparse` argument.
#'
#' Both parameters are estimated freely, with no prior on them. A single
#' realisation of a spatial field identifies them only loosely, and a field
#' whose range approaches the size of the domain is close to improper: the
#' range then runs off and the marginal standard deviation follows it, with
#' the fitted surface barely changing. If that happens, the fit is still
#' usable but the two parameters are not, and a shorter mesh `max.edge` will
#' not help -- it is the data that are uninformative. INLA's PC priors exist
#' for this.
#'
#' The field is proper for any positive \eqn{\kappa}, so no identifiability
#' constraint is imposed, as in INLA. It is not centred, and its overall level
#' is only weakly separated from the intercept when the range is large.
#'
#' @param object,data,knots As [mgcv::smooth.construct()].
#' @return A `smoothCon` object of class `spde.smooth`, with sparse `X` and
#'   sparse penalties.
#' @references
#' Lindgren, F., Rue, H. and Lindstrom, J. (2011). An explicit link between
#' Gaussian fields and Gaussian Markov random fields. \emph{JRSS-B} 73, 423-498.
#' @examples
#' \donttest{
#' if (requireNamespace("fmesher", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(x = runif(300), y = runif(300))
#'   d$z <- rnorm(300, sin(6 * d$x) * cos(6 * d$y), 0.3)
#'   mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.15, 0.4),
#'                               cutoff = 0.06, offset = c(0.1, 0.3))
#'   fit <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde",
#'                                      xt = list(mesh = mesh)), sd = ~ 1),
#'                  data = d)
#'   edf(fit)          # the sp column reports the field's sd and range
#' }
#' }
#' @exportS3Method mgcv::smooth.construct
smooth.construct.spde.smooth.spec <- function(object, data, knots) {
  if (!requireNamespace("fmesher", quietly = TRUE))
    stop("bs = \"spde\" needs the fmesher package", call. = FALSE)
  dim <- length(object$term)
  if (dim > 2L || dim < 1L)
    stop("an SPDE smooth is one- or two-dimensional, but this one has ",
         dim, " terms", call. = FALSE)

  mesh <- object$xt$mesh
  if (is.null(mesh)) {
    if (dim == 2L)
      stop("a two-dimensional SPDE smooth needs a mesh: build one with ",
           "fmesher::fm_mesh_2d() and pass it as ",
           "s(x, y, bs = \"spde\", xt = list(mesh = mesh)). The right mesh ",
           "depends on the domain and on the range you expect, so there is ",
           "no useful default.", call. = FALSE)
    k <- if (object$bs.dim < 0L) 20L else object$bs.dim
    x <- data[[object$term]]
    mesh <- fmesher::fm_mesh_1d(seq(min(x), max(x), length.out = k),
                                degree = 2, boundary = "free")
  }
  if (!inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d", "inla.mesh", "inla.mesh.1d")))
    stop("xt$mesh must be an fmesher mesh, from fm_mesh_2d() or fm_mesh_1d()",
         call. = FALSE)

  object$X <- fmesher::fm_basis(mesh, .spde_loc(object$term, data))
  fem <- fmesher::fm_fem(mesh)
  ## The mass matrix is the lumped one, as in INLA: `fm_fem` already forms
  ## g2 = g1 C0^-1 g1 with it, and using the consistent c1 in the kappa^4 term
  ## instead would not change the sparsity pattern -- both sit inside g2's --
  ## but would make the three matrices inconsistent with each other.
  ## lambda = exp(L theta) with theta = (log sigma, log rho). Substituting
  ## kappa = 2 sqrt(2) / rho and tau = 1 / (sigma kappa sqrt(4 pi)) into
  ## tau^2 (kappa^4 C + 2 kappa^2 G1 + G2) leaves the constants below in front
  ## of the three matrices and these powers of sigma and rho. Checked against
  ## the tau/kappa form to 5e-16 relative.
  object$S <- list((2 / pi) * fem$c0, fem$g1 / (2 * pi), fem$g2 / (32 * pi))
  object$L <- matrix(c(-2, -2, -2, -2, 0, 2), ncol = 2L)
  object$theta.names <- c("sd", "range")
  ## A range of a fifth of the domain and unit marginal standard deviation:
  ## smooth without being flat, and on the scale of a standardised predictor.
  object$theta.start <- c(0, log(.spde_extent(mesh, object$X) / 5))

  object$rank <- rep(ncol(object$X), length(object$S))
  object$null.space.dim <- 0L      # proper for any positive kappa
  ## Two flags that keep `smoothCon` out of the way. `no.rescale` suppresses
  ## its penalty rescaling, which would divide the finite element matrices by
  ## a norm of X and so break the meaning of tau and kappa -- and which is
  ## also the only place it takes a matrix norm of X, so the design matrix
  ## can stay the sparse object `fm_basis()` returned. A zero-row `C` says
  ## there is no centring constraint, which is right for a proper field and
  ## is what INLA does.
  object$no.rescale <- TRUE
  object$C <- matrix(0, 0L, ncol(object$X))
  object$df <- ncol(object$X)
  object$mesh <- mesh
  object$te.ok <- 0L
  class(object) <- "spde.smooth"
  object
}

#' @param object,data As [mgcv::Predict.matrix()].
#' @rdname smooth.construct.spde.smooth.spec
#' @exportS3Method mgcv::Predict.matrix
Predict.matrix.spde.smooth <- function(object, data)
  fmesher::fm_basis(object$mesh, .spde_loc(object$term, data))

#' Observation locations for an SPDE smooth
#' @keywords internal
.spde_loc <- function(term, data) {
  if (length(term) == 1L) return(as.numeric(data[[term]]))
  cbind(as.numeric(data[[term[1L]]]), as.numeric(data[[term[2L]]]))
}

#' A length scale for the meshed domain, for starting values
#'
#' The larger side of the bounding box of the mesh nodes that actually carry
#' basis weight, which is the region the data occupy rather than the outer
#' extension. A one-dimensional mesh stores its knots as a plain vector, and
#' its degree-2 basis has one more function than it has knots, so both the
#' shape and the length of `loc` have to be taken as they come.
#'
#' @keywords internal
.spde_extent <- function(mesh, X) {
  lc <- if (!is.null(mesh$loc)) mesh$loc else mesh$mid
  lc <- as.matrix(lc)
  lc <- lc[, seq_len(min(2L, ncol(lc))), drop = FALSE]
  used <- which(Matrix::colSums(abs(X)) > 0)
  if (length(used) > 1L && max(used) <= nrow(lc)) lc <- lc[used, , drop = FALSE]
  max(apply(lc, 2L, function(z) diff(range(z))), .Machine$double.eps)
}
