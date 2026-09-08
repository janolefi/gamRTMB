## Family objects. A family declares what can be modelled, with what link,
## and how to evaluate the log density; it knows nothing about smooths or
## fitting.

.links <- list(
  identity = list(linkfun = function(x) x,  linkinv = function(x) x),
  log      = list(linkfun = log,            linkinv = exp),
  logit    = list(linkfun = stats::qlogis,  linkinv = RTMB::plogis)
)

#' Default links by native parameter name
#'
#' RTMBdist ships 83 densities and exposes no metadata, so nothing reports a
#' parameter's support. Names, however, repeat across the library and mostly
#' imply one. Only unambiguous names belong here; the rest are resolved per
#' distribution in [.family_overrides].
#'
#' @keywords internal
.link_dict <- c(
  ## strictly positive
  sigma = "log", sd = "log", sdlog = "log", scale = "log", rate = "log",
  shape = "log", shape1 = "log", shape2 = "log", powershape = "log",
  phi = "log", kappa = "log", lambda = "log", df = "log", tau = "log",
  omega = "log", eta = "log", a = "log", b = "log", beta = "log",
  alpha = "log", mu1 = "log", mu2 = "log",
  ## `size` here is the nbinom-type overdispersion parameter. In binomial-type
  ## densities `size` is the known number of trials, and .family_overrides
  ## marks it fixed there, which removes it before this lookup happens.
  size = "log",
  ## unconstrained
  mu = "identity", mean = "identity", location = "identity",
  meanlog = "identity", skew = "identity", xi = "identity",
  ## (0, 1)
  prob = "logit", zeroprob = "logit", oneprob = "logit", rho = "logit"
)

#' Per-distribution corrections to the name-based defaults
#'
#' `links` overrides the dictionary, `fixed` marks arguments that are known
#' data or constants rather than parameters, `modelled` pins the subset.
#'
#' The ambiguous names, and why each needs resolving per distribution:
#' \describe{
#'   \item{alpha}{a positive shape in frechet/llogis/kumar, an unconstrained
#'     skewness in skewnorm/skewnorm2.}
#'   \item{nu}{an unconstrained Box-Cox power in bccg/bcpe/bct/gengamma, a
#'     positive power in powerexp, positive in combinom.}
#'   \item{theta}{a positive rate in bell, a direction vector in vmf2.}
#'   \item{size}{the overdispersion \emph{parameter} in nbinom2-type
#'     densities, the known number of trials in binomial-type ones. The most
#'     consequential distinction in the table.}
#'   \item{eps, min, max}{numerical guards and truncation bounds, never
#'     parameters.}
#' }
#'
#' @keywords internal
.family_overrides <- list(
  skewnorm    = list(links = c(alpha = "identity")),
  skewnorm2   = list(links = c(alpha = "identity")),
  bccg        = list(links = c(nu = "identity")),
  bcpe        = list(links = c(nu = "identity")),
  bct         = list(links = c(nu = "identity")),
  gengamma    = list(links = c(nu = "identity")),
  powerexp    = list(links = c(nu = "log")),
  powerexp2   = list(links = c(nu = "log")),
  jsu         = list(links = c(nu = "identity")),
  jsu2        = list(links = c(nu = "identity")),
  bell        = list(links = c(theta = "log")),
  bell2       = list(links = c(mu = "log")),
  pareto      = list(links = c(mu = "log")),
  beta        = list(fixed = "eps"),
  beta2       = list(fixed = "eps"),
  laplace     = list(fixed = "eps"),
  ## size = known number of binomial trials
  betabinom   = list(fixed = "size"),
  zibinom     = list(fixed = "size"),
  zibetabinom = list(fixed = "size"),
  hbetabinom  = list(fixed = "size"),
  hbinom      = list(fixed = "size"),
  ztbinom     = list(fixed = "size"),
  ztbetabinom = list(fixed = "size"),
  combinom    = list(fixed = "size", links = c(nu = "log")),
  dirmult     = list(fixed = "size"),
  ## truncation bounds
  truncnorm   = list(fixed = c("min", "max")),
  trunct      = list(fixed = c("min", "max")),
  trunct2     = list(fixed = c("min", "max"))
)

## Location parameters whose support the NAME does not reveal. `mu`/`mean` is
## unconstrained in a Gaussian but strictly positive in dgamma2, dnbinom2 and
## friends, and confined to (0, 1) in beta regression. An identity link there
## lets the linear predictor leave the support and the density returns NaN, so
## these are corrections rather than preferences -- found by fitting the
## families, not by reading their argument names.
.pos_location <- c("gamma2", "invgauss", "nbinom2", "hnbinom2", "zinbinom2",
                   "ztnbinom2", "zigamma2", "ziinvgauss", "bccg", "bcpe",
                   "bct", "gengamma")
.unit_location <- c("beta2", "zibeta2", "oibeta2", "zoibeta2")

local({
  for (d in .pos_location)
    .family_overrides[[d]]$links <<-
      c(.family_overrides[[d]]$links, c(mu = "log", mean = "log"))
  for (d in .unit_location)
    .family_overrides[[d]]$links <<-
      c(.family_overrides[[d]]$links, c(mu = "logit"))
})

## Not univariate regression families: vector- or matrix-valued responses, or
## copula constructions taking other densities as arguments. Rejected up front
## rather than failing somewhere inside the AD tape.
.not_families <- c("copula", "dcopula", "mvcopula", "mvt", "wishart",
                   "dirichlet", "dirmult", "vmf", "vmf2", "cmvgauss", "gmrf")

.link_neutral <- c(identity = 0, log = 0, logit = -2.2)

## Magnitude for shape parameters that must not start at zero; see the
## zero-score discussion in rtmbdist_family(). Verified to work from 0.1 up.
.shape_start_mag <- 0.25

.sample_skew <- function(y) {
  y <- as.numeric(y)
  mean(((y - mean(y)) / max(stats::sd(y), 1e-8))^3)
}

#' Build a family object from an RTMBdist density
#'
#' Reads the density's `formals()` and derives everything needed to model it:
#' the native parameter names in the density's own order, a link per
#' parameter, which arguments are data rather than parameters, and starting
#' values. 66 of RTMBdist's 83 densities resolve with no hand-written spec.
#'
#' @section What is read off the density:
#' Three things are available from `formals()` and all three are used. The
#' \strong{native names}, kept as-is (no mu/sigma/nu/tau renaming). The
#' author's own \strong{default values}, which are neutral in-support points
#' and so make good starting values. And \strong{redundant arguments}, visible
#' as a default that is an expression over other arguments
#' (`dinvgamma`'s `scale = 1/rate`, `dinvchisq`'s `scale = 1/df`): these are
#' dropped so the density derives them itself. Note that having a default is
#' \emph{not} a fixed-versus-modelled signal — `dgamma2` defaults
#' `mean = 1, sd = 1` and both are modelled.
#'
#' @section Modelled versus fixed arguments:
#' Modelled parameters get a formula, a link and smooths. Fixed arguments are
#' known data or constants — the number of binomial trials, truncation bounds,
#' a numerical `eps` — and are passed straight to the density, unreachable
#' from the formula interface. Supply one with `fixed = list(size = "trials")`
#' naming a column of the data, or a constant.
#'
#' @section Starting values:
#' Location parameters get data-driven starts and the rest keep the density's
#' default, with one exception. A shape parameter entering an already mean/sd
#' standardised density can have an \strong{identically zero score} at the
#' symmetric point: in `dskewnorm2` the direct effect of `alpha` on the log
#' density cancels exactly against the shift in the internal location needed
#' to hold `mean` and `sd` fixed, so the derivative is zero for every
#' observation at `alpha = 0` (measured at ~1e-15, i.e. exactly zero). The
#' default `alpha = 0` is therefore a perfectly neutral value and a useless
#' starting point: the optimiser has no descent direction, and because that
#' coefficient block has no curvature either, the inner Newton solve is
#' singular and the Laplace approximation is undefined. Such parameters start
#' off the symmetric point instead, with the sign taken from the sample
#' skewness — starting at `-0.5` on right-skewed data gets stuck just as
#' badly as starting at `0`.
#'
#' @param dist Density name, with or without the leading `d`
#'   (`"skewnorm2"` or `"dskewnorm2"`).
#' @param links Named character vector overriding the derived links.
#' @param fixed Named list of values for fixed arguments; each is a constant
#'   or the name of a column of the data.
#' @param start,eta_scale Optional replacements for the starting-value and
#'   linear-predictor-scale heuristics.
#' @return An object of class `gamRTMB_family`.
#' @examples
#' rtmbdist_family("gamma2")
#' rtmbdist_family("skewnorm2")
#' rtmbdist_family("betabinom", fixed = list(size = "trials"))
#' @export
rtmbdist_family <- function(dist, links = NULL, fixed = NULL, start = NULL,
                            eta_scale = NULL) {
  ns <- asNamespace("RTMBdist")
  nm <- as.character(dist)[1L]
  dfun_name <- if (startsWith(nm, "d") && exists(nm, ns, inherits = FALSE)) nm
               else paste0("d", nm)
  dist <- sub("^d", "", dfun_name)
  if (dist %in% .not_families)
    stop("'", dist, "' is not a univariate regression family (vector- or ",
         "matrix-valued response, or a copula construction); gamRTMB models ",
         "one scalar response per observation")
  if (!exists(dfun_name, ns, inherits = FALSE))
    stop("RTMBdist has no density '", dfun_name, "'")
  dfun <- get(dfun_name, envir = ns)

  fo <- formals(dfun)
  fo <- fo[setdiff(names(fo), c("x", "log"))]
  ov <- .family_overrides[[dist]]

  ## classify the arguments -------------------------------------------------
  has_default <- vapply(fo, function(z) !identical(z, quote(expr = )), TRUE)
  ## Evaluate each default in an empty environment. A literal or self-contained
  ## expression (0, -Inf, NULL) evaluates; one referring to other arguments
  ## (scale = 1/rate) does not, which is exactly the signal that the argument
  ## is a redundant reparameterisation the density derives for itself.
  dval <- lapply(names(fo), function(a) if (!has_default[[a]]) NULL else
    tryCatch(list(v = eval(fo[[a]], baseenv())), error = function(e) NULL))
  names(dval) <- names(fo)
  is_num <- vapply(names(fo), function(a)
    !is.null(dval[[a]]) && is.numeric(dval[[a]]$v) && length(dval[[a]]$v) == 1L, TRUE)
  derived <- names(fo)[has_default &
                         vapply(names(fo), function(a) is.null(dval[[a]]), TRUE)]
  fixed_names <- unique(c(ov$fixed, names(fixed)))
  modelled <- setdiff(names(fo), c(derived, fixed_names))
  if (!is.null(ov$modelled)) modelled <- ov$modelled
  if (!length(modelled)) stop("no modelled parameters left for ", dfun_name)

  ## links ------------------------------------------------------------------
  lk <- .link_dict[modelled]; names(lk) <- modelled
  ovl <- ov$links[intersect(names(ov$links), modelled)]
  if (length(ovl)) lk[names(ovl)] <- ovl
  if (length(links)) {
    unk <- setdiff(names(links), modelled)
    if (length(unk))
      stop("links given for non-modelled parameter(s) ", paste(unk, collapse = ", "),
           "; ", dfun_name, " models ", paste(modelled, collapse = ", "))
    lk[names(links)] <- links
  }
  if (anyNA(lk))
    stop("no default link known for parameter(s) ",
         paste(modelled[is.na(lk)], collapse = ", "), " of ", dfun_name,
         ". Pass links = c(", modelled[is.na(lk)][1L], " = \"log\"), or add ",
         "an entry to .family_overrides.")
  bad <- setdiff(lk, names(.links))
  if (length(bad)) stop("unsupported link(s): ", paste(bad, collapse = ", "))

  def_start <- vapply(modelled, function(nm) {
    v <- if (is_num[[nm]]) as.numeric(dval[[nm]]$v) else NA_real_
    z <- if (is.na(v)) NA_real_ else .links[[lk[[nm]]]]$linkfun(v)
    if (is.na(z) || !is.finite(z)) .link_neutral[[lk[[nm]]]] else z
  }, numeric(1))

  fixed_vals <- list()
  for (nm in fixed_names) {
    if (nm %in% names(fixed)) fixed_vals[[nm]] <- fixed[[nm]]
    else if (is_num[[nm]] && is.finite(dval[[nm]]$v))
      fixed_vals[[nm]] <- as.numeric(dval[[nm]]$v)
    else if (has_default[[nm]]) next          # density supplies its own default
    else stop("family '", dist, "' needs the fixed argument '", nm,
              "' supplied from the data, e.g. rtmbdist_family(\"", dist,
              "\", fixed = list(", nm, " = \"", nm, "_column\"))")
  }

  logdens <- function(y, theta, fx = list()) {
    args <- c(list(y), theta[modelled], fx, list(log = TRUE))
    names(args)[1L] <- "x"
    do.call(dfun, args)
  }

  zero_score <- modelled[lk[modelled] == "identity" &
                           !modelled %in% c("mu", "mean", "location", "meanlog") &
                           def_start[modelled] == 0]

  start_fun <- if (!is.null(start)) start else function(y) {
    s <- def_start
    if (length(zero_score)) {
      sgn <- if (.sample_skew(y) < 0) -1 else 1
      for (nm in zero_score) s[[nm]] <- sgn * .shape_start_mag
    }
    for (nm in intersect(modelled, c("mu", "mean", "location")))
      s[[nm]] <- switch(lk[[nm]], identity = mean(y),
                        log = log(max(mean(y), 1e-3)),
                        logit = stats::qlogis(min(max(mean(y), .01), .99)), s[[nm]])
    if ("meanlog" %in% modelled) s[["meanlog"]] <- mean(log(pmax(y, 1e-8)))
    for (nm in intersect(modelled, c("sigma", "sd", "scale")))
      if (lk[[nm]] == "log") s[[nm]] <- log(max(stats::sd(y), 1e-3))
    if ("zeroprob" %in% modelled)
      s[["zeroprob"]] <- stats::qlogis(min(max(mean(y == 0), 0.02), 0.5))
    s
  }

  scale_fun <- if (!is.null(eta_scale)) eta_scale else function(y) {
    stats::setNames(vapply(modelled, function(nm)
      if (lk[[nm]] == "identity" && nm %in% c("mu", "mean", "location"))
        max(stats::sd(y), 1e-3) else 0.5, numeric(1)), modelled)
  }

  structure(list(family = dist, dist = dfun_name, parnames = modelled,
                 links = stats::setNames(as.character(lk), modelled),
                 fixed = fixed_vals, derived = derived,
                 logdens = logdens, start = start_fun, eta_scale = scale_fun),
            class = "gamRTMB_family")
}

#' Gaussian location-scale family
#'
#' Written out rather than derived, because RTMB's own `dnorm` is not part of
#' RTMBdist and because this is the family the pipeline is checked against
#' [mgcv::gaulss()] with.
#'
#' @param link_mu,link_sigma Link names.
#' @return An object of class `gamRTMB_family`.
#' @examples
#' gaussian_ls()
#' @export
gaussian_ls <- function(link_mu = "identity", link_sigma = "log") {
  structure(list(
    family = "gaussian_ls", dist = "dnorm",
    parnames = c("mu", "sigma"),
    links = c(mu = link_mu, sigma = link_sigma),
    fixed = list(), derived = character(0),
    logdens = function(y, p, fx = list()) dnorm(y, p$mu, p$sigma, log = TRUE),
    start = function(y) c(mu = mean(y), sigma = log(stats::sd(y))),
    eta_scale = function(y) c(mu = stats::sd(y), sigma = 0.5)
  ), class = "gamRTMB_family")
}

#' @param x A `gamRTMB_family`.
#' @param ... Ignored.
#' @rdname rtmbdist_family
#' @export
print.gamRTMB_family <- function(x, ...) {
  cat("gamRTMB family: ", x$family, "  [", x$dist, "]\n", sep = "")
  cat("  modelled: ",
      paste(sprintf("%s (%s)", x$parnames, x$links), collapse = ", "), "\n", sep = "")
  if (length(x$fixed)) cat("  fixed:    ", paste(names(x$fixed), collapse = ", "), "\n", sep = "")
  if (length(x$derived))
    cat("  derived:  ", paste(x$derived, collapse = ", "),
        " (computed by the density)\n", sep = "")
  invisible(x)
}

#' Resolve a family's fixed arguments against the data
#'
#' A numeric constant is recycled; a character string or one-sided formula
#' names a column of `data`.
#'
#' @keywords internal
.resolve_fixed <- function(family, data, n) {
  fx <- family$fixed
  if (!length(fx)) return(list())
  stats::setNames(lapply(names(fx), function(nm) {
    v <- fx[[nm]]
    if (inherits(v, "formula")) v <- eval(v[[2L]], data, environment(v))
    else if (is.character(v) && length(v) == 1L && v %in% names(data)) v <- data[[v]]
    if (!is.numeric(v))
      stop("fixed argument '", nm, "' must be numeric, or name a numeric ",
           "column of data")
    if (length(v) == 1L || length(v) == n) v
    else stop("fixed argument '", nm, "' has length ", length(v),
              ", need 1 or ", n)
  }), names(fx))
}
