## What can be fitted. A family declares which parameters are modelled, with
## what link, and how to evaluate the log density. It knows nothing about
## smooths or fitting.

.links <- list(
  identity = list(linkfun = function(x) x,  linkinv = function(x) x),
  log      = list(linkfun = log,            linkinv = exp),
  logit    = list(linkfun = stats::qlogis,  linkinv = RTMB::plogis)
)

#' Default links, by native parameter name
#'
#' RTMBdist exposes no metadata, so nothing reports a parameter's support.
#' Names, however, repeat across the library and mostly imply one. Only
#' unambiguous names belong here; the rest are resolved per distribution in
#' [.family_overrides].
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
  ## densities it is the known number of trials, and .family_overrides marks
  ## it fixed there, which removes it before this lookup happens.
  size = "log",
  ## unconstrained
  mu = "identity", mean = "identity", location = "identity",
  meanlog = "identity", skew = "identity", xi = "identity",
  ## (0, 1)
  prob = "logit", zeroprob = "logit", oneprob = "logit", rho = "logit"
)

## Argument names that are never regression parameters wherever they appear:
## numerical guards, non-centrality, truncation bounds. Removing them by name
## saves an override entry for every density that has one.
.never_modelled <- c("eps", "ncp", "min", "max")

#' Per-distribution corrections to the name-based defaults
#'
#' `links` overrides the dictionary, `fixed` marks arguments that are known
#' data rather than parameters, `modelled` pins the subset.
#'
#' The ambiguous names, and why each must be resolved per distribution:
#' \describe{
#'   \item{alpha}{a positive shape in frechet/llogis/kumar, an unconstrained
#'     skewness in skewnorm/skewnorm2/sn.}
#'   \item{nu}{an unconstrained Box-Cox power in bccg/bcpe/bct/gengamma, a
#'     positive power in powerexp, positive in combinom.}
#'   \item{theta}{a positive rate in bell, a direction vector in vmf2.}
#'   \item{size}{the overdispersion \emph{parameter} in nbinom2-type
#'     densities, the known number of trials in binomial-type ones. The most
#'     consequential distinction in the table.}
#' }
#'
#' @keywords internal
.family_overrides <- list(
  skewnorm    = list(links = c(alpha = "identity")),
  skewnorm2   = list(links = c(alpha = "identity")),
  sn          = list(links = c(alpha = "identity")),
  bccg        = list(links = c(nu = "identity")),
  bcpe        = list(links = c(nu = "identity")),
  bct         = list(links = c(nu = "identity")),
  gengamma    = list(links = c(nu = "identity")),
  powerexp    = list(links = c(nu = "log")),
  powerexp2   = list(links = c(nu = "log")),
  jsu         = list(links = c(nu = "identity")),
  jsu2        = list(links = c(nu = "identity")),
  SHASHo      = list(links = c(nu = "identity")),
  bell        = list(links = c(theta = "log")),
  bell2       = list(links = c(mu = "log")),
  pareto      = list(links = c(mu = "log")),
  ## size = known number of trials
  binom       = list(fixed = "size"),
  betabinom   = list(fixed = "size"),
  zibinom     = list(fixed = "size"),
  zibetabinom = list(fixed = "size"),
  hbetabinom  = list(fixed = "size"),
  hbinom      = list(fixed = "size"),
  ztbinom     = list(fixed = "size"),
  ztbetabinom = list(fixed = "size"),
  combinom    = list(fixed = "size", links = c(nu = "log"))
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

## Densities that take a vector- or matrix-valued response: rejected up front
## rather than failing inside the AD tape. Only these six need naming -- the
## copula constructions take other densities rather than an `x` argument, so
## the shape check in .find_density() already excludes them.
.not_families <- c("mvt", "wishart", "dirichlet", "dirmult", "vmf", "vmf2")

## RTMBdist is searched first; these standard densities, which RTMB makes
## AD-aware, fill in the families it does not carry. The rest of RTMB is
## multivariate, a robust reparameterisation, or has no sensible default link.
.rtmb_families <- c("norm", "pois", "binom", "gamma", "exp", "lnorm",
                    "weibull", "cauchy", "logis", "t", "chisq")

.link_neutral <- c(identity = 0, log = 0, logit = -2.2)

## Magnitude for shape parameters that must not start at zero; see the
## zero-score discussion in fam(). Verified to work from 0.1 upward.
.shape_start_mag <- 0.25

.sample_skew <- function(y) {
  y <- as.numeric(y)
  mean(((y - mean(y)) / max(stats::sd(y), 1e-8))^3)
}

#' Locate a density by name
#'
#' RTMBdist first, then the curated RTMB list. Accepts either `"gamma2"` or
#' `"dgamma2"`.
#'
#' @keywords internal
.find_density <- function(dist) {
  nm <- as.character(dist)[1L]
  ## Try the name as given before prefixing, so that a density whose own name
  ## begins with d ("ddirichlet", "dcopula") is not mistaken for an already
  ## prefixed one and stripped to nonsense.
  for (src in c("RTMBdist", "RTMB")) {
    ns <- asNamespace(src)
    for (cand in unique(c(nm, paste0("d", nm)))) {
      if (!startsWith(cand, "d")) next
      if (src == "RTMB" && !sub("^d", "", cand) %in% .rtmb_families) next
      if (!exists(cand, ns, inherits = FALSE)) next
      f <- get(cand, envir = ns)
      if (is.function(f) && all(c("x", "log") %in% names(formals(f))))
        return(list(dfun = f, name = cand, dist = sub("^d", "", cand),
                    source = src))
    }
  }
  NULL
}

#' Derive a family's structure from its density
#'
#' Everything that can be read off `formals()` plus the correction tables, and
#' nothing that needs the data. Shared by [fam()] and [families()], so that
#' listing families does not mean catching errors from constructing them.
#'
#' @return `NULL` if no such density; otherwise a list with the modelled
#'   parameters, their links, dropped (derived) arguments, the fixed arguments
#'   and which of them the user must still supply, and default starts.
#' @keywords internal
.classify <- function(dist, fixed = NULL) {
  d <- .find_density(dist)
  if (is.null(d)) return(NULL)
  if (d$dist %in% .not_families)
    return(structure(list(dist = d$dist), class = "gamRTMB_notfamily"))

  fo <- formals(d$dfun)
  fo <- fo[setdiff(names(fo), c("x", "log"))]
  ov <- .family_overrides[[d$dist]]

  has_default <- vapply(fo, function(z) !identical(z, quote(expr = )), TRUE)
  ## Evaluate each default in an empty environment. A literal or a
  ## self-contained expression (0, -Inf, NULL) evaluates; one referring to
  ## other arguments (scale = 1/rate) does not, which is exactly the signal
  ## that the argument is a redundant reparameterisation the density derives
  ## for itself. Note that HAVING a default is not a fixed-versus-modelled
  ## signal: dgamma2 defaults mean = 1, sd = 1, and both are modelled.
  dval <- lapply(names(fo), function(a) if (!has_default[[a]]) NULL else
    tryCatch(list(v = eval(fo[[a]], baseenv())), error = function(e) NULL))
  names(dval) <- names(fo)
  is_num <- vapply(names(fo), function(a) !is.null(dval[[a]]) &&
                     is.numeric(dval[[a]]$v) && length(dval[[a]]$v) == 1L, TRUE)
  derived <- names(fo)[has_default &
                         vapply(names(fo), function(a) is.null(dval[[a]]), TRUE)]

  ## Two kinds of non-modelled argument. A guard named in .never_modelled is
  ## never required: if it has no default (dt's `ncp`), it is simply left out
  ## and the density does whatever it does without it. One named in the
  ## override table is real data (`size` = number of trials) and, lacking a
  ## default, must be supplied.
  guards <- intersect(names(fo), .never_modelled)
  fixed_names <- unique(c(ov$fixed, guards, names(fixed)))
  modelled <- if (!is.null(ov$modelled)) ov$modelled else
    setdiff(names(fo), c(derived, fixed_names))

  lk <- .link_dict[modelled]; names(lk) <- modelled
  ## support corrections the parameter name cannot reveal
  if (d$dist %in% .pos_location)
    lk[names(lk) %in% c("mu", "mean")] <- "log"
  if (d$dist %in% .unit_location)
    lk[names(lk) == "mu"] <- "logit"
  ovl <- ov$links[intersect(names(ov$links), modelled)]
  if (length(ovl)) lk[names(ovl)] <- ovl

  ## fixed arguments: a usable default needs nothing from the user; one with a
  ## non-numeric default (eps = NULL) is left to the density; otherwise it is
  ## genuinely data and must be supplied
  vals <- list(); needs <- character(0)
  for (nm in fixed_names) {
    if (nm %in% names(fixed)) vals[[nm]] <- fixed[[nm]]   # the user's value wins
    else if (nm %in% guards) next                         # leave it to the density
    else if (is_num[[nm]] && is.finite(dval[[nm]]$v))
      vals[[nm]] <- as.numeric(dval[[nm]]$v)
    else if (has_default[[nm]]) next
    else needs <- c(needs, nm)
  }

  def_start <- vapply(modelled, function(nm) {
    v <- if (is_num[[nm]]) as.numeric(dval[[nm]]$v) else NA_real_
    z <- if (is.na(v) || is.na(lk[[nm]])) NA_real_ else .links[[lk[[nm]]]]$linkfun(v)
    if (is.na(z) || !is.finite(z)) .link_neutral[[if (is.na(lk[[nm]])) "identity"
                                                  else lk[[nm]]]] else z
  }, numeric(1))

  list(dfun = d$dfun, name = d$name, dist = d$dist, source = d$source,
       modelled = modelled, links = lk, derived = derived,
       fixed = vals, needs = needs, def_start = def_start)
}

#' Build a family object
#'
#' Reads a density's `formals()` and derives everything needed to model it:
#' the native parameter names in the density's own order, a link per
#' parameter, which arguments are data rather than parameters, and starting
#' values. Most of RTMBdist works with no hand-written family; use
#' [families()] to see what is available.
#'
#' @section Where densities come from:
#' \pkg{RTMBdist} is searched first, then a curated list of the standard
#' densities that \pkg{RTMB} makes AD-aware (`norm`, `pois`, `binom`, `gamma`,
#' `exp`, `lnorm`, `weibull`, `cauchy`, `logis`, `t`, `chisq`). Parameter
#' names are always the density's own: a skew normal is `xi`, `omega`,
#' `alpha`, and a Gaussian is `mean`, `sd`.
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
#' own default, with one exception. A shape parameter entering an already
#' mean/sd standardised density can have an \strong{identically zero score} at
#' the symmetric point: in `dskewnorm2` the direct effect of `alpha` on the log
#' density cancels exactly against the shift in the internal location needed
#' to hold `mean` and `sd` fixed, so the derivative is zero for every
#' observation at `alpha = 0` (measured at ~1e-15, i.e. exactly zero). The
#' default `alpha = 0` is therefore a perfectly neutral value and a useless
#' starting point: the optimiser has no descent direction, and because that
#' coefficient block has no curvature either, the inner Newton solve is
#' singular and the Laplace approximation is undefined. Such parameters start
#' off the symmetric point instead, with the sign taken from the sample
#' skewness — starting at `-0.5` on right-skewed data gets stuck just as badly
#' as starting at `0`.
#'
#' @param dist Density name, with or without the leading `d`
#'   (`"skewnorm2"` or `"dskewnorm2"`).
#' @param links Named character vector overriding the derived links.
#' @param fixed Named list of values for fixed arguments; each is a constant
#'   or the name of a column of the data.
#' @param start,eta_scale Optional replacements for the starting-value and
#'   linear-predictor-scale heuristics.
#' @return An object of class `gamRTMB_family`.
#' @seealso [families()] for what is available, [gamRTMB()] to fit.
#' @examples
#' fam("norm")
#' fam("gamma2")
#' fam("skewnorm2")
#' fam("betabinom", fixed = list(size = "trials"))
#' @export
fam <- function(dist, links = NULL, fixed = NULL, start = NULL,
                eta_scale = NULL) {
  sp <- .classify(dist, fixed)
  if (is.null(sp))
    stop("no density for '", as.character(dist)[1L],
         "' in RTMBdist, and it is not one of the standard RTMB densities ",
         "offered (", paste(.rtmb_families, collapse = ", "),
         "). See families().")
  if (inherits(sp, "gamRTMB_notfamily"))
    stop("'", sp$dist, "' is not a univariate regression family (vector- or ",
         "matrix-valued response, or a copula construction); gamRTMB models ",
         "one scalar response per observation")
  if (length(sp$needs))
    stop("family '", sp$dist, "' needs the fixed argument(s) ",
         paste0("'", sp$needs, "'", collapse = ", "), " supplied from the ",
         "data, e.g. fam(\"", sp$dist, "\", fixed = list(", sp$needs[1L],
         " = \"", sp$needs[1L], "_column\"))")

  modelled <- sp$modelled
  lk <- sp$links
  if (length(links)) {
    unk <- setdiff(names(links), modelled)
    if (length(unk))
      stop("links given for non-modelled parameter(s) ",
           paste(unk, collapse = ", "), "; ", sp$name, " models ",
           paste(modelled, collapse = ", "))
    lk[names(links)] <- links
  }
  if (anyNA(lk))
    stop("no default link known for parameter(s) ",
         paste(modelled[is.na(lk)], collapse = ", "), " of ", sp$name,
         ". Pass links = c(", modelled[is.na(lk)][1L], " = \"log\"), or add ",
         "an entry to .family_overrides.")

  dfun <- sp$dfun
  logdens <- function(y, theta, fx = list()) {
    args <- c(list(y), theta[modelled], fx, list(log = TRUE))
    names(args)[1L] <- "x"
    do.call(dfun, args)
  }

  def_start <- sp$def_start
  zero_score <- modelled[lk == "identity" &
                           !modelled %in% c("mu", "mean", "location", "meanlog") &
                           def_start == 0]

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

  structure(list(family = sp$dist, dist = sp$name, source = sp$source,
                 parnames = modelled,
                 links = stats::setNames(as.character(lk), modelled),
                 fixed = sp$fixed, derived = sp$derived,
                 logdens = logdens, start = start_fun, eta_scale = scale_fun),
            class = "gamRTMB_family")
}

#' Available families
#'
#' Every density [fam()] can turn into a family: all of \pkg{RTMBdist} that is
#' a univariate regression family, plus the standard \pkg{RTMB} densities.
#'
#' @param pattern Optional regular expression to filter family names.
#' @return A data frame with one row per family: its name, the modelled
#'   parameters with their links, any fixed arguments that must be supplied
#'   from the data, and which package the density comes from.
#' @seealso [fam()]
#' @examples
#' head(families(), 10)
#' families("beta")
#' @export
families <- function(pattern = NULL) {
  cand <- unique(c(
    grep("^d", ls(asNamespace("RTMBdist")), value = TRUE),
    paste0("d", .rtmb_families)))
  rows <- lapply(cand, function(nm) {
    sp <- tryCatch(.classify(nm), error = function(e) NULL)
    if (is.null(sp) || inherits(sp, "gamRTMB_notfamily") || anyNA(sp$links))
      return(NULL)
    data.frame(family = sp$dist,
               parameters = paste(sprintf("%s/%s", sp$modelled, sp$links),
                                  collapse = ", "),
               needs = paste(sp$needs, collapse = ", "),
               source = sp$source, row.names = NULL)
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$family), ]
  if (!is.null(pattern)) out <- out[grepl(pattern, out$family), ]
  rownames(out) <- NULL
  out
}

#' @param x A `gamRTMB_family`.
#' @param ... Ignored.
#' @rdname fam
#' @export
print.gamRTMB_family <- function(x, ...) {
  cat("gamRTMB family: ", x$family, "  [", x$dist, ", ", x$source, "]\n", sep = "")
  cat("  modelled: ",
      paste(sprintf("%s (%s)", x$parnames, x$links), collapse = ", "), "\n", sep = "")
  if (length(x$fixed))
    cat("  fixed:    ", paste(names(x$fixed), collapse = ", "), "\n", sep = "")
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
