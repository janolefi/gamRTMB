## Getting a fit started at all: the cheap probe, the sigma_frac ladder it
## drives, and the curvature diagnosis for the models no ladder can rescue.
## The Box-Cox power exponential is the worked example throughout -- see
## .inner_indefinite() for why it is the family that breaks this.

sim_ls <- function(n = 200, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x1), exp(-1 + d$x2))
  d
}

## the minimal model that breaks: one smooth, on mu, everything else constant
bcpe_obj <- function(frac = 0.05, method = "REML") {
  skip_if_not_installed("gamlss.data")
  film90 <- gamlss.data::film90
  fo <- lborev1 ~ list(mu = ~ s(lboopen, bs = "ps"), sigma = ~ 1,
                       nu = ~ 1, tau = ~ 1)
  family <- fam("bcpe")
  pf <- .parse_formula(fo, family$parnames)
  md <- .model_data(pf$response, pf$par_formulas, film90, NULL,
                    stats::na.omit, environment())
  des <- .build_design(pf$par_formulas, md$data, family$parnames,
                       knots = NULL, sparse = "auto")
  fx <- .resolve_fixed(family, md$data, length(md$y))
  obj <- RTMB::MakeADFun(.make_nll(des, family, md$y, fx, NULL),
                         .init_pars(des, family, md$y, frac, NULL),
                         random = if (method == "REML") c("beta", "b") else "b",
                         silent = TRUE, inner.control = list(maxit = 1000))
  list(obj = obj, design = des, method = method)
}

test_that("the probe answers, and leaves the inner cap as it found it", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  obj <- fit$obj
  expect_true(.probe_finite(obj))
  expect_equal(obj$env$inner.control$maxit, 1000)
})

test_that("the probe rejects an inner problem the full cap also rejects", {
  z <- bcpe_obj()
  expect_false(.probe_finite(z$obj))
  ## and the full cap agrees, which is the claim the probe is standing in for
  expect_false(is.finite(z$obj$fn(z$obj$par)))
})

test_that("inner_control reaches MakeADFun's inner solver", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ s(x1, k = 8), data = d, inner_control = list(maxit = 7))
  expect_equal(fit$obj$env$inner.control$maxit, 7)
  ## and the default is left alone
  expect_equal(gamRTMB(y ~ s(x1, k = 8), data = d)$obj$env$inner.control$maxit,
               1000)
})

test_that("random coefficients are labelled in MakeADFun's own order", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  des <- fit$design
  lab <- .random_labels(des, "REML")
  expect_length(lab, des$nbeta + des$nb)
  expect_identical(lab[seq_len(des$nbeta)], .beta_labels(des))
  ## under ML only the penalized coefficients are random
  expect_length(.random_labels(des, "ML"), des$nb)
})

test_that("a positive definite inner Hessian is not diagnosed as anything", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  expect_null(.inner_indefinite(fit$obj, fit$design, "REML", "norm"))
})

test_that("an indefinite inner Hessian is named as such, with the direction", {
  z <- bcpe_obj()
  invisible(suppressWarnings(tryCatch(z$obj$fn(z$obj$par), error = function(e) NULL)))
  msg <- .inner_indefinite(z$obj, z$design, z$method, "bcpe")
  expect_type(msg, "character")
  expect_match(msg, "not positive definite")
  ## the flattest negative direction implicates the smooth's unpenalized
  ## null-space column, which is what makes bs = "cs" the fix
  expect_match(msg, "s\\(lboopen\\)\\.null1")
  expect_match(msg, "bs = \"cs\"")
})

test_that("the ladder reports the start it settled on, and is off by default", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  ## a model that starts cleanly never walks the ladder
  expect_null(fit$sigma_frac_used)
  expect_true(all(.sigma_frac_ladder > 0))
})

test_that("an explicit start is never overwritten by the ladder", {
  skip_if_not_installed("gamlss.data")
  film90 <- gamlss.data::film90
  ## the ladder cannot rescue this model at any frac, so the fit fails either
  ## way; what is being checked is that it fails *without* retrying, which a
  ## user-supplied start is entitled to
  expect_error(
    expect_no_message(
      gamRTMB(lborev1 ~ list(mu = ~ s(lboopen, bs = "ps"), sigma = ~ 1,
                             nu = ~ 1, tau = ~ 1),
              family = fam("bcpe"), data = film90,
              start = list(beta = c(2.6, 0, -1.5, 1, 0.69)))),
    "not positive definite")
})

## ---------------------------------------------------------------------------
## EDF at a point the optimiser never reached

test_that("healthy EDF sit inside the [0, k] the guard checks", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  ph <- .penalized_hessian(fit)
  e <- 1 - Matrix::diag(Matrix::solve(ph$H,
                                      .penalty_matrix(fit$design, fit$log_sigma, ph)))
  for (bl in fit$design$blocks) {
    expect_gte(sum(e[ph$i_b[bl$idx]]), -.edf_tol * bl$q)
    expect_lte(sum(e[ph$i_b[bl$idx]]), bl$q * (1 + .edf_tol))
  }
  ## an `iid` block has a diagonal penalty, so its coefficients are confined
  ## one by one as well -- which is exactly what a non-diagonal one is not
  expect_gte(min(e), -.edf_tol)
  expect_lte(max(e), 1 + .edf_tol)
  ## and the guard therefore stays out of the way
  expect_silent(edf(fit))
})

test_that("a block EDF is confined where its coefficients are not", {
  set.seed(4); n <- 400
  d <- data.frame(x = stats::runif(n), z = stats::runif(n))
  d$y <- stats::rnorm(n, sin(5 * d$x) * cos(4 * d$z), 0.3)
  fit <- gamRTMB(y ~ list(mean = ~ te(x, z), sd = ~ 1), data = d)
  ph <- .penalized_hessian(fit)
  e <- 1 - Matrix::diag(Matrix::solve(ph$H,
                                      .penalty_matrix(fit$design, fit$log_sigma, ph)))
  bl <- fit$design$blocks[[1L]]
  ## the block total is inside [0, q] and is what gets reported ...
  expect_gte(sum(e[ph$i_b[bl$idx]]), 0)
  expect_lte(sum(e[ph$i_b[bl$idx]]), bl$q)
  ## ... while single coefficients are not confined to [0, 1], because the
  ## penalty is not diagonal and H^-1 S is not symmetric
  expect_gt(max(e[ph$i_b[bl$idx]]), 1 + .edf_tol)
  expect_silent(edf(fit))
})

test_that("EDF outside [0, 1] are reported as NA rather than as numbers", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  ## Shrink the variance components without moving the mode, which is what a
  ## fit stopped at a bad point amounts to as far as `1 - diag(H^-1 S)` is
  ## concerned: S grows, the "projection" leaves [0, 1], and the EDF stop
  ## meaning anything. Cheaper and more exact than fitting a family that does
  ## this for real -- see .inner_indefinite() for one that does.
  bad <- fit
  bad$log_sigma <- fit$log_sigma - 10
  expect_warning(e <- edf(bad), "not defined at these values")
  expect_true(all(is.na(e$edf)))
  expect_true(is.na(attr(e, "edf.total")))
  ## the rest of the table still describes the model, so it is kept
  expect_equal(e$term, edf(fit)$term)
  expect_equal(e$k, edf(fit)$k)
})

test_that("summary says the EDF are undefined rather than printing a total", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  bad <- fit
  bad$log_sigma <- fit$log_sigma - 10
  s <- suppressWarnings(summary(bad))
  expect_false(s$edf_defined)
  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "undefined at these values")
  expect_no_match(out, "Total EDF")
  ## AIC counts the EDF as its parameter count, so it goes too
  expect_no_match(out, "AIC")
  ## a healthy fit is untouched by all of this
  ok <- paste(utils::capture.output(print(summary(fit))), collapse = "\n")
  expect_match(ok, "Total EDF")
  expect_match(ok, "AIC")
})

## ---------------------------------------------------------------------------
## Finishing: the convergence test, and the ladder walked a second time when
## the outer optimisation stalls rather than when the start is unusable.

test_that("a flat minimum nlminb calls false convergence is still converged", {
  ## nlminb answers a flat minimum with a nonzero code however small the
  ## gradient, which used to be reported as a failed fit.
  expect_true(.outer_ok(list(convergence = 1L), 1.86e-04))
  expect_false(.outer_ok(list(convergence = 1L), 5.95))
  ## the threshold is the EFS engine's, so the two mean the same thing
  expect_true(.outer_ok(list(convergence = 1L), .efs_defaults$gtol))
  expect_false(.outer_ok(list(convergence = 1L), .efs_defaults$gtol * 1.01))
  ## a zero code is never overruled by a large gradient, and NA is not small
  expect_true(.outer_ok(list(convergence = 0L), 1e6))
  expect_false(.outer_ok(list(convergence = 1L), NA_real_))
})

test_that(".outer_grad reports rather than throws", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8)), data = d)
  expect_equal(.outer_grad(fit$obj, fit$opt$par), fit$max_grad)
  ## no smoothing parameters to be stationary in
  expect_equal(.outer_grad(fit$obj, numeric(0)), 0)
  ## a gradient that cannot be evaluated is a diagnostic, not a lost fit
  expect_true(is.na(.outer_grad(list(gr = function(p) stop("nope")), 1)))
})

test_that("a stalled fit is retried over the ladder, and improves", {
  ## film90 with a JSU: the default sigma_frac is the only value in the set
  ## that fails, and the first rung reaches a converged optimum. The start is
  ## perfectly finite here, so the probe never fires -- this is the second
  ## walk, driven by the outer optimisation rather than by the start.
  skip_on_cran()
  skip_if_not_installed("gamlss.data")
  film90 <- gamlss.data::film90
  fo <- lborev1 ~ list(mu = ~ s(lboopen, k = 20), sigma = ~ s(lboopen, k = 20),
                       nu = ~ 1, tau = ~ 1)
  said <- character(0)
  fit <- withCallingHandlers(
    gamRTMB(fo, family = fam("jsu2"), data = film90),
    message = function(m) {
      said <<- c(said, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  ## the first message only fires when the fit as it stood had not converged,
  ## and the second only when a rung beat it -- which together are the claim
  expect_match(said[1], "the outer optimisation did not converge; retrying")
  expect_match(said[2], "reached a lower criterion and converged")
  expect_true(fit$convergence)
  expect_equal(fit$sigma_frac_used, 0.005)
  expect_lt(fit$max_grad, .efs_defaults$gtol)
})

test_that("the retry selects on the criterion, not on which rung converged", {
  ## MASS::mcycle with a t response is the case that decides the rule: the two
  ## rungs that converge sit at 685.4 and 685.9 where the unconverged
  ## incumbent sits at 579.9, so a retry that took a converged candidate on
  ## sight would trade a good fit for a badly over-smoothed one.
  skip_on_cran()
  skip_if_not_installed("MASS")
  mcycle <- MASS::mcycle
  fo <- accel ~ list(mu = ~ s(times, k = 20), sigma = ~ s(times, k = 10),
                     df = ~ 1)
  fit <- suppressMessages(gamRTMB(fo, family = fam("t2"), data = mcycle))
  ## it did not take either converged-but-worse rung
  expect_false(fit$sigma_frac_used %in% c(0.005, 0.001))
  expect_lt(fit$objective, 600)
  ## an honest flag: the best point found is still not a stationary one
  expect_false(fit$convergence)
})

test_that("a start the user supplied switches the retry off too", {
  ## Same guard as the probe's ladder: `repars` is NULL, so there is nothing
  ## to walk to and the fit is returned as it came out.
  skip_on_cran()
  skip_if_not_installed("MASS")
  mcycle <- MASS::mcycle
  fo <- accel ~ list(mu = ~ s(times, k = 20), sigma = ~ s(times, k = 10),
                     df = ~ 1)
  free <- suppressMessages(gamRTMB(fo, family = fam("t2"), data = mcycle))
  pinned <- suppressWarnings(gamRTMB(fo, family = fam("t2"), data = mcycle,
                                     start = list(log_sigma = free$obj$env$par[
                                       names(free$obj$env$par) == "log_sigma"])))
  expect_null(pinned$sigma_frac_used)
})
