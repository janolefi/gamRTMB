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

test_that("healthy EDF sit inside the [0, 1] the guard checks", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  ph <- .penalized_hessian(fit)
  e <- 1 - Matrix::diag(Matrix::solve(ph$H,
                                      .penalty_matrix(fit$design, fit$log_sigma, ph)))
  expect_gte(min(e), -.edf_tol)
  expect_lte(max(e), 1 + .edf_tol)
  ## and the guard therefore stays out of the way
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
