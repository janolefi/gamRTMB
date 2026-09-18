## Smooths that penalise one coefficient vector several times over: te(), ti()
## and bs = "ad". They share a route with the GMRF blocks -- reduce by the
## common null space, keep the penalties, hand the block to dgmrf() -- and the
## fact that makes it legitimate is that null(sum lambda_i S_i) is the
## intersection of the null(S_i) and so does not move with lambda.

sim2 <- function(n = 600, seed = 11) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n), z = stats::runif(n), w = stats::runif(n))
  d$y  <- stats::rnorm(n, sin(5 * d$x) * cos(4 * d$z), 0.3)      # interaction
  d$y2 <- stats::rnorm(n, sin(5 * d$x) + 2 * d$z, 0.3)           # additive
  d$y3 <- stats::rnorm(n, ifelse(d$x < 0.5, 0, 1) + sin(3 * d$x), 0.15)
  d
}

specs <- function(d)
  list(te = mgcv::te(x, z), ti = mgcv::ti(x, z),
       ad = mgcv::s(x, bs = "ad", k = 25))

con <- function(spec, d)
  mgcv::smoothCon(spec, data = d, absorb.cons = FALSE,
                  null.space.penalty = FALSE)[[1L]]


## ---------------------------------------------------------------------------
## Which route a smooth takes, and why

test_that("overlap is what sends a smooth down the sparse route", {
  d <- sim2(200)
  ## t2 penalises disjoint sets of coefficients, which is what lets
  ## smooth2random split it into one iid block per penalty
  t2sm <- con(mgcv::t2(x, z, k = 4), d)
  expect_gt(length(t2sm$S), 1L)
  expect_false(.penalties_overlap(t2sm$S))
  expect_false(.use_sparse(t2sm, "auto"))
  expect_false(.use_sparse(t2sm, "always"))

  for (sm in lapply(specs(d), con, d = d)) {
    expect_gt(length(sm$S), 1L)
    expect_true(.penalties_overlap(sm$S))
    expect_true(.use_sparse(sm, "auto"))
    expect_false(.use_sparse(sm, "never"))
  }
})

test_that("mgcv itself declines these smooths, which is the reason for the route", {
  d <- sim2(200)
  for (sm in lapply(specs(d), function(s)
        mgcv::smoothCon(s, data = d, absorb.cons = TRUE)[[1L]]))
    expect_error(mgcv::smooth2random(sm, "", type = 2))
})


## ---------------------------------------------------------------------------
## The reduction

test_that("the common null space is the intersection and does not move with lambda", {
  d <- sim2(200)
  for (sm in lapply(specs(d), con, d = d)) {
    ns <- .null_space(lapply(sm$S, .spm), as.integer(sm$null.space.dim), sm$label)
    expect_identical(ncol(ns$N), as.integer(sm$null.space.dim))
    expect_identical(length(ns$drop), as.integer(sm$null.space.dim))
    ## S_i N = 0 for every penalty separately, which is what lets one
    ## reparameterisation serve all of them
    for (S in sm$S) expect_lt(max(abs(S %*% ns$N)), 1e-10 * max(abs(S)))
    ## and the null space really is lambda-free: the rank of the sum is the
    ## same at wildly different weights
    q <- ncol(sm$X)
    for (pw in c(0, 6, 12)) {
      lam <- 10^(pw * (seq_along(sm$S) - 1) / max(1, length(sm$S) - 1))
      M <- Reduce(`+`, Map(`*`, sm$S, lam))
      expect_lt(max(abs(M %*% ns$N)), 1e-6 * max(abs(M)))
    }
  }
})

test_that("the reduced penalty is positive definite at any weights, and stays sparse", {
  d <- sim2(200)
  for (sm in lapply(specs(d), con, d = d)) {
    B <- .gmrf_block(sm)
    sp <- B$spec[[1L]]
    expect_identical(sp$kind, "multi")
    expect_identical(sp$ntheta, length(sm$S))
    expect_identical(sp$L, diag(length(sm$S)))         # one free log lambda each
    expect_equal(sp$sp_pow, rep(1, length(sm$S)))
    for (pw in c(-8, 0, 8)) {
      Q <- .block_prec(sp, pw * (seq_len(sp$ntheta) - 1) / max(1, sp$ntheta - 1))
      expect_true(min(eigen(as.matrix(Q), symmetric = TRUE)$values) >
                    1e-12 * max(abs(Q)))
    }
  }
  ## absorbing mgcv's constraint instead is what would destroy the sparsity
  smc <- mgcv::smoothCon(mgcv::te(x, z, k = c(10, 10)), data = sim2(200),
                         absorb.cons = TRUE)[[1L]]
  expect_equal(mean(smc$S[[1L]] != 0), 1)
  smu <- con(mgcv::te(x, z, k = c(10, 10)), sim2(200))
  Q <- .gmrf_block(smu)$spec[[1L]]$Smats[[1L]]
  expect_lt(Matrix::nnzero(as(Q, "generalMatrix")) / prod(dim(Q)), 0.25)
})

test_that("a null space mgcv has mis-stated is refused rather than fitted", {
  d <- sim2(200)
  sm <- con(mgcv::te(x, z), d)
  expect_identical(as.integer(sm$null.space.dim), 4L)
  expect_error(.null_space(lapply(sm$S, .spm), 6L, sm$label), "null space")
  expect_error(.null_space(lapply(sm$S, .spm), 2L, sm$label), "null space")
})

test_that("a te() hands its main effects to beta and its constant to the intercept", {
  d <- sim2(200)
  B <- .gmrf_block(con(mgcv::te(x, z), d))
  ## the null space of a tensor product is {1, x, z, xz}; the intercept takes
  ## the constant and the other three become unpenalized columns
  expect_identical(ncol(B$Xf), 3L)
  expect_true(B$intrinsic)
  ## which is why a te() plus a main effect for one of its margins is caught
  expect_error(.build_design(list(mu = ~ te(x, z, k = 3) + s(z)), d, "mu"),
               "rank deficient")
})


## ---------------------------------------------------------------------------
## Against mgcv

test_that("te, ti and adaptive fits agree with mgcv", {
  skip_on_cran()
  d <- sim2()
  agrees <- function(rhs, gterm, resp, tol_edf = 0.01, tol_fit = 0.01) {
    f <- gamRTMB(stats::reformulate(rhs, response = resp), data = d)
    g <- mgcv::gam(list(stats::reformulate(gterm, response = resp), ~ 1),
                   family = mgcv::gaulss(), data = d)
    ge <- as.numeric(summary(g)$edf); ge <- ge[ge > 1e-6]
    expect_true(f$convergence)
    expect_equal(edf(f)$edf, ge, tolerance = tol_edf)
    expect_equal(unname(stats::predict(f, type = "response")$mean),
                 unname(stats::predict(g, type = "response")[, 1L]),
                 tolerance = tol_fit)
  }
  agrees("te(x, z)", "te(x, z)", "y")
  agrees("te(x, z)", "te(x, z)", "y2")            # one margin driven flat
  agrees("s(x) + s(z) + ti(x, z)", "s(x) + s(z) + ti(x, z)", "y2")
  agrees("s(x, bs = \"ad\", k = 25)", "s(x, bs = \"ad\", k = 25)", "y3")
})

test_that("qREML reaches the same place as REML on a te()", {
  skip_on_cran()
  d <- sim2()
  a <- gamRTMB(y ~ list(mean = ~ te(x, z), sd = ~ 1), data = d)
  b <- gamRTMB(y ~ list(mean = ~ te(x, z), sd = ~ 1), data = d, method = "qREML")
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-3)
  expect_equal(sum(edf(a)$edf), sum(edf(b)$edf), tolerance = 0.01)
})


## ---------------------------------------------------------------------------
## The Fellner-Schall step, and the bounds

test_that("free log lambdas get the exact multiplicative update, an SPDE does not", {
  d <- sim2(200)
  D <- .build_design(list(mu = ~ te(x, z, k = 4)), d, "mu")
  pm <- .efs_pairs(D)
  ## L = I, so each weight moves on its own and Wood & Fasiolo's update is
  ## exact -- sharing a block with other penalties does not tie a parameter
  expect_identical(pm$nfree, 2L)
  expect_false(any(pm$tied))
  expect_true(all(pm$multi))

  ## an L that is not diagonal does tie them
  sm <- con(mgcv::te(x, z, k = 4), d)
  sm$L <- matrix(c(-2, -2, 0, 2), 2L, 2L)
  sm$theta.names <- c("a", "b"); sm$theta.start <- c(0, 0)
  sp <- .penalty_spec(sm, lapply(sm$S, .spm))
  expect_null(sp$sp_pow)
})

test_that("a penalty weight cannot walk off to nothing", {
  skip_on_cran()
  ## mcycle's adaptive smooth has a weight the criterion is exactly flat in:
  ## unbounded, nlminb walks it to log lambda = -466 and reports failure
  data(mcycle, package = "MASS")
  f <- gamRTMB(accel ~ s(times, bs = "ad", k = 20), data = mcycle)
  bl <- f$design$blocks[[1L]]
  expect_true(f$convergence)
  expect_lte(diff(range(f$log_sigma[bl$theta_idx])), 2 * .theta_halfwidth + 1e-6)
  ## the box is what stops it: unbounded, one weight reaches log lambda = -466
  b <- .theta_bounds(f$design, .init_pars(f$design, f$family,
                                          mcycle$accel)$log_sigma)
  expect_true(all(is.finite(b$lower)))
  expect_gte(min(f$log_sigma[bl$theta_idx]), min(b$lower) - 1e-6)
  ## and the fit itself is the one mgcv finds
  g <- mgcv::gam(list(accel ~ s(times, bs = "ad", k = 20), ~ 1),
                 family = mgcv::gaulss(), data = mcycle)
  expect_equal(sum(edf(f)$edf), sum(as.numeric(summary(g)$edf)), tolerance = 0.01)
})

test_that("an id ties a multi-penalty smooth entry by entry", {
  d <- sim2(300)
  D <- .build_design(list(mu = ~ te(x, z, k = 4, id = 1),
                          sg = ~ te(x, z, k = 4, id = 1)),
                     d, c("mu", "sg"))
  expect_identical(D$nsigma, 4L)          # two smooths, two penalties each
  expect_identical(D$nsigma_free, 2L)     # tied entry by entry, not as a pair
})


## ---------------------------------------------------------------------------
## Reporting

test_that("the smoothing parameters are reported on the right scale", {
  d <- sim2(300)
  ## s() on a third covariate, since a te()'s null space already carries the
  ## main effects of its own margins
  f <- gamRTMB(y ~ list(mean = ~ te(x, z, k = 4) + s(w, k = 6), sd = ~ 1), data = d)
  bl <- f$design$blocks
  te_bl <- bl[[which(vapply(bl, function(z) z$kind == "multi", NA))]]
  s_bl  <- bl[[which(vapply(bl, function(z) z$kind == "iid", NA))]]
  ## theta is log lambda for the tensor product ...
  expect_equal(as.numeric(sub(",.*", "", .sp_show(te_bl, f$log_sigma)[1L])),
               exp(f$log_sigma[te_bl$theta_idx[1L]]), tolerance = 1e-3)
  ## ... and log sigma for the ordinary smooth, reported as sigma^-2
  expect_equal(as.numeric(.sp_show(s_bl, f$log_sigma)),
               exp(-2 * f$log_sigma[s_bl$theta_idx]), tolerance = 1e-3)
  e <- edf(f)
  expect_identical(nrow(e), 2L)
  expect_true(grepl(",", e$sp[e$term == "te(x,z)"]))   # two weights, one row
  expect_identical(length(coef(f)$log_sigma), 3L)
})
