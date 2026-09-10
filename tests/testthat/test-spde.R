skip_if_not_installed("fmesher")

spde_data <- function(n = 600, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n), y = stats::runif(n))
  d$truth <- sin(5 * d$x) * cos(5 * d$y)
  d$z <- stats::rnorm(n, d$truth, 0.25)
  d
}

spde_mesh <- function(d)
  fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.14, 0.4),
                      cutoff = 0.06, offset = c(0.1, 0.3))


test_that("the smooth is built from fmesher and stays sparse", {
  d <- spde_data(); mesh <- spde_mesh(d)
  sm <- mgcv::smoothCon(mgcv::s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        data = d, absorb.cons = FALSE)[[1]]
  expect_s3_class(sm, "spde.smooth")
  expect_true(inherits(sm$X, "Matrix"))            # not densified
  expect_lt(Matrix::nnzero(sm$X) / length(sm$X), 0.02)
  expect_length(sm$S, 3L)
  expect_equal(dim(sm$L), c(3L, 2L))
  expect_equal(sm$null.space.dim, 0L)              # proper for any kappa > 0
  expect_equal(nrow(sm$C), 0L)                     # and so left unconstrained
})

test_that("the precision is exactly tau^2 (kappa^4 C + 2 kappa^2 G1 + G2)", {
  d <- spde_data(); mesh <- spde_mesh(d)
  sm <- mgcv::smoothCon(mgcv::s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        data = d, absorb.cons = FALSE)[[1]]
  B <- .gmrf_block(sm)
  bl <- c(B$spec[[1]], list(q = ncol(sm$X)))
  ## The block is parameterised by (log sd, log range); this checks that the
  ## constants folded into S and the powers in L put it back exactly onto the
  ## closed form in (tau, kappa).
  tau <- 0.35; kappa <- 6; k2 <- kappa^2
  fem <- fmesher::fm_fem(mesh)
  want <- tau^2 * (k2 * k2 * fem$c0 + 2 * k2 * fem$g1 + fem$g2)
  got <- .block_prec(bl, log(c(1 / (tau * kappa * sqrt(4 * pi)),
                               2 * sqrt(2) / kappa)))
  expect_lt(max(abs(as.matrix(got - want))) / max(abs(want)), 1e-12)
  expect_lt(Matrix::nnzero(got) / length(got), 0.1)
})

test_that("an SPDE term gets two named smoothing parameters", {
  d <- spde_data(); mesh <- spde_mesh(d)
  f <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        sd = ~ 1), data = d)
  expect_true(f$convergence)
  expect_equal(f$design$nsigma, 2L)
  expect_equal(f$design$blocks[[1]]$theta_names, c("sd", "range"))
  expect_match(edf(f)$sp[1], "sd=.*range=")
  expect_match(names(coef(f)$log_sigma)[2], "range")
})

test_that("a simulated Matern field is recovered", {
  set.seed(42)
  n <- 1200
  d <- data.frame(x = stats::runif(n), y = stats::runif(n))
  mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.08, 0.25),
                              cutoff = 0.03, offset = c(0.15, 0.4))
  fem <- fmesher::fm_fem(mesh)
  tau <- 0.35; kappa <- 6; k2 <- kappa^2
  Q <- methods::as(methods::as(tau^2 * (k2 * k2 * fem$c0 + 2 * k2 * fem$g1 + fem$g2),
                               "symmetricMatrix"), "CsparseMatrix")
  L <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)
  set.seed(7)
  b <- as.vector(Matrix::solve(L, Matrix::solve(L, stats::rnorm(ncol(Q)),
                                                system = "Lt"), system = "Pt"))
  fld <- as.vector(fmesher::fm_basis(mesh, cbind(d$x, d$y)) %*% b)
  d$z <- fld + stats::rnorm(n, 0, 0.05)            # field must beat the noise
  f <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        sd = ~ 1), data = d)
  ## One realisation of a spatial field identifies its parameters only
  ## loosely, so this asks for a factor of two, which is still tight enough to
  ## catch a wrong power in L or a missing constant. The field itself is
  ## recovered much more sharply than its parameters are.
  th <- exp(coef(f)$log_sigma)
  true_sd <- 1 / (tau * kappa * sqrt(4 * pi))
  true_range <- 2 * sqrt(2) / kappa
  expect_gt(unname(th[1]) / true_sd, 0.5);    expect_lt(unname(th[1]) / true_sd, 2)
  expect_gt(unname(th[2]) / true_range, 0.5); expect_lt(unname(th[2]) / true_range, 2)
  expect_gt(stats::cor(stats::predict(f, type = "response")$mean, fld), 0.95)
})

test_that("a field can go on a parameter other than the mean", {
  set.seed(4)
  d <- data.frame(x = stats::runif(600), y = stats::runif(600))
  lsd <- -1 + 0.8 * sin(4 * d$x)
  d$z <- stats::rnorm(600, 0, exp(lsd))
  mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.2, 0.5),
                              cutoff = 0.08, offset = c(0.1, 0.3))
  f <- gamRTMB(z ~ list(mean = ~ 1,
                        sd = ~ s(x, y, bs = "spde", xt = list(mesh = mesh))),
               data = d)
  expect_true(f$convergence)
  expect_gt(stats::cor(log(stats::predict(f, type = "response")$sd), lsd), 0.9)
})

test_that("a field works for a non-Gaussian family", {
  set.seed(5)
  d <- data.frame(x = stats::runif(700), y = stats::runif(700))
  lam <- exp(1 + sin(5 * d$x) * cos(5 * d$y))
  d$cnt <- stats::rpois(700, lam)
  mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y), max.edge = c(0.15, 0.4),
                              cutoff = 0.06, offset = c(0.1, 0.3))
  f <- gamRTMB(cnt ~ list(lambda = ~ s(x, y, bs = "spde", xt = list(mesh = mesh))),
               family = fam("pois"), data = d)
  expect_true(f$convergence)
  expect_gt(stats::cor(log(stats::predict(f, type = "response")$lambda),
                       log(lam)), 0.9)
})

test_that("predictions at new locations work and carry standard errors", {
  d <- spde_data(); mesh <- spde_mesh(d)
  f <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        sd = ~ 1), data = d)
  nd <- data.frame(x = c(0.2, 0.5, 0.8), y = c(0.3, 0.5, 0.7))
  p <- stats::predict(f, newdata = nd, type = "link", se.fit = TRUE)
  expect_length(p$fit$mean, 3L)
  expect_true(all(p$se.fit$mean > 0))
  expect_lt(max(abs(p$fit$mean - sin(5 * nd$x) * cos(5 * nd$y))), 0.35)
})

test_that("a one-dimensional field builds its own mesh", {
  set.seed(6)
  d <- data.frame(t = stats::runif(400))
  d$y <- stats::rnorm(400, sin(2 * pi * d$t), 0.3)
  f <- gamRTMB(y ~ list(mean = ~ s(t, bs = "spde", k = 40), sd = ~ 1), data = d)
  expect_true(f$convergence)
  expect_gt(stats::cor(stats::predict(f, type = "response")$mean,
                       sin(2 * pi * d$t)), 0.95)
})

test_that("a field sits alongside parametric and ordinary smooth terms", {
  d <- spde_data(); mesh <- spde_mesh(d)
  set.seed(8)
  d$w <- stats::runif(nrow(d))
  d$g <- factor(sample(c("a", "b"), nrow(d), TRUE))
  d$z <- d$truth + 0.5 * d$w + ifelse(d$g == "b", 0.4, 0) +
    stats::rnorm(nrow(d), 0, 0.25)
  f <- gamRTMB(z ~ list(mean = ~ g + s(w, k = 8) +
                          s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        sd = ~ 1), data = d)
  expect_true(f$convergence)
  expect_equal(nrow(edf(f)), 2L)
  expect_equal(f$design$nsigma, 3L)               # one for s(w), two for the field
})

test_that("the constructor refuses what it cannot do", {
  d <- spde_data()
  expect_error(gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde"), sd = ~ 1),
                       data = d), "needs a mesh")
  expect_error(gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde",
                                           xt = list(mesh = "nope")), sd = ~ 1),
                       data = d), "fmesher mesh")
  mesh <- spde_mesh(d)
  expect_error(
    gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh)),
                     sd = ~ 1), data = d, sparse = "never"),
    "smooth2random")
})


test_that("the field's standard deviation and range converge with more data", {
  skip_on_cran()
  true_sd <- 1; true_range <- 0.3
  one <- function(n, max_edge) {
    set.seed(1)
    d <- data.frame(x = stats::runif(n), y = stats::runif(n))
    mesh <- fmesher::fm_mesh_2d(loc = cbind(d$x, d$y),
                                max.edge = c(max_edge, 3 * max_edge),
                                cutoff = max_edge / 2, offset = c(0.2, 0.5))
    fem <- fmesher::fm_fem(mesh)
    kap <- 2 * sqrt(2) / true_range
    tau <- 1 / (true_sd * kap * sqrt(4 * pi)); k2 <- kap^2
    Q <- methods::as(methods::as(tau^2 * (k2 * k2 * fem$c0 + 2 * k2 * fem$g1 +
                                            fem$g2), "symmetricMatrix"),
                     "CsparseMatrix")
    L <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)
    set.seed(9)
    b <- as.vector(Matrix::solve(L, Matrix::solve(L, stats::rnorm(ncol(Q)),
                                                  system = "Lt"), system = "Pt"))
    fld <- as.vector(fmesher::fm_basis(mesh, cbind(d$x, d$y)) %*% b)
    d$z <- fld + stats::rnorm(n, 0, 0.3)
    f <- gamRTMB(z ~ list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh)),
                          sd = ~ 1), data = d)
    list(fit = f, th = exp(coef(f)$log_sigma), fld = fld)
  }
  small <- one(2000, 0.10)
  big <- one(10000, 0.06)
  expect_true(small$fit$convergence)
  expect_true(big$fit$convergence)
  ## closer to the truth with more data, on both parameters
  expect_lt(abs(big$th[1] - true_sd), abs(small$th[1] - true_sd))
  expect_lt(abs(big$th[2] - true_range), abs(small$th[2] - true_range))
  expect_gt(stats::cor(stats::predict(big$fit, type = "response")$mean,
                       big$fld), 0.99)
})
