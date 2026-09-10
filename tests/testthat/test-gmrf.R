## Helpers: a rook-adjacency lattice, optionally with one region cut loose.
lattice_nb <- function(g, island = FALSE) {
  ids <- as.vector(outer(seq_len(g), seq_len(g), function(i, j) paste0(i, "_", j)))
  nb <- stats::setNames(vector("list", g * g), ids)
  for (i in seq_len(g)) for (j in seq_len(g)) {
    nm <- paste0(i, "_", j); n <- character(0)
    if (i > 1) n <- c(n, paste0(i - 1, "_", j))
    if (i < g) n <- c(n, paste0(i + 1, "_", j))
    if (j > 1) n <- c(n, paste0(i, "_", j - 1))
    if (j < g) n <- c(n, paste0(i, "_", j + 1))
    nb[[nm]] <- n
  }
  if (island) {
    for (nm in names(nb)) nb[[nm]] <- setdiff(nb[[nm]], ids[1])
    nb[[ids[1]]] <- character(0)
  }
  nb
}

lattice_data <- function(g, rep = 4, seed = 1) {
  nb <- lattice_nb(g); ids <- names(nb)
  ij <- do.call(rbind, lapply(strsplit(ids, "_"), as.numeric))
  set.seed(seed)
  d <- data.frame(reg = factor(rep(ids, each = rep), levels = ids))
  k <- as.integer(factor(as.character(d$reg), levels = ids))
  d$y <- stats::rnorm(nrow(d), sin(ij[k, 1] / 2) + cos(ij[k, 2] / 2), 0.4)
  d$x <- stats::runif(nrow(d))
  d
}

mrf_smooth <- function(nb, absorb = FALSE, g = 6) {
  d <- data.frame(reg = factor(rep(names(nb), each = 2), levels = names(nb)))
  mgcv::smoothCon(mgcv::s(reg, bs = "mrf", xt = list(nb = nb)), data = d,
                  absorb.cons = absorb)[[1]]
}

is_pd <- function(Q)
  min(eigen(as.matrix(Q), symmetric = TRUE)$values) > 1e-10 * max(abs(Q))


test_that("a corner constraint leaves the penalty sparse and positive definite", {
  sm <- mrf_smooth(lattice_nb(8))
  expect_equal(sm$null.space.dim, 1)
  B <- .gmrf_block(sm)
  Q <- B$Q[[1]]
  expect_equal(ncol(Q), 63L)                 # one coefficient pinned
  expect_true(is_pd(Q))
  expect_lt(Matrix::nnzero(Q) / length(Q), 0.1)
  ## the whole point: absorbing the constraint instead would fill it in
  expect_equal(mean(mrf_smooth(lattice_nb(8), absorb = TRUE)$S[[1]] != 0), 1)
  ## a connected field's null space is exactly the constant, which the
  ## intercept carries, so nothing is handed to beta
  expect_equal(ncol(B$Xf), 0L)
  expect_true(B$intrinsic)
})

test_that("a disconnected graph keeps one free level per extra component", {
  sm <- mrf_smooth(lattice_nb(8, island = TRUE))
  expect_equal(sm$null.space.dim, 2)
  B <- .gmrf_block(sm)
  expect_equal(ncol(B$Q[[1]]), 62L)          # one per component
  expect_true(is_pd(B$Q[[1]]))
  expect_equal(ncol(B$Xf), 1L)               # one contrast between components
})

test_that("a proper precision is left alone", {
  nb <- lattice_nb(8); ids <- names(nb)
  A <- matrix(0, 64, 64, dimnames = list(ids, ids))
  for (nm in ids) A[nm, nb[[nm]]] <- 1
  sm <- mrf_smooth(nb)                        # replaced below
  d <- data.frame(reg = factor(rep(ids, each = 2), levels = ids))
  sm <- mgcv::smoothCon(mgcv::s(reg, bs = "mrf",
                                xt = list(penalty = diag(rowSums(A)) - 0.9 * A)),
                        data = d, absorb.cons = FALSE)[[1]]
  expect_equal(sm$null.space.dim, 0)
  B <- .gmrf_block(sm)
  expect_equal(ncol(B$Q[[1]]), 64L)          # nothing dropped
  expect_false(B$intrinsic)
  expect_equal(ncol(B$Xf), 0L)
})

test_that("the null space is found even when it is not just the constant", {
  d <- data.frame(x = stats::runif(200))
  sm <- mgcv::smoothCon(mgcv::s(x, bs = "ps", k = 30, m = c(2, 2)), data = d,
                        absorb.cons = FALSE)[[1]]
  expect_equal(sm$null.space.dim, 2)         # constant and linear
  B <- .gmrf_block(sm)
  expect_equal(ncol(B$Q[[1]]), 28L)
  expect_true(is_pd(B$Q[[1]]))
  expect_equal(ncol(B$Xf), 1L)               # linear survives, constant does not
})

test_that("the automatic rule fires on a field and not on an ordinary smooth", {
  expect_true(.use_sparse(mrf_smooth(lattice_nb(10))$S[[1]], "auto"))
  expect_false(.use_sparse(mrf_smooth(lattice_nb(6))$S[[1]], "auto"))  # too small
  d <- data.frame(x = stats::runif(200))
  S <- mgcv::smoothCon(mgcv::s(x, k = 20), data = d, absorb.cons = FALSE)[[1]]$S[[1]]
  expect_false(.use_sparse(S, "auto"))
  expect_true(.use_sparse(S, "always"))
  expect_false(.use_sparse(S, "never"))
})


## The contract that matters: the sparse route is a reparameterisation, not an
## approximation, so everything except the individual coefficients must agree.
expect_same_fit <- function(a, b, tol = 1e-5) {
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = tol)
  expect_equal(sum(edf(a)$edf), sum(edf(b)$edf), tolerance = tol)
  expect_equal(AIC(a), AIC(b), tolerance = tol)
  pa <- stats::predict(a, type = "response")
  pb <- stats::predict(b, type = "response")
  for (nm in names(pa)) expect_equal(pa[[nm]], pb[[nm]], tolerance = tol)
}

test_that("sparse and dense give the same Markov random field fit", {
  nb <- lattice_nb(7); d <- lattice_data(7)
  f <- y ~ list(mean = ~ s(reg, bs = "mrf", xt = list(nb = nb)), sd = ~ 1)
  expect_same_fit(gamRTMB(f, data = d, sparse = "never"),
                  gamRTMB(f, data = d, sparse = "always"))
})

test_that("sparse and dense agree with a field on more than one parameter", {
  nb <- lattice_nb(7); d <- lattice_data(7)
  f <- y ~ list(mean = ~ s(reg, bs = "mrf", xt = list(nb = nb)),
                sd   = ~ s(reg, bs = "mrf", xt = list(nb = nb)))
  expect_same_fit(gamRTMB(f, data = d, sparse = "never"),
                  gamRTMB(f, data = d, sparse = "always"))
})

test_that("sparse and dense agree on ordinary smooths too", {
  d <- lattice_data(7)
  for (f in list(y ~ list(mean = ~ s(x, k = 10), sd = ~ 1),
                 y ~ list(mean = ~ s(x, bs = "ps", k = 30, m = c(2, 2)), sd = ~ 1),
                 y ~ list(mean = ~ s(reg, bs = "re"), sd = ~ 1)))
    expect_same_fit(gamRTMB(f, data = d, sparse = "never"),
                    gamRTMB(f, data = d, sparse = "always"))
})

test_that("sparse and dense agree with an island", {
  nb <- lattice_nb(7, island = TRUE); d <- lattice_data(7)
  f <- y ~ list(mean = ~ s(reg, bs = "mrf", xt = list(nb = nb)), sd = ~ 1)
  expect_same_fit(gamRTMB(f, data = d, sparse = "never"),
                  gamRTMB(f, data = d, sparse = "always"))
})

test_that("predictions and their standard errors agree", {
  nb <- lattice_nb(7); d <- lattice_data(7)
  f <- y ~ list(mean = ~ s(reg, bs = "mrf", xt = list(nb = nb)), sd = ~ 1)
  a <- gamRTMB(f, data = d, sparse = "never")
  b <- gamRTMB(f, data = d, sparse = "always")
  nd <- data.frame(reg = factor(names(nb)[c(1, 5, 20, 49)], levels = names(nb)))
  pa <- stats::predict(a, newdata = nd, type = "link", se.fit = TRUE)
  pb <- stats::predict(b, newdata = nd, type = "link", se.fit = TRUE)
  expect_equal(pa$fit$mean, pb$fit$mean, tolerance = 1e-6)
  expect_equal(pa$se.fit$mean, pb$se.fit$mean, tolerance = 1e-6)
})

test_that("an intrinsic field on a parameter with no intercept is refused", {
  nb <- lattice_nb(7); d <- lattice_data(7)
  expect_error(
    gamRTMB(y ~ list(mean = ~ 0 + s(reg, bs = "mrf", xt = list(nb = nb)),
                     sd = ~ 1), data = d, sparse = "always"),
    "no intercept")
})

test_that("a sparse field works for a non-Gaussian family with an offset", {
  nb <- lattice_nb(8); ids <- names(nb)
  ij <- do.call(rbind, lapply(strsplit(ids, "_"), as.numeric))
  set.seed(3)
  sp <- sin(ij[, 1] / 2) + cos(ij[, 2] / 2)
  d <- data.frame(reg = factor(ids, levels = ids), E = stats::runif(64, 50, 300))
  d$y <- stats::rpois(64, d$E * exp(0.3 * sp))
  f <- y ~ list(lambda = ~ offset(log(E)) + s(reg, bs = "mrf", xt = list(nb = nb)))
  a <- gamRTMB(f, family = fam("pois"), data = d, sparse = "never")
  b <- gamRTMB(f, family = fam("pois"), data = d, sparse = "always")
  expect_true(b$convergence)
  expect_same_fit(a, b, tol = 1e-4)
  expect_gt(stats::cor(log(stats::predict(b, type = "response")$lambda) -
                         log(d$E), 0.3 * sp), 0.9)
})
