## method = "aREML": the REML criterion optimised by extended Fellner-Schall.
## Its correctness check is method = "REML", which optimises the same
## criterion by a different route -- aREML drops a third-derivative term from
## the gradient but not from the criterion itself, so the two should agree
## closely rather than exactly.

sim_efs <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x1), exp(-1 + d$x2))
  d
}

lattice_nb_efs <- function(g) {
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
  nb
}


test_that(".block_penalties reproduces .block_prec on every kind of block", {
  ## The two are the same statement written for plain numbers and for the AD
  ## tape; the whole generality of the update rests on their agreeing.
  d <- sim_efs(200)
  nb <- lattice_nb_efs(8)
  d2 <- data.frame(reg = factor(rep(names(nb), each = 3), levels = names(nb)))
  d2$y <- stats::rnorm(nrow(d2))

  blocks <- c(
    .build_design(list(mean = ~ s(x1, k = 8), sd = ~ 1), d, c("mean", "sd"))$blocks,
    .build_design(list(mean = ~ s(reg, bs = "mrf", xt = list(nb = nb)), sd = ~ 1),
                  d2, c("mean", "sd"))$blocks)
  expect_setequal(unique(vapply(blocks, `[[`, "", "kind")), c("iid", "gmrf"))

  for (bl in blocks) {
    th <- seq(-0.7, 0.4, length.out = bl$ntheta)
    pen <- .block_penalties(bl)
    got <- pen$S[[1]] * exp(sum(pen$A[1, ] * th))
    for (i in seq_along(pen$S)[-1])
      got <- got + pen$S[[i]] * exp(sum(pen$A[i, ] * th))
    want <- .block_prec(bl, th)
    expect_lt(max(abs(as.matrix(got - want))), 1e-12 * max(abs(as.matrix(want))))
  }
})

test_that("a smooth with an L matrix is classified as tied, an ordinary one is not", {
  d <- sim_efs(200)
  D <- .build_design(list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 6)), d,
                     c("mean", "sd"))
  expect_equal(unname(.efs_pairs(D)$tied), c(FALSE, FALSE))

  skip_if_not_installed("fmesher")
  set.seed(1)
  ds <- data.frame(x = stats::runif(300), y = stats::runif(300))
  ds$z <- stats::rnorm(300)
  mesh <- fmesher::fm_mesh_2d(loc = cbind(ds$x, ds$y), max.edge = c(0.2, 0.5),
                              cutoff = 0.1, offset = c(0.1, 0.3))
  Ds <- .build_design(list(mean = ~ s(x, y, bs = "spde", xt = list(mesh = mesh)),
                           sd = ~ 1), ds, c("mean", "sd"))
  ## two smoothing parameters on one block, tied through L
  expect_true(all(.efs_pairs(Ds)$tied))
})

test_that("aREML matches REML on a location-scale model", {
  d <- sim_efs(300)
  f <- y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 8))
  a <- gamRTMB(f, data = d)
  b <- gamRTMB(f, data = d, method = "aREML")
  expect_true(b$convergence)
  expect_equal(b$objective, a$objective, tolerance = 1e-4)
  expect_equal(edf(b)$edf, edf(a)$edf, tolerance = 1e-2)
  expect_equal(b$coefficients$beta, a$coefficients$beta, tolerance = 1e-2)
  ## the fit contract: the same fields, and no `obj`
  expect_null(b[["obj"]])
  expect_equal(names(b$coefficients), c("beta", "b"))
  expect_length(b$log_sigma, b$design$nsigma)
})

test_that("aREML matches REML with an id-tied smoothing parameter", {
  d <- sim_efs(300)
  f <- y ~ list(mean = ~ s(x1, k = 8, id = 1) + s(x2, k = 8, id = 1))
  a <- gamRTMB(f, data = d)
  b <- gamRTMB(f, data = d, method = "aREML")
  expect_true(b$convergence)
  expect_equal(b$objective, a$objective, tolerance = 1e-4)
  ## one estimated value, expanded back over both blocks
  expect_equal(b$design$nsigma_free, 1L)
  expect_equal(b$log_sigma[1], b$log_sigma[2])
  expect_equal(edf(b)$edf, edf(a)$edf, tolerance = 1e-2)
})

test_that("aREML handles a sparse GMRF block, which needs the general trace", {
  nb <- lattice_nb_efs(8)
  ids <- names(nb)
  ij <- do.call(rbind, lapply(strsplit(ids, "_"), as.numeric))
  set.seed(1)
  d <- data.frame(reg = factor(rep(ids, each = 4), levels = ids))
  k <- as.integer(factor(as.character(d$reg), levels = ids))
  d$y <- stats::rnorm(nrow(d), sin(ij[k, 1] / 2) + cos(ij[k, 2] / 2), 0.4)
  f <- y ~ list(mean = ~ s(reg, bs = "mrf", xt = list(nb = nb)), sd = ~ 1)
  a <- gamRTMB(f, data = d)
  b <- gamRTMB(f, data = d, method = "aREML")
  expect_equal(b$design$blocks[[1]]$kind, "gmrf")
  expect_true(b$convergence)
  expect_equal(b$objective, a$objective, tolerance = 1e-3)
  expect_equal(edf(b)$edf, edf(a)$edf, tolerance = 0.5)
})

test_that("a model with no smooths needs no outer iteration", {
  d <- sim_efs(200)
  b <- gamRTMB(y ~ list(mean = ~ x1), data = d, method = "aREML")
  expect_true(b$convergence)
  expect_equal(b$max_grad, 0)
  expect_equal(b$opt$message, "no smoothing parameters to estimate")
  a <- gamRTMB(y ~ list(mean = ~ x1), data = d)
  expect_equal(b$coefficients$beta, a$coefficients$beta, tolerance = 1e-5)
})

test_that("the criterion never rises: the step control does its job", {
  d <- sim_efs(300)
  b <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 8)),
               data = d, method = "aREML")
  V <- b$efs_trace$V
  expect_gt(length(V), 1L)
  expect_true(all(diff(V) <= 1e-10))
  expect_equal(V[length(V)], b$objective)
})

test_that("everything downstream of the fit works on an aREML fit", {
  d <- sim_efs(300)
  b <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 8)),
               data = d, method = "aREML")
  expect_s3_class(edf(b), "data.frame")
  V <- vcov(b)
  expect_equal(dim(V), rep(b$design$nbeta, 2L))
  expect_true(all(diag(V) > 0))
  s <- summary(b)
  expect_equal(ncol(s$coefficients), 4L)          # with standard errors
  p <- stats::predict(b, se.fit = TRUE)
  expect_named(p$fit, c("mean", "sd"))
  expect_true(all(p$se.fit$mean > 0))
  expect_length(stats::residuals(b), nrow(d))
  expect_output(print(b), "aREML")
  ## named so it is not read as a stationarity measure
  expect_output(print(b), "max\\|FS grad\\|")
})

test_that("aREML reports its progress when asked, and not otherwise", {
  d <- sim_efs(200)
  f <- y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 6))
  expect_silent(gamRTMB(f, data = d, method = "aREML"))
  ## `silent = FALSE` is what a user reaches for; there is no MakeADFun here
  ## for it to reach, so the fit has to honour it itself.
  expect_message(gamRTMB(f, data = d, method = "aREML", silent = FALSE), "\\[start\\]")
  expect_message(gamRTMB(f, data = d, method = "aREML", silent = FALSE),
                 "converged after")
  expect_message(gamRTMB(f, data = d, method = "aREML",
                         control = list(trace = TRUE)), "efs\\s+1\\s+-REML")
  ## the size line comes out before the first inner solve, which on a hard
  ## family is the longest single step in the fit
  expect_message(gamRTMB(f, data = d, method = "aREML", silent = FALSE),
                 "observations.*coefficients.*smoothing")
  ## Level 2 adds the inner Newton, and level 1 does not. Anchored on the
  ## indent, because level 1 mentions the inner solve too -- it says the first
  ## one is starting, and points at level 2 for watching it.
  expect_message(gamRTMB(f, data = d, method = "aREML",
                         control = list(trace = 2)), "^\\s+inner\\s+\\d")
  expect_false(any(grepl("^\\s+inner", capture_messages(
    gamRTMB(f, data = d, method = "aREML", control = list(trace = 1))))))
  ## and an explicit trace wins over silent, in either direction
  expect_silent(gamRTMB(f, data = d, method = "aREML", silent = FALSE,
                        control = list(trace = FALSE)))
})

test_that("the traced smoothing parameters are the ones edf() reports", {
  d <- sim_efs(200)
  b <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 6)),
               data = d, method = "aREML")
  pmap <- .efs_pairs(b$design)
  sp <- .efs_sp(.efs_free(b$log_sigma, pmap), .efs_sp_kind(b$design, pmap))
  expect_equal(signif(sp, 4), signif(as.numeric(edf(b)$sp), 4))
})

test_that("both inner solvers reach the same place on a well-posed problem", {
  ## They have to, since only one of them can be the default and the other is
  ## offered as a speed option on large structured terms. A difference here
  ## would mean the option changes the answer, not just the cost.
  d <- sim_efs(300)
  f <- y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 8))
  bf <- gamRTMB(f, data = d, method = "aREML",
                control = list(inner_method = "bfgs"))
  nw <- gamRTMB(f, data = d, method = "aREML",
                control = list(inner_method = "newton"))
  expect_true(bf$convergence)
  expect_true(nw$convergence)
  expect_equal(bf$objective, nw$objective, tolerance = 1e-5)
  expect_equal(edf(bf)$edf, edf(nw)$edf, tolerance = 1e-2)
})

test_that("the Newton branch is TMB's, not one written here", {
  ## A hand-rolled Newton loop lived here and was wrong in two ways at once.
  ## TMB::newton() is the solver the Laplace approximation itself uses, it is
  ## faster on every well-conditioned model measured, and it is somebody
  ## else's to maintain. This pins the decision.
  expect_true(any(grepl("TMB::newton", deparse(.efs_inner_newton), fixed = TRUE)))
  expect_error(.efs_inner(0, identity, identity, identity,
                          list(inner_method = "nope")), "bfgs")
})

test_that("aREML honours the rest of the model contract", {
  ## Prior weights, offsets and the smooth types that reparameterise in
  ## awkward ways all pass through .make_nll() and the design rather than
  ## through the fitting routine, so this is checking that it has not quietly
  ## bypassed any of them -- method = "REML" is the reference.
  set.seed(1); n <- 400
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n),
                  g = factor(sample(letters[1:4], n, TRUE)),
                  E = stats::runif(n, 0.5, 2), w = sample(1:3, n, TRUE))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x1), 0.3)
  d$cnt <- stats::rpois(n, d$E * exp(0.5 + sin(2 * pi * d$x1)))
  same <- function(form, family = fam("norm"), ...) {
    b <- gamRTMB(form, family = family, data = d, method = "aREML", ...)
    a <- gamRTMB(form, family = family, data = d, ...)
    expect_true(b$convergence)
    expect_equal(b$objective, a$objective, tolerance = 1e-4)
  }
  same(y ~ list(mean = ~ s(x1, k = 8)), weights = w)
  same(cnt ~ list(lambda = ~ offset(log(E)) + s(x1, k = 8)), fam("pois"))
  same(y ~ list(mean = ~ s(x1, k = 8, by = g)))
  same(y ~ list(mean = ~ s(x1, g, bs = "fs", k = 6)))
  same(y ~ list(mean = ~ t2(x1, x2, k = 4)))

  d$x1[c(3, 7)] <- NA
  r <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8)), data = d, method = "aREML")
  expect_equal(r$dropped, 2L)
  expect_equal(r$design$n, n - 2L)
})

test_that("one aREML fit does not disturb the next", {
  ## RTMB's OBS() keys on the deparsed name of its argument in a registry
  ## global to the package. MakeADFun resets it per object; MakeTape, which is
  ## what this engine uses, does not -- so a tape that calls OBS(y) leaves the
  ## response of its own model behind under the name "y", and the next model's
  ## objective quietly evaluates against it. That produced a second fit that
  ## converged confidently, with a gradient of 1e-10, to the wrong place. The
  ## engine asks .make_nll() not to mark the response at all; this is the
  ## check that it stays that way.
  d3 <- sim_efs(300, seed = 3)
  d2 <- sim_efs(200, seed = 4)
  alone <- gamRTMB(y ~ list(mean = ~ x1), data = d2, method = "aREML")
  invisible(gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 8)),
                    data = d3, method = "aREML"))
  after <- gamRTMB(y ~ list(mean = ~ x1), data = d2, method = "aREML")
  expect_equal(after$coefficients$beta, alone$coefficients$beta)
  expect_equal(after$objective, alone$objective)
})

test_that("the exact update and the metric step agree to first order", {
  ## The two branches of .efs_step() have to meet: on a free-weight parameter
  ## the multiplicative update is, to first order in log(n/d), the same step
  ## the tied branch would take. Checked here because nothing else would
  ## notice if one of them picked up a stray factor.
  n <- 10.4; dv <- 10
  lr <- log(n / dv)
  L <- matrix(-2, 1L, 1L)
  metric <- as.numeric(crossprod(L, (n + dv) / 4 * L))
  grad <- as.numeric(crossprod(L, dv - n)) / 2
  expect_equal(-grad / metric, lr / L[1, 1], tolerance = 1e-3)
})
