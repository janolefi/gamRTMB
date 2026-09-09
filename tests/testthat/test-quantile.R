test_that("fitted quantiles recover the truth and are calibrated", {
  set.seed(1); n <- 900L
  d <- data.frame(x = runif(n))
  mu <- sin(2 * pi * d$x); sg <- exp(-1 + 0.8 * cos(2 * pi * d$x))
  d$y <- stats::rnorm(n, mu, sg)
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 10), sd = ~ s(x, k = 10)), data = d)
  pr <- c(0.05, 0.5, 0.95)
  q <- stats::predict(f, type = "quantile", prob = pr)
  expect_identical(dim(q), c(n, 3L))
  expect_identical(colnames(q), c("q0.05", "q0.5", "q0.95"))
  ## quantiles must be ordered
  expect_true(all(q[, 1] < q[, 2] & q[, 2] < q[, 3]))
  ## and calibrated against the data that generated them
  for (k in seq_along(pr)) expect_equal(mean(d$y < q[, k]), pr[k], tolerance = 0.03)
  ## close to the true quantiles
  tq <- vapply(pr, function(p) stats::qnorm(p, mu, sg), numeric(n))
  expect_lt(max(abs(q - tq)), 0.5 * stats::sd(d$y))
})

test_that("quantile standard errors are right, by two independent checks", {
  set.seed(1); n <- 800L
  d <- data.frame(x = runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x), exp(-1 + 0.8 * cos(2 * pi * d$x)))
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 10), sd = ~ s(x, k = 10)), data = d)
  qs <- stats::predict(f, type = "quantile", prob = c(0.1, 0.5), se.fit = TRUE)
  expect_identical(dim(qs$se.fit), c(n, 2L))
  expect_true(all(qs$se.fit > 0))

  ## (1) for a symmetric family the median IS the mean parameter, so the
  ## delta-method SE must reproduce that parameter's SE exactly
  expect_equal(as.numeric(qs$se.fit[, 2]),
               as.numeric(stats::predict(f, se.fit = TRUE)$se.fit$mean),
               tolerance = 1e-8)

  ## (2) against simulation from the joint posterior of the coefficients
  Vj <- gamRTMB:::.joint_cov(f)
  set.seed(2); B <- 1500
  L <- chol(Vj$V + diag(1e-12, nrow(Vj$V)))
  i <- 40L
  draws <- vapply(seq_len(B), function(b) {
    z <- as.vector(crossprod(L, stats::rnorm(nrow(Vj$V))))
    f2 <- f
    f2$coefficients$beta <- f$coefficients$beta + z[Vj$ib]
    f2$coefficients$b <- f$coefficients$b + z[Vj$ir]
    stats::predict(f2, newdata = d[i, , drop = FALSE], type = "quantile",
                   prob = 0.1)[1]
  }, 0)
  expect_equal(as.numeric(qs$se.fit[i, 1]), stats::sd(draws), tolerance = 0.1)
})

test_that("an asymmetric family needs the cross-parameter covariance", {
  set.seed(3); n <- 900L
  d <- data.frame(x = runif(n))
  mn <- exp(1 + 0.7 * sin(2 * pi * d$x)); sv <- exp(0.2 + 0.6 * cos(2 * pi * d$x))
  d$y <- RTMBdist::rgamma2(n, mn, sv)
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 10), sd = ~ s(x, k = 10)),
               family = fam("gamma2"), data = d)
  qs <- stats::predict(f, type = "quantile", prob = 0.5, se.fit = TRUE)
  expect_equal(mean(d$y < qs$fit[, 1]), 0.5, tolerance = 0.03)
  ## the median of a gamma is not its mean parameter, so the SEs must differ
  expect_gt(max(abs(qs$se.fit[, 1] - stats::predict(f, se.fit = TRUE)$se.fit$mean)),
            1e-3)
})

test_that("quantiles are refused where they are not available", {
  set.seed(1); n <- 400
  d <- data.frame(x = runif(n))
  d$y <- stats::rpois(n, exp(1.2 + sin(2 * pi * d$x)))
  ## no quantile function
  fz <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8), zeroprob = ~ 1),
                family = fam("zipois"), data = d)
  expect_null(fz$family$qf)
  expect_error(stats::predict(fz, type = "quantile"), "no quantile function")
  ## a lattice family has quantiles but no meaningful derivative for them
  fp <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8)), family = fam("pois"), data = d)
  q <- stats::predict(fp, type = "quantile", prob = c(0.1, 0.9))
  expect_true(all(q == round(q)))                # integer-valued
  expect_true(all(q[, 1] <= q[, 2]))
  expect_error(stats::predict(fp, type = "quantile", se.fit = TRUE),
               "continuous")
  ## families() says which support what
  fs <- families()
  expect_true(fs$quantiles[fs$family == "gamma2"])
  expect_false(fs$quantiles[fs$family == "zipois"])
})

test_that("a single row keeps its matrix shape", {
  set.seed(1); n <- 300
  d <- data.frame(x = runif(n)); d$y <- stats::rnorm(n, sin(2 * pi * d$x), 0.3)
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d)
  expect_identical(dim(stats::predict(f, newdata = d[1, ], type = "quantile",
                                      prob = 0.5)), c(1L, 1L))
  expect_identical(dim(stats::predict(f, newdata = d[1:2, ], type = "quantile",
                                      prob = c(0.25, 0.75))), c(2L, 2L))
})

test_that("the quantile plot draws and returns its curves", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  data(mcycle, package = "MASS", envir = environment())
  f <- gamRTMB(accel ~ list(mean = ~ s(times, k = 20), sd = ~ s(times, k = 12)),
               data = mcycle)
  expect_silent(r <- plot(f, type = "quantile"))
  expect_s3_class(r, "data.frame")
  expect_identical(names(r)[1], "times")
  expect_true("se.mid" %in% names(r))
  expect_identical(nrow(r), 200L)
  ## symmetric probabilities pair up in the line coding
  expect_silent(plot(f, type = "quantile", prob = c(0.05, 0.5, 0.95)))
  expect_silent(plot(f, type = "quantile", se = FALSE, n = 50))
  expect_error(plot(f, type = "quantile", xvar = "nope"), "must name a column")
  ## a lattice family gets step curves and no band
  set.seed(1); n <- 400
  d <- data.frame(x = runif(n)); d$y <- stats::rpois(n, exp(1 + sin(2 * pi * d$x)))
  fp <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8)), family = fam("pois"), data = d)
  rp <- plot(fp, type = "quantile")
  expect_false("se.mid" %in% names(rp))
})

test_that("shrinkage bases resolve overlapping null spaces and select terms", {
  set.seed(2); n <- 500
  d <- data.frame(x = runif(n), z = runif(n), noise = runif(n),
                  g = factor(sample(3, n, TRUE)))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x), 0.4)
  ## bs = "ts" penalizes the null space too, so there is none left to duplicate
  expect_error(gamRTMB(y ~ list(mean = ~ x + s(x, k = 8)), data = d),
               "rank deficient")
  expect_no_error(gamRTMB(y ~ list(mean = ~ x + s(x, bs = "ts", k = 8)), data = d))
  expect_no_error(gamRTMB(y ~ list(mean = ~ t2(x, z, k = 4) +
                                     s(z, bs = "ts", k = 8)), data = d))
  ## and with no unpenalized part a whole term can shrink away
  f <- gamRTMB(y ~ list(mean = ~ s(x, bs = "ts", k = 10) +
                          s(noise, bs = "ts", k = 10)), data = d)
  e <- edf(f)
  expect_lt(e$edf[e$term == "s(noise)"], 1)
  expect_gt(e$edf[e$term == "s(x)"], 4)
})

test_that("the quantile fan uses one hue, fading outwards from a black median", {
  pr <- seq(0.05, 0.95, by = 0.05)
  cl <- gamRTMB:::.fan_cols(pr)
  expect_length(cl, length(pr))
  ## the median is black and everything else is the same hue
  expect_identical(cl[pr == 0.5], "black")
  hues <- unique(substr(cl[pr != 0.5], 1, 7))
  expect_length(hues, 1L)
  ## opacity falls away from the median, and pairs match
  alpha <- strtoi(substr(cl[pr != 0.5], 8, 9), 16L)
  d <- abs(pr[pr != 0.5] - 0.5)
  expect_true(all(diff(alpha[order(d)]) <= 0))
  expect_equal(alpha[which.min(pr[pr != 0.5])],
               alpha[which.max(pr[pr != 0.5])])
})

test_that("a fixed argument is taken from newdata, not from the fitting data", {
  set.seed(7); n <- 400
  d <- data.frame(x = runif(n), trials = sample(10:40, n, TRUE))
  d$y <- stats::rbinom(n, d$trials, stats::plogis(-0.5 + sin(2 * pi * d$x)))
  f <- gamRTMB(y ~ list(prob = ~ s(x, k = 8)),
               family = fam("binom", fixed = list(size = "trials")), data = d)
  ## the same covariate, different numbers of trials: quantiles must scale
  nd <- data.frame(x = c(0.5, 0.5), trials = c(10, 40))
  q <- stats::predict(f, newdata = nd, type = "quantile", prob = 0.5)
  expect_identical(dim(q), c(2L, 1L))
  expect_gt(q[2, 1], q[1, 1])
  expect_equal(as.numeric(q[2, 1] / q[1, 1]), 4, tolerance = 0.35)
  ## and the resolver reads the new rows
  expect_identical(gamRTMB:::.resolve_fixed(f$family, nd, 2L)$size, c(10, 40))
})

test_that("conditional density plots work wherever the density does", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  data(mcycle, package = "MASS", envir = environment())
  f <- gamRTMB(accel ~ list(mean = ~ s(times, k = 20), sd = ~ s(times, k = 12)),
               data = mcycle)
  r <- plot(f, type = "density")
  expect_named(r, c("y", paste0("d", 1:5)))
  expect_identical(length(r$y), 200L)
  ## each density is scaled to a common width
  expect_true(all(vapply(r[-1], max, 0) == 1))
  expect_true(all(vapply(r[-1], function(z) all(z >= 0), TRUE)))
  ## `at` chooses the positions
  expect_named(plot(f, type = "density", at = c(15, 30)), c("y", "d1", "d2"))

  ## the point of it: this works for a family with no quantile function
  set.seed(1); n <- 500
  d <- data.frame(x = runif(n))
  d$y <- ifelse(runif(n) < 0.3, 0, stats::rpois(n, exp(1.5 + sin(2 * pi * d$x))))
  fz <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8), zeroprob = ~ 1),
                family = fam("zipois"), data = d)
  expect_null(fz$family$qf)
  expect_error(plot(fz, type = "quantile"), "no quantile function")
  expect_silent(rz <- plot(fz, type = "density"))
  expect_named(rz, c("y", paste0("d", 1:5)))
  ## a lattice response is evaluated on the integers
  expect_true(all(rz$y == round(rz$y)))
  expect_error(plot(f, type = "density", xvar = "nope"), "must name a column")
})
