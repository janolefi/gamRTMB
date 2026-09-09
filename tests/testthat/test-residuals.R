## A correct model must give uniform PIT values; the check is run over several
## randomisation draws because discrete residuals are not deterministic.
ks_range <- function(fit, m = 5) {
  p <- vapply(seq_len(m), function(i) {
    set.seed(500 + i)
    suppressWarnings(stats::ks.test(stats::residuals(fit, type = "uniform"),
                                    "punif")$p.value)
  }, 0)
  range(p)
}

test_that("continuous responses give standard normal residuals", {
  set.seed(1); n <- 1200
  d <- data.frame(x = runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x), exp(-1 + d$x))
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 8), sd = ~ s(x, k = 8)), data = d)
  r <- stats::residuals(f)
  expect_length(r, n)
  expect_lt(abs(mean(r)), 0.1)
  expect_lt(abs(stats::sd(r) - 1), 0.1)
  expect_gt(min(ks_range(f)), 0.01)
  ## no randomisation for a continuous response, so draws agree exactly
  expect_identical(diff(ks_range(f)), 0)
  ## the uniform scale is the PIT itself
  u <- stats::residuals(f, type = "uniform")
  expect_true(all(u > 0 & u < 1))
  expect_equal(sort(stats::qnorm(u)), sort(as.numeric(r)), tolerance = 1e-10)
})

test_that("lattice responses use F(y-1) and come out uniform", {
  set.seed(3); n <- 1500
  d <- data.frame(x = runif(n))
  d$y <- stats::rpois(n, exp(1 + sin(2 * pi * d$x)))
  f <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8)), family = fam("pois"), data = d)
  expect_identical(f$family$support, "lattice")
  expect_gt(min(ks_range(f)), 0.01)
  set.seed(1); r <- stats::residuals(f)
  expect_lt(abs(stats::sd(r) - 1), 0.12)
  ## the CDF must never be asked for a value below the support
  expect_true(any(d$y == 0))
  expect_true(all(is.finite(r)))
})

test_that("a lattice family with structural zeros is handled", {
  set.seed(4); n <- 1500
  d <- data.frame(x = runif(n), x2 = runif(n))
  d$y <- ifelse(runif(n) < stats::plogis(-1 + 0.8 * cos(2 * pi * d$x2)), 0,
                stats::rpois(n, exp(1.5 + sin(2 * pi * d$x))))
  f <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8), zeroprob = ~ s(x2, k = 8)),
               family = fam("zipois"), data = d)
  expect_gt(min(ks_range(f)), 0.01)
})

test_that("mixed responses take the atom mass from the density", {
  set.seed(5); n <- 1500
  d <- data.frame(x = runif(n), x2 = runif(n))
  d$y <- ifelse(runif(n) < stats::plogis(-1 + 0.8 * cos(2 * pi * d$x2)), 0,
                stats::rlnorm(n, sin(2 * pi * d$x), 0.6))
  f <- gamRTMB(y ~ list(meanlog = ~ s(x, k = 8), zeroprob = ~ s(x2, k = 8)),
               family = fam("zilnorm"), data = d)
  expect_identical(f$family$support, "mixed")
  expect_identical(f$family$atoms, 0)
  expect_gt(min(ks_range(f)), 0.01)
})

test_that("residuals detect a wrong distribution", {
  set.seed(6); n <- 1500
  d <- data.frame(x = runif(n))
  d$y <- ifelse(runif(n) < 0.35, 0, stats::rpois(n, exp(1.5 + sin(2 * pi * d$x))))
  bad  <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8)), family = fam("pois"), data = d)
  good <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8), zeroprob = ~ 1),
                  family = fam("zipois"), data = d)
  expect_lt(max(ks_range(bad)), 0.001)     # rejected
  expect_gt(min(ks_range(good)), 0.01)     # not rejected
  set.seed(1)
  expect_gt(stats::sd(stats::residuals(bad)), 1.4)
  expect_lt(stats::AIC(good), stats::AIC(bad))
})

test_that("randomise = FALSE returns the bounding interval", {
  set.seed(3); n <- 500
  d <- data.frame(x = runif(n))
  d$y <- stats::rpois(n, exp(1 + sin(2 * pi * d$x)))
  f <- gamRTMB(y ~ list(lambda = ~ s(x, k = 6)), family = fam("pois"), data = d)
  b <- stats::residuals(f, randomise = FALSE)
  expect_s3_class(b, "data.frame")
  expect_named(b, c("lower", "upper", "mid"))
  expect_true(all(b$lower <= b$mid & b$mid <= b$upper))
  ## a continuous response has no interval to return
  d2 <- d; d2$y <- stats::rnorm(n)
  f2 <- gamRTMB(y ~ list(mean = ~ s(x, k = 6)), data = d2)
  expect_type(stats::residuals(f2, randomise = FALSE), "double")
})

test_that("families without a CDF refuse rather than guess", {
  set.seed(7); n <- 400
  d <- data.frame(x = runif(n), n = sample(10:30, n, TRUE))
  d$y <- RTMBdist::rbetabinom(n, d$n, 2, 3)
  f <- gamRTMB(y ~ list(shape1 = ~ s(x, k = 6)),
               family = fam("betabinom", fixed = list(size = "n")), data = d)
  expect_null(f$family$cdf)
  expect_error(stats::residuals(f), "no CDF")
  ## and families() says so up front
  fs <- families()
  expect_false(fs$residuals[fs$family == "betabinom"])
  expect_true(fs$residuals[fs$family == "pois"])
})

test_that("support is derived and can be overridden", {
  expect_identical(fam("norm")$support, "continuous")
  expect_identical(fam("pois")$support, "lattice")
  expect_identical(fam("zipois")$support, "lattice")   # lattice beats the atom
  expect_identical(fam("zilnorm")$support, "mixed")
  expect_identical(fam("zoibeta2")$atoms, c(0, 1))
  expect_identical(fam("norm", support = "lattice")$support, "lattice")
  expect_identical(fam("zilnorm", support = "continuous")$atoms, numeric(0))
})
