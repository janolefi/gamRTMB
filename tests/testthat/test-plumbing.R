sim1 <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x), 0.3)
  d
}

test_that("prior weights match replicating rows", {
  d <- sim1(200)
  w <- rep(c(1, 3), 100)
  fw <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d, weights = w)
  ## weighting an observation must be the same as repeating it
  fr <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)),
                data = d[rep(seq_len(200), w), ])
  expect_equal(edf(fw)$edf, edf(fr)$edf, tolerance = 1e-3)
  expect_true(fw$convergence)
  expect_identical(length(fw$weights), 200L)
})

test_that("weights are validated", {
  d <- sim1(100)
  expect_error(gamRTMB(y ~ list(mean = ~ s(x, k = 5)), data = d,
                       weights = rep(1, 7)), "length")
  expect_error(gamRTMB(y ~ list(mean = ~ s(x, k = 5)), data = d,
                       weights = c(-1, rep(1, 99))), "non-negative")
  expect_error(gamRTMB(y ~ list(mean = ~ s(x, k = 5)), data = d,
                       weights = as.character(rep(1, 100))), "numeric")
  ## a scalar weight is recycled and changes nothing but the scale
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 5)), data = d, weights = 1)
  expect_true(f$convergence)
})

test_that("an offset enters one parameter and is recomputed on newdata", {
  set.seed(2); n <- 400
  d <- data.frame(x = runif(n), E = runif(n, 1, 5))
  d$y <- rpois(n, d$E * exp(sin(2 * pi * d$x)))
  f <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8) + offset(log(E))),
               family = fam("pois"), data = d)
  expect_true(f$convergence)
  ## the offset is in the linear predictor, so removing it recovers the signal
  expect_gt(stats::cor(stats::predict(f)$lambda - log(d$E),
                       sin(2 * pi * d$x)), 0.7)
  ## and it must be rebuilt from newdata, not carried over
  expect_equal(stats::predict(f, newdata = d), stats::predict(f),
               tolerance = 1e-10)
  half <- d; half$E <- d$E / 2
  expect_equal(stats::predict(f, newdata = half)$lambda,
               stats::predict(f)$lambda - log(2), tolerance = 1e-10)
})

test_that("na.action drops incomplete rows and reports them", {
  d <- sim1(200)
  d$y[3] <- NA; d$x[7] <- NA
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d)
  expect_identical(f$dropped, 2L)
  expect_identical(f$design$n, 198L)
  expect_identical(length(f$y), 198L)
  expect_error(gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d,
                       na.action = stats::na.fail))
  ## weights are dropped alongside the rows they belong to
  fw <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d,
                weights = rep(1, 200))
  expect_identical(length(fw$weights), 198L)
})

test_that("model variables must live in the data", {
  d <- sim1(100)
  z <- stats::runif(100)              # only in the calling environment
  expect_error(gamRTMB(y ~ list(mean = ~ s(z, k = 5)), data = d),
               "must be a column of `data`")
})

test_that("logLik is the likelihood, not the fitting criterion", {
  d <- sim1(300)
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d)
  ll <- stats::logLik(f)
  ## the REML criterion is a different number, and must not be used for AIC
  expect_false(isTRUE(all.equal(as.numeric(ll), -f$objective)))
  ## it is the plain log-density sum at the fitted parameters
  th <- stats::predict(f, type = "response")
  expect_equal(as.numeric(ll),
               sum(stats::dnorm(f$y, th$mean, th$sd, log = TRUE)),
               tolerance = 1e-8)
  ## df is the total effective degrees of freedom (the mgcv/GAIC convention)
  expect_equal(attr(ll, "df"), attr(edf(f), "edf.total"), tolerance = 1e-8)
  expect_identical(attr(ll, "nobs"), 300L)
  expect_identical(stats::nobs(f), 300L)
  expect_equal(stats::AIC(f), -2 * as.numeric(ll) + 2 * attr(ll, "df"),
               tolerance = 1e-8)
  ## and it prefers a real signal over an intercept-only mean
  f0 <- gamRTMB(y ~ list(mean = ~ 1), data = d)
  expect_lt(stats::AIC(f), stats::AIC(f0))
})

test_that("weighted logLik uses the weights", {
  d <- sim1(200)
  w <- rep(c(1, 2), 100)
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 8)), data = d, weights = w)
  th <- stats::predict(f, type = "response")
  expect_equal(as.numeric(stats::logLik(f)),
               sum(w * stats::dnorm(f$y, th$mean, th$sd, log = TRUE)),
               tolerance = 1e-8)
})

test_that("summary and vcov report per-parameter coefficients", {
  set.seed(3); n <- 300
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- rnorm(n, sin(2 * pi * d$x1), exp(-1 + d$x2))
  f <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
  s <- summary(f)
  expect_s3_class(s, "summary.gamRTMB")
  expect_identical(colnames(s$coefficients),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  ## one row per parametric coefficient; smooth null-space bases excluded
  expect_identical(rownames(s$coefficients),
                   c("mean:(Intercept)", "sd:(Intercept)"))
  expect_true(all(s$coefficients[, "Std. Error"] > 0))
  expect_identical(nrow(s$smooth), 2L)
  expect_length(s$formulas, 2L)
  expect_output(print(s), "Parametric coefficients")
  expect_output(print(s), "Smooth terms")

  V <- vcov(f)
  expect_identical(dim(V), c(f$design$nbeta, f$design$nbeta))
  expect_identical(rownames(V)[1:2], c("mean:(Intercept)", "mean:s(x1).null1"))
  expect_true(all(diag(V) > 0))
})

test_that("summary degrades gracefully without the joint precision", {
  d <- sim1(150)
  f <- gamRTMB(y ~ list(mean = ~ s(x, k = 6)), data = d,
               joint_precision = FALSE)
  s <- summary(f)
  expect_identical(colnames(s$coefficients), "Estimate")
  expect_output(print(s), "joint_precision")
})
