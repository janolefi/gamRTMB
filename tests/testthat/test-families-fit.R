## Simulate from the same RTMBdist density that is then fitted, so the
## data-generating process and the likelihood come from identical code.
## NOTE: RTMBdist's rzipois (and eight other zero/one-inflated r* functions)
## does not subset its shape parameters when drawing the non-inflated part, so
## a per-observation lambda is misaligned with its observation. Simulated
## locally here.

skip_slow <- function() testthat::skip_on_cran()

test_that("gamma2: two positive parameters are recovered", {
  skip_slow()
  set.seed(1); n <- 600
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  eta <- list(mean = 1.2 + 0.8 * sin(2 * pi * d$x1),
              sd = 0.2 + 0.6 * cos(2 * pi * d$x2))
  d$y <- RTMBdist::rgamma2(n, exp(eta$mean), exp(eta$sd))
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 10)),
                 family = fam("gamma2"), data = d)
  expect_true(fit$convergence)
  p <- stats::predict(fit)
  expect_lt(sqrt(mean((p$mean - eta$mean)^2)), 0.2 * stats::sd(eta$mean))
  expect_lt(sqrt(mean((p$sd - eta$sd)^2)), 0.5 * stats::sd(eta$sd))
})

test_that("skewnorm2: smooths on three parameters at once", {
  skip_slow()
  set.seed(2); n <- 1200
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  eta <- list(mean = sin(2 * pi * d$x1),
              sd = -0.2 + 0.5 * cos(2 * pi * d$x2),
              alpha = 4 * (d$x3 - 0.5))
  d$y <- RTMBdist::rskewnorm2(n, eta$mean, exp(eta$sd), eta$alpha)
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 10),
                          alpha = ~ s(x3, k = 10)),
                 family = fam("skewnorm2"), data = d)
  expect_true(fit$convergence)
  expect_identical(nrow(edf(fit)), 3L)
  p <- stats::predict(fit)
  expect_lt(sqrt(mean((p$mean - eta$mean)^2)), 0.2 * stats::sd(eta$mean))
  expect_lt(sqrt(mean((p$alpha - eta$alpha)^2)), stats::sd(eta$alpha))
})

test_that("zipois: discrete, with a logit-linked parameter", {
  skip_slow()
  set.seed(3); n <- 1000
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  eta <- list(lambda = 1.5 + sin(2 * pi * d$x1),
              zeroprob = -1 + 0.8 * cos(2 * pi * d$x2))
  d$y <- ifelse(runif(n) < stats::plogis(eta$zeroprob), 0,
                stats::rpois(n, exp(eta$lambda)))
  fit <- gamRTMB(y ~ list(lambda = ~ s(x1, k = 10), zeroprob = ~ s(x2, k = 10)),
                 family = fam("zipois"), data = d)
  expect_true(fit$convergence)
  p <- stats::predict(fit)
  expect_lt(sqrt(mean((p$lambda - eta$lambda)^2)), 0.3 * stats::sd(eta$lambda))
})

test_that("betabinom: a fixed argument taken from the data", {
  skip_slow()
  set.seed(5); n <- 800
  d <- data.frame(x1 = runif(n), x2 = runif(n), trials = sample(10:30, n, TRUE))
  eta <- list(shape1 = 1 + 0.8 * sin(2 * pi * d$x1),
              shape2 = 1 + 0.8 * cos(2 * pi * d$x2))
  d$y <- RTMBdist::rbetabinom(n, d$trials, exp(eta$shape1), exp(eta$shape2))
  fam <- fam("betabinom", fixed = list(size = "trials"))
  fit <- gamRTMB(y ~ list(shape1 = ~ s(x1, k = 10), shape2 = ~ s(x2, k = 10)),
                 family = fam, data = d)
  expect_true(fit$convergence)
  ## the fixed argument must be unreachable from the formula
  expect_error(gamRTMB(y ~ list(size = ~ s(x1)), family = fam, data = d),
               "unknown distributional parameter")
})

test_that("a spurious smooth on a flat parameter collapses to its null space", {
  skip_slow()
  set.seed(11); n <- 600
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x1) + d$x2, exp(-1))   # sigma constant
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10) + s(x2, k = 10),
                          sd = ~ s(x2, k = 10)),
                 data = d)
  e <- edf(fit)
  expect_lt(e$edf[e$parameter == "sd"], 2)
})
