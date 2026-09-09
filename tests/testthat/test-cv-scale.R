test_that("a relative scale parameter starts on the right scale", {
  ## The Box-Cox families parameterise sigma as a coefficient of variation,
  ## not a standard deviation. Starting it at log(sd(y)) is wrong by orders of
  ## magnitude and the Laplace approximation diverges from there.
  skip_if_not_installed("gamlss.data")
  data(abdom, package = "gamlss.data", envir = environment())
  cv <- stats::sd(abdom$y) / mean(abdom$y)

  s_bct <- fam("bct")$start(abdom$y)
  expect_equal(exp(s_bct[["sigma"]]), cv, tolerance = 1e-8)
  expect_lt(exp(s_bct[["sigma"]]), 1)                 # relative, so order 0.1

  ## and a family whose scale really is a standard deviation is untouched
  expect_equal(exp(fam("gamma2")$start(abdom$y)[["sd"]]), stats::sd(abdom$y),
               tolerance = 1e-8)
  expect_equal(gaussian <- exp(fam("norm")$start(abdom$y)[["sd"]]),
               stats::sd(abdom$y), tolerance = 1e-8)

  ## the four-parameter Box-Cox t then fits
  f <- gamRTMB(y ~ list(mu = ~ s(x), sigma = ~ s(x), nu = ~ 1, tau = ~ 1),
               family = fam("bct"), data = abdom)
  expect_true(f$convergence)
  expect_identical(length(f$family$parnames), 4L)
})
