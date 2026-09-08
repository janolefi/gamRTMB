## mgcv::gaulss fits the same model by a different route, so its EDF are the
## sharpest available check on the whole pipeline at once: basis construction,
## the smooth2random reparameterisation, the Laplace fit and the EDF formula
## all have to be right simultaneously for these to agree.

sim_ls <- function(n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x1) + d$x2^2,
                      exp(-1 + 0.8 * cos(2 * pi * d$x2)))
  d
}

test_that("location-scale EDF agree with mgcv::gaulss", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mu = ~ s(x1, k = 10) + s(x2, k = 10),
                          sigma = ~ s(x2, k = 10)),
                 family = gaussian_ls(), data = d)
  g <- mgcv::gam(list(y ~ s(x1, k = 10) + s(x2, k = 10), ~ s(x2, k = 10)),
                 family = mgcv::gaulss(), data = d)
  expect_true(fit$convergence)
  ## two different optimisers land on slightly different smoothing
  ## parameters, so agreement is checked at the 1% level, not to the last digit
  expect_equal(edf(fit)$edf, as.numeric(summary(g)$edf), tolerance = 0.01)
})

test_that("id-tied smoothing parameters agree with mgcv, and reduce the count", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mu = ~ s(x1, k = 10, id = 1) + s(x2, k = 10, id = 1)),
                 family = gaussian_ls(), data = d)
  g <- mgcv::gam(list(y ~ s(x1, k = 10, id = 1) + s(x2, k = 10, id = 1), ~ 1),
                 family = mgcv::gaulss(), data = d)
  expect_identical(fit$design$nsigma_free, 1L)
  expect_identical(length(g$sp), 1L)
  e <- edf(fit)
  expect_equal(e$edf, as.numeric(summary(g)$edf[summary(g)$edf > 0.01]),
               tolerance = 0.01)
  expect_identical(unique(e$sp), e$sp[1])        # one smoothing parameter
})

test_that("predict reuses the fitted bases exactly", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mu = ~ s(x1, k = 10), sigma = ~ s(x2, k = 10)),
                 family = gaussian_ls(), data = d)
  ins <- stats::predict(fit)
  ## PredictMat on the fitting data must reproduce the in-sample predictors
  expect_equal(stats::predict(fit, newdata = d), ins, tolerance = 1e-10)
  ## a subset predicts as the corresponding rows, i.e. no refitting of knots
  sub <- stats::predict(fit, newdata = d[10:20, ])
  expect_equal(sub$mu, ins$mu[10:20], tolerance = 1e-10)
  ## response scale applies the link inverse
  expect_equal(stats::predict(fit, type = "response")$sigma, exp(ins$sigma),
               tolerance = 1e-10)
  ## per-term contributions sum to the predictor minus the intercept
  tm <- stats::predict(fit, type = "terms")
  expect_equal(as.numeric(tm$mu$fit[, 1]),
               ins$mu - stats::coef(fit)$beta[1], tolerance = 1e-10)
})

test_that("standard errors need the joint precision, and are finite with it", {
  d <- sim_ls(200)
  fit <- gamRTMB(y ~ list(mu = ~ s(x1, k = 8)), family = gaussian_ls(), data = d)
  expect_error(stats::predict(fit, se.fit = TRUE), "joint_precision")
  fit2 <- gamRTMB(y ~ list(mu = ~ s(x1, k = 8)), family = gaussian_ls(),
                  data = d, joint_precision = TRUE)
  p <- stats::predict(fit2, newdata = d[1:5, ], se.fit = TRUE)
  expect_true(all(is.finite(p$se.fit$mu)))
  expect_true(all(p$se.fit$mu > 0))
})

test_that("ML runs but cannot report EDF", {
  d <- sim_ls(200)
  fit <- gamRTMB(y ~ list(mu = ~ s(x1, k = 8)), family = gaussian_ls(),
                 data = d, method = "ML")
  expect_true(fit$convergence)
  expect_error(edf(fit), "REML")
})

test_that("the efs engine is declared but not implemented", {
  d <- sim_ls(100)
  expect_error(gamRTMB(y ~ list(mu = ~ s(x1, k = 5)), family = gaussian_ls(),
                       data = d, engine = "efs"),
               "not implemented")
})

test_that("bad input is caught before fitting", {
  d <- sim_ls(100)
  expect_error(gamRTMB(y ~ list(mu = ~ s(x1)), family = "gaussian", data = d),
               "gamRTMB_family")
  expect_error(gamRTMB(y ~ list(mu = ~ s(x1)), family = gaussian_ls(), data = 1),
               "data frame")
})
