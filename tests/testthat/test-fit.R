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
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10) + s(x2, k = 10),
                          sd = ~ s(x2, k = 10)),
                 data = d)
  g <- mgcv::gam(list(y ~ s(x1, k = 10) + s(x2, k = 10), ~ s(x2, k = 10)),
                 family = mgcv::gaulss(), data = d)
  expect_true(fit$convergence)
  ## two different optimisers land on slightly different smoothing
  ## parameters, so agreement is checked at the 1% level, not to the last digit
  expect_equal(edf(fit)$edf, as.numeric(summary(g)$edf), tolerance = 0.01)
})

test_that("id-tied smoothing parameters agree with mgcv, and reduce the count", {
  d <- sim_ls()
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10, id = 1) + s(x2, k = 10, id = 1)),
                 data = d)
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
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 10)),
                 data = d)
  ins <- stats::predict(fit)
  ## PredictMat on the fitting data must reproduce the in-sample predictors
  expect_equal(stats::predict(fit, newdata = d), ins, tolerance = 1e-10)
  ## a subset predicts as the corresponding rows, i.e. no refitting of knots
  sub <- stats::predict(fit, newdata = d[10:20, ])
  expect_equal(sub$mean, ins$mean[10:20], tolerance = 1e-10)
  ## response scale applies the link inverse
  expect_equal(stats::predict(fit, type = "response")$sd, exp(ins$sd),
               tolerance = 1e-10)
  ## per-term contributions sum to the predictor minus the intercept
  tm <- stats::predict(fit, type = "terms")
  expect_equal(as.numeric(tm$mean$fit[, 1]),
               ins$mean - stats::coef(fit)$beta[1], tolerance = 1e-10)
})

test_that("standard errors need the joint precision, and are finite with it", {
  d <- sim_ls(200)
  ## on by default, so the error path needs it switched off explicitly
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8)), data = d,
                 joint_precision = FALSE)
  expect_error(stats::predict(fit, se.fit = TRUE), "joint_precision")
  fit2 <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8)), data = d)
  p <- stats::predict(fit2, newdata = d[1:5, ], se.fit = TRUE)
  expect_true(all(is.finite(p$se.fit$mean)))
  expect_true(all(p$se.fit$mean > 0))
})

test_that("ML runs but cannot report EDF", {
  d <- sim_ls(200)
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 8)), 
                 data = d, method = "ML")
  expect_true(fit$convergence)
  expect_error(edf(fit), "REML")
})

test_that("the efs engine is declared but not implemented", {
  d <- sim_ls(100)
  expect_error(gamRTMB(y ~ list(mean = ~ s(x1, k = 5)), 
                       data = d, engine = "efs"),
               "not implemented")
})

test_that("bad input is caught before fitting", {
  d <- sim_ls(100)
  expect_error(gamRTMB(y ~ list(mean = ~ s(x1)), family = "gaussian", data = d),
               "gamRTMB_family")
  expect_error(gamRTMB(y ~ list(mean = ~ s(x1)), data = 1),
               "data frame")
})

test_that("a non-finite gradient is diagnosed, not passed to the optimiser", {
  ## .flat_start must survive a gradient it cannot use
  expect_null(gamRTMB:::.flat_start(c(NA_real_, NaN), function(p) 0,
                                    c(beta = 0), TRUE,
                                    list(parnames = character(0))))

  ## RTMBdist's dtruncnorm has a finite value but a NaN derivative w.r.t. the
  ## scale when a bound is infinite. Whether that is still true upstream or
  ## not, the fit must either work or say precisely why -- never hand NaN to
  ## nlminb.
  set.seed(1); n <- 300
  d <- data.frame(x = stats::runif(n))
  d$y <- RTMBdist::rtruncnorm(n, 1 + sin(2 * pi * d$x), 1, min = 0, max = Inf)
  res <- tryCatch(gamRTMB(y ~ list(mean = ~ s(x, k = 6)),
                          family = fam("truncnorm", fixed = list(min = 0)),
                          data = d),
                  error = function(e) conditionMessage(e))
  expect_true(inherits(res, "gamRTMB") || grepl("gradient", res))

  ## with finite bounds it fits either way
  fit <- gamRTMB(y ~ list(mean = ~ s(x, k = 6)),
                 family = fam("truncnorm", fixed = list(min = 0, max = 1e6)),
                 data = d)
  expect_true(fit$convergence)
  expect_gt(stats::cor(stats::predict(fit)$mean, 1 + sin(2 * pi * d$x)), 0.9)
})

test_that("the inner method can be switched, and says it is experimental", {
  set.seed(1)
  d <- data.frame(x = stats::runif(300))
  d$y <- stats::rnorm(300, sin(2 * pi * d$x), 0.3)
  f <- y ~ list(mean = ~ s(x, k = 10), sd = ~ 1)
  ## BFGS also leaves sdreport unable to report the random effects' standard
  ## deviations, which is one more symptom of the same warm-start problem.
  withr::local_options(warn = -1)
  expect_warning(b <- gamRTMB(f, data = d, inner_method = "BFGS"),
                 "experimental")
  expect_s3_class(b, "gamRTMB")
  expect_equal(b$obj$env$inner.method, "BFGS")
  expect_equal(b$obj$env$inner.control$maxit, 1e4)
  ## overrides reach optim, and the default route is untouched
  expect_warning(b2 <- gamRTMB(f, data = d, inner_method = "BFGS",
                               inner_control = list(reltol = 1e-14)))
  expect_equal(b2$obj$env$inner.control$reltol, 1e-14)
  n <- gamRTMB(f, data = d)
  expect_false(identical(n$obj$env$inner.method, "BFGS"))
})
