## Plots are checked by rendering to a null device and inspecting the data
## they return, rather than by comparing images.
local_null_device <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  withr_defer <- function(expr) do.call(on.exit, list(substitute(expr), add = TRUE),
                                        envir = env)
  withr_defer(grDevices::dev.off())
}

fit_two <- function(n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x1), exp(-1 + d$x2))
  gamRTMB(y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)), data = d)
}

test_that("term plots return the plotted curves", {
  local_null_device()
  f <- fit_two()
  r <- plot(f)
  expect_length(r, 2L)
  expect_named(r, c("mean: s(x1)", "sd: s(x2)"))
  expect_named(r[[1]], c("x", "fit", "se"))
  expect_identical(nrow(r[[1]]), 200L)
  expect_true(all(is.finite(r[[1]]$se)) && all(r[[1]]$se > 0))
  ## the grid spans the observed covariate range, read from the stored data
  expect_equal(range(r[[1]]$x), range(f$data$x1), tolerance = 1e-12)
  expect_identical(nrow(plot(f, n = 50)[[1]]), 50L)
})

test_that("select picks single smooths, by index or pattern", {
  local_null_device()
  f <- fit_two()
  expect_named(plot(f, select = 1), "mean: s(x1)")
  expect_named(plot(f, select = "sd"), "sd: s(x2)")
  expect_named(plot(f, select = "s\\(x1\\)"), "mean: s(x1)")
  expect_error(plot(f, select = "nope"), "matched none")
  expect_error(plot(f, select = 99), "must index")
})

test_that("caller's graphical arguments override the defaults", {
  local_null_device()
  f <- fit_two()
  ## bty = "n" is a default, not a fixture: it can be overridden, and any
  ## plot() argument passes through without a duplicate-match error
  expect_silent(plot(f, select = 1, bty = "o", main = "mine", col = "red",
                     lwd = 4, xlim = c(0, 1)))
  expect_silent(plot(f, select = 1, rug = FALSE, se = FALSE))
})

test_that("non-1D smooths are reported and skipped, not drawn", {
  local_null_device()
  set.seed(2); n <- 500
  d <- data.frame(x = runif(n), z = runif(n), g = factor(sample(3, n, TRUE)))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x), 0.4)
  ## note: t2(x, z) must not be paired with s(x) or s(z) -- overlapping null
  ## spaces are rank deficient and rejected, so use a third covariate
  d$w <- stats::runif(n)
  f <- gamRTMB(y ~ list(mean = ~ t2(x, z, k = 4) + s(w, k = 8) + s(g, bs = "re")),
               data = d)
  expect_message(r <- plot(f), "not drawable")
  expect_named(r, "mean: s(w)")
  ## a model with nothing drawable says so
  f2 <- gamRTMB(y ~ list(mean = ~ s(g, bs = "re")), data = d)
  expect_error(suppressMessages(plot(f2)), "no one-dimensional")
})

test_that("by= factor smooths are drawn per level", {
  local_null_device()
  set.seed(3); n <- 600L
  d <- data.frame(x = runif(n), g = factor(sample(3, n, TRUE)))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x) * as.numeric(d$g), 0.4)
  f <- gamRTMB(y ~ list(mean = ~ g + s(x, by = g, k = 8)), data = d)
  r <- plot(f, all.terms = FALSE)
  expect_length(r, 3L)
  expect_true(all(grepl("s\\(x\\):g", names(r))))
  ## the three level curves must differ
  expect_gt(max(abs(r[[1]]$fit - r[[3]]$fit)), 0.1)
  ## by default the g main effect gets a panel of its own, ahead of the smooths
  expect_named(plot(f)[1], "mean: g")
})

test_that("parametric terms get panels, centred as predict.lm does", {
  local_null_device()
  set.seed(4); n <- 400L
  d <- data.frame(z = stats::rnorm(n), g = factor(sample(letters[1:3], n, TRUE)))
  d$y <- stats::rnorm(n, 1 + 2 * d$z + as.numeric(d$g), 0.5)
  f <- gamRTMB(y ~ list(mean = ~ z + g, sd = ~ 1), data = d)

  r <- plot(f)
  expect_named(r, c("mean: z", "mean: g"))
  expect_true(all(is.finite(r[["mean: z"]]$se)))

  ## a factor term is one row per level, in level order; a numeric one a grid
  expect_identical(as.character(r[["mean: g"]]$x), levels(d$g))
  expect_identical(nrow(r[["mean: z"]]), 100L)
  expect_equal(range(r[["mean: z"]]$x), range(d$z), tolerance = 1e-12)

  ## the centring convention is stats', so the same terms fitted by lm agree
  pt <- stats::predict(stats::lm(y ~ z + g, data = d), type = "terms")
  expect_equal(r[["mean: g"]]$fit, as.vector(tapply(pt[, "g"], d$g, mean)),
               tolerance = 1e-8)
  fz <- stats::approxfun(r[["mean: z"]]$x, r[["mean: z"]]$fit)
  expect_equal(fz(d$z), unname(pt[, "z"]), tolerance = 1e-6)

  ## and the standard errors are those of the same linear form
  V <- stats::vcov(f)
  expect_equal(r[["mean: z"]]$se,
               abs(r[["mean: z"]]$x - mean(d$z)) * sqrt(V["mean:z", "mean:z"]),
               tolerance = 1e-10)
})

test_that("a model without smooths still plots its terms", {
  local_null_device()
  set.seed(5); n <- 300L
  d <- data.frame(z = stats::rnorm(n), x = runif(n))
  d$y <- stats::rnorm(n, d$z, exp(-1 + d$x))
  f <- gamRTMB(y ~ list(mean = ~ z, sd = ~ x), data = d)
  expect_named(plot(f), c("mean: z", "sd: x"))
  expect_error(plot(f, all.terms = FALSE), "no smooth terms")
  ## an intercept-only model has nothing to draw either way
  f0 <- gamRTMB(y ~ list(mean = ~ 1), data = d)
  expect_error(plot(f0), "every distributional parameter is an intercept")
})

test_that("parametric interactions are reported and skipped", {
  local_null_device()
  set.seed(6); n <- 400L
  d <- data.frame(z = stats::rnorm(n), g = factor(sample(letters[1:2], n, TRUE)))
  d$y <- stats::rnorm(n, d$z * as.numeric(d$g), 0.5)
  f <- gamRTMB(y ~ list(mean = ~ z * g), data = d)
  expect_message(r <- plot(f), "not drawable")
  expect_named(r, c("mean: z", "mean: g"))
})

test_that("layout is restored after plotting", {
  local_null_device()
  f <- fit_two()
  before <- graphics::par("mfrow")
  plot(f)
  expect_identical(graphics::par("mfrow"), before)
  plot(f, ask = TRUE)
  expect_identical(graphics::par("mfrow"), before)
})

test_that("qq and worm plots return their coordinates", {
  local_null_device()
  set.seed(4); n <- 600L
  d <- data.frame(x = runif(n))
  d$y <- stats::rpois(n, exp(1 + sin(2 * pi * d$x)))
  f <- gamRTMB(y ~ list(lambda = ~ s(x, k = 8)), family = fam("pois"), data = d)
  q <- plot(f, type = "qq")
  expect_identical(nrow(q), n)
  expect_true(all(diff(q$theoretical) >= 0))
  w <- plot(f, type = "worm")
  expect_identical(nrow(w), n)
  ## a worm plot is the detrended QQ plot
  expect_equal(w$deviation + w$theoretical, sort(w$deviation + w$theoretical),
               tolerance = 1e-8)
  ## nsim overlays draws, which differ for a discrete response
  q5 <- plot(f, type = "qq", nsim = 5)
  expect_identical(ncol(q5), 6L)
  expect_false(isTRUE(all.equal(q5[[2]], q5[[3]])))
})

test_that("nsim draws coincide for a continuous response", {
  local_null_device()
  f <- fit_two()
  q <- plot(f, type = "qq", nsim = 3)
  expect_equal(q[[2]], q[[3]], tolerance = 1e-12)
})

test_that("a boundary smoothing parameter does not break the covariance", {
  ## A term shrunk onto its null space leaves the criterion flat in its own
  ## log_sigma, so the joint precision is singular; the singularity is
  ## confined to that block and must not take the coefficient covariance
  ## with it.
  local_null_device()
  set.seed(1); n <- 600
  d <- data.frame(x1 = runif(n), x2 = runif(n), z = runif(n),
                  g = factor(sample(3, n, TRUE)))
  d$y <- stats::rnorm(n, sin(2 * pi * d$x1) + d$x2^2, 0.4)
  f <- gamRTMB(y ~ list(mean = ~ t2(x1, z, k = 4) + s(x2, k = 8) + s(g, bs = "re")),
               data = d)
  expect_true(any(exp(-2 * f$log_sigma) > 1e6))          # at least one at the boundary
  expect_error(solve(as.matrix(f$sdr$jointPrecision)), "singular")
  V <- vcov(f)                                            # must still work
  expect_true(all(is.finite(diag(V))) && all(diag(V) > 0))
  expect_true("Std. Error" %in% colnames(summary(f)$coefficients))
  r <- suppressMessages(plot(f))
  expect_true(all(is.finite(r[[1]]$se)))
})

test_that("the Schur route equals the full inverse when it is well conditioned", {
  f <- fit_two()
  Q <- as.matrix(f$sdr$jointPrecision)
  ic <- which(colnames(Q) %in% c("beta", "b"))
  expect_equal(gamRTMB:::.joint_cov(f)$V, solve(Q)[ic, ic],
               tolerance = 1e-10, ignore_attr = TRUE)
})
