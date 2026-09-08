test_that("links are derived from native parameter names", {
  expect_identical(rtmbdist_family("skewnorm2")$links,
                   c(mean = "identity", sd = "log", alpha = "identity"))
  expect_identical(rtmbdist_family("zipois")$links,
                   c(lambda = "log", zeroprob = "logit"))
})

test_that("positive- and unit-support locations override the name default", {
  # `mean`/`mu` reads as unconstrained but is not in these families
  expect_identical(rtmbdist_family("gamma2")$links[["mean"]], "log")
  expect_identical(rtmbdist_family("nbinom2")$links[["mu"]], "log")
  expect_identical(rtmbdist_family("beta2")$links[["mu"]], "logit")
})

test_that("redundant reparameterisations are dropped, not modelled", {
  f <- rtmbdist_family("invgamma")   # scale = 1/rate
  expect_identical(f$parnames, c("shape", "rate"))
  expect_identical(f$derived, "scale")
})

test_that("`size` is a parameter or data depending on the density", {
  expect_true("size" %in% rtmbdist_family("nbinom2")$parnames)
  expect_error(rtmbdist_family("betabinom"), "fixed argument 'size'")
  f <- rtmbdist_family("betabinom", fixed = list(size = "trials"))
  expect_identical(f$parnames, c("shape1", "shape2"))
  expect_identical(names(f$fixed), "size")
})

test_that("non-families are rejected up front", {
  for (d in c("mvt", "wishart", "dirichlet", "copula"))
    expect_error(rtmbdist_family(d), "not a univariate regression family")
  expect_error(rtmbdist_family("nosuchthing"), "no density")
})

test_that("most of RTMBdist resolves with no hand-written spec", {
  dens <- grep("^d", ls(asNamespace("RTMBdist")), value = TRUE)
  dens <- dens[vapply(dens, function(x) is.function(get(x, asNamespace("RTMBdist"))), TRUE)]
  ok <- vapply(dens, function(x)
    !is.null(tryCatch(rtmbdist_family(x), error = function(e) NULL)), TRUE)
  expect_gt(sum(ok), 60)
})

test_that("a shape parameter with a zero score does not start at zero", {
  f <- rtmbdist_family("skewnorm2")
  set.seed(1)
  right <- RTMBdist::rskewnorm2(400, 0, 1,  3)
  left  <- RTMBdist::rskewnorm2(400, 0, 1, -3)
  expect_gt(f$start(right)[["alpha"]], 0)   # sign follows the sample skewness
  expect_lt(f$start(left)[["alpha"]],  0)
  expect_identical(gaussian_ls()$start(right)[["mu"]], mean(right))
})

test_that("fixed arguments are resolved against the data", {
  f <- rtmbdist_family("betabinom", fixed = list(size = "trials"))
  d <- data.frame(trials = rep(10, 5))
  expect_identical(gamRTMB:::.resolve_fixed(f, d, 5)$size, d$trials)
  f2 <- rtmbdist_family("betabinom", fixed = list(size = 10))
  expect_identical(gamRTMB:::.resolve_fixed(f2, d, 5)$size, 10)
  f3 <- rtmbdist_family("betabinom", fixed = list(size = "nope"))
  expect_error(gamRTMB:::.resolve_fixed(f3, d, 5), "must be numeric")
})
