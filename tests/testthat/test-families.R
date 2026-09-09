test_that("links are derived from native parameter names", {
  expect_identical(fam("skewnorm2")$links,
                   c(mean = "identity", sd = "log", alpha = "identity"))
  expect_identical(fam("zipois")$links,
                   c(lambda = "log", zeroprob = "logit"))
})

test_that("positive- and unit-support locations override the name default", {
  # `mean`/`mu` reads as unconstrained but is not in these families
  expect_identical(fam("gamma2")$links[["mean"]], "log")
  expect_identical(fam("nbinom2")$links[["mu"]], "log")
  expect_identical(fam("beta2")$links[["mu"]], "logit")
})

test_that("redundant reparameterisations are dropped, not modelled", {
  f <- fam("invgamma")   # scale = 1/rate
  expect_identical(f$parnames, c("shape", "rate"))
  expect_identical(f$derived, "scale")
})

test_that("`size` is a parameter or data depending on the density", {
  expect_true("size" %in% fam("nbinom2")$parnames)
  expect_error(fam("betabinom"), "fixed argument")
  f <- fam("betabinom", fixed = list(size = "trials"))
  expect_identical(f$parnames, c("shape1", "shape2"))
  expect_identical(names(f$fixed), "size")
})

test_that("non-families are rejected up front", {
  ## vector- or matrix-valued responses are named explicitly
  for (d in c("mvt", "wishart", "dirichlet", "dirmult", "vmf"))
    expect_error(fam(d), "not a univariate regression family")
  ## copulas take other densities rather than a response, so they never even
  ## look like a density
  for (d in c("copula", "mvcopula", "nosuchthing"))
    expect_error(fam(d), "no density for")
})

test_that("most of RTMBdist resolves with no hand-written spec", {
  dens <- grep("^d", ls(asNamespace("RTMBdist")), value = TRUE)
  dens <- dens[vapply(dens, function(x) is.function(get(x, asNamespace("RTMBdist"))), TRUE)]
  ok <- vapply(dens, function(x)
    !is.null(tryCatch(fam(x), error = function(e) NULL)), TRUE)
  expect_gt(sum(ok), 60)
})

test_that("a shape parameter with a zero score does not start at zero", {
  f <- fam("skewnorm2")
  set.seed(1)
  right <- RTMBdist::rskewnorm2(400, 0, 1,  3)
  left  <- RTMBdist::rskewnorm2(400, 0, 1, -3)
  expect_gt(f$start(right)[["alpha"]], 0)   # sign follows the sample skewness
  expect_lt(f$start(left)[["alpha"]],  0)
  expect_identical(fam("norm")$start(right)[["mean"]], mean(right))
})

test_that("fixed arguments are resolved against the data", {
  f <- fam("betabinom", fixed = list(size = "trials"))
  d <- data.frame(trials = rep(10, 5))
  expect_identical(gamRTMB:::.resolve_fixed(f, d, 5)$size, d$trials)
  f2 <- fam("betabinom", fixed = list(size = 10))
  expect_identical(gamRTMB:::.resolve_fixed(f2, d, 5)$size, 10)
  f3 <- fam("betabinom", fixed = list(size = "nope"))
  expect_error(gamRTMB:::.resolve_fixed(f3, d, 5), "must be numeric")
})

test_that("families() lists what fam() can actually build", {
  f <- families()
  expect_named(f, c("family", "parameters", "needs", "support", "residuals",
                    "quantiles", "source"))
  expect_true(all(c("norm", "gamma2", "skewnorm2", "zipois") %in% f$family))
  expect_gt(nrow(f), 60)
  ## nothing listed as ready to use may fail to build
  ready <- f$family[!nzchar(f$needs)]
  bad <- ready[!vapply(ready, function(x)
    !is.null(tryCatch(fam(x), error = function(e) NULL)), TRUE)]
  expect_identical(bad, character(0))
  ## and the ones that need something say what it is
  expect_true(all(nzchar(f$needs[f$family %in% c("betabinom", "binom")])))
})

test_that("families() filters by pattern", {
  expect_true(all(grepl("beta", families("beta")$family)))
  expect_identical(nrow(families("^norm$")), 1L)
})

test_that("standard RTMB densities fill in what RTMBdist lacks", {
  expect_identical(fam("norm")$links, c(mean = "identity", sd = "log"))
  expect_identical(fam("pois")$parnames, "lambda")
  expect_identical(fam("norm")$source, "RTMB")
  expect_identical(fam("gamma2")$source, "RTMBdist")
  ## dgamma's scale = 1/rate is redundant and must be dropped
  expect_identical(fam("gamma")$derived, "scale")
  ## non-centrality is never a regression parameter
  expect_identical(fam("t")$parnames, "df")
})
