## The reparameterisation map is the one place a silent error would corrupt
## every coefficient, so check it against the identity it must satisfy:
##   cbind(rand..., Xf) %*% c(br, bf)  ==  X %*% Tmap %*% c(br, bf)
reconstructs <- function(sm) {
  re <- mgcv::smooth2random(sm, "", type = 2)
  Tm <- gamRTMB:::.reconstruct_map(re)
  nb <- sum(vapply(re$rand, ncol, 1L)); nf <- ncol(re$Xf)
  set.seed(7); cf <- stats::rnorm(nb + nf)
  Xall <- do.call(cbind, c(lapply(re$rand, as.matrix), list(re$Xf)))
  max(abs(as.vector(Xall %*% cf) - as.vector(sm$X %*% (Tm %*% cf))))
}

test_that("smooth2random reconstruction is exact for every supported basis", {
  set.seed(1); n <- 200
  d <- data.frame(x = runif(n), z = runif(n), g = factor(sample(3, n, TRUE)))
  expect_lt(reconstructs(mgcv::smoothCon(mgcv::s(x, k = 10), d, absorb.cons = TRUE)[[1]]), 1e-12)
  expect_lt(reconstructs(mgcv::smoothCon(mgcv::s(x, bs = "cr", k = 8), d, absorb.cons = TRUE)[[1]]), 1e-12)
  expect_lt(reconstructs(mgcv::smoothCon(mgcv::t2(x, z, k = 4), d, absorb.cons = TRUE)[[1]]), 1e-12)
  expect_lt(reconstructs(mgcv::smoothCon(mgcv::s(g, bs = "re"), d, absorb.cons = TRUE)[[1]]), 1e-12)
  ## bs = "fs" is the case where rind is a real permutation
  expect_lt(reconstructs(mgcv::smoothCon(mgcv::s(x, g, bs = "fs", k = 5), d, absorb.cons = TRUE)[[1]]), 1e-12)
  ## by = returns one smooth per level
  for (sm in mgcv::smoothCon(mgcv::s(x, by = g, k = 6), d, absorb.cons = TRUE))
    expect_lt(reconstructs(sm), 1e-12)
})

test_that("id links the bases so one shared variance is meaningful", {
  set.seed(1); n <- 300
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  D <- gamRTMB:::.build_design(list(mu = ~ s(x1, k = 8, id = 1) + s(x2, k = 8, id = 1)),
                               d, "mu")
  ss <- vapply(D$parts$mu$smooths, function(z) z$sm$S.scale, numeric(1))
  expect_equal(ss[[1]], ss[[2]])                 # pooled penalty scaling
  expect_identical(D$nsigma_free, 1L)
  expect_identical(D$nsigma, 2L)

  Du <- gamRTMB:::.build_design(list(mu = ~ s(x1, k = 8) + s(x2, k = 8)), d, "mu")
  ssu <- vapply(Du$parts$mu$smooths, function(z) z$sm$S.scale, numeric(1))
  expect_false(isTRUE(all.equal(ssu[[1]], ssu[[2]])))
  expect_identical(Du$nsigma_free, 2L)
})

test_that("ids are resolved across distributional parameters", {
  set.seed(1); n <- 200
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  D <- gamRTMB:::.build_design(list(mu = ~ s(x1, k = 8, id = "sh"),
                                    sigma = ~ s(x2, k = 8, id = "sh")),
                               d, c("mu", "sigma"))
  expect_identical(D$nsigma, 2L)
  expect_identical(D$nsigma_free, 1L)
})

test_that("unsupported smooths are refused", {
  set.seed(1); d <- data.frame(x = runif(100), z = runif(100))
  expect_error(gamRTMB:::.build_design(list(mu = ~ mgcv::te(x, z, k = 3)), d, "mu"))
  expect_error(gamRTMB:::.build_design(list(mu = ~ s(x, k = 5, fx = TRUE)), d, "mu"),
               "fx = TRUE")
})
