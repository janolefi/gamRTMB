test_that("missing parameters are filled with ~1", {
  pf <- gamRTMB:::.parse_formula(y ~ list(mu = ~ s(x)), c("mu", "sigma"))
  expect_identical(names(pf$par_formulas), c("mu", "sigma"))
  expect_equal(pf$par_formulas$sigma, ~1, ignore_attr = TRUE)
  expect_identical(pf$response, as.name("y"))
})

test_that("a plain right-hand side models the first parameter", {
  pf <- gamRTMB:::.parse_formula(y ~ s(x) + z, c("mu", "sigma"))
  expect_equal(pf$par_formulas$mu, ~ s(x) + z, ignore_attr = TRUE)
  expect_equal(pf$par_formulas$sigma, ~1, ignore_attr = TRUE)
  expect_identical(pf$response, as.name("y"))
  ## and it is the same object as the long form
  long <- gamRTMB:::.parse_formula(y ~ list(mu = ~ s(x) + z), c("mu", "sigma"))
  expect_equal(pf$par_formulas, long$par_formulas, ignore_attr = TRUE)
})

test_that("formula misuse is rejected clearly", {
  expect_error(gamRTMB:::.parse_formula(y ~ list(~ s(x)), "mu"), "named")
  expect_error(gamRTMB:::.parse_formula(y ~ list(nope = ~ s(x)), "mu"),
               "unknown distributional parameter")
  expect_error(gamRTMB:::.parse_formula(y ~ list(mu = y ~ s(x)), "mu"),
               "one-sided")
  expect_error(gamRTMB:::.parse_formula(~ list(mu = ~ s(x)), "mu"), "two-sided")
})
