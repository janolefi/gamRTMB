## ---------------------------------------------------------------------------
## Validation of the RTMBdist family layer.
## Run with:  Rscript dev/demo-families.R
##
## Each case simulates from the same RTMBdist density that is then fitted, so
## the data-generating process and the likelihood come from identical code and
## any disagreement is the pipeline's fault, not a parameterisation mismatch.
##
## Four families, chosen to stress different axes rather than to cover ground:
##   gamma2     two strictly positive parameters, both smoothed
##   skewnorm2  THREE parameters smoothed at once -- the novelty claim
##   zipois     discrete response, log x logit links
##   nbinom2    discrete, `size` as a modelled overdispersion parameter
## plus betabinom, where `size` is instead known data -- the fixed-argument
## path, which is the part of the family schema most likely to be wrong.
## ---------------------------------------------------------------------------

.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- "dev"
if (requireNamespace("pkgload", quietly = TRUE))
  suppressMessages(pkgload::load_all(dirname(.here), quiet = TRUE)) else
  library(gamRTMB)
suppressPackageStartupMessages(library(RTMBdist))

## report a fit against the true linear predictors ---------------------------
## Smooths are sum-to-zero constrained and the intercept carries the level, so
## the fitted linear predictor is directly comparable to the truth.

report <- function(fit, eta_true, label) {
  cat("\n---", label, "-------------------------------------------\n")
  print(fit$family)
  cat("  converged:", fit$convergence,
      "| max|grad|:", sprintf("%.2g", fit$max_grad),
      "| -logLik:", sprintf("%.2f", fit$objective), "\n")
  eta <- predict(fit)                       # link scale
  e <- edf(fit)
  for (p in fit$family$parnames) {
    et <- eta_true[[p]]
    ei <- e[e$parameter == p, , drop = FALSE]
    cat(sprintf("  %-9s link=%-8s rmse(eta)=%.4f  sd(eta_true)=%.3f  edf: %s\n",
                p, fit$family$links[[p]],
                sqrt(mean((eta[[p]] - et)^2)), stats::sd(et),
                if (nrow(ei)) paste(sprintf("%s=%.2f", ei$term, ei$edf), collapse = " ")
                else "(no smooth)"))
  }
  invisible(fit)
}

f1 <- function(x) sin(2 * pi * x)
f2 <- function(x) 0.8 * cos(2 * pi * x)
f3 <- function(x) 2 * (x - 0.5)

## --- 1. gamma2: two positive parameters ------------------------------------
demo_gamma2 <- function(n = 800, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  eta <- list(mean = 1.2 + 0.8 * f1(d$x1), sd = 0.2 + 0.6 * f2(d$x2))
  d$y <- rgamma2(n, exp(eta$mean), exp(eta$sd))
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 10)),
                     family = rtmbdist_family("gamma2"), data = d)
  report(fit, eta, "gamma2(mean, sd)")
}

## --- 2. skewnorm2: smooths on THREE parameters simultaneously --------------
## Distinct covariates per parameter, so this measures the machinery rather
## than cross-parameter identifiability (probed separately in
## demo-gaussian-ls.R, which shares a covariate on purpose).
demo_skewnorm2 <- function(n = 2000, seed = 2) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  eta <- list(mean = f1(d$x1), sd = -0.2 + 0.5 * f2(d$x2), alpha = 2 * f3(d$x3))
  d$y <- rskewnorm2(n, eta$mean, exp(eta$sd), eta$alpha)
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 10),
                              alpha = ~ s(x3, k = 10)),
                     family = rtmbdist_family("skewnorm2"), data = d)
  report(fit, eta, "skewnorm2(mean, sd, alpha)")
}

## --- 3. zipois: discrete, log x logit --------------------------------------
## NOTE: simulated locally rather than with RTMBdist::rzipois, which draws
##   res[!is_zero] <- rpois(sum(!is_zero), lambda)
## Without subsetting `lambda`, a vector lambda is silently truncated to its
## first sum(!is_zero) elements and those draws land on the !is_zero
## positions, so each observation receives another observation's lambda. The
## covariate association is destroyed whenever zeroprob > 0 and lambda varies
## per observation -- i.e. exactly in distributional regression. Confirmed
## against mgcv::gam(family = ziplss()), which recovers the same non-signal
## from the same data (cor 0.0799 vs 0.0803), so the fitter is fine and the
## simulator is not. Nine RTMBdist r* functions share the pattern: roibeta,
## rzibeta, rzibinom, rzigamma, rziinvgauss, rzilnorm, rzipois, rziweibull,
## rzoibeta.
.rzipois_aligned <- function(n, lambda, zeroprob) {
  ifelse(stats::runif(n) < zeroprob, 0, stats::rpois(n, lambda))
}

demo_zipois <- function(n = 1500, seed = 3) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  eta <- list(lambda = 1.5 + f1(d$x1), zeroprob = -1 + f2(d$x2))
  d$y <- .rzipois_aligned(n, exp(eta$lambda), plogis(eta$zeroprob))
  fit <- gamRTMB(y ~ list(lambda = ~ s(x1, k = 10), zeroprob = ~ s(x2, k = 10)),
                     family = rtmbdist_family("zipois"), data = d)
  report(fit, eta, "zipois(lambda, zeroprob)")
}

## --- 4. nbinom2: `size` as a modelled overdispersion parameter -------------
demo_nbinom2 <- function(n = 2000, seed = 4) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n))
  eta <- list(mu = 2 + 0.8 * f1(d$x1), size = 1 + f2(d$x2))
  d$y <- rnbinom2(n, exp(eta$mu), exp(eta$size))
  fit <- gamRTMB(y ~ list(mu = ~ s(x1, k = 10), size = ~ s(x2, k = 10)),
                     family = rtmbdist_family("nbinom2"), data = d)
  report(fit, eta, "nbinom2(mu, size)")
}

## --- 5. betabinom: `size` as known data (fixed-argument path) --------------
demo_betabinom <- function(n = 1500, seed = 5) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n), x2 = runif(n), trials = sample(10:30, n, TRUE))
  eta <- list(shape1 = 1 + 0.8 * f1(d$x1), shape2 = 1 + 0.8 * f2(d$x2))
  d$y <- rbetabinom(n, d$trials, exp(eta$shape1), exp(eta$shape2))
  fam <- rtmbdist_family("betabinom", fixed = list(size = "trials"))
  fit <- gamRTMB(y ~ list(shape1 = ~ s(x1, k = 10), shape2 = ~ s(x2, k = 10)),
                     family = fam, data = d)
  report(fit, eta, "betabinom(size = trials, shape1, shape2)")
  ## the fixed argument must not be reachable from the formula interface
  bad <- tryCatch(gamRTMB(y ~ list(size = ~ s(x1)), family = fam, data = d),
                  error = conditionMessage)
  cat("  formula on a fixed argument is rejected:",
      is.character(bad) && grepl("unknown distributional parameter", bad), "\n")
  invisible(fit)
}

if (!interactive()) {
  cat("family layer: ")
  nres <- length(grep("^d", ls("package:RTMBdist"), value = TRUE))
  auto <- sum(vapply(grep("^d", ls("package:RTMBdist"), value = TRUE),
                     function(x) !is.null(tryCatch(rtmbdist_family(x),
                                                   error = function(e) NULL)), TRUE))
  cat(auto, "of", nres, "RTMBdist densities resolve with no hand-written spec\n")
  demo_gamma2(); demo_skewnorm2(); demo_zipois(); demo_nbinom2(); demo_betabinom()
  cat("\n")
}
