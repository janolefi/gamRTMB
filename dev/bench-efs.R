## method = "aREML" against method = "REML": criterion, effective degrees of
## freedom and time, on models where both run, plus one where only one does.
##
## Usage: Rscript dev/bench-efs.R    (a few minutes)
##
## The point is not to declare a winner. The two optimise the same criterion,
## and aREML drops a third-derivative term from the *gradient* but not from
## the criterion itself, so the question is how far short of the REML optimum
## it stops and what it buys for that.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

quiet <- function(expr)
  withCallingHandlers(tryCatch(expr, error = function(e) structure(
    conditionMessage(e), class = "bench_error")),
    warning = function(w) invokeRestart("muffleWarning"))

run <- function(label, f, family, data) {
  out <- list()
  for (eng in c("REML", "aREML")) {
    tt <- system.time(fit <- quiet(gamRTMB(f, family = family, data = data,
                                           method = eng)))
    out[[eng]] <- if (inherits(fit, "bench_error"))
      list(ok = FALSE, msg = substr(fit, 1, 70))
    else list(ok = TRUE, V = fit$objective, conv = fit$convergence,
              g = fit$max_grad, edf = quiet(edf(fit)),
              rep = .or_else(fit$psd_repairs, 0L),
              rep_mode = isTRUE(fit$psd_repaired_at_mode), t = tt[[3]])
  }
  cat(sprintf("\n--- %s ---\n", label))
  for (eng in names(out)) {
    r <- out[[eng]]
    if (!r$ok) { cat(sprintf("  %-8s FAILED: %s\n", eng, r$msg)); next }
    e <- if (inherits(r$edf, "bench_error")) NA else sum(r$edf$edf)
    ## A repair at an intermediate iteration is routine -- the starting values
    ## often have negative curvature and the iteration walks away from them.
    ## Only a repair at the final point changes what is being reported.
    cat(sprintf("  %-8s -REML %12.4f  conv %-5s  max|g| %8.2g  edf %7.2f  %6.2fs%s\n",
                eng, r$V, r$conv, r$g, e, r$t,
                if (r$rep > 0) sprintf("  [%d repairs%s]", r$rep,
                  if (r$rep_mode) ", incl. the last" else "") else ""))
  }
  if (all(vapply(out, `[[`, NA, "ok"))) {
    both_edf <- !vapply(out, function(r) inherits(r$edf, "bench_error"), NA)
    cat(sprintf("  %-8s dV %+.4f   %s   speed %.2fx\n", "delta",
                out$aREML$V - out$REML$V,
                if (all(both_edf))
                  sprintf("dEDF %+.3f", sum(out$aREML$edf$edf) - sum(out$REML$edf$edf))
                else "dEDF     n/a",
                out$REML$t / out$aREML$t))
  }
  invisible(out)
}


set.seed(1)
n <- 500
d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n),
                g = factor(sample(letters[1:5], n, TRUE)))
d$ynorm <- rnorm(n, sin(2 * pi * d$x1) + d$x2^2, exp(-1 + 0.8 * cos(2 * pi * d$x2)))
d$ygam  <- rgamma(n, 6, 6 / exp(1 + sin(2 * pi * d$x1)))
d$ypois <- rpois(n, exp(1 + sin(2 * pi * d$x1)))
d$ybeta <- rbeta(n, 4, 4 / plogis(sin(2 * pi * d$x1)) - 4)

run("gaussian location-scale, two smooths",
    ynorm ~ list(mean = ~ s(x1, k = 10) + s(x2, k = 10), sd = ~ s(x2, k = 10)),
    fam("norm"), d)

run("gaussian, three smooths and a random effect",
    ynorm ~ list(mean = ~ s(x1, k = 10) + s(x2, k = 10) + s(g, bs = "re"),
                 sd = ~ s(x2, k = 8)),
    fam("norm"), d)

run("gaussian, id-tied smoothing parameter",
    ynorm ~ list(mean = ~ s(x1, k = 10, id = 1) + s(x2, k = 10, id = 1),
                 sd = ~ s(x2, k = 8)),
    fam("norm"), d)

run("gamma, mean and sd smoothed",
    ygam ~ list(mean = ~ s(x1, k = 10), sd = ~ s(x2, k = 8)),
    fam("gamma2"), d)

run("poisson, one smooth plus a spurious one",
    ypois ~ list(lambda = ~ s(x1, k = 10) + s(x3, k = 10)),
    fam("pois"), d)

## A markov random field: the block keeps its own sparse penalty, so this is
## the route where the EFS numerator needs the general trace rather than the
## iid shortcut.
g <- 8
ids <- as.vector(outer(seq_len(g), seq_len(g), function(i, j) paste0(i, "_", j)))
nb <- setNames(vector("list", g * g), ids)
for (i in seq_len(g)) for (j in seq_len(g)) {
  nm <- paste0(i, "_", j); v <- character(0)
  if (i > 1) v <- c(v, paste0(i - 1, "_", j)); if (i < g) v <- c(v, paste0(i + 1, "_", j))
  if (j > 1) v <- c(v, paste0(i, "_", j - 1)); if (j < g) v <- c(v, paste0(i, "_", j + 1))
  nb[[nm]] <- v
}
ij <- do.call(rbind, lapply(strsplit(ids, "_"), as.numeric))
set.seed(1)
dm <- data.frame(reg = factor(rep(ids, each = 6), levels = ids))
k <- as.integer(factor(as.character(dm$reg), levels = ids))
dm$y <- rnorm(nrow(dm), sin(ij[k, 1] / 2) + cos(ij[k, 2] / 2), 0.4)
run("markov random field, 64 regions (sparse penalty)",
    y ~ list(mean = ~ s(reg, bs = "mrf", xt = list(nb = nb)), sd = ~ 1),
    fam("norm"), dm)

## Where the Laplace engine cannot start at all: the Box-Cox power
## exponential's inner Hessian is indefinite at the starting values, and no
## starting value repairs it -- see .inner_indefinite(). EFS does not need a
## positive definite inner Hessian; it repairs the data Hessian instead, and
## says how often it had to.
if (requireNamespace("gamlss.data", quietly = TRUE)) {
  data(film90, package = "gamlss.data")
  set.seed(1)
  df <- film90[sample(nrow(film90), 600), ]
  run("bcpe on film90 (subsample), four smooths",
      lborev1 ~ list(mu = ~ s(lboopen, k = 8), sigma = ~ s(lboopen, k = 8),
                     nu = ~ s(lboopen, k = 8), tau = ~ s(lboopen, k = 8)),
      fam("bcpe"), df)
}
