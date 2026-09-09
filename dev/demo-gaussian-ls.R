## ---------------------------------------------------------------------------
## Demo: Gaussian location-scale on simulated data.
## Run with:  Rscript dev/demo-gaussian-ls.R
## ---------------------------------------------------------------------------

.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- "dev"
if (requireNamespace("pkgload", quietly = TRUE))
  suppressMessages(pkgload::load_all(dirname(.here), quiet = TRUE)) else
  library(gamRTMB)

demo_gaussian_ls <- function(n = 600, seed = 42, plot_file = NULL) {
  set.seed(seed)
  f1 <- function(x) sin(2 * pi * x)
  f2 <- function(x) 0.6 * (exp(2 * x) - 3) / 3
  f3 <- function(x) -1 + 0.9 * cos(2 * pi * x)          # log sd
  dat <- data.frame(x1 = runif(n), x2 = runif(n))
  dat$mean_true <- f1(dat$x1) + f2(dat$x2)
  dat$ls_true <- f3(dat$x1)
  dat$y <- rnorm(n, dat$mean_true, exp(dat$ls_true))

  cat("\n== 1. fit ================================================\n")
  t0 <- proc.time()[3]
  fit <- gamRTMB(y ~ list(mean = ~ s(x1, k = 12) + s(x2, k = 12),
                              sd = ~ s(x1, k = 12)),
                     data = dat,
                     method = "REML", joint_precision = TRUE)
  cat("  fitted in", sprintf("%.2f", proc.time()[3] - t0), "s\n")
  print(fit)

  cat("\n== 2. sdreport ===========================================\n")
  print(summary(fit$sdr, "fixed"))
  cat("  implied smoothing parameters (1/sigma^2):\n")
  cs <- coef(fit)
  print(setNames(round(exp(-2 * cs$log_sigma), 4),
                 vapply(fit$design$blocks, `[[`, "", "label")))

  cat("\n== 3. EDF per smooth =====================================\n")
  e <- edf(fit); print(e)
  cat("  total EDF (incl. fixed effects):", sprintf("%.2f", attr(e, "edf.total")), "\n")

  cat("\n== 4. predict on new data ================================\n")
  grid <- data.frame(x1 = seq(0, 1, length.out = 101), x2 = 0.5)
  pr <- predict(fit, newdata = grid, se.fit = TRUE)
  cat("  link-scale prediction + se on a fresh grid:\n")
  print(utils::head(data.frame(x1 = grid$x1,
                               mean = pr$fit$mean, mean.se = pr$se$mean,
                               log.sigma = pr$fit$sd,
                               sigma = exp(pr$fit$sd)), 3))
  ## in-sample accuracy against the truth
  ins <- predict(fit)
  cat(sprintf("  RMSE(mean_hat, mean_true)        = %.4f  (sd(y) = %.3f)\n",
              sqrt(mean((ins$mean - dat$mean_true)^2)), sd(dat$y)))
  cat(sprintf("  RMSE(log sd_hat, truth)   = %.4f\n",
              sqrt(mean((ins$sd - dat$ls_true)^2))))
  ## PredictMat really is reused: predicting at the original data must
  ## reproduce the in-sample linear predictors exactly
  chk <- predict(fit, newdata = dat)
  cat(sprintf("  max |predict(newdata = data) - predict()| = %.2e\n",
              max(abs(unlist(chk) - unlist(ins)))))

  cat("\n== 5. cross-check against mgcv::gam(gaulss) ==============\n")
  g <- mgcv::gam(list(y ~ s(x1, k = 12) + s(x2, k = 12), ~ s(x1, k = 12)),
                 family = mgcv::gaulss(), data = dat)
  gm <- stats::predict(g, type = "response")
  cat(sprintf("  cor(mu_gamRTMB, mu_mgcv)             = %.5f\n", cor(ins$mean, gm[, 1])))
  cat(sprintf("  cor(log sd_gamRTMB, log sd_mgcv)= %.5f\n",
              cor(ins$sd, log(1 / gm[, 2]))))
  cat("  mgcv EDF: ", paste(sprintf("%.2f", summary(g)$edf), collapse = ", "), "\n")

  cat("\n== 6. cross-parameter identifiability probe (note 11) ====\n")
  cat("  s(x2) added to sigma, where the truth is flat, while mu also uses x2.\n")
  fit2 <- gamRTMB(y ~ list(mean = ~ s(x1, k = 12) + s(x2, k = 12),
                               sd = ~ s(x1, k = 12) + s(x2, k = 12)),
                      data = dat, method = "REML")
  print(fit2)
  print(edf(fit2))
  cat("  -> the spurious sigma ~ s(x2) term should collapse toward ~1 EDF;\n",
      "     an EDF well above that on a flat truth is the failure mode to watch.\n")

  cat("\n== 7. shared smoothing parameters via s(..., id = ) ============\n")
  ## (a) matches mgcv exactly, including the linked-basis rescaling
  g2 <- mgcv::gam(list(y ~ s(x1, k = 12, id = 1) + s(x2, k = 12, id = 1), ~ 1),
                  family = mgcv::gaulss(), data = dat)
  fid <- gamRTMB(y ~ list(mean = ~ s(x1, k = 12, id = 1) + s(x2, k = 12, id = 1)),
                     data = dat)
  cat(sprintf("  mean ~ s(x1, id=1) + s(x2, id=1):  gamRTMB edf %s | mgcv gaulss edf %s\n",
              paste(sprintf("%.3f", edf(fid)$edf), collapse = ", "),
              paste(sprintf("%.3f", summary(g2)$edf[summary(g2)$edf > 0.01]), collapse = ", ")))
  cat(sprintf("  free smoothing parameters: %d of %d (gamRTMB) vs %d (mgcv)\n",
              fid$design$nsigma_free, fid$design$nsigma, length(g2$sp)))
  cat(sprintf("  S.scale equal across the id group: %s\n",
              isTRUE(all.equal(fid$design$parts$mean$smooths[[1]]$sm$S.scale,
                               fid$design$parts$mean$smooths[[2]]$sm$S.scale))))
  ## (b) an id tying smooths on DIFFERENT distributional parameters, which a
  ## gam formula cannot express at all
  fx <- gamRTMB(y ~ list(mean = ~ s(x1, k = 12, id = "sh") + s(x2, k = 12),
                             sd = ~ s(x1, k = 12, id = "sh")),
                    data = dat)
  cat("  one smoothness shared between mu and sd:\n")
  print(edf(fx))

  if (!is.null(plot_file)) {
    grDevices::png(plot_file, width = 1100, height = 450, res = 110)
    graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    g1 <- data.frame(x1 = seq(0, 1, length.out = 200), x2 = 0.5)
    p1 <- predict(fit, newdata = g1, se.fit = TRUE)
    tr <- predict(fit, newdata = g1, type = "terms", se.fit = TRUE)
    plot(g1$x1, tr$mean$fit[, 1], type = "l", lwd = 2, ylim = c(-2.2, 2.2),
         xlab = "x1", ylab = "s(x1)", main = "mean: s(x1)")
    graphics::polygon(c(g1$x1, rev(g1$x1)),
                      c(tr$mean$fit[, 1] + 2 * tr$mean$se[, 1],
                        rev(tr$mean$fit[, 1] - 2 * tr$mean$se[, 1])),
                      col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
    graphics::lines(g1$x1, f1(g1$x1) - mean(f1(dat$x1)), col = 2, lty = 2, lwd = 2)
    graphics::legend("topright", c("fitted", "truth"), col = c(1, 2), lty = c(1, 2), bty = "n")

    g2 <- data.frame(x1 = 0.5, x2 = seq(0, 1, length.out = 200))
    tr2 <- predict(fit, newdata = g2, type = "terms", se.fit = TRUE)
    plot(g2$x2, tr2$mean$fit[, 2], type = "l", lwd = 2, ylim = c(-2, 2),
         xlab = "x2", ylab = "s(x2)", main = "mean: s(x2)")
    graphics::polygon(c(g2$x2, rev(g2$x2)),
                      c(tr2$mean$fit[, 2] + 2 * tr2$mean$se[, 2],
                        rev(tr2$mean$fit[, 2] - 2 * tr2$mean$se[, 2])),
                      col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
    graphics::lines(g2$x2, f2(g2$x2) - mean(f2(dat$x2)), col = 2, lty = 2, lwd = 2)

    plot(g1$x1, p1$fit$sd, type = "l", lwd = 2, ylim = c(-2.5, 0.5),
         xlab = "x1", ylab = "log sd", main = "sd: log-scale linear predictor")
    graphics::lines(g1$x1, f3(g1$x1), col = 2, lty = 2, lwd = 2)
    grDevices::dev.off()
    cat("\n  plot written to", plot_file, "\n")
  }

  invisible(list(fit = fit, fit2 = fit2, data = dat))
}

if (!interactive())
  .res <- demo_gaussian_ls(plot_file = file.path(.here, "prototype-fit.png"))
