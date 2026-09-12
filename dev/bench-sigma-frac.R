## How much does the variance-component starting value matter?
##
## Produces the table in dev/NOTES-sigma-frac.md. Run from the package root:
##
##   Rscript dev/bench-sigma-frac.R
##
## Takes a few minutes. Everything is simulated, so it needs no data package.

suppressMessages({library(gamRTMB); library(RTMBdist)})

n <- 300
mkdat <- function(seed) { set.seed(seed); data.frame(x1 = runif(n), x2 = runif(n)) }

## Two-parameter families first, then the four-parameter ones the choice is
## actually supposed to matter for. Each `sim` draws from the model its
## formula then fits, so a failure is the fitter's and not a misspecification.
specs <- list(
  norm = list(
    f = fam("norm"), fo = y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)),
    sim = function(d) rnorm(n, sin(2 * pi * d$x1), exp(-1 + d$x2))),
  gamma2 = list(
    f = fam("gamma2"), fo = y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8)),
    sim = function(d) rgamma2(n, exp(1 + sin(2 * pi * d$x1)), exp(d$x2))),
  beta2 = list(
    f = fam("beta2"), fo = y ~ list(mu = ~ s(x1, k = 8), phi = ~ s(x2, k = 8)),
    sim = function(d) { m <- plogis(sin(2 * pi * d$x1)); p <- exp(2 + d$x2)
                        rbeta(n, m * p, (1 - m) * p) }),
  skewnorm2 = list(
    f = fam("skewnorm2"),
    fo = y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8), alpha = ~ s(x1, k = 8)),
    sim = function(d) rskewnorm2(n, sin(2 * pi * d$x1), exp(-1 + d$x2), 2)),
  bct = list(
    f = fam("bct"),
    fo = y ~ list(mu = ~ s(x1, k = 8), sigma = ~ s(x2, k = 8),
                  nu = ~ s(x1, k = 8), tau = ~ s(x2, k = 8)),
    sim = function(d) rbct(n, exp(1 + 0.5 * sin(2 * pi * d$x1)),
                           exp(-1.5 + 0.3 * d$x2), 1, 6)),
  bcpe = list(
    f = fam("bcpe"),
    fo = y ~ list(mu = ~ s(x1, k = 8), sigma = ~ s(x2, k = 8),
                  nu = ~ s(x1, k = 8), tau = ~ s(x2, k = 8)),
    sim = function(d) rbcpe(n, exp(1 + 0.5 * sin(2 * pi * d$x1)),
                            exp(-1.5 + 0.3 * d$x2), 1, 2)),
  jsu2 = list(
    f = fam("jsu2"),
    fo = y ~ list(mu = ~ s(x1, k = 8), sigma = ~ s(x2, k = 8),
                  nu = ~ s(x1, k = 8), tau = ~ s(x2, k = 8)),
    sim = function(d) rjsu2(n, sin(2 * pi * d$x1), exp(-0.5 + 0.3 * d$x2), -0.5, 2)),
  skewt2 = list(
    f = fam("skewt2"),
    fo = y ~ list(mean = ~ s(x1, k = 8), sd = ~ s(x2, k = 8),
                  skew = ~ s(x1, k = 8), df = ~ s(x2, k = 8)),
    sim = function(d) rskewt2(n, sin(2 * pi * d$x1), exp(-0.5 + 0.3 * d$x2), 1, 8))
)

fracs <- c(0.2, 0.05, 0.02, 0.01, 0.005, 0.001)
seeds <- 1:5

res <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]
  for (seed in seeds) {
    d <- mkdat(seed); d$y <- sp$sim(d)
    for (fr in fracs) {
      ## The sigma_frac ladder is left on, and `used` records the rung the
      ## fit actually settled on. A row the ladder had to rescue is a row
      ## where `fr` did not work, so `ok` below requires used == frac; that
      ## attributes the outcome to the value being measured rather than to
      ## whatever the ladder fell back to.
      tt <- system.time(fit <- try(suppressWarnings(suppressMessages(
        gamRTMB(sp$fo, family = sp$f, data = d, sigma_frac = fr,
                joint_precision = FALSE))), silent = TRUE))
      bad <- inherits(fit, "try-error")
      res[[length(res) + 1L]] <- data.frame(
        fam = nm, seed = seed, frac = fr, err = bad,
        conv = if (bad) NA else isTRUE(fit$convergence),
        grad = if (bad) NA_real_ else fit$max_grad,
        obj  = if (bad) NA_real_ else fit$objective,
        used = if (bad) NA_real_ else if (is.null(fit$sigma_frac_used)) fr else fit$sigma_frac_used,
        secs = unname(tt[3]))
      cat(sprintf("%-10s s%d frac%-7g %s\n", nm, seed, fr,
                  if (bad) "ERR" else if (isTRUE(fit$convergence)) "OK" else "NOCV"))
    }
  }
}
res <- do.call(rbind, res)
saveRDS(res, "dev/bench-sigma-frac.rds")

## A fit counts only if the outer optimiser both reported convergence and
## left a gradient small enough to believe it.
res$ok <- !res$err & res$conv %in% TRUE & !is.na(res$grad) &
  res$grad < 1e-2 & res$used == res$frac

cat("\n== converged fits (of", length(seeds), "seeds) ==\n")
print(xtabs(ok ~ fam + frac, res))

## Converging to the wrong place is the other failure mode, and the one the
## small fracs are prone to: the smooth stays pinned near its null space.
b <- do.call(rbind, lapply(split(res, list(res$fam, res$seed), drop = TRUE),
  function(z) { best <- suppressWarnings(min(z$obj[z$ok], na.rm = TRUE))
                z$worse <- z$ok & is.finite(best) & z$obj - best > 0.01; z }))
cat("\n== converged, but to a worse optimum than another frac found ==\n")
print(xtabs(worse ~ fam + frac, b))
