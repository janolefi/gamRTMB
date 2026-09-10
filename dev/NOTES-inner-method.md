# TMB's inner method: BFGS was tried, and does not work

Not an open question. Measured, negative, removed. Written down so the
experiment is not repeated.

RTMB solves the inner problem -- maximising over the coefficients at fixed
smoothing parameters -- with its own Newton iteration, forming and factorising
the Hessian at every step. On a dense Hessian each step costs `O(p^3)`, so a
quasi-Newton method needing only the gradient looks like it should win. TMB
exposes the switch:

```r
obj$env$inner.method  <- "BFGS"
obj$env$inner.control <- list(maxit = 1e4, reltol = 1e-10)
```

## It is not faster

Dense inner problems (the `smooth2random` route on an `mrf`), cold start, where
BFGS does reach the same mode as Newton -- the objectives agree to 6e-09 and
their gradients to six significant figures:

| dimension | newton | BFGS | ratio |
|---|---|---|---|
| 65 | 0.155 s | 0.157 s | 0.99x |
| 145 | 1.645 s | 1.613 s | 1.02x |
| 257 | 10.982 s | 10.630 s | 1.03x |

The per-step saving is real and the step count rises to cancel it. The Laplace
approximation needs the Hessian's log determinant and the outer gradient needs
third derivatives, so the dense factorisation happens once an outer iteration
either way -- skipping it inside the inner loop buys less than it looks like.

## It is not stable

TMB warm-starts each inner solve from the previous one, which is what makes the
outer loop affordable. Newton re-converges from anywhere; `optim` does not.
Walking a path of eight smoothing parameters and comparing each marginal
objective against a freshly started Newton solve at the same point:

| theta | newton | BFGS | error |
|---|---|---|---|
| -1.50 | 708.03 | 708.03 | 0.00 |
| -1.17 | 606.89 | 616.21 | 9.32 |
| -0.84 | 530.45 | 549.10 | 18.65 |
| -0.51 | 484.19 | 511.91 | 27.72 |
| -0.19 | 458.80 | 495.12 | 36.32 |
| 0.14 | 446.17 | 490.50 | 44.33 |
| 0.47 | 441.16 | 492.82 | 51.65 |
| 0.80 | 440.75 | 498.94 | 58.19 |

Monotone, compounding, because each warm start begins from the previous drifted
point. Calling the same point six times running returns the same wrong value to
six decimals, so `optim` is stopping at once rather than iterating and falling
short. No setting of `reltol` (down to 1e-16), `abstol` or `maxit` changes any
number in that table.

The consequence is not slowness but a silently wrong fit: the outer optimiser
minimises a criterion wrong by tens of nats and reports convergence at a point
where the true REML gradient is `(-9.4, -13.7)`. On a two-smooth Gaussian model
it landed at 14.20 effective degrees of freedom where Newton gives 17.90.
`sdreport()` also cannot report standard deviations for the random effects.

## Why there is no fix worth having

Forcing a cold start on every inner call would restore correctness. But
cold-start BFGS is already no faster than warm-started Newton, so that trade is
strictly worse. There is no version of this that wins.

One thing worth keeping from the exercise: the suspicion that a quasi-Newton
inner mode makes the *outer* criterion noisy, and so needs a looser outer
tolerance, is wrong. From a cold start the two objectives agree to 6e-09. The
problem is entirely the warm start.
