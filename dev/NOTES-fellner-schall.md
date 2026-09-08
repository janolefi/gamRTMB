# Holding a path open for extended Fellner–Schall

Not implemented. This note records what EFS needs so the package structure
does not have to change when it is added, and one result that makes the
implementation much smaller than it looks.

## Why it is a separate engine, not an option

The current engine hands the smoothing parameters to `nlminb` and lets
`RTMB::MakeADFun(random = c("beta", "b"))` supply the REML criterion and its
gradient. That gradient differentiates `log|H|` with respect to the smoothing
parameters, which needs third derivatives of the log-likelihood in the
coefficients — the expensive term, in both time and tape memory.

EFS (Wood & Fasiolo 2017) drops exactly that term and replaces the gradient
step with a multiplicative update that is still guaranteed to increase the
REML criterion. For that to pay off, **nothing may be declared `random`**: the
Laplace machinery must not be built at all, or its tape cost is paid anyway.
So the engine has to own both halves of the fit:

1. an inner Newton loop maximising the penalized log-likelihood over
   `c = (beta, b)` at fixed smoothing parameters, and
2. an outer multiplicative update of the smoothing parameters.

That is a genuinely different control flow, not a different optimiser setting.

## The update is nearly free in this parameterisation

The general EFS step for penalty `S_k` with weight `lambda_k` is

    lambda_k <- lambda_k * [ tr(S_lambda^- S_k) - tr(H^-1 S_k) ] / (c' S_k c)

with `H` the penalized Hessian at the inner mode and `S_lambda = sum_k
lambda_k S_k`. In a general basis each of those traces is real work.

After `smooth2random`, the penalty on block `k` is the **identity** on that
block's coefficients and zero elsewhere, with `lambda_k = 1 / sigma_k^2`. So,
writing `q_k` for the block's dimension and `t_k = sum_{j in k} (H^-1)_jj`:

    tr(S_lambda^- S_k) = q_k / lambda_k
    tr(H^-1 S_k)       = t_k
    c' S_k c           = ||b_k||^2

    lambda_k <- [ q_k - lambda_k * t_k ] / ||b_k||^2

and `q_k - lambda_k t_k = sum_{j in k} (1 - lambda_k (H^-1)_jj)` is precisely
the penalized part of the block's EDF, which `edf()` already computes. The
whole update collapses to

    lambda_k <- edf_k / ||b_k||^2      i.e.   sigma_k^2 <- ||b_k||^2 / edf_k

the classical Fellner–Schall variance update. Nothing new has to be derived:
the two ingredients are `||b_k||^2` and `diag(H^-1)` restricted to each block,
and the second is the same quantity the EDF code needs.

## What the engine still has to supply

- **Inner Newton loop.** Joint gradient and sparse Hessian of the penalized
  negative log-likelihood in `c`, at fixed smoothing parameters, plus step
  halving. Without `random =`, `spHess` is unavailable, so this needs its own
  sparse-Hessian route (`RTMB::MakeTape` / a Hessian tape over `c`).
- **`diag(H^-1)` for a sparse `H`.** Takahashi's recursion off the Cholesky
  factor gives just the needed entries; a dense `solve()` would defeat the
  purpose on larger models.
- **REML criterion value** for monitoring and convergence:
  `-l(c) + 0.5 c' S_lambda c - 0.5 log|S_lambda|_+ + 0.5 log|H|`, minus the
  usual constant. Cheap once the Cholesky of `H` is in hand.
- **Step control.** Wood & Fasiolo's safeguard: halve the step in `log lambda`
  until the criterion does not decrease.

## Seams this requires, and how the package provides them

| Requirement | Provided by |
|---|---|
| The objective must be buildable with or without random effects | `.make_nll()` returns a plain closure over the design; it never mentions `random =`. Only the engine decides that. |
| Penalty structure available as explicit blocks | `design$blocks[[k]]` carries `idx` (coefficient positions), `q`, and `sig_key`; `design$sig_group` carries the id-tying. EFS needs no other penalty representation. |
| `edf()` / `vcov()` must not reach into TMB internals | Both go through `penalized_hessian(fit)`, which returns `H` in the design's own `(beta, b)` ordering plus the index vectors. One implementation per engine; everything downstream is shared. |
| Fitted objects must be interchangeable | The engine fills a fixed contract: `coefficients` (`beta`, `b`), `log_sigma` (always full length, never the mapped-down vector), `design`, `family`, `method`, `engine`. `predict()`, `edf()`, `summary()` use only that. |
| Smoothing parameters tied by `id` | The engine reads `design$sig_group`; the Laplace engine passes it to TMB's `map`, and EFS would instead apply one update per group, using the group's pooled `||b_k||^2` and `edf_k`. |

The one thing the contract deliberately does **not** promise is
`fit$obj`: it is the Laplace engine's RTMB object and may be absent. Anything
needing it (currently `sdreport`-based standard errors) checks first and says
so, rather than assuming.
