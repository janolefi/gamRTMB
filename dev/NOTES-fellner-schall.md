# Extended Fellner-Schall: what it is, and what it cost to build

Implemented, in `R/efs.R`, reachable as `method = "aREML"`. This note keeps
the reasoning that led there, corrects the things it got wrong before the code
existed, and records the sharp edges found on the way.

## It is a criterion, not a second axis

The interface was `method = c("REML", "ML")` crossed with
`engine = c("laplace", "efs")` for a while. That was wrong, and the package
author said so: the way anyone actually thinks about this is three paths --
REML, ML, approximate REML -- not a two-by-two whose fourth cell nobody wants.

The fourth cell was `engine = "efs"` with `method = "ML"`, and it deserved to
go on its own merits. It was approximate twice over: the Fellner-Schall
gradient drops the third-derivative term, *and* the unpenalized coefficients
were profiled at the penalized likelihood's mode rather than at the maximiser
of the criterion. It measured +0.07 nats against Laplace ML, ten times the
~0.002 that aREML costs against REML.

So there is one argument, `method`, with three values, and `"aREML"` means
"the REML criterion, optimised by extended Fellner-Schall". The
approximation is in the *gradient*, not in the criterion, which is why the
reported value is labelled `-REML` and is directly comparable with a
`method = "REML"` fit's.

## Why it is a separate fitting routine, not an optimiser setting

The Laplace engine hands the smoothing parameters to `nlminb` and lets
`RTMB::MakeADFun(random = c("beta", "b"))` supply the REML criterion and its
gradient. That gradient differentiates `log|H|` with respect to the smoothing
parameters, which needs third derivatives of the log-likelihood in the
coefficients — the expensive term, in both time and tape memory.

EFS (Wood & Fasiolo 2017) drops exactly that term and replaces the gradient
step with a multiplicative update. For that to pay off, **nothing may be
declared `random`**: the Laplace machinery must not be built at all, or its
tape cost is paid anyway. So `"aREML"` owns both halves of the fit:

1. an inner Newton loop maximising the penalized log-likelihood over
   `c = (beta, b)` at fixed smoothing parameters, and
2. an outer multiplicative update of the smoothing parameters.

That is a genuinely different control flow, not a different optimiser setting.

## The update is nearly free — but only for iid blocks

The earlier version of this note said the general EFS step

    lambda_k <- lambda_k * [ tr(S_lambda^- S_k) - tr(H^-1 S_k) ] / (c' S_k c)

collapses, after `smooth2random`, to `lambda_k <- edf_k / ||b_k||^2`, and that
"nothing new has to be derived". Both halves are true of an `"iid"` block and
neither is true of the package as it now stands. A `"gmrf"` block keeps its
own sparse precision, and a `"multi"` block combines several penalty matrices
through mgcv's `L` convention. The engine has to carry the general form from
the start, which it does, through one representation:

    Q_k(theta) = sum_i exp(A_i . theta) S_ki

with `A` the `L` matrix generalised to every block — `-2` for the two kinds
that have a single penalty and a log standard deviation. `.block_penalties()`
supplies it, `.block_prec()` is the same statement on the AD tape, and a test
holds the two together. The payoff is that `te()`, `ti()` and adaptive
smooths, when they arrive, are `"multi"` blocks and need nothing added here.

The simplification does survive where it applies: for an `"iid"` block
`tr(Q^-1 lambda S) = q_k` and the numerator really is the block's penalized
EDF. The general case pays for one extra sparse solve per penalty matrix, and
only for blocks that have more than one.

## The tied-weight case needs a different step, and this is the thing to know

Wood & Fasiolo's update is multiplicative in the penalty weights, and that is
valid only where the weights are free to move one at a time. A block with an
`L` matrix ties them: `lambda_i = exp(L_i . theta)`.

Applying the update weight by weight and projecting the result back onto
`theta` by least squares looks reasonable and is wrong. It has a fixed point
where the projected log-ratios vanish, which is not where the gradient
vanishes. On the SPDE test field it converged smoothly over 25 iterations to
a point **1.3 nats short** of the Laplace optimum with `|dV/dtheta| = 3.9`,
and every diagnostic said converged. Nothing about the trajectory looked
wrong; the criterion was monotone and the steps shrank geometrically.

What works is a step in the metric the multiplicative update induces.
Linearising `log(n/d)` about `n = d` gives `2(n-d)/(n+d)`, so the exact update
is, to first order,

    delta theta = -M^-1 g,    M = L' diag((n_i + d_i)/4) L.

`M` is positive definite, so the step is always downhill; it reproduces the
multiplicative update exactly on the free-weight case (there is a unit test
for that meeting point); and its fixed point is `g = 0`. Same SPDE field:
0.08 nats from the Laplace optimum, gradient driven to 1e-3.

The free-weight case keeps the logarithmic form rather than using the metric
throughout, because far from the optimum `log(n/d)` takes larger and better
steps than its linearisation.

`.efs_pairs()` decides which parameters are which, and tying spreads: an `id`
can put a tied parameter and a free one on the same penalty matrix.

## What the engine supplies, and how

- **Inner Newton loop** over `c`, with step halving, warm-started from the
  previous outer iteration. Newton rather than quasi-Newton for the reason in
  `NOTES-inner-method.md`: warm-started quasi-Newton inner solves drift, and
  the outer criterion ends up wrong while every flag says success.
- **Block inverses of `H`**, one block at a time — `q_k` right-hand sides
  against the factor already in hand. Solving for all of `H^-1` at once would
  be one call instead of several but would hold an `n` by `nb` dense matrix,
  which is the thing this engine exists to avoid. A Takahashi recursion would
  be better still and is the obvious next improvement.
- **One factorisation per criterion evaluation**, reused for the log
  determinant, the block inverses and the outer step.
- **Step control.** Halve the step in the log parameters until the criterion
  falls. The criterion is the real one — re-solved and re-factorised at the
  candidate — so a step that looks good to the dropped-third-derivative
  gradient and is not gets caught.
- **PSD repair of `H - S`.** Wood & Fasiolo's Theorem 1 needs the data
  Hessian positive semi-definite, and so do the EDF. Eigenvalue clamping
  below 1000 coefficients, a ridge above it. Repairs at intermediate
  iterations are routine — the starting values are often a place with
  negative curvature and the iteration walks away from it — and are counted
  and printed; a repair at the *final* point changes what is reported and
  warns.

## Two things that cost a day

**`RTMB::OBS()` keys on the deparsed name of its argument in a registry global
to RTMB.** `MakeADFun` resets that registry per object; `MakeTape` does not.
So an EFS tape calling `OBS(y)` left `"y"` behind pointing at its own
response, and the *next* model's `OBS(y)` — on a tape or in plain R — returned
it. The second fit in a session converged confidently, gradient 1e-10, to the
minimum of the previous model's likelihood. `.make_nll(obs = FALSE)` is the
fix: `OBS` exists for `obj$simulate()` and one-step-ahead residuals, which
only the Laplace engine builds. `test-efs.R` has the regression test.

**RTMB will not take a sparse second-order Jacobian through `dgmrf`.**
`tape$jacfun(sparse = TRUE)$jacfun(sparse = TRUE)` over a parameter vector
that includes `log_sigma` fails with "Inverse subset: order 2 not yet
implemented": `log|Q(theta)|` already has an inverse subset in its first
derivative. The fix is also the better design, and is what the smoothing
parameters wanted anyway — keep them out of the tape, as `DataEval` data that
can be moved between outer iterations without re-taping. The tape is then over
the coefficients alone, which is the only thing the engine ever differentiates
twice, and the objective really is built once.

While checking that: `DataEval` does re-read on every pass, despite what the
caching discussion in `?Tape` suggests — that only starts after `reorder()`.
What looks like a stale `DataEval` is usually TMB caching `obj$fn(p)` for a
repeated `p`. Nudge the parameter vector when testing this.

## What it is worth: measured

`dev/bench-efs.R`, 500 observations, one run each. `dV` is aREML minus REML on
the criterion, so positive is worse; `dEDF` likewise.

```
                                            dV       dEDF   speed
  gaussian location-scale, two smooths   +0.0022    +0.062   0.79x
  three smooths and a random effect      +0.0045    +0.132   1.06x
  id-tied smoothing parameter            +0.0020    +0.061   1.24x
  gamma, mean and sd smoothed            +0.0001    +0.005   1.97x
  poisson, one smooth plus a spurious    +0.0007    +0.010   0.56x
  markov random field, 64 regions        +0.0014    +0.096   0.35x
  bcpe on film90 subsample               REML does not start
```

And on the full `film90`, against `LaMa::qreml()` rather than against REML,
which cannot start: 42s for aREML and 56s for `qreml`, logLik -5882.8 against
-5883.2, total EDF 33.7 against 33.8.

Two things to read off this. **The approximation is small**: a few
thousandths of a nat on the criterion and about a tenth of an effective
degree of freedom, on models where both converge.

**The speed is not the reason to use it.** Between 0.35x and 2.0x, with no
clear pattern, at this problem size. The third-derivative term REML pays for
is not yet the dominant cost on a few hundred observations and a few dozen
coefficients; where aREML should win is where that term's tape gets large,
which these models do not reach. What it buys at this size is the last row: a
fit that exists at all.

On that row aREML returns a fit where REML cannot start, with effective
degrees of freedom in `[0, 1]` per coefficient, and reports honestly that it
did not get there cleanly on the 600-row subsample: `converged: FALSE`, and
the data Hessian repaired at three iterations including the last, so the
criterion and the EDF are those of the repaired problem. The full 4031-row
fit does converge, in 12 outer iterations. That the smaller, noisier problem
is the harder one is worth remembering when reading the row.

## The inner solve: BFGS, not Newton, and why that was got wrong

The first version of this engine solved the inner problem with a hand-rolled
Newton loop. On `bcpe`/`film90` — four smooths, the model this engine exists
for — it ran for more than twenty minutes without converging, while
`LaMa::qreml()`, which uses `optim(method = "BFGS")`, fits the same model in
56 seconds. Replacing the inner solve with BFGS: **54 seconds, converged in 12
outer iterations**, matching `qreml` to 0.4 on the log-likelihood and 0.2 on
the total EDF. The `film90` subsample in `dev/bench-efs.R` went from 105
seconds to 4.4.

The Newton loop had two defects and one bad argument behind it.

**The bad argument.** "Newton, because the Hessian is needed at the end
anyway." It is needed once per *outer* iteration, at the converged inner
point, for the update's traces and the criterion's log determinant. A Newton
inner loop forms and factorises it once per inner *step* — on this model a
hundred per outer iteration instead of one. The premise was true and the
conclusion did not follow from it.

**No sufficient-decrease condition.** The backtracking accepted any
`fn_new <= f`. That is the textbook way to converge to a point that is not
stationary, and it duly did: with `control = list(trace = 2)` the inner
iterations read

```
   inner 11.1  25  f = 673.474666  max|g| = 899  step = 0.5
   inner 11.1  26  f = 673.473740  max|g| = 510  step = 0.5
   inner 11.1  27  f = 673.473438  max|g| = 899  step = 0.5
   inner 11.1  28  f = 673.472512  max|g| = 510  step = 0.5
```

a two-cycle, taking 0.001 off the objective per iteration until it hit
`inner_maxit = 100`. Adding an Armijo condition was not enough on its own —
still unconverged after twenty minutes.

**An arbitrary step scale.** `bcpe`'s inner Hessian is indefinite by the
family's parameterisation, not by the starting values (that is what
`.inner_indefinite()` establishes). So `.efs_chol()` escalates a ridge by
factors of four until a Cholesky succeeds, and the step that comes back has
whatever length that correction happened to give. BFGS sidesteps the whole
question: it never forms the Hessian during the solve, its direction comes
from a positive definite approximation, and `optim`'s line search enforces a
curvature condition as well as sufficient decrease.

**And the loop itself is gone.** `inner_method = "newton"` is now
[TMB::newton()], which is the solver the Laplace approximation itself uses and
which the hand-rolled one was a poor imitation of: it creates the Cholesky
symbolically once and updates the factor in place rather than refactorising
from scratch, and its `smartsearch` carries an adaptive interpolation between
the Newton step and a gradient step instead of an escalating ridge. It is
faster than the hand-rolled version on every well-conditioned model measured,
and it is not this package's to maintain.

It is the option, not the default, and the reason is the same model. On
`bcpe` it returns in 138 seconds **reporting success**, at a criterion 4500
nats worse, with every smooth collapsed onto its null space (EDF 1, 1, 1, 1).
Not a tolerance artefact: 1e-8 and 1e-10 agree on the wrong answer.
`smartsearch` regularises toward a gradient step when the Hessian is not
positive definite, which is the right response to negative curvature *at a
point*, and `bcpe`'s inner Hessian is indefinite by the family's
parameterisation everywhere. A wrong answer delivered confidently is worse
than a slow one.

Where it does belong, to the same criterion and EDF in every row:

```
                     newton    bfgs
  gaussian, 2 smooths  0.51s   1.26s
  gamma                0.43s   0.62s
  mrf, 196 regions     0.13s   0.34s
  spde, 267 nodes      1.46s   4.06s
```

So: worth reaching for when the coefficient vector is large and the family is
ordinary, a trap when it is not.

## A repaired data Hessian: a caveat, not the observed failure

Worth stating precisely, because an earlier version of this note stated it
wrongly. When `.psd_repair()` fires, `log|H|` and therefore the criterion
belong to the repaired problem, and the repair is a function of the current
point — so successive iterations are comparing criteria of slightly different
problems, and the step control's monotonicity does not formally mean what it
says.

That is a real caveat and it is worth knowing. It is **not** what caused the
96-iteration wander this note previously blamed on it. That fit repaired at
every iteration *and* had a broken inner solve, and the inner solve was the
cause: with BFGS the full-data fit converges in 12 outer iterations despite
repairing at 7 of them. A repair being active is not on its own a reason to
distrust the iteration; the fit still warns when the *final* point is
repaired, which is when the reported criterion and EDF are affected.

## What it does not do

- Under `"ML"` the unpenalized coefficients are left at the penalized
  likelihood's mode rather than at the maximiser of the criterion, which
  ignores the dependence of `log|H_bb|` on them. The Laplace engine does not
  make that approximation. `"REML"` is the better-founded criterion here.
- `vcov()` is conditional on the fitted smoothing parameters — mgcv's
  `unconditional = FALSE`. The correction needs the second derivative of the
  criterion in the smoothing parameters, which is the term EFS drops.
- The metric step is first order, so the tied case converges linearly: 43
  iterations on the SPDE field, where the free-weight cases take 5 to 25. A
  BFGS correction on the tied block would fix that and is not implemented.

## The engine contract, as it turned out

Unchanged from the plan, and it held. `edf()` and `vcov()` go through
`.penalized_hessian()` and `.joint_cov()`, each of which grew one branch.
`fit$obj` is still not promised — EFS has tapes, not an RTMB object — and
`.joint_cov()` checks the engine rather than assuming. One thing the plan did
not anticipate: because EFS forms the penalized Hessian over *every*
coefficient itself, `edf()` now works under `"ML"` as well, which it never
could on the Laplace engine.
