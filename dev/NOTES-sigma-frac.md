# sigma_frac: measured, and it is not a dial with a good setting

The starting value for the variance components is `sigma_frac`, and the
docstring on `.init_pars()` used to justify the default of 0.05 with a study
of 180 fits that is not in this repository. This is that study redone, with a
script that is: `dev/bench-sigma-frac.R`, about ten minutes, no data package
needed. It disagrees with the old claim on one point and confirms it on
another, and the conclusion is that `frac` is the wrong thing to be tuning.

Eight families, five seeds each, `n = 300`, `k = 8`, every response simulated
from the model that is then fitted. A fit counts only if the outer optimiser
reported convergence *and* left `max_grad < 1e-2`; the `sigma_frac` ladder is
on, and a row it had to rescue counts against the value being measured rather
than the rung it fell back to.

## Converged fits, of 5 seeds

```
             frac 0.001  0.005  0.01  0.02  0.05  0.2
  norm              5      5     5     5     5     5
  gamma2            5      5     5     5     5     5
  beta2             5      5     5     5     5     5
  ---
  bcpe              4      3     1     3     1     4
  bct               3      2     1     3     3     3
  skewnorm2         5      5     2     0     1     2
  skewt2            3      2     1     1     2     2
  jsu2              1      0     0     0     0     0
  ---
  2-param total    15     15    15    15    15    15   (of 15)
  4-param total    16     12     5     7     7    11   (of 25)
```

**Confirmed: the well-behaved families do not care.** Every two-parameter
family is 5/5 across the whole range, and their objectives agree to seven
significant figures from 0.2 down to 0.001. The optimum really is wide and
flat, and for these models `sigma_frac` is not worth a thought.

**Not reproduced: "collapsing below 0.02".** That was the load-bearing claim
behind the default, and the four-parameter block does not show it. Below 0.02
is if anything the better half of the range. Either the original designs
covered something these eight families do not, or the claim drifted; it should
not be repeated without a script behind it.

**What is true instead: there is no ordering.** 0.01 is the *worst* value
here, worse than both 0.2 and 0.001 on either side of it. A curve with a peak
can have its default moved onto the peak; this is noise, and moving the
default just buys a different lottery ticket. Note also how much larger the
spread down a column is than across a row: `jsu2` never converges at any
`frac`, `beta2` always does. The family dominates the starting value.

## The other failure mode, which the small values buy into

Converging is not the same as converging to the right place. Counting fits
that converged to an objective more than 0.01 worse than the best any `frac`
found for that family and seed:

```
             frac 0.001  0.005  0.01  0.02  0.05  0.2
  worse optimum     8      2     1     0     0     1
```

This is the half of the old docstring that survives: too clamped a start
leaves the smooth pinned near its null space, where the REML gradient in
`log_sigma` is nearly zero, and the outer optimiser stalls and then reports
success. `frac = 0.001` converges most often and lands in the wrong place most
often. So the extra convergences at the small end are not free.

## What was done about it

The default stayed at 0.05. On this evidence the honest alternative would be
0.2 — 11/25 with a single bad optimum — but the margin over 0.05 is four fits
out of twenty-five across five families, which is not enough to move a default
that every existing fit depends on.

Instead `.sigma_frac_ladder` retries `c(0.005, 0.2, 0.001)` when the first
marginal evaluation is not finite, and reports the rung it settled on. The
rungs are ordered by the table above: the two that did best on count, then the
one that converges most but misses the optimum most. Because a failed start is
detectable on the very first evaluation (see `.probe_finite()`), the ladder
costs nothing on a model that was going to work, so the choice stops being
something a user has to guess right in advance.

The ladder is coverage, not a better default. Nothing here says these three
values are good; what the table says is that no single value is.

## What this does not address

None of it helps `jsu2`, which converges once in thirty. And on the model that
started this — `bcpe` on `film90` with four smooths — the ladder gets the fit
to run but it still reports `max_grad = 7.8e8` and EDF that do not exist. That
is not a starting-value problem at all; see `.inner_indefinite()` for what it
is, and `dev/NOTES-fellner-schall.md` for the engine that would not have it.

## Retrying when the fit stalls, not only when it will not start

Added after fitting 33 textbook GAMLSS models both ways
(`inst/examples/gamlss.R`). Four of them needed hand-tuning, and the pattern
was the same every time: the objective was perfectly finite at the start, so
`.probe_finite()` never fired, and the fit stalled later. The ladder existed
and could not see the failure it was built for.

Measured across the ladder's rungs plus the default, `-REML` at the point each
one reaches:

```
                 sf=0.05      sf=0.005     sf=0.2       sf=0.001
mcycle-TF     579.85 (g 6.0)  685.42 (ok)  579.29 (0.8) 685.86 (ok)
CD4-BCT      4278.49 (NaN)   4274.36 (NaN) 4273.10 (4)  4271.16 (726)
film90-JSU   6080.42 (129)   6065.35 (ok)  6065.35 (ok) 6065.35 (ok)
```

Three things fall out of that table.

**The default is not a good value on `film90`.** It is the only rung that
fails, and every other one reaches the same optimum 15 units of criterion
better. That is the case the retry is for, and it now converges untouched.

**Selection has to be on the criterion.** On `mcycle-TF` the two rungs that
*converge* sit 105 units above the incumbent that does not: a retry that took
a converged candidate on sight would trade a good fit for a badly over-smoothed
one, which is a worse failure than the one being fixed, and a silent one.
Keeping the lowest criterion with the incumbent in the running cannot return a
worse point than it was handed. `mcycle-TF` accordingly moves 579.85 -> 579.29
and still reports `convergence = FALSE`, which is the honest answer.

**`nlminb`'s code is not the convergence test.** Every row above comes back
with a nonzero code and `false convergence (8)`, including `mcycle-TF` at
`sigma_frac = 0.001` sitting on `max|g| = 1.9e-04`. That is what
`.outer_ok()` is for; the threshold is `.efs_defaults$gtol` so that the two
engines mean the same thing by the word.

### What it costs, and what it does not fix

Each rung is a full refit, taken only by a fit that has already failed. On the
four-parameter `dbbmi` BCPE (n = 7294, k = 25) that is 265s becoming 726s, and
the walk finds nothing: `sigma_frac = 0.15` converges there and is not on the
ladder. Whether to widen the ladder is a separate question from this one and
wants its own measurements -- the rungs are shared with the start-time probe,
which has its own tuning above. For now `dbbmi` is pinned in the example suite
and the retry announces itself before spending the time.

`CD4-BCT` and `mcycle-TF` are not fixed by any rung under REML; both converge
under `method = "qREML"`. A criterion fallback is the obvious next thing to
try and is deliberately not done here: it would have to compare two different
criteria to choose between them, which the objective-selection rule above
cannot do.
