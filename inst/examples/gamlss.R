### ---------------------------------------------------------------------------
### Thirty-three standard GAMLSS examples, fitted twice
###
### Each example below is a textbook `gamlss` fit -- most of them on Rigby and
### Stasinopoulos's own datasets -- paired with the `gamRTMB` model that is
### meant to be the same model. Every dataset comes from a CRAN package:
### gamlss.data, MASS, or base R's datasets.
###
### Run it:
###
###   Rscript inst/examples/gamlss.R              # everything, about four minutes
###   Rscript inst/examples/gamlss.R abdom rent   # only names matching these
###
### or `source()` it and inspect `results`.
###
###
### What is being compared
### ----------------------
### The two packages are not doing the same arithmetic. `gamlss` fits by
### backfitting (the RS algorithm) with P-splines whose smoothing parameter
### comes from a local maximum-likelihood criterion; `gamRTMB` integrates the
### spline coefficients out with a Laplace approximation and takes the
### smoothing parameters from the resulting REML criterion. So the fitted
### curves should agree closely, the effective degrees of freedom roughly, and
### nothing to the last digit.
###
### Where an example has no smooth at all, both are plain maximum likelihood on
### the same likelihood. All nine such rows -- `usair-GA`, `polio-PO`,
### `fabric-PO`, `species-NBI`, `quine-NBI`, `glass-WEI`, `lice-NBI`,
### `cysts-ZIP`, `stylo-ZAP` -- come out with the same AIC to a tenth. They are
### the anchors: they say the densities, the offsets and the prior weights are
### the same on both sides, so whatever the smooth rows disagree about, it is
### the smoothing and not the likelihood.
###
### `diff%` is the largest gap between the two fits' fitted parameter vectors,
### as a percentage of that parameter's own spread, maximised over parameters;
### `worst` names the parameter it came from. Single digits mean the curves lie
### on top of each other. Read a large `diff%` on a constant parameter with
### care: when both fits put an intercept near zero -- a Box-Cox `nu`, a hurdle
### `zeroprob` on data with no zeros -- a ratio to its own magnitude is large
### however negligible the absolute difference.
###
###
### Three things that have to be got right for the models to match
### -------------------------------------------------------------
### **Links.** gamlss and gamRTMB agree on every default link used here except
### one: the Box-Cox families BCCG/BCT/BCPE take `mu.link = "identity"` in
### gamlss and a log link in gamRTMB. Every Box-Cox example below therefore
### passes `mu.link = "log"` explicitly. (On `CD4` the identity link does not
### merely differ, it fails: the backfitted mu goes negative and gamlss stops.)
###
### **Basis size.** mgcv's `s()` defaults to `k = 10`, which caps a smooth at 9
### effective degrees of freedom; gamlss's `pb()` lays down 20 inner knots and
### can spend about 20. On small data neither cap binds and the default is left
### alone. Where it does bind -- `dbbmi`, `film90`, `mcycle`, `rent`,
### `brownfat`, `VictimsOfCrime` -- `k` is raised to 20 or 25, and stated
### explicitly in the call so it is visible that it was raised. Leave it at the
### default on `dbbmi` and gamRTMB reports 25 edf against gamlss's 35, which
### looks like a disagreement about smoothing and is nothing of the sort: the
### mean alone wants 21.8 df and the basis cannot hold them. With `k = 25` it
### spends 21.8 too.
###
### **Parameterisation.** gamlss renames every density's parameters to
### mu/sigma/nu/tau; gamRTMB keeps the density's own names. Usually this is
### only a renaming, but twice below it changes which models are reachable:
###
###   * A gamma with *constant dispersion* is a constant CV in gamlss's GA, and
###     a constant CV is not a constant `sd`. `gamma2(mean, sd)` with `sd = ~ 1`
###     is a different model, and a worse-fitting one. The equivalent is
###     RTMBdist's `gamma(shape, rate)`: hold `shape` constant and put the
###     covariates on `rate`, since log(mu) = log(shape) - log(rate). That
###     reproduces gamlss's GA fit to the digit -- see `usair-GA`, where both
###     packages report a deviance of 303.1602.
###   * A beta-binomial in gamlss is a probability and a dispersion; in
###     RTMBdist it is the two shape parameters. A formula on `shape1`/`shape2`
###     is not a formula on `mu`/`sigma`, so `aep-BB` is deliberately *not* the
###     same model in the two packages. Its likelihoods are still comparable --
###     same family, same terms, same number of parameters -- but its edf are
###     not, and its `diff%` is measured on the fitted mean probability rather
###     than on the parameters.
###
###
### The mapping used in each example, all verified against gamlss.dist to 1e-6
### ---------------------------------------------------------------------------
###   gamlss  gamRTMB      mu ->       sigma ->          nu ->        tau ->
###   NO      norm         mean        sd
###   TF      t2           mu          sigma             df
###   LOGNO   lnorm        meanlog     sdlog
###   GA      gamma2       mean        sd = mu * sigma
###   GA      gamma        shape = 1/sigma^2, rate = shape/mu
###   IG      invgauss     mean        shape = 1/sigma^2
###   WEI     weibull      scale       shape
###   JSU     jsu2         mu          sigma             nu           tau
###   BCCG    bccg         mu          sigma             nu
###   BCT     bct          mu          sigma             nu           tau
###   BCPE    bcpe         mu          sigma             nu           tau
###   PO      pois         lambda
###   NBI     nbinom2      mu          size = 1/sigma
###   ZIP     zipois       lambda      zeroprob
###   ZAP     hpois        lambda      zeroprob
###   BI      binom        prob                                 (size fixed)
###   BB      betabinom    shape1 = mu/sigma, shape2 = (1-mu)/sigma
###
### Others that hold but are not used below: LO/logis, EXP/exp (rate = 1/mu),
### SN1/skewnorm, ST2/skewt, PE/powerexp, PE2/powerexp2, GG/gengamma,
### exGAUS/exgauss (lambda = 1/nu), ZINBI/zinbinom2, ZAGA/zigamma2,
### BE/beta2 (phi = 1/sigma^2 - 1), GEOM/geom.ad (prob = 1/(1 + mu)).
### ---------------------------------------------------------------------------

## ---- packages ---------------------------------------------------------------

need <- c("gamlss", "gamlss.data", "gamlss.dist", "MASS")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss))
  stop("this script needs: ", paste(miss, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(gamRTMB)
  library(gamlss)
  library(gamlss.data)
})

## gamlss masks edf(), which gamRTMB also exports
edf <- gamRTMB::edf

options(width = 170)

## Every gamlss fit gets the same control: silent, and a generous iteration
## budget. Twenty backfitting cycles, the package default, is not enough for
## the four-parameter fits here -- `film90-JSU` needs thirty-six -- and a model
## that has not converged is not worth comparing to one that has.
glc <- gamlss.control(n.cyc = 200, trace = FALSE)

## ---- the comparison machinery -----------------------------------------------

## A `map` entry is a function of the two fits returning the pair of vectors
## whose agreement is to be measured. `same()` covers the usual case, where one
## gamlss parameter corresponds to one gamRTMB parameter after an optional
## transform of the gamlss side; anything else is written out longhand.
same <- function(gpar, rpar, f = NULL) {
  g <- function(m, p) {
    a <- fitted(m, gpar)
    list(if (is.null(f)) a else f(a, m), p[[rpar]])
  }
  attr(g, "par") <- rpar
  g
}

ex <- function(name, data, family, note, gl, rt, map = list()) {
  if (length(map) && is.null(names(map)))
    names(map) <- vapply(map, function(g) {
      nm <- attr(g, "par"); if (is.null(nm)) "?" else nm
    }, "")
  list(name = name, data = data, family = family, note = note,
       gl = gl, rt = rt, map = map)
}

## Fit, timing it, and decide whether it counts as converged. gamlss signals
## non-convergence with a warning rather than a flag, so the warning is caught.
fit_gamlss <- function(f, d) {
  m <- NULL; conv <- TRUE; msg <- NA_character_
  sec <- system.time(
    m <- withCallingHandlers(
      tryCatch(f(d), error = function(e) { msg <<- conditionMessage(e); NULL }),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w), ignore.case = TRUE))
          conv <<- FALSE
        invokeRestart("muffleWarning")
      })
  )[["elapsed"]]
  list(fit = m, sec = sec, msg = msg,
       ok = !is.null(m) && conv && isTRUE(m$converged))
}

fit_gamRTMB <- function(f, d) {
  m <- NULL; msg <- NA_character_
  sec <- system.time(
    m <- suppressWarnings(
      tryCatch(f(d), error = function(e) { msg <<- conditionMessage(e); NULL }))
  )[["elapsed"]]
  list(fit = m, sec = sec, msg = msg,
       ok = !is.null(m) && isTRUE(m$convergence))
}

## Largest absolute discrepancy, as a fraction of the reference's own spread.
## The spread is the sd, unless the reference is near constant, when its level
## is the only scale there is.
discrep <- function(a, b) {
  s <- max(stats::sd(b), mean(abs(b)), 1e-8)
  max(abs(a - b)) / s
}

agreement <- function(e, m, p) {
  if (!length(e$map)) return(list(value = NA_real_, par = ""))
  v <- vapply(e$map, function(g) {
    ab <- tryCatch(g(m, p), error = function(err) NULL)
    if (is.null(ab)) return(NA_real_)
    discrep(as.numeric(ab[[1]]), as.numeric(ab[[2]]))
  }, numeric(1))
  if (all(is.na(v))) return(list(value = NA_real_, par = ""))
  i <- which.max(replace(v, is.na(v), -Inf))
  list(value = v[[i]], par = names(e$map)[i])
}

run_one <- function(e) {
  d <- e$data()
  g <- fit_gamlss(e$gl, d)
  r <- fit_gamRTMB(e$rt, d)
  a <- if (g$ok && r$ok)
    agreement(e, g$fit, predict(r$fit, type = "response")) else
    list(value = NA_real_, par = "")
  status <- paste(c(if (!g$ok) "gamlss failed", if (!r$ok) "gamRTMB failed"),
                  collapse = "; ")
  data.frame(
    example = e$name,
    n       = nrow(d),
    family  = e$family,
    edf.gl  = if (g$ok) round(g$fit$df.fit, 2) else NA_real_,
    edf.rt  = if (r$ok) round(attr(edf(r$fit), "edf.total"), 2) else NA_real_,
    AIC.gl  = if (g$ok) round(AIC(g$fit), 1) else NA_real_,
    AIC.rt  = if (r$ok) round(AIC(r$fit), 1) else NA_real_,
    s.gl    = round(g$sec, 2),
    s.rt    = round(r$sec, 2),
    `diff%` = if (is.na(a$value)) NA_real_ else round(100 * a$value, 1),
    worst   = a$par,
    status  = if (nzchar(status)) status else "",
    stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
}

## ---- the examples -----------------------------------------------------------

examples <- list(

  ## == Location and scale, and the LMS centile models =========================

  ex("abdom-NO", function() { data(abdom); abdom }, "NO / norm",
     "Foetal abdominal circumference against gestational age: the introductory
      GAMLSS example, a normal with both parameters smooth.",
     function(d) gamlss(y ~ pb(x), sigma.fo = ~ pb(x), data = d, control = glc),
     function(d) gamRTMB(y ~ list(mean = ~ s(x), sd = ~ s(x)), data = d),
     list(same("mu", "mean"), same("sigma", "sd"))),

  ex("abdom-TF", function() { data(abdom); abdom }, "TF / t2",
     "The same data with a t: a third parameter for the tails, held constant.",
     function(d) gamlss(y ~ pb(x), sigma.fo = ~ pb(x), nu.fo = ~ 1,
                        family = TF, data = d, control = glc),
     function(d) gamRTMB(y ~ list(mu = ~ s(x), sigma = ~ s(x), df = ~ 1),
                         family = fam("t2"), data = d),
     list(same("mu", "mu"), same("sigma", "sigma"), same("nu", "df"))),

  ex("abdom-BCCG", function() { data(abdom); abdom }, "BCCG / bccg",
     "LMS: Box-Cox Cole-Green, the original centile method. sigma is a
      coefficient of variation here, not a standard deviation, and nu is a
      Box-Cox power -- which both fits put near zero, so its diff% is a ratio
      of two negligible numbers and means nothing.",
     function(d) gamlss(y ~ pb(x), sigma.fo = ~ pb(x), nu.fo = ~ 1,
                        family = BCCG(mu.link = "log"), data = d,
                        control = glc),
     function(d) gamRTMB(y ~ list(mu = ~ s(x), sigma = ~ s(x), nu = ~ 1),
                         family = fam("bccg"), data = d),
     list(same("mu", "mu"), same("sigma", "sigma"), same("nu", "nu"))),

  ex("abdom-BCT", function() { data(abdom); abdom }, "BCT / bct",
     "LMS with a fourth parameter for kurtosis: GAMLSS's flagship model, and
      the one place on easy data where the two genuinely part company. Neither
      the link nor the basis explains it: on 610 points the two criteria simply
      choose different amounts of smoothing, and here gamlss chooses better --
      more edf and a lower AIC. Worth knowing rather than glossing over.",
     function(d) gamlss(y ~ pb(x), sigma.fo = ~ pb(x), nu.fo = ~ 1, tau.fo = ~ 1,
                        family = BCT(mu.link = "log"), data = d, control = glc),
     function(d) gamRTMB(y ~ list(mu = ~ s(x), sigma = ~ s(x), nu = ~ 1,
                                  tau = ~ 1),
                         family = fam("bct"), data = d),
     list(same("mu", "mu"), same("sigma", "sigma"), same("nu", "nu"),
          same("tau", "tau"))),

  ex("dbbmi-BCPE", function() { data(dbbmi); dbbmi }, "BCPE / bcpe",
     "Dutch boys' BMI against age, all four parameters smooth: the LMS-P
      reference fit, and by a wide margin the largest and slowest model in this
      file -- about two minutes per side. k = 25 on the mean because gamlss
      spends 21.8 df there and mgcv's default basis cannot reach it; with k = 25
      gamRTMB spends 21.8 as well. sigma_frac is raised for the same reason as
      in film90-JSU, and without it this fit stops short of converging.",
     function(d) gamlss(bmi ~ pb(age), sigma.fo = ~ pb(age), nu.fo = ~ pb(age),
                        tau.fo = ~ pb(age), family = BCPE(mu.link = "log"),
                        data = d, control = glc),
     function(d) gamRTMB(bmi ~ list(mu = ~ s(age, k = 25),
                                    sigma = ~ s(age, k = 20),
                                    nu = ~ s(age, k = 20),
                                    tau = ~ s(age, k = 20)),
                         family = fam("bcpe"), data = d, sigma_frac = 0.15),
     list(same("mu", "mu"), same("sigma", "sigma"), same("nu", "nu"),
          same("tau", "tau"))),

  ex("CD4-BCT", function() { data(CD4); CD4[CD4$cd4 > 0, ] }, "BCT / bct",
     "CD4 counts in uninfected children against age: skew, heteroscedastic, and
      the hardest likelihood here for both packages. The one response that has
      to be filtered -- the Box-Cox families have no mass at zero and the data
      have one zero -- and the one gamRTMB fit that needs both a different
      smoothness criterion and a smaller sigma_frac to converge. BCCG on these
      data does not converge in gamlss at any number of cycles; the extra
      kurtosis parameter fixes that.",
     function(d) gamlss(cd4 ~ pb(age), sigma.fo = ~ pb(age), nu.fo = ~ 1,
                        tau.fo = ~ 1, family = BCT(mu.link = "log"), data = d,
                        control = glc),
     function(d) gamRTMB(cd4 ~ list(mu = ~ s(age), sigma = ~ s(age), nu = ~ 1,
                                    tau = ~ 1),
                         family = fam("bct"), data = d,
                         method = "aREML", sigma_frac = 0.02),
     list(same("mu", "mu"), same("sigma", "sigma"))),

  ex("mcycle-NO", function() { data(mcycle, package = "MASS"); mcycle },
     "NO / norm",
     "The motorcycle crash data: a wiggly mean and a variance that moves with
      it, the standard argument for modelling the scale at all. Both packages
      spend 25-ish effective degrees of freedom on 133 points.",
     function(d) gamlss(accel ~ pb(times), sigma.fo = ~ pb(times), data = d,
                        control = glc),
     function(d) gamRTMB(accel ~ list(mean = ~ s(times, k = 20),
                                      sd = ~ s(times, k = 10)), data = d),
     list(same("mu", "mean"), same("sigma", "sd"))),

  ex("mcycle-TF", function() { data(mcycle, package = "MASS"); mcycle },
     "TF / t2",
     "The same, with a t response: does the apparent heteroscedasticity survive
      letting the tails be heavy? One of two fits here needing
      method = \"aREML\" -- the default REML optimisation of the smoothing
      parameters stalls on this likelihood with a gradient of order 1.",
     function(d) gamlss(accel ~ pb(times), sigma.fo = ~ pb(times), nu.fo = ~ 1,
                        family = TF, data = d, control = glc),
     function(d) gamRTMB(accel ~ list(mu = ~ s(times, k = 20),
                                      sigma = ~ s(times, k = 10), df = ~ 1),
                         family = fam("t2"), data = d, method = "aREML"),
     list(same("mu", "mu"), same("sigma", "sigma"))),

  ex("hodges-RE", function() { data(hodges); hodges }, "NO / norm",
     "Prescription costs by US state: a random intercept, from gamlss's
      random() against mgcv's bs = \"re\". Two different estimators of the same
      variance component, and they land in the same place.",
     function(d) gamlss(prind ~ random(state), data = d, control = glc),
     function(d) gamRTMB(prind ~ list(mean = ~ s(state, bs = "re"), sd = ~ 1),
                         data = d),
     list(same("mu", "mean"))),

  ## == Positive continuous responses ==========================================

  ex("rent-GA", function() { data(rent); rent }, "GA / gamma2",
     "Munich rents against floor space, both parameters smooth. GA's second
      parameter is a coefficient of variation, so it is multiplied by the mean
      to be compared with a standard deviation.",
     function(d) gamlss(R ~ pb(Fl), sigma.fo = ~ pb(Fl), family = GA, data = d,
                        control = glc),
     function(d) gamRTMB(R ~ list(mean = ~ s(Fl, k = 20), sd = ~ s(Fl, k = 20)),
                         family = fam("gamma2"), data = d),
     list(same("mu", "mean"),
          same("sigma", "sd", function(a, m) a * fitted(m, "mu")))),

  ex("rent-GA-multi", function() { data(rent); rent }, "GA / gamma2",
     "The same rents with the covariates the textbook uses: two smooths and
      three factors in the mean, one smooth in the scale.",
     function(d) gamlss(R ~ pb(Fl) + pb(A) + H + loc + B, sigma.fo = ~ pb(Fl),
                        family = GA, data = d, control = glc),
     function(d) gamRTMB(R ~ list(mean = ~ s(Fl, k = 20) + s(A, k = 20) + H +
                                    loc + B, sd = ~ s(Fl, k = 20)),
                         family = fam("gamma2"), data = d),
     list(same("mu", "mean"),
          same("sigma", "sd", function(a, m) a * fitted(m, "mu")))),

  ex("rent-IG", function() { data(rent); rent }, "IG / invgauss",
     "Inverse Gaussian on the same response: a heavier right tail than the
      gamma, and a shape rather than a dispersion in RTMBdist's version.",
     function(d) gamlss(R ~ pb(Fl), sigma.fo = ~ pb(Fl), family = IG, data = d,
                        control = glc),
     function(d) gamRTMB(R ~ list(mean = ~ s(Fl, k = 20),
                                  shape = ~ s(Fl, k = 20)),
                         family = fam("invgauss"), data = d),
     list(same("mu", "mean"),
          same("sigma", "shape", function(a, m) 1 / a^2))),

  ex("usair-GA", function() { data(usair); usair }, "GA / gamma",
     "US air pollution against six city characteristics, purely parametric. The
      sharpest anchor in the file: constant dispersion, no smooth anywhere, and
      the shape/rate parameterisation, so the two packages maximise the same
      eight-dimensional likelihood and report the same deviance to six figures.",
     function(d) gamlss(y ~ x1 + x2 + x3 + x4 + x5 + x6, family = GA, data = d,
                        control = glc),
     function(d) gamRTMB(y ~ list(shape = ~ 1,
                                  rate = ~ x1 + x2 + x3 + x4 + x5 + x6),
                         family = fam("gamma"), data = d),
     list(mean = function(m, p) list(fitted(m, "mu"), p$shape / p$rate),
          sigma = function(m, p) list(fitted(m, "sigma"), 1 / sqrt(p$shape)))),

  ex("airquality-GA",
     function() stats::na.omit(datasets::airquality[, c("Ozone", "Solar.R",
                                                        "Wind", "Temp")]),
     "GA / gamma",
     "New York ozone against three weather variables: additive smooths in the
      mean, constant dispersion. The gamma GAM everyone fits first, and again
      through shape/rate so that 'constant dispersion' means the same thing on
      both sides.",
     function(d) gamlss(Ozone ~ pb(Solar.R) + pb(Wind) + pb(Temp), family = GA,
                        data = d, control = glc),
     function(d) gamRTMB(Ozone ~ list(shape = ~ 1,
                                      rate = ~ s(Solar.R) + s(Wind) + s(Temp)),
                         family = fam("gamma"), data = d),
     list(mean = function(m, p) list(fitted(m, "mu"), p$shape / p$rate),
          sigma = function(m, p) list(fitted(m, "sigma"), 1 / sqrt(p$shape)))),

  ex("plasma-GA", function() { data(plasma); plasma[plasma$betaplasma > 0, ] },
     "GA / gamma",
     "Plasma beta-carotene against age, BMI, supplement use and smoking status:
      smooths and factors together. One patient records a zero, which a gamma
      cannot accommodate.",
     function(d) gamlss(betaplasma ~ pb(age) + pb(bmi) + vituse + smokstat,
                        family = GA, data = d, control = glc),
     function(d) gamRTMB(betaplasma ~ list(shape = ~ 1,
                                           rate = ~ s(age) + s(bmi) + vituse +
                                             smokstat),
                         family = fam("gamma"), data = d),
     list(mean = function(m, p) list(fitted(m, "mu"), p$shape / p$rate),
          sigma = function(m, p) list(fitted(m, "sigma"), 1 / sqrt(p$shape)))),

  ex("film90-LOGNO", function() { data(film90); film90 }, "LOGNO / lnorm",
     "Hollywood film revenue against opening-week takings, on the log scale.
      Four thousand points, both parameters smooth, k = 20 to match pb()'s
      basis.",
     function(d) gamlss(lborev1 ~ pb(lboopen), sigma.fo = ~ pb(lboopen),
                        family = LOGNO, data = d, control = glc),
     function(d) gamRTMB(lborev1 ~ list(meanlog = ~ s(lboopen, k = 20),
                                        sdlog = ~ s(lboopen, k = 20)),
                         family = fam("lnorm"), data = d),
     list(same("mu", "meanlog"), same("sigma", "sdlog"))),

  ex("film90-JSU", function() { data(film90); film90 }, "JSU / jsu2",
     "The same revenues with Johnson's SU: skewness and kurtosis on top of the
      location and scale. The slowest gamlss fit here -- thirty-six backfitting
      cycles, well past the package default of twenty -- and the one gamRTMB
      fit that needs sigma_frac raised, without which REML settles on a much
      too smooth solution.",
     function(d) gamlss(lborev1 ~ pb(lboopen), sigma.fo = ~ pb(lboopen),
                        nu.fo = ~ 1, tau.fo = ~ 1, family = JSU, data = d,
                        control = glc),
     function(d) gamRTMB(lborev1 ~ list(mu = ~ s(lboopen, k = 20),
                                        sigma = ~ s(lboopen, k = 20),
                                        nu = ~ 1, tau = ~ 1),
                         family = fam("jsu2"), data = d, sigma_frac = 0.25),
     list(same("mu", "mu"), same("sigma", "sigma"), same("nu", "nu"),
          same("tau", "tau"))),

  ex("glass-WEI", function() { data(glass); glass }, "WEI / weibull",
     "Breaking strength of glass fibres, no covariates: the textbook 'fit a
      distribution' exercise. Both parameters are intercepts, so this is two
      optimisers on the same two-dimensional surface.",
     function(d) gamlss(strength ~ 1, sigma.fo = ~ 1, family = WEI, data = d,
                        control = glc),
     function(d) gamRTMB(strength ~ list(shape = ~ 1, scale = ~ 1),
                         family = fam("weibull"), data = d),
     list(same("mu", "scale"), same("sigma", "shape"))),

  ## == Counts =================================================================

  ex("aids-PO", function() { data(aids); aids }, "PO / pois",
     "Quarterly AIDS registrations in England and Wales: a smooth trend and a
      quarter effect. The Poisson fit is the one that shows the overdispersion
      -- compare its AIC with the next row's.",
     function(d) gamlss(y ~ pb(x) + qrt, family = PO, data = d, control = glc),
     function(d) gamRTMB(y ~ list(lambda = ~ s(x) + qrt), family = fam("pois"),
                         data = d),
     list(same("mu", "lambda"))),

  ex("aids-NBI", function() { data(aids); aids }, "NBI / nbinom2",
     "The same counts with a negative binomial, which is what the AIC asks for.
      NBI's sigma is the reciprocal of RTMBdist's size.",
     function(d) gamlss(y ~ pb(x) + qrt, sigma.fo = ~ 1, family = NBI, data = d,
                        control = glc),
     function(d) gamRTMB(y ~ list(mu = ~ s(x) + qrt, size = ~ 1),
                         family = fam("nbinom2"), data = d),
     list(same("mu", "mu"), same("sigma", "size", function(a, m) 1 / a))),

  ex("polio-PO", function() {
       data(polio)
       tt <- seq_along(polio) / 12
       data.frame(y = as.numeric(polio), t = tt,
                  cos1 = cos(2 * pi * tt), sin1 = sin(2 * pi * tt),
                  cos2 = cos(4 * pi * tt), sin2 = sin(4 * pi * tt))
     }, "PO / pois",
     "Monthly US polio incidence, Zeger's trend-plus-harmonics model. Purely
      parametric, so another exact-agreement anchor.",
     function(d) gamlss(y ~ t + cos1 + sin1 + cos2 + sin2, family = PO, data = d,
                        control = glc),
     function(d) gamRTMB(y ~ list(lambda = ~ t + cos1 + sin1 + cos2 + sin2),
                         family = fam("pois"), data = d),
     list(same("mu", "lambda"))),

  ex("fabric-PO", function() { data(fabric); fabric }, "PO / pois",
     "Faults in rolls of fabric, with log(length) as an offset: a rate model,
      and a check that offset() survives inside the formula list.",
     function(d) gamlss(y ~ x + offset(log(leng)), family = PO, data = d,
                        control = glc),
     function(d) gamRTMB(y ~ list(lambda = ~ x + offset(log(leng))),
                         family = fam("pois"), data = d),
     list(same("mu", "lambda"))),

  ex("species-NBI", function() { data(species); species }, "NBI / nbinom2",
     "Fish species against lake area, on the log scale, with the dispersion
      modelled too. Parametric in both parameters.",
     function(d) gamlss(fish ~ log(lake), sigma.fo = ~ log(lake), family = NBI,
                        data = d, control = glc),
     function(d) gamRTMB(fish ~ list(mu = ~ log(lake), size = ~ log(lake)),
                         family = fam("nbinom2"), data = d),
     list(same("mu", "mu"), same("sigma", "size", function(a, m) 1 / a))),

  ex("quine-NBI", function() { data(quine, package = "MASS"); quine },
     "NBI / nbinom2",
     "Days absent from school in rural New South Wales: the standard
      overdispersed-count dataset, four factors, no smooths.",
     function(d) gamlss(Days ~ Eth + Sex + Age + Lrn, sigma.fo = ~ 1,
                        family = NBI, data = d, control = glc),
     function(d) gamRTMB(Days ~ list(mu = ~ Eth + Sex + Age + Lrn, size = ~ 1),
                         family = fam("nbinom2"), data = d),
     list(same("mu", "mu"), same("sigma", "size", function(a, m) 1 / a))),

  ex("LGAclaims-NBI", function() { data(LGAclaims); LGAclaims },
     "NBI / nbinom2",
     "Third-party motor claims by local government area, population as an
      offset and two smooth covariates.",
     function(d) gamlss(Claims ~ pb(L_Popdensity) + pb(L_Accidents) +
                          offset(L_Population), sigma.fo = ~ 1, family = NBI,
                        data = d, control = glc),
     function(d) gamRTMB(Claims ~ list(mu = ~ s(L_Popdensity) + s(L_Accidents) +
                                         offset(L_Population), size = ~ 1),
                         family = fam("nbinom2"), data = d),
     list(same("mu", "mu"), same("sigma", "size", function(a, m) 1 / a))),

  ex("tidal-NBI", function() { data(tidal); tidal }, "NBI / nbinom2",
     "Counts of tidal events against height of tide, with a three-level factor.
      Fitting this as a Poisson instead sends both packages chasing the
      overdispersion with the smooth, and they chase it by different amounts;
      with the negative binomial they agree.",
     function(d) gamlss(number ~ pb(vertht) + factor(ht), sigma.fo = ~ 1,
                        family = NBI, data = d, control = glc),
     function(d) gamRTMB(number ~ list(mu = ~ s(vertht) + factor(ht),
                                       size = ~ 1),
                         family = fam("nbinom2"), data = d),
     list(same("mu", "mu"), same("sigma", "size", function(a, m) 1 / a))),

  ex("lice-NBI", function() { data(lice); lice }, "NBI / nbinom2",
     "Head lice counts on schoolchildren, supplied as a frequency table: the
      example for prior weights. 71 distinct counts standing for 1083 children.",
     function(d) gamlss(head ~ 1, sigma.fo = ~ 1, family = NBI, data = d,
                        weights = d$freq, control = glc),
     function(d) gamRTMB(head ~ list(mu = ~ 1, size = ~ 1),
                         family = fam("nbinom2"), data = d, weights = freq),
     list(same("mu", "mu"), same("sigma", "size", function(a, m) 1 / a))),

  ex("cysts-ZIP", function() { data(cysts); cysts }, "ZIP / zipois",
     "Kidney cysts in embryonic mice, again a frequency table. Zero-inflated:
      most of the mass at zero is structural rather than Poisson.",
     function(d) gamlss(y ~ 1, sigma.fo = ~ 1, family = ZIP, data = d,
                        weights = d$f, control = glc),
     function(d) gamRTMB(y ~ list(lambda = ~ 1, zeroprob = ~ 1),
                         family = fam("zipois"), data = d, weights = f),
     list(same("mu", "lambda"), same("sigma", "zeroprob"))),

  ex("stylo-ZAP", function() { data(stylo); stylo }, "ZAP / hpois",
     "Word-frequency counts from a stylometric study, weighted by frequency. A
      hurdle Poisson: the zero probability is free of the Poisson mean. These
      data have no zeros at all, so both fits drive that probability to the
      boundary -- which is why the diff% is meaningless and the AICs identical.",
     function(d) gamlss(word ~ 1, sigma.fo = ~ 1, family = ZAP, data = d,
                        weights = d$freq, control = glc),
     function(d) gamRTMB(word ~ list(lambda = ~ 1, zeroprob = ~ 1),
                         family = fam("hpois"), data = d, weights = freq),
     list(same("mu", "lambda"), same("sigma", "zeroprob"))),

  ## == Binary and binomial ====================================================

  ex("brownfat-BI", function() { data(brownfat); brownfat }, "BI / binom",
     "Whether brown adipose tissue is detected, against age, BMI and outside
      temperature: logistic regression with three smooths on 4842 patients.",
     function(d) gamlss(brownfat ~ pb(age) + pb(BMI) + pb(exttemp) +
                          factor(sex), family = BI, data = d, control = glc),
     function(d) gamRTMB(brownfat ~ list(prob = ~ s(age, k = 20) +
                                           s(BMI, k = 20) +
                                           s(exttemp, k = 20) + factor(sex)),
                         family = fam("binom", fixed = list(size = 1)),
                         data = d),
     list(same("mu", "prob"))),

  ex("VictimsOfCrime-BI", function() { data(VictimsOfCrime); VictimsOfCrime },
     "BI / binom",
     "Whether a crime was reported, against the victim's age. One smooth,
      10590 observations: the largest n here, and a smooth gamlss spends 13 df
      on, so k = 20.",
     function(d) gamlss(reported ~ pb(age), family = BI, data = d,
                        control = glc),
     function(d) gamRTMB(reported ~ list(prob = ~ s(age, k = 20)),
                         family = fam("binom", fixed = list(size = 1)),
                         data = d),
     list(same("mu", "prob"))),

  ex("InfMort-BI", function() { data(InfMort); InfMort }, "BI / binom",
     "Infant deaths out of live births by region, against illiteracy and log
      GDP. Binomial counts with a genuine denominator, so gamRTMB takes the
      number of trials from a column with fixed = list(size = \"bornalive\"),
      which is what gamlss's cbind() response says. Wildly overdispersed for a
      binomial, and the two smoothing criteria disagree about what to do with
      that: gamlss shrinks both smooths away entirely, gamRTMB keeps 3.5 df.",
     function(d) gamlss(cbind(dead, bornalive - dead) ~ pb(illit) + pb(lGDP),
                        family = BI, data = d, control = glc),
     function(d) gamRTMB(dead ~ list(prob = ~ s(illit) + s(lGDP)),
                         family = fam("binom",
                                      fixed = list(size = "bornalive")),
                         data = d),
     list(same("mu", "prob"))),

  ex("aep-BB", function() { data(aep); aep }, "BB / betabinom",
     "Inappropriate days of hospital stay out of total days, the standard
      beta-binomial example. NOT the same model in the two packages -- gamlss
      puts the terms on a probability and a dispersion, gamRTMB on the two
      shape parameters -- so read the likelihoods, not the edf. The diff% is on
      the fitted mean probability.",
     function(d) gamlss(cbind(noinap, los - noinap) ~ ward + year + pb(age),
                        sigma.fo = ~ 1, family = BB, data = d, control = glc),
     function(d) gamRTMB(noinap ~ list(shape1 = ~ ward + year + s(age),
                                       shape2 = ~ ward + year + s(age)),
                         family = fam("betabinom", fixed = list(size = "los")),
                         data = d),
     list(prob = function(m, p)
       list(fitted(m, "mu"), p$shape1 / (p$shape1 + p$shape2))))
)

## ---- run --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
nms <- vapply(examples, `[[`, "", "name")
keep <- if (length(args))
  Reduce(`|`, lapply(args, grepl, x = nms, fixed = TRUE)) else
  rep(TRUE, length(examples))
chosen <- examples[keep]
if (!length(chosen)) stop("no example matched ", paste(args, collapse = " "))

cat(sprintf("Fitting %d example%s with gamlss %s and gamRTMB %s\n\n",
            length(chosen), if (length(chosen) == 1) "" else "s",
            packageVersion("gamlss"), packageVersion("gamRTMB")))

results <- do.call(rbind, lapply(chosen, function(e) {
  cat(sprintf("  %-20s ", e$name)); utils::flush.console()
  row <- run_one(e)
  cat(sprintf("%7.1fs / %7.1fs%s\n", row$s.gl, row$s.rt,
              if (nzchar(row$status)) paste0("   [", row$status, "]") else ""))
  utils::flush.console()
  row
}))

cat("\n")
print(results, row.names = FALSE)

cat("\n.gl = gamlss, .rt = gamRTMB.",
    "\ndiff%: largest gap between the two fitted parameter vectors, as a",
    "\n       percentage of that parameter's own spread; `worst` names it.",
    "\n       Large values on a near-zero intercept mean nothing -- see the",
    "\n       header, and the abdom-BCCG and stylo-ZAP notes.\n")

if (any(nzchar(results$status)))
  cat("\nFits that did not converge are left blank above.\n")

invisible(results)
