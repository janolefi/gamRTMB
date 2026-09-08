#' gamRTMB: distributional regression with mgcv smooths and RTMB
#'
#' Fits GAMLSS-style models in which every parameter of a distribution may
#' carry smooth terms. The pieces are borrowed rather than rebuilt:
#' \pkg{mgcv} constructs the bases and penalties, [mgcv::smooth2random()]
#' turns the penalized coefficients into iid Gaussian random effects,
#' \pkg{RTMB} supplies automatic differentiation and the Laplace
#' approximation, and \pkg{RTMBdist} supplies the log-densities in each
#' distribution's own native parameterisation.
#'
#' @section Entry points:
#' [gamRTMB()] fits a model, [rtmbdist_family()] builds a family object from
#' any suitable RTMBdist density, and [edf()] reports effective degrees of
#' freedom per smooth.
#'
#' @references
#' Wood, S. N. (2011) Fast stable restricted maximum likelihood and marginal
#' likelihood estimation of semiparametric generalized linear models.
#' \emph{JRSS-B} 73(1), 3--36.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73(4), 1071--1081.
#'
#' @keywords internal
#' @import RTMB
#' @importFrom RTMBdist dgamma2
"_PACKAGE"
NULL

## RTMBdist densities are resolved dynamically by name in rtmbdist_family(),
## since the point is to support the whole library rather than a fixed list.
## The importFrom above therefore looks unused: it is there to declare the
## hard dependency that the dynamic lookup relies on.

## `getAll()` creates the parameter objects in the objective's evaluation
## frame, so static analysis cannot see where they come from. (`beta` is not
## listed because base R already has a function of that name.)
utils::globalVariables(c("b", "log_sigma"))
