#' Parse a distributional model formula
#'
#' The response is taken from the left-hand side of the outer formula and each
#' distributional parameter gets a one-sided formula from the `list()` on the
#' right, in the style of `brms::bf()`. Any parameter the family declares but
#' the user does not mention is given `~1`, which is why the family spec has
#' to exist before formula processing.
#'
#' @param formula e.g. `y ~ list(mu = ~ s(x1) + s(x2), sigma = ~ s(x1))`.
#' @param parnames Character vector of the family's modelled parameters.
#' @return A list with the response expression and one formula per parameter,
#'   ordered as `parnames`.
#' @keywords internal
.parse_formula <- function(formula, parnames) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("`formula` must be two-sided, e.g. y ~ list(mu = ~ s(x))")
  rhs <- formula[[3L]]
  if (!(is.call(rhs) && identical(rhs[[1L]], as.name("list"))))
    stop("the right-hand side must be a list() of per-parameter formulas, ",
         "e.g. y ~ list(mu = ~ s(x), sigma = ~ 1)")
  args <- as.list(rhs)[-1L]
  if (length(args) && (is.null(names(args)) || any(!nzchar(names(args)))))
    stop("every element of the list() must be named after a distributional parameter")
  unknown <- setdiff(names(args), parnames)
  if (length(unknown))
    stop("unknown distributional parameter(s) ", paste(unknown, collapse = ", "),
         "; this family models ", paste(parnames, collapse = ", "))
  env <- environment(formula)
  out <- lapply(args, function(a) {
    f <- eval(a, env)
    if (!inherits(f, "formula")) stop("each list() element must be a formula")
    if (length(f) == 3L) stop("per-parameter formulas must be one-sided")
    environment(f) <- env
    f
  })
  for (p in setdiff(parnames, names(out))) out[[p]] <- stats::as.formula("~1", env)
  list(response = formula[[2L]], par_formulas = out[parnames])
}
