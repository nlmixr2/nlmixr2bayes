# Translate `ini({})` population-parameter priors into NIMBLE BUGS syntax
# (Spec/phase 2). Reuses stanPriors.R's parsing and bound machinery
# (.stanPriorLookup, .stanArgIsLiteral, .stanNum) rather than re-deriving
# it -- only the SYNTAX layer (which NIMBLE distribution, which argument
# names, in what order) is new. Deliberately narrower than Stan's
# catalogue: a fixed, individually-verified subset of univariate
# distributions (see .nimblePriorCatalog), refusing anything else with a
# clear message rather than silently approximating it.
#
# Bounds use NIMBLE's native T(dist, lower, upper) truncation, which -- unlike
# Stan's approach in stanPriors.R -- always supplies the correct normalizing
# constant regardless of whether the prior's own arguments are literals or
# reference another parameter, so there is no literal/truncated distinction
# to make here (verified empirically: T(dnorm(0,sd=1), 0, ) at theta=0.5
# gives the exact truncated-normal log-density, dnorm(0.5,0,1,log=TRUE) -
# log(1-pnorm(0,0,1))).
#
# Argument order in each template matches Stan's own canonical order for
# that distribution (the same order stanPriors()/.stanPriorLookup() already
# assumes user-written prior() calls follow), translated to NIMBLE's named
# arguments -- named throughout, never positional, because several of
# NIMBLE's BUGS-registered distributions default a bare second positional
# argument to a PRECISION, not the scale Stan's convention (and this
# package's `ini()` prior catalogue) uses: dnorm(mean, tau, sd, var),
# dlnorm(meanlog, taulog, sdlog, varlog), dlogis(location, rate, scale),
# dt(mu, tau, df, sigma, sigma2) -- confirmed against NIMBLE's own BUGS
# distribution registrations (nimble:::distributionsInputList), not assumed
# from the R-level mirrors of these functions (which use different defaults
# again).

#' NIMBLE template per supported Stan distribution name
#'
#' `template` args are in STAN's canonical order for that distribution
#' (matching what a `prior(x) ~ <distName>(...)` statement is written in);
#' `nargs` is the expected argument count, checked before templating.
#'
#' No `cauchy` entry: confirmed against `nimble:::distributionsInputList`
#' (the authoritative list of natively BUGS-usable distributions, not the
#' R-level mirror functions nimble also exports) that `dcauchy` is NOT
#' among them -- using it in generated model code fails at `nimbleModel()`
#' build time with "density function for dcauchy is not available. It must
#' be a nimbleFunction." Supporting it would need a custom
#' nimbleFunction/registerDistributions() pair (the same pattern
#' dFoceiCondLik uses) and its own validation pass; left for later rather
#' than risk an under-verified addition here.
#' @noRd
.nimblePriorCatalog <- list(
  normal      = list(template = "dnorm(%s, sd = %s)",                  nargs = 2L),
  std_normal  = list(template = "dnorm(0, sd = 1)",                    nargs = 0L),
  lognormal   = list(template = "dlnorm(meanlog = %s, sdlog = %s)",    nargs = 2L),
  gamma       = list(template = "dgamma(shape = %s, rate = %s)",       nargs = 2L),
  beta        = list(template = "dbeta(shape1 = %s, shape2 = %s)",     nargs = 2L),
  exponential = list(template = "dexp(rate = %s)",                     nargs = 1L),
  uniform     = list(template = "dunif(min = %s, max = %s)",           nargs = 2L),
  student_t   = list(template = "dt(df = %s, mu = %s, sigma = %s)",    nargs = 3L),
  logistic    = list(template = "dlogis(location = %s, scale = %s)",   nargs = 2L)
)

#' Translate one prior-call argument to NIMBLE source text
#'
#' Mirrors `.stanArgDeparse()` (R/stanPriors.R): a numeric literal formats
#' via `.stanNum()`; a bare name must reference another population
#' parameter this model map knows about, translated to its `theta[]` slot
#' (NIMBLE has no per-parameter named scalar the way Stan's generated
#' program does -- theta is one array, referenced positionally).
#' @noRd
.nimbleArgText <- function(arg, thetaNames) {
  if (.stanArgIsLiteral(arg)) return(.stanNum(eval(arg, envir = baseenv())))
  if (is.name(arg)) {
    .nm <- as.character(arg)
    .idx <- match(.nm, thetaNames)
    if (is.na(.idx)) {
      stop("prior argument '", .nm, "' does not reference a population ",
           "parameter nimble linked sampling has a slot for", call. = FALSE)
    }
    return(sprintf("theta[%d]", .idx))
  }
  stop("cannot translate prior argument '", deparse1(arg), "' to nimble: ",
       "arguments must be numeric literals or parameter names", call. = FALSE)
}

#' The NIMBLE distribution text for one free theta's prior, truncated to its
#' effective bounds if finite
#'
#' @param name the theta's nlmixr2 name (`map$theta$name[i]`)
#' @param defaultEst,defaultSd used when no `ini({})` prior matches (falls
#'   back to a weakly-informative normal at the theta's own declared bounds)
#' @param priRow the matching row of `stanPriors(ui)$pop` (`kind !=
#'   "multivariate"`), or `NULL`
#' @param lowerFallback,upperFallback the theta's own declared bounds,
#'   used when `priRow` is `NULL`
#' @param thetaNames `map$theta$name`, for resolving parameter-name prior
#'   arguments to `theta[]` slots
#' @noRd
.nimbleThetaPriorText <- function(name, defaultEst, defaultSd, priRow,
                                  lowerFallback, upperFallback, thetaNames) {
  if (is.null(priRow)) {
    .dist <- sprintf("dnorm(%s, sd = %s)", .stanNum(defaultEst), .stanNum(defaultSd))
    .lo <- lowerFallback
    .hi <- upperFallback
  } else {
    .lk <- .stanPriorLookup(priRow$prior)
    .cat <- .nimblePriorCatalog[[.lk$stanName]]
    if (is.null(.cat)) {
      stop("prior distribution '", .lk$name, "' on '", name, "' is not yet ",
           "supported for nimble linked sampling", call. = FALSE)
    }
    if (length(.lk$args) != .cat$nargs) {
      stop("prior '", .lk$name, "' on '", name, "' has ", length(.lk$args), # nocov
           " argument(s), expected ", .cat$nargs, call. = FALSE)            # nocov
    }
    .dist <- if (.cat$nargs == 0L) {
      .cat$template
    } else {
      .args <- lapply(.lk$args, .nimbleArgText, thetaNames = thetaNames)
      do.call(sprintf, c(list(.cat$template), .args))
    }
    .lo <- priRow$lower
    .hi <- priRow$upper
  }
  .loFin <- is.finite(.lo)
  .hiFin <- is.finite(.hi)
  if (.loFin || .hiFin) {
    .dist <- sprintf("T(%s, %s, %s)", .dist,
                     if (.loFin) .stanNum(.lo) else "",
                     if (.hiFin) .stanNum(.hi) else "")
  }
  .dist
}
