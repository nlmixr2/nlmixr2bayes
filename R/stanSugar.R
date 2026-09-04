# Sugar estimation methods over the est="stan" engine, the way "foce"/
# "focei" are members of nlmixr2est's focei family: est="nuts", est="advi",
# and est="pathfinder" run the SAME generated program and linked likelihood
# with the matching stanControl(algorithm=) forced, and the fit records
# the sugar name as its $est -- letting a user pick the ALGORITHM
# ("nuts"/"advi"/"pathfinder") instead of the software tool ("stan").
#
# Each sugar method carries its own control constructor -- nutsControl(),
# adviControl(), pathfinderControl() -- built the way nlmixr2est's
# foceControl() is built from foceiControl(): the SAME list contents as
# stanControl(), the distinguishing option forced, re-classed so
# class(control)[1] is "<est>Control".  That length-1 class is what lets
# `nlmixr2(model, data, nutsControl())` -- and the adviControl() and
# pathfinderControl() forms -- infer est= from the control with no est=
# argument at all (nlmixr2est's .nlmixr2inferEst regexes
# `^(.*?)Control$` off `class(control)[1]`).
# The engine downstream still runs on a plain stanControl (the class swap
# is `.stanSugarAsStanControl()`, the twin of .foceControlToFoceiControl);
# the sugar control rides along onto the fit so `nmObjGetControl` and
# `rxUiDeparse` report the control the user actually wrote.

#' Algorithm a sugar est= accepts, and the one it forces
#'
#' Contradictions are errors, not silent overrides -- the one exception
#' is `"NUTS"` under `est="advi"`, which is only ever `stanControl()`'s
#' default leaking through, so it selects the ADVI default instead.
#' @param algorithm requested algorithm
#' @param est sugar estimation name
#' @return the forced algorithm
#' @noRd
.stanSugarAlgorithm <- function(algorithm, est) {
  if (length(algorithm) != 1L) algorithm <- algorithm[1]
  .bad <- function(what) {
    stop("est=\"", est, "\" runs ", what, "; use stanControl(algorithm=\"",
         algorithm, "\") with est=\"stan\" instead", call. = FALSE)
  }
  switch(est,
         nuts = {
           if (!identical(algorithm, "NUTS")) .bad("Stan's NUTS sampler")
           "NUTS"
         },
         advi = {
           # "NUTS" here is stanControl()'s default; est="advi" selects
           # the ADVI default rather than refusing it
           if (identical(algorithm, "NUTS")) algorithm <- "meanfield"
           if (!(algorithm %in% c("meanfield", "fullrank"))) .bad("ADVI")
           algorithm
         },
         pathfinder = {
           if (!(algorithm %in% c("NUTS", "pathfinder"))) .bad("Pathfinder")
           "pathfinder"
         },
         stop("unknown stan sugar est \"", est, "\"", call. = FALSE)) # nocov
}

#' Build a sugar control from stanControl()
#' @param est sugar estimation name
#' @param algorithm requested algorithm
#' @param dots the stanControl() arguments
#' @return an `<est>Control`
#' @noRd
.stanSugarControl <- function(est, algorithm, dots) {
  .ctl <- do.call(stanControl,
                  c(dots, list(algorithm = .stanSugarAlgorithm(algorithm, est))))
  class(.ctl) <- paste0(est, "Control")
  .ctl
}

#' Coerce whatever the user handed the sugar method to its own control
#' @param control a stanControl or any sugar control
#' @param est sugar estimation name
#' @return an `<est>Control`
#' @noRd
.stanSugarAsControl <- function(control, est) {
  .cls <- paste0(est, "Control")
  if (inherits(control, .cls)) return(control)
  control$algorithm <- .stanSugarAlgorithm(control$algorithm, est)
  class(control) <- .cls
  control
}

#' The plain stanControl the engine runs (the .foceControlToFoceiControl twin)
#' @param control a sugar control
#' @return the same list classed `stanControl`
#' @noRd
.stanSugarAsStanControl <- function(control) {
  class(control) <- "stanControl"
  control
}

#' Shared `getValidNlmixrCtl` body for the sugar methods
#' @param control length-1 list holding the control
#' @param est sugar estimation name
#' @return an `<est>Control`
#' @noRd
.stanSugarValidCtl <- function(control, est) {
  .ctl <- control[[1]]
  .cls <- paste0(est, "Control")
  .fn <- switch(est, nuts = nutsControl, advi = adviControl,
                pathfinder = pathfinderControl)
  if (is.null(.ctl)) return(.fn())
  if (inherits(.ctl, .cls)) return(do.call(.fn, .ctl))
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    return(do.call(.fn, .ctl))
  }
  .known <- c("stanControl", "nutsControl", "adviControl", "pathfinderControl")
  if (inherits(.ctl, .known)) {
    if (!inherits(.ctl, "stanControl")) {
      cli::cli_inform(paste0("converting ", class(.ctl)[1], " to ", .cls))
    }
    return(do.call(.fn, .ctl))
  }
  cli::cli_inform(paste0("invalid control for est=\"", est, "\", using default"))
  .fn()
}

#' Shared `nmObjGetControl` body for the sugar methods
#' @param x the fit (a length-1 list holding its env)
#' @param est sugar estimation name
#' @return the control found on the fit
#' @noRd
.stanSugarGetControl <- function(x, est) {
  .env <- x[[1]]
  .cls <- paste0(est, "Control")
  for (.nm in c(.cls, "control")) {
    if (exists(.nm, .env)) {
      .control <- get(.nm, .env)
      if (inherits(.control, .cls)) return(.control)
    }
  }
  # a fit run as est="<est>" with a plain stanControl() still has one
  for (.nm in c("stanControl", "control")) {
    if (exists(.nm, .env)) {
      .control <- get(.nm, .env)
      if (inherits(.control, "stanControl")) return(.control)
    }
  }
  stop("cannot find ", est, " related control object", call. = FALSE)
}

#' Shared `nlmixr2Est` body for the sugar methods
#'
#' Normalizes `env$control` (a sugar control, a stanControl, or nothing)
#' into the pair the engine and the fit each need.
#' @param env nlmixr2 estimation environment
#' @param est sugar estimation name
#' @return `env`, invisibly
#' @noRd
.stanSugarPrepEnv <- function(env, est) {
  .control <- env$control
  if (!inherits(.control, c(paste0(est, "Control"), "stanControl"))) {
    .control <- switch(est, nuts = nutsControl(), advi = adviControl(),
                       pathfinder = pathfinderControl()) # nocov
  }
  .control <- .stanSugarAsControl(.control, est)
  assign(paste0(est, "Control"), .control, envir = env)
  env$control <- .stanSugarAsStanControl(.control)
  env$stanEstName <- est
  invisible(env)
}

#' Control options for `est="nuts"` (Stan's NUTS sampler)
#'
#' `nutsControl()` is [stanControl()] with `algorithm = "NUTS"` forced --
#' the same options, named after the SAMPLING ALGORITHM instead of the
#' software tool, the way nlmixr2est's `foceControl()` is `foceiControl()`
#' with `interaction = FALSE` forced.  Because `class(nutsControl())` is
#' `"nutsControl"`, `nlmixr2(model, data, nutsControl())` infers
#' `est="nuts"` with no `est=` argument.
#'
#' `"NUTS"` is already [stanControl()]'s default, so this is mostly a
#' naming convenience; asking for a variational algorithm here is an error
#' rather than a silent override.
#'
#' @param ... [stanControl()] options
#' @param algorithm must be `"NUTS"` (present so a [stanControl()] can be
#'   converted to a `nutsControl`; a variational algorithm errors)
#' @return a `nutsControl` list
#' @seealso [stanControl()], [adviControl()], [pathfinderControl()]
#' @export
#' @author Matthew L Fidler
nutsControl <- function(..., algorithm = "NUTS") {
  .stanSugarControl("nuts", algorithm, list(...))
}

#' Control options for `est="advi"` (Stan's ADVI variational Bayes)
#'
#' `adviControl()` is [stanControl()] with an ADVI `algorithm` forced --
#' the same options, named after the INFERENCE ALGORITHM instead of the
#' software tool.  Because `class(adviControl())` is `"adviControl"`,
#' `nlmixr2(model, data, adviControl())` infers `est="advi"` with no
#' `est=` argument.
#'
#' ADVI is typically 10-100x faster than NUTS and is a fair first look,
#' but it is an APPROXIMATION -- it understates tails and correlations;
#' the Pareto-k diagnostic (`khat`) replaces Rhat/divergences and
#' `khat > 0.7` means NUTS should be used instead.
#'
#' @param ... [stanControl()] options
#' @param algorithm the ADVI variant: `"meanfield"` (default; independent
#'   Gaussians on the unconstrained scale) or `"fullrank"` (one
#'   full-covariance Gaussian).  `"NUTS"` is accepted only because it is
#'   [stanControl()]'s default and selects `"meanfield"`; `"pathfinder"`
#'   errors
#' @return an `adviControl` list
#' @seealso [stanControl()], [nutsControl()], [pathfinderControl()]
#' @export
#' @author Matthew L Fidler
adviControl <- function(..., algorithm = c("meanfield", "fullrank")) {
  if (length(algorithm) > 1L) algorithm <- match.arg(algorithm)
  .stanSugarControl("advi", algorithm, list(...))
}

#' Control options for `est="pathfinder"` (Stan's multi-path Pathfinder)
#'
#' `pathfinderControl()` is [stanControl()] with
#' `algorithm = "pathfinder"` forced -- the same options, named after the
#' INFERENCE ALGORITHM instead of the software tool.  Because
#' `class(pathfinderControl())` is `"pathfinderControl"`,
#' `nlmixr2(model, data, pathfinderControl())` infers `est="pathfinder"`
#' with no `est=` argument.
#'
#' Pathfinder needs the `StanEstimators` package installed (rstan does not
#' expose Pathfinder yet); see [stanControl()] for `pathfinderPaths` and
#' the `khat` diagnostic.
#'
#' @param ... [stanControl()] options
#' @param algorithm must be `"pathfinder"` (`"NUTS"` is accepted only
#'   because it is [stanControl()]'s default; an ADVI variant errors)
#' @return a `pathfinderControl` list
#' @seealso [stanControl()], [nutsControl()], [adviControl()]
#' @export
#' @author Matthew L Fidler
pathfinderControl <- function(..., algorithm = c("pathfinder", "NUTS")) {
  if (length(algorithm) > 1L) algorithm <- match.arg(algorithm)
  .stanSugarControl("pathfinder", algorithm, list(...))
}

#' @author Matthew L Fidler
#' @export
rxUiDeparse.nutsControl <- function(object, var) {
  .default <- nutsControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

#' @author Matthew L Fidler
#' @export
rxUiDeparse.adviControl <- function(object, var) {
  .default <- adviControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

#' @author Matthew L Fidler
#' @export
rxUiDeparse.pathfinderControl <- function(object, var) {
  .default <- pathfinderControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

#' @author Matthew L Fidler
#' @export
nmObjHandleControlObject.nutsControl <- function(control, env) {
  assign("nutsControl", control, envir = env)
}

#' @author Matthew L Fidler
#' @export
nmObjHandleControlObject.adviControl <- function(control, env) {
  assign("adviControl", control, envir = env)
}

#' @author Matthew L Fidler
#' @export
nmObjHandleControlObject.pathfinderControl <- function(control, env) {
  assign("pathfinderControl", control, envir = env)
}

#' @author Matthew L Fidler
#' @export
getValidNlmixrCtl.nuts <- function(control) {
  .stanSugarValidCtl(control, "nuts")
}

#' @author Matthew L Fidler
#' @export
getValidNlmixrCtl.advi <- function(control) {
  .stanSugarValidCtl(control, "advi")
}

#' @author Matthew L Fidler
#' @export
getValidNlmixrCtl.pathfinder <- function(control) {
  .stanSugarValidCtl(control, "pathfinder")
}

#' @author Matthew L Fidler
#' @export
nmObjGetControl.nuts <- function(x, ...) {
  .stanSugarGetControl(x, "nuts")
}

#' @author Matthew L Fidler
#' @export
nmObjGetControl.advi <- function(x, ...) {
  .stanSugarGetControl(x, "advi")
}

#' @author Matthew L Fidler
#' @export
nmObjGetControl.pathfinder <- function(x, ...) {
  .stanSugarGetControl(x, "pathfinder")
}

#' NUTS estimation through Stan (sugar for est="stan" + algorithm)
#'
#' `est="nuts"` is `est="stan"` with `stanControl(algorithm="NUTS")`
#' forced -- naming the SAMPLING ALGORITHM (Stan's No-U-Turn full-Bayes
#' HMC sampler) the way `est="advi"` and `est="pathfinder"` name theirs,
#' instead of naming the software tool (`est="stan"`).  `"NUTS"` is
#' already `stanControl()`'s default, so this is mostly a naming
#' convenience; asking for a variational algorithm under `est="nuts"` is
#' an error rather than a silent override.  [nutsControl()] is its
#' control, and `nlmixr2(model, data, nutsControl())` infers `est="nuts"`
#' from it.
#'
#' @param env nlmixr2 estimation environment
#' @param ... passed through
#' @return nlmixr2 fit
#' @keywords internal
#' @author Matthew L Fidler
#' @export
nlmixr2Est.nuts <- function(env, ...) {
  .stanSugarPrepEnv(env, "nuts")
  nlmixr2Est.stan(env, ...)
}
attr(nlmixr2Est.nuts, "nlmixr2Priors") <- "all"
attr(nlmixr2Est.nuts, "type") <- "Markov chain Monte Carlo"
attr(nlmixr2Est.nuts, "description") <-
  "Stan NUTS (Hamiltonian Monte Carlo) linked to the rxode2/nlmixr2est likelihood"
attr(nlmixr2Est.nuts, "covPresent") <- TRUE
# a declared non-normal random effect arrives here already expanded by
# nlmixr2est's pre-processing hook: the latent eta is a FIXED unit
# variance (which the generator already supports), the transform lives
# inside the linked nlmixr2_cond_all2() likelihood, and the copula
# correlation is an ordinary unbounded theta -- so nothing in the
# generated Stan program is special cased for it
attr(nlmixr2Est.nuts, "etaDist") <- TRUE
attr(nlmixr2Est.nuts, "mu") <- FALSE
attr(nlmixr2Est.nuts, "unbounded") <- FALSE
attr(nlmixr2Est.nuts, "iov") <- function(control) .stanHasIovSens()

#' ADVI estimation through Stan (sugar for est="stan" + algorithm)
#'
#' `est="advi"` is `est="stan"` with `stanControl(algorithm="meanfield")`
#' unless the control already names an ADVI variant (`"fullrank"` is
#' honored; asking for `"NUTS"`/`"pathfinder"` under `est="advi"` is an
#' error rather than a silent override).  [adviControl()] is its control,
#' and `nlmixr2(model, data, adviControl())` infers `est="advi"` from it.
#'
#' @param env nlmixr2 estimation environment
#' @param ... passed through
#' @return nlmixr2 fit
#' @keywords internal
#' @author Matthew L Fidler
#' @export
nlmixr2Est.advi <- function(env, ...) {
  .stanSugarPrepEnv(env, "advi")
  nlmixr2Est.stan(env, ...)
}
attr(nlmixr2Est.advi, "nlmixr2Priors") <- "all"
attr(nlmixr2Est.advi, "type") <- "Variational inference"
attr(nlmixr2Est.advi, "description") <-
  "Stan ADVI (variational Bayes) linked to the rxode2/nlmixr2est likelihood"
attr(nlmixr2Est.advi, "covPresent") <- TRUE
# a declared non-normal random effect arrives here already expanded by
# nlmixr2est's pre-processing hook: the latent eta is a FIXED unit
# variance (which the generator already supports), the transform lives
# inside the linked nlmixr2_cond_all2() likelihood, and the copula
# correlation is an ordinary unbounded theta -- so nothing in the
# generated Stan program is special cased for it
attr(nlmixr2Est.advi, "etaDist") <- TRUE
attr(nlmixr2Est.advi, "mu") <- FALSE
attr(nlmixr2Est.advi, "unbounded") <- FALSE
attr(nlmixr2Est.advi, "iov") <- function(control) .stanHasIovSens()

#' Pathfinder estimation through Stan (sugar for est="stan" + algorithm)
#'
#' `est="pathfinder"` is `est="stan"` with
#' `stanControl(algorithm="pathfinder")` forced.
#' [pathfinderControl()] is its control, and
#' `nlmixr2(model, data, pathfinderControl())` infers `est="pathfinder"`
#' from it.
#'
#' @inheritParams nlmixr2Est.advi
#' @return nlmixr2 fit
#' @keywords internal
#' @author Matthew L Fidler
#' @export
nlmixr2Est.pathfinder <- function(env, ...) {
  .stanSugarPrepEnv(env, "pathfinder")
  nlmixr2Est.stan(env, ...)
}
attr(nlmixr2Est.pathfinder, "nlmixr2Priors") <- "all"
attr(nlmixr2Est.pathfinder, "type") <- "Variational inference"
attr(nlmixr2Est.pathfinder, "description") <-
  "Stan Pathfinder linked to the rxode2/nlmixr2est likelihood"
attr(nlmixr2Est.pathfinder, "covPresent") <- TRUE
# a declared non-normal random effect arrives here already expanded by
# nlmixr2est's pre-processing hook: the latent eta is a FIXED unit
# variance (which the generator already supports), the transform lives
# inside the linked nlmixr2_cond_all2() likelihood, and the copula
# correlation is an ordinary unbounded theta -- so nothing in the
# generated Stan program is special cased for it
attr(nlmixr2Est.pathfinder, "etaDist") <- TRUE
attr(nlmixr2Est.pathfinder, "mu") <- FALSE
attr(nlmixr2Est.pathfinder, "unbounded") <- FALSE
attr(nlmixr2Est.pathfinder, "iov") <- function(control) .stanHasIovSens()
