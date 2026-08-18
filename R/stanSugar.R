# Sugar estimation methods over the est="stan" engine, the way "foce"/
# "focei" are members of nlmixr2est's focei family: est="nuts", est="advi",
# and est="pathfinder" run the SAME generated program and linked likelihood
# with the matching stanControl(algorithm=) forced, and the fit records
# the sugar name as its $est -- letting a user pick the ALGORITHM
# ("nuts"/"advi"/"pathfinder") instead of the software tool ("stan").

#' @author Matthew L Fidler
#' @export
getValidNlmixrCtl.nuts <- function(control) {
  getValidNlmixrCtl.stan(control)
}

#' @author Matthew L Fidler
#' @export
getValidNlmixrCtl.advi <- function(control) {
  getValidNlmixrCtl.stan(control)
}

#' @author Matthew L Fidler
#' @export
getValidNlmixrCtl.pathfinder <- function(control) {
  getValidNlmixrCtl.stan(control)
}

#' @author Matthew L Fidler
#' @export
nmObjGetControl.nuts <- function(x, ...) {
  nmObjGetControl.stan(x, ...)
}

#' @author Matthew L Fidler
#' @export
nmObjGetControl.advi <- function(x, ...) {
  nmObjGetControl.stan(x, ...)
}

#' @author Matthew L Fidler
#' @export
nmObjGetControl.pathfinder <- function(x, ...) {
  nmObjGetControl.stan(x, ...)
}

#' NUTS estimation through Stan (sugar for est="stan" + algorithm)
#'
#' `est="nuts"` is `est="stan"` with `stanControl(algorithm="NUTS")`
#' forced -- naming the SAMPLING ALGORITHM (Stan's No-U-Turn full-Bayes
#' HMC sampler) the way `est="advi"` and `est="pathfinder"` name theirs,
#' instead of naming the software tool (`est="stan"`).  `"NUTS"` is
#' already `stanControl()`'s default, so this is mostly a naming
#' convenience; asking for a variational algorithm under `est="nuts"` is
#' an error rather than a silent override.
#'
#' @param env nlmixr2 estimation environment
#' @param ... passed through
#' @return nlmixr2 fit
#' @keywords internal
#' @author Matthew L Fidler
#' @export
nlmixr2Est.nuts <- function(env, ...) {
  .control <- env$control
  if (!inherits(.control, "stanControl")) .control <- stanControl() # nocov
  if (!identical(.control$algorithm, "NUTS")) {
    stop("est=\"nuts\" runs Stan's NUTS sampler; use stanControl(algorithm=\"",
         .control$algorithm, "\") with est=\"stan\" instead", call. = FALSE)
  }
  env$control <- .control
  env$stanEstName <- "nuts"
  nlmixr2Est.stan(env, ...)
}
attr(nlmixr2Est.nuts, "nlmixr2Priors") <- "all"
attr(nlmixr2Est.nuts, "type") <- "Markov chain Monte Carlo"
attr(nlmixr2Est.nuts, "description") <-
  "Stan NUTS (Hamiltonian Monte Carlo) linked to the rxode2/nlmixr2est likelihood"
attr(nlmixr2Est.nuts, "covPresent") <- TRUE
attr(nlmixr2Est.nuts, "mu") <- FALSE
attr(nlmixr2Est.nuts, "unbounded") <- FALSE
attr(nlmixr2Est.nuts, "iov") <- function(control) .stanHasIovSens()

#' ADVI estimation through Stan (sugar for est="stan" + algorithm)
#'
#' `est="advi"` is `est="stan"` with `stanControl(algorithm="meanfield")`
#' unless the control already names an ADVI variant (`"fullrank"` is
#' honored; asking for `"NUTS"`/`"pathfinder"` under `est="advi"` is an
#' error rather than a silent override).
#'
#' @param env nlmixr2 estimation environment
#' @param ... passed through
#' @return nlmixr2 fit
#' @keywords internal
#' @author Matthew L Fidler
#' @export
nlmixr2Est.advi <- function(env, ...) {
  .control <- env$control
  if (!inherits(.control, "stanControl")) .control <- stanControl() # nocov
  if (identical(.control$algorithm, "NUTS")) {
    # the stanControl() default; est="advi" selects the ADVI default
    .control$algorithm <- "meanfield"
  }
  if (!(.control$algorithm %in% c("meanfield", "fullrank"))) {
    stop("est=\"advi\" runs ADVI; use stanControl(algorithm=\"",
         .control$algorithm, "\") with est=\"stan\" instead", call. = FALSE)
  }
  env$control <- .control
  env$stanEstName <- "advi"
  nlmixr2Est.stan(env, ...)
}
attr(nlmixr2Est.advi, "nlmixr2Priors") <- "all"
attr(nlmixr2Est.advi, "type") <- "Variational inference"
attr(nlmixr2Est.advi, "description") <-
  "Stan ADVI (variational Bayes) linked to the rxode2/nlmixr2est likelihood"
attr(nlmixr2Est.advi, "covPresent") <- TRUE
attr(nlmixr2Est.advi, "mu") <- FALSE
attr(nlmixr2Est.advi, "unbounded") <- FALSE
attr(nlmixr2Est.advi, "iov") <- function(control) .stanHasIovSens()

#' Pathfinder estimation through Stan (sugar for est="stan" + algorithm)
#'
#' `est="pathfinder"` is `est="stan"` with
#' `stanControl(algorithm="pathfinder")` forced.
#'
#' @inheritParams nlmixr2Est.advi
#' @return nlmixr2 fit
#' @keywords internal
#' @author Matthew L Fidler
#' @export
nlmixr2Est.pathfinder <- function(env, ...) {
  .control <- env$control
  if (!inherits(.control, "stanControl")) .control <- stanControl() # nocov
  if (!(.control$algorithm %in% c("NUTS", "pathfinder"))) {
    stop("est=\"pathfinder\" runs Pathfinder; use stanControl(algorithm=\"",
         .control$algorithm, "\") with est=\"stan\" instead", call. = FALSE)
  }
  .control$algorithm <- "pathfinder"
  env$control <- .control
  env$stanEstName <- "pathfinder"
  nlmixr2Est.stan(env, ...)
}
attr(nlmixr2Est.pathfinder, "nlmixr2Priors") <- "all"
attr(nlmixr2Est.pathfinder, "type") <- "Variational inference"
attr(nlmixr2Est.pathfinder, "description") <-
  "Stan Pathfinder linked to the rxode2/nlmixr2est likelihood"
attr(nlmixr2Est.pathfinder, "covPresent") <- TRUE
attr(nlmixr2Est.pathfinder, "mu") <- FALSE
attr(nlmixr2Est.pathfinder, "unbounded") <- FALSE
attr(nlmixr2Est.pathfinder, "iov") <- function(control) .stanHasIovSens()
