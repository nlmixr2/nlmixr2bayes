# The link lifecycle: load the FOCEi conditional-likelihood problem
# (nlmixr2est::foceiLikLoad) and take the package-level re-entrancy lock.
# op_focei is a process-wide global on the nlmixr2est side -- one loaded
# problem at a time, main R thread only -- so a second est="stan" run starting
# while one is active must error, not silently corrupt state.

.stanLinkEnv <- new.env(parent = emptyenv())
.stanLinkEnv$handle <- NULL
.stanLinkEnv$hash <- NULL

#' Load the linked likelihood for a model/data pair
#'
#' Wraps [nlmixr2est::foceiLikLoad()] with `scale="natural"` (so the theta
#' vector is directly comparable with `ui$iniDf$est`), records the handle and
#' a setup hash, and takes the package-level re-entrancy lock.
#'
#' @param ui rxode2 ui (or model function)
#' @param data estimation data
#' @param likelihood individual likelihood type; only `"focei"` and `"foce"`
#'   are self-consistent value/gradient pairs (`"focep"` is refused by the
#'   capability flags)
#' @param rxControl solving options
#' @param thetaSens also wire the theta-sensitivity model (`d/d(theta)` of the
#'   conditional at fixed eta)
#' @param literalFix inline fixed population parameters as model literals (the
#'   nlmixr2est literal-fix hook); keep this matched with the est method's
#'   `stanControl(literalFix=)` so the linked problem and the generated Stan
#'   program agree on the parameter vector
#' @param cores subject-parallel thread count for each likelihood evaluation
#' @return invisibly, the link handle (the `foceiLikLoad` handle plus
#'   `flags`, `setupHash`)
#' @export
#' @author Matthew L. Fidler
stanLinkSetup <- function(ui, data, likelihood = c("focei", "foce"),
                          rxControl = rxode2::rxControl(),
                          thetaSens = FALSE, literalFix = TRUE,
                          cores = rxode2::getRxThreads()) {
  likelihood <- match.arg(likelihood)
  if (!is.null(.stanLinkEnv$handle)) {
    stop("an nlmixr2stan likelihood is already linked; call stanLinkFree() first\n",
         "(op_focei is a process-wide global: one linked problem at a time)",
         call. = FALSE)
  }
  .h <- nlmixr2est::foceiLikLoad(ui, data, likelihood, rxControl = rxControl,
                                 scale = "natural", thetaSens = thetaSens,
                                 est = "stan", literalFix = literalFix)
  .d <- .Call(`_nlmixr2stan_dims`)
  if (.d[["status"]] != 0L) {
    nlmixr2est::foceiLikUnload()
    stop("the linked problem did not come up (dims status ", .d[["status"]], ")",
         call. = FALSE)
  }
  .flags <- .d[["flags"]]
  # refuse-at-load hazards (see inst/include/nlmixr2estFoceiPtr.h in
  # nlmixr2est): focep/fo have inconsistent or absent gradients; FD etas add
  # ~1e-6 gradient noise a gradient-based sampler sees as value/gradient
  # mismatch; mixtures re-index subjects
  if (bitwAnd(.flags, 0x02L) != 0L || bitwAnd(.flags, 0x04L) != 0L) {
    nlmixr2est::foceiLikUnload()
    stop("the loaded likelihood's value and gradient are not a consistent pair",
         " (focep/fo); use likelihood=\"focei\" or \"foce\"", call. = FALSE)
  }
  if (bitwAnd(.flags, 0x08L) != 0L) {
    nlmixr2est::foceiLikUnload()
    stop("some eta uses finite-difference event sensitivities; the gradient",
         " noise breaks a gradient-based sampler", call. = FALSE)
  }
  if (bitwAnd(.flags, 0x10L) != 0L) {
    nlmixr2est::foceiLikUnload()
    stop("mixture models are not supported by the Stan linkage", call. = FALSE)
  }
  stanSetCores(cores)
  .h$flags <- .flags
  .h$setupHash <- digest::digest(list(.h$thetaNames, .h$etaNames, .h$idLvl,
                                      .h$nid, .h$neta, likelihood))
  .stanLinkEnv$handle <- .h
  .stanLinkEnv$hash <- .h$setupHash
  invisible(.h)
}

#' Free the linked likelihood
#'
#' @return invisibly `TRUE` if a link was freed, `FALSE` if none was active
#' @export
#' @author Matthew L. Fidler
stanLinkFree <- function() {
  if (is.null(.stanLinkEnv$handle)) return(invisible(FALSE))
  if (isTRUE(.stanLinkEnv$handle$pop)) {
    nlmixr2est::.nlmFreeEnv()
  } else {
    nlmixr2est::foceiLikUnload()
  }
  .stanLinkEnv$handle <- NULL
  .stanLinkEnv$hash <- NULL
  invisible(TRUE)
}

#' The active link handle (or NULL)
#' @noRd
.stanLinkHandle <- function() .stanLinkEnv$handle

#' value + d/d(eta) of the conditional log-likelihood, through the exact C
#' path the compiled Stan model uses (R_RegisterCCallable ->
#' nlmixr2stan_cond_batch -> nlmixr2est's pointer table)
#'
#' @param eta nid x neta matrix
#' @return list(value=, grad=, nBad=)
#' @noRd
.condBatch <- function(eta) {
  .Call(`_nlmixr2stan_condBatch`, as.matrix(eta))
}

#' Write the natural-scale parameter vector into the linked problem
#' @noRd
.linkSetTheta <- function(theta) {
  .rc <- .Call(`_nlmixr2stan_setTheta`, as.double(theta))
  if (.rc != 0L) {
    stop("setting theta failed with status ", .rc, call. = FALSE)
  }
  invisible(.rc)
}

#' Install a caller-side Omega^-1 for gradient conditioning
#' @noRd
.linkSetOmegaInv <- function(omegaInv) {
  .rc <- .Call(`_nlmixr2stan_setOmegaInv`, as.matrix(omegaInv))
  if (.rc != 0L) {
    stop("setting Omega^-1 failed with status ", .rc, call. = FALSE)
  }
  invisible(.rc)
}

# ---- tier 0: population-only (no-eta) link --------------------------------

#' Load the tier-0 (population-only) likelihood for a model/data pair
#'
#' Wraps [nlmixr2est::nlmObjectiveSetup()] with `gradient=TRUE` and
#' `scale="natural"` (the nlm C API, nlmixr2/nlmixr2est#953), so repeated
#' evaluations return minus the log-likelihood and its analytic theta
#' gradient on the model's own scale.  Shares the package-level one-linked-
#' problem-at-a-time lock with [stanLinkSetup()].
#'
#' @param ui rxode2 ui (or model function) WITHOUT etas
#' @param data estimation data
#' @param rxControl solving options
#' @return invisibly, a handle (ntheta, nobs, flags, initPar, setupHash)
#' @export
#' @author Matthew L. Fidler
stanPopLinkSetup <- function(ui, data, rxControl = rxode2::rxControl()) {
  if (!is.null(.stanLinkEnv$handle)) {
    stop("an nlmixr2stan likelihood is already linked; call stanLinkFree() ",
         "first", call. = FALSE)
  }
  if (!.stanHasNlmApi()) {
    stop("this nlmixr2est does not provide the nlm population-likelihood C ",
         "API; update nlmixr2est (nlmixr2/nlmixr2est#953)", call. = FALSE)
  }
  .ini <- nlmixr2est::nlmObjectiveSetup(
    ui, data, control = nlmixr2est::nlmControl(rxControl = rxControl),
    gradient = TRUE, scale = "natural")
  .d <- .Call(`_nlmixr2stan_nlmDims`)
  if (.d[["status"]] != 0L) {
    nlmixr2est::.nlmFreeEnv() # nocov
    stop("the tier-0 problem did not come up (dims status ", .d[["status"]],
         ")", call. = FALSE) # nocov
  }
  if (bitwAnd(.d[["flags"]], 0x01L) == 0L ||
        bitwAnd(.d[["flags"]], 0x02L) == 0L) {
    nlmixr2est::.nlmFreeEnv() # nocov
    stop("the tier-0 load must carry the gradient model on the natural ",
         "scale", call. = FALSE) # nocov
  }
  if (bitwAnd(.d[["flags"]], 0x04L) != 0L) {
    nlmixr2est::.nlmFreeEnv()
    stop("some theta's sensitivity is finite-differenced; the gradient ",
         "noise breaks a gradient-based sampler", call. = FALSE)
  }
  .h <- list(pop = TRUE, ntheta = .d[["ntheta"]], nobs = .d[["nobs"]],
             flags = .d[["flags"]], initPar = as.numeric(.ini))
  .h$setupHash <- digest::digest(list("pop", .h$ntheta, .h$nobs, .h$initPar))
  .stanLinkEnv$handle <- .h
  .stanLinkEnv$hash <- .h$setupHash
  invisible(.h)
}

#' Free the tier-0 link (also freed by [stanLinkFree()])
#' @noRd
.stanPopLinkFree <- function() {
  if (is.null(.stanLinkEnv$handle)) return(invisible(FALSE))
  if (isTRUE(.stanLinkEnv$handle$pop)) {
    nlmixr2est::.nlmFreeEnv()
    .stanLinkEnv$handle <- NULL
    .stanLinkEnv$hash <- NULL
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

#' value + gradient of the tier-0 population log-likelihood (log p = -value
#' of the raw nlm convention; this mirror returns the RAW convention, the
#' injected header negates)
#' @noRd
.popEval <- function(theta) {
  .Call(`_nlmixr2stan_popEval`, as.double(theta))
}
