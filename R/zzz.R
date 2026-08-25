#' @useDynLib nlmixr2bayes, .registration=TRUE
#' @importFrom stats setNames
#' @importFrom methods is
#' @importFrom nlmixr2est nlmixr2Est getValidNlmixrCtl nmObjGetControl nmObjHandleControlObject
#' @importFrom rxode2 rxUiDeparse
NULL

.rxsEnv <- new.env(parent = emptyenv())
.rxsEnv$handles <- list()
.rxsEnv$nextHandle <- 1L

.onLoad <- function(libname, pkgname) {
  # Install nlmixr2est's FOCEi conditional-likelihood entry points on EVERY
  # load, not just the first: a reloaded nlmixr2est hands back new addresses,
  # and cached ones would point into an unmapped DLL.
  .iniFoceiPtrs()
  # The nlm (tier-0, population-only) table is OPTIONAL: an older nlmixr2est
  # without nlmixr2/nlmixr2est#953 just leaves tier 0 unavailable.
  .iniNlmPtrs()
  # Install the rxode2 function-pointer table for the rxstan bridge.
  .Call(`_nlmixr2bayes_iniRxodePtrs`, rxode2::.rxode2ptrs())
  rxode2::.s3register("nlmixr2est::nmObjGet", "stanfit")
  rxode2::.s3register("nlmixr2est::nmObjGet", "posteriorSummary")
  # nlmixr2save is a Suggests, so its saveFitItem() method is registered here
  # rather than in NAMESPACE (see R/stanSave.R)
  rxode2::.s3register("nlmixr2save::saveFitItem", "stanfit")
  rxode2::.s3register("nlmixr2save::saveFitItem", "nlmixr2bayesPathfinder")
  invisible()
}

#' Report whether rxode2's C function-pointer table was installed
#'
#' The bridge cannot work without it, so this is the first thing to check if
#' solves fail in a way that looks like rxode2 is not there.
#'
#' @return named logical vector, one entry per probed rxode2 entry point
#' @author Lukas A. Widmer
#' @export
rxsProbeRxode2 <- function() {
  .Call(C_rxstanProbeRxode2)
  invisible()
}

#' Install (or refresh) the nlmixr2est nlm pointer table (tier 0, #953);
#' returns FALSE (tier 0 unavailable) on an older nlmixr2est
#' @noRd
.iniNlmPtrs <- function() {
  .fun <- tryCatch(getExportedValue("nlmixr2est", ".nlmixr2estNlmPtrs"),
                   error = function(e) NULL)
  if (is.null(.fun)) return(invisible(FALSE))
  .Call(`_nlmixr2bayes_iniNlmPtrs`, .fun())
  .v <- .Call(`_nlmixr2bayes_nlmApiVersion`)
  if (!identical(.v, 1L)) {
    stop("nlmixr2bayes was built against nlm C API version 1, but the loaded ",
         "nlmixr2est provides version ", .v,
         "; reinstall nlmixr2bayes against this nlmixr2est", call. = FALSE)
  }
  invisible(TRUE)
}

#' Does the loaded nlmixr2est support the combined eta+theta sensitivity
#' build (nlmixr2/nlmixr2est#958)?  Probed structurally: the feature added
#' the combSens argument to foceiLikLoad.
#' @noRd
.stanHasCombSens <- function() {
  "combSens" %in% names(formals(nlmixr2est::foceiLikLoad))
}

#' Is the tier-0 (population-only) nlm C API available?
#' @noRd
.stanHasNlmApi <- function() {
  identical(.Call(`_nlmixr2bayes_nlmApiVersion`), 1L)
}

#' Install (or refresh) the nlmixr2est FOCEi pointer table
#'
#' @return invisibly, the API version
#' @noRd
.iniFoceiPtrs <- function() {
  .fun <- tryCatch(getExportedValue("nlmixr2est", ".nlmixr2estFoceiPtrs"),
                   error = function(e) NULL)
  if (is.null(.fun)) {
    warning("this nlmixr2est does not provide the FOCEi conditional-likelihood ",
            "C API (nlmixr2est::.nlmixr2estFoceiPtrs); update nlmixr2est ",
            "(nlmixr2/nlmixr2est#937); Stan-based estimation will not work",
            call. = FALSE)
    return(invisible(FALSE))
  }
  .Call(`_nlmixr2bayes_iniFoceiPtrs`, .fun())
  .v <- .Call(`_nlmixr2bayes_apiVersion`)
  if (!identical(.v, 1L)) {
    stop("nlmixr2bayes was built against FOCEi C API version 1, but the ",
         "loaded nlmixr2est provides version ", .v,
         "; reinstall nlmixr2bayes against this nlmixr2est",
         call. = FALSE)
  }
  invisible(.v)
}

#' Set the subject-parallel thread count for the linked likelihood
#'
#' Stan itself runs single-threaded in nlmixr2bayes; the parallelism lives in
#' the rxode2/nlmixr2est subject loop inside each likelihood evaluation.  This
#' sets how many threads that loop uses (clamped to what the loaded solve pool
#' supports, and to 1 when the ODE method is not thread safe).
#'
#' @param cores number of threads (default [rxode2::getRxThreads()])
#' @return invisibly, the thread count set
#' @author Matthew L Fidler
#' @export
stanSetCores <- function(cores = rxode2::getRxThreads()) {
  checkmate::assertIntegerish(cores, lower = 1, len = 1, any.missing = FALSE)
  invisible(.Call(`_nlmixr2bayes_setCores`, as.integer(cores)))
}
