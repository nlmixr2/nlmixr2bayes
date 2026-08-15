#' @useDynLib nlmixr2stan, .registration=TRUE
#' @importFrom stats setNames
NULL

.onLoad <- function(libname, pkgname) {
  # Install nlmixr2est's FOCEi conditional-likelihood entry points on EVERY
  # load, not just the first: a reloaded nlmixr2est hands back new addresses,
  # and cached ones would point into an unmapped DLL.
  .iniFoceiPtrs()
}

#' Install (or refresh) the nlmixr2est FOCEi pointer table
#'
#' @return invisibly, the API version
#' @noRd
.iniFoceiPtrs <- function() {
  .fun <- tryCatch(getExportedValue("nlmixr2est", ".nlmixr2estFoceiPtrs"),
                   error = function(e) NULL)
  if (is.null(.fun)) {
    stop("this nlmixr2est does not provide the FOCEi conditional-likelihood ",
         "C API (nlmixr2est::.nlmixr2estFoceiPtrs); update nlmixr2est ",
         "(nlmixr2/nlmixr2est#937)",
         call. = FALSE)
  }
  .Call(`_nlmixr2stan_iniFoceiPtrs`, .fun())
  .v <- .Call(`_nlmixr2stan_apiVersion`)
  if (!identical(.v, 1L)) {
    stop("nlmixr2stan was built against FOCEi C API version 1, but the ",
         "loaded nlmixr2est provides version ", .v,
         "; reinstall nlmixr2stan against this nlmixr2est",
         call. = FALSE)
  }
  invisible(.v)
}

#' Set the subject-parallel thread count for the linked likelihood
#'
#' Stan itself runs single-threaded in nlmixr2stan; the parallelism lives in
#' the rxode2/nlmixr2est subject loop inside each likelihood evaluation.  This
#' sets how many threads that loop uses (clamped to what the loaded solve pool
#' supports, and to 1 when the ODE method is not thread safe).
#'
#' @param cores number of threads (default [rxode2::getRxThreads()])
#' @return invisibly, the thread count set
#' @export
stanSetCores <- function(cores = rxode2::getRxThreads()) {
  checkmate::assertIntegerish(cores, lower = 1, len = 1, any.missing = FALSE)
  invisible(.Call(`_nlmixr2stan_setCores`, as.integer(cores)))
}
