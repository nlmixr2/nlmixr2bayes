## nlmixr2save integration.
##
## The vignettes and README run their fits for real and cache them with
## nlmixr2save's `:=` operator, so the rendered documents carry live output
## without refitting on every build.  A cached est="stan" fit is only useful
## if it is small enough to commit, and by default it is not: rstan keeps the
## compiled model's shared object INSIDE the stanfit
## (`@stanmodel@dso@.CXXDSOMISC$dso_bin`), which serializes to tens of
## megabytes -- ~99% of a saved fit -- while the draws themselves are a
## fraction of that.

#' Drop the compiled shared object from a `stanfit`
#'
#' rstan carries the compiled model's DSO (its `.so` bytes, the Rcpp module and
#' the `CFunc`) inside every `stanfit`, in the `.CXXDSOMISC` environment of the
#' `dso` slot.  For a *saved* fit that payload is dead weight twice over: it is
#' the bulk of the serialized size, and a DSO restored into another session
#' cannot be re-entered anyway.
#'
#' Everything that reads a finished posterior keeps working without it --
#' `rstan::summary()`, `extract()`, `traceplot()`, `get_sampler_params()`,
#' `print()`, `posterior::as_draws_df()` -- because those read the `sim` slot.
#' Only calls that need to re-enter the compiled model (`log_prob()`,
#' `grad_log_prob()`, more `sampling()` from the same object, and so
#' [rxsCheckGradient()]) require the live model, and those belong to the
#' session that compiled it.
#'
#' The slot is replaced rather than emptied in place, so the caller's own
#' `stanfit` is left fully functional.
#'
#' @param x a `stanfit`; anything else is returned unchanged
#' @return `x` with an empty `.CXXDSOMISC` environment
#' @noRd
#' @author Matthew L. Fidler
.stanTrimStanfit <- function(x) {
  if (inherits(x, "nlmixr2bayesPathfinder")) {
    # a Pathfinder result carries the chains=0 stanfit it optimized against;
    # its constrained draws and generated quantities were already forced into
    # $cache before the linked likelihood was torn down, so the compiled model
    # is not needed to read the result back either
    x$fit0 <- .stanTrimStanfit(x$fit0)
    return(x)
  }
  if (!methods::is(x, "stanfit")) return(x)
  .dso <- x@stanmodel@dso
  .dso@.CXXDSOMISC <- new.env(parent = emptyenv())
  # dso_saved = TRUE promises the bytes are there to reload from; they are not
  .dso@dso_saved <- FALSE
  x@stanmodel@dso <- .dso
  x
}

#' nlmixr2save `saveFitItem()` method for the `stanfit` on an est="stan" fit
#'
#' Without it nlmixr2save falls back to `saveRDS()` on the whole `stanfit`
#' (with a "could not determine how to save" warning), which writes ~15 MB per
#' fit.  Trimming the DSO first takes the same fit to well under a megabyte,
#' which is what makes a committed `inst/cache` practical.
#'
#' Registered against `nlmixr2save::saveFitItem` in [.onLoad()] rather than
#' through `NAMESPACE`, because nlmixr2save is only a `Suggests`.
#'
#' @param item the `stanfit` being saved
#' @param name the fit-environment item name (`"stanfit"`)
#' @param file the saved fit's base path
#' @return `TRUE` if the item was written
#' @noRd
#' @author Matthew L. Fidler
saveFitItem.stanfit <- function(item, name, file) {
  .v <- try(saveRDS(.stanTrimStanfit(item), paste0(file, "-", name, ".rds")),
            silent = TRUE)
  !inherits(.v, "try-error")
}

#' @rdname saveFitItem.stanfit
#' @noRd
saveFitItem.nlmixr2bayesPathfinder <- saveFitItem.stanfit
