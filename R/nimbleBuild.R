# Build-environment probes for the NIMBLE backend, mirroring stanBuild.R's
# philosophy (issue #1, Spec 3): each probe costs milliseconds (no
# compilation), and turns a class of silent future breakage into a
# load-time refusal with an actionable message, rather than a cryptic
# downstream NIMBLE compiler/linker error.
#
# What could silently drift under us, and what each probe catches:
#   * a NIMBLE release renames/removes a BUGS-registered distribution
#     argument R/nimblePriors.R's .nimblePriorCatalog templates assume is
#     valid (e.g. dnorm's "sd") -- probed against
#     nimble:::distributionsInputList directly, not the R-level mirror
#     functions of the same names (which use different, non-BUGS
#     defaults again -- see the design note atop R/nimblePriors.R).
#   * nimble::nimbleExternalCall is renamed or dropped -- probed by
#     formals() existing with the expected argument names.

.nimbleBuildEnv <- new.env(parent = emptyenv())

#' BUGS-registered argument names R/nimblePriors.R's `.nimblePriorCatalog`
#' templates assume are valid for each distribution. Every template uses
#' these as NAMED arguments (never positional) specifically because
#' several of them default an unnamed second positional argument to
#' something other than what the template needs (dnorm's raw 2nd
#' positional is precision, not sd; similarly dlnorm/dlogis/dt) -- so the
#' named argument itself existing and meaning what it's assumed to mean is
#' the load-bearing fact this probes.
#' @noRd
.nimbleExpectedDistArgs <- list(
  dnorm  = c("mean", "sd"),
  dlnorm = c("meanlog", "sdlog"),
  dgamma = c("shape", "rate"),
  dbeta  = c("shape1", "shape2"),
  dexp   = "rate",
  dunif  = c("min", "max"),
  dt     = c("mu", "df", "sigma"),
  dlogis = c("location", "scale")
)

#' Check one distribution's BUGS-registered argument list contains every
#' name a `.nimblePriorCatalog` template assumes
#' @noRd
.nimbleCheckDistArgs <- function(distList, distName, expectedArgs) {
  .bugs <- vapply(distList, function(x) {
    .b <- x$BUGSdist
    if (is.null(.b) || length(.b) != 1L) NA_character_ else .b
  }, character(1))
  .w <- which(grepl(paste0("^", distName, "\\("), .bugs))
  if (length(.w) != 1L) return(FALSE)
  all(vapply(expectedArgs, function(a) grepl(a, .bugs[.w], fixed = TRUE),
             logical(1)))
}

#' Probe the NIMBLE build environment
#'
#' Inspects `nimble`'s own distribution registry and the `nimbleExternalCall`
#' formals (cached per session): no NIMBLE compilation happens here, so this
#' is cheap enough to call before every `nimbleLinkedSample()`.
#'
#' @param force re-run the probe even if cached
#' @return a list with `ok`, `distArgsOk` (named per distribution in
#'   [.nimblePriorCatalog]), `externalCallOk`, `nimbleVersion`; `ok=NA` when
#'   nimble is not installed
#' @export
#' @author Matthew L Fidler
nimbleBuildInfo <- function(force = FALSE) {
  if (!force && !is.null(.nimbleBuildEnv$info)) return(.nimbleBuildEnv$info)
  if (!requireNamespace("nimble", quietly = TRUE)) {
    .ret <- list(ok = NA, distArgsOk = NA, externalCallOk = NA,
                 nimbleVersion = NA_character_)
    .nimbleBuildEnv$info <- .ret
    return(.ret)
  }
  .dl <- tryCatch(get("distributionsInputList", envir = asNamespace("nimble")),
                  error = function(e) NULL) # nocov
  .distArgsOk <- if (is.null(.dl)) {
    stats::setNames(rep(NA, length(.nimbleExpectedDistArgs)),
                    names(.nimbleExpectedDistArgs)) # nocov
  } else {
    vapply(names(.nimbleExpectedDistArgs), function(.nm) {
      .nimbleCheckDistArgs(.dl, .nm, .nimbleExpectedDistArgs[[.nm]])
    }, logical(1))
  }
  .ecFormals <- tryCatch(names(formals(nimble::nimbleExternalCall)),
                         error = function(e) character(0)) # nocov
  .externalCallOk <- all(c("prototype", "returnType", "Cfun", "headerFile",
                          "oFile") %in% .ecFormals)
  .ret <- list(ok = isTRUE(all(.distArgsOk)) && isTRUE(.externalCallOk),
               distArgsOk = .distArgsOk,
               externalCallOk = .externalCallOk,
               nimbleVersion = as.character(utils::packageVersion("nimble")))
  .nimbleBuildEnv$info <- .ret
  .ret
}

#' Error unless the NIMBLE build environment probes pass
#' @noRd
.nimbleAssertBuildOk <- function() {
  .i <- nimbleBuildInfo()
  if (isTRUE(is.na(.i$ok))) {
    stop("nimble is not installed; install nimble to use nimbleLinkedSample()",
         call. = FALSE)
  }
  .badDist <- names(.i$distArgsOk)[!.i$distArgsOk]
  if (length(.badDist) > 0L) {
    stop("this nimble (", .i$nimbleVersion, ") no longer registers the ",
         "expected argument names for: ", paste(.badDist, collapse = ", "),
         "; nlmixr2bayes's nimble prior catalog (R/nimblePriors.R) needs ",
         "updating for this nimble version", call. = FALSE)
  }
  if (!isTRUE(.i$externalCallOk)) {
    stop("this nimble (", .i$nimbleVersion, ") no longer provides ",
         "nimbleExternalCall() with the expected arguments; nlmixr2bayes's ",
         "nimble shim linkage (R/nimbleGen.R) needs updating for this ",
         "nimble version", call. = FALSE)
  }
  invisible(TRUE)
}
