# IOV prior repair (issue #15).
#
# nlmixr2est's IOV preprocessing hook (`.uiApplyIov`) rewrites an occasion
# random effect `iov.x ~ v | OCC` into an estimated MAGNITUDE THETA plus
# per-occasion fixed unit-variance etas.  The new theta row is built by
# copying the FIRST theta row of the model as a template and overwriting a
# few fields -- `prior` is not one of them.  Two consequences, both silent:
#
#  1. the magnitude theta inherits theta #1's prior, so `iov_cl` samples
#     against a prior belonging to an unrelated parameter, and
#  2. the prior the user actually declared with `prior(iov.x) ~ ...` lives on
#     the ORIGINAL eta row, which the rewrite deletes -- so it never reaches
#     the generated program at all.
#
# The original ui is not reachable from an estimation method (the hooks
# replace it in place), so the declared prior is captured by a pre-processing
# hook of our own, registered ahead of `.uiApplyIov`, and re-attached to the
# magnitude theta here.  The capture is checked against the rewritten ui
# before it is used: a mismatch is an error, never a guess.
#
# NOT repaired here, because the template copy breaks them INSIDE the rewrite,
# before an estimation method is reached -- both fail loudly, upstream, and
# want the same one-line nlmixr2est fix (clear `prior` on the copied template
# rows):
#   * `iov.x ~ fix(v) | OCC` with a prior on theta #1 -- the fixed magnitude
#     theta inherits it and rxode2 refuses "a prior given for fixed
#     parameter(s): 'iov.x'";
#   * a prior on the model's FIRST eta -- `.eta1 <- .etas[1, ]` is the
#     template for the per-occasion etas too, so every `rx.<iov>.<occ>` row
#     inherits it and rxode2 refuses one of those instead.
#
# Scale: est="stan" never sets `iovXform`, so the hook always uses its "sd"
# default -- the magnitude theta IS the occasion standard deviation, and a
# declared `prior(iov.x)` is emitted on that SD scale.  Without a declared
# prior the theta gets the same default half-Cauchy an ordinary omega
# diagonal SD gets (`stanControl(diagOmegaSdPrior=)`), so the rewrite does
# not change which prior a random effect ends up with.

.stanIovEnv <- new.env(parent = emptyenv())
.stanIovEnv$priors <- character(0)

# `backTransform` of an IOV magnitude theta: one of nlmixr2est's
# nlmixr2iov{Sd,Var,Logsd,Logvar}{Cv,Sd} back-transformations, which nothing
# else in an iniDf carries
.stanIovBackTransform <- "^nlmixr2iov(Sd|Var|Logsd|Logvar)(Cv|Sd)$"

#' IOV magnitude theta names in a ui rewritten by nlmixr2est's IOV hook
#' @noRd
.stanIovMagnitude <- function(iniDf) {
  if (!is.data.frame(iniDf) || !("backTransform" %in% names(iniDf))) {
    return(character(0))
  }
  .w <- which(!is.na(iniDf$ntheta) & !is.na(iniDf$backTransform) &
                grepl(.stanIovBackTransform, iniDf$backTransform))
  iniDf$name[.w]
}

#' Occasion (IOV) eta rows of an ORIGINAL, not-yet-rewritten ui
#'
#' The same rows `.uiApplyIov()` keys on: a `condition` that is neither
#' `NA` nor `"id"`, and no residual-error model (which is how an error
#' parameter's endpoint `condition` is told apart from an occasion).
#' @noRd
.stanIovEtaRows <- function(iniDf) {
  if (!is.data.frame(iniDf) || !("condition" %in% names(iniDf))) {
    return(integer(0))
  }
  which(!is.na(iniDf$neta1) & iniDf$neta1 == iniDf$neta2 &
          !is.na(iniDf$condition) & iniDf$condition != "id" &
          is.na(iniDf$err))
}

#' Pre-processing hook: remember the priors declared on occasion etas
#'
#' Registered with [nlmixr2est::preProcessHooksAdd()] so it runs before
#' `.uiApplyIov` deletes those rows.  Pure capture -- it never overrides the
#' ui, est, data or control, and never errors.
#' @noRd
.stanCaptureIovPriors <- function(ui, est, data, control) {
  .stanIovEnv$priors <- character(0)
  .iniDf <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(.iniDf) || !is.data.frame(.iniDf) ||
        !("prior" %in% names(.iniDf))) {
    return(NULL)
  }
  .w <- .stanIovEtaRows(.iniDf)
  if (length(.w) == 0L) return(NULL)
  .stanIovEnv$priors <- stats::setNames(as.character(.iniDf$prior[.w]),
                                        .iniDf$name[.w])
  NULL
}

#' Register (or re-register) the capture hook with nlmixr2est
#'
#' Ordering is load-bearing: the hook must see the ui before `.uiApplyIov`
#' rewrites it.  nlmixr2est runs pre-processing hooks in `ls()` order, and
#' `".stanCaptureIovPriors"` sorts before `".uiApplyIov"`; the capture is
#' re-checked against the rewritten ui in [.stanIovRepairPriors()] anyway, so
#' a future ordering change becomes an error rather than a wrong prior.
#' @noRd
.stanIovRegisterHook <- function() {
  .add <- tryCatch(getExportedValue("nlmixr2est", "preProcessHooksAdd"),
                   error = function(e) NULL)
  .lst <- tryCatch(getExportedValue("nlmixr2est", "preProcessHooks"),
                   error = function(e) NULL)
  if (is.null(.add) || is.null(.lst)) return(invisible(FALSE)) # nocov
  if (".stanCaptureIovPriors" %in% .lst()) {
    .rm <- tryCatch(getExportedValue("nlmixr2est", "preProcessHooksRm"),
                    error = function(e) NULL)
    if (is.null(.rm)) return(invisible(FALSE)) # nocov
    .rm(".stanCaptureIovPriors")
  }
  .add(".stanCaptureIovPriors", .stanCaptureIovPriors)
  invisible(TRUE)
}

#' Default Stan prior statements for the IOV magnitude thetas without one
#'
#' A magnitude theta is the occasion SD, so it gets the same
#' `stanControl(diagOmegaSdPrior=)` half-Cauchy an ordinary omega diagonal
#' SD gets, scaled off the same `2.5 * sd0`.
#'
#' @param names IOV magnitude thetas that carry no declared prior
#' @param free the generator's free-theta data frame (`name`, `est`)
#' @return list(statements=, notes=), both possibly empty
#' @noRd
.stanIovDefaultPriors <- function(names, free, ctl) {
  .statements <- character(0)
  .notes <- character(0)
  for (.n in names) {
    .p <- sprintf(ctl$diagOmegaSdPrior,
                  .stanNum(2.5 * abs(free$est[match(.n, free$name)])))
    .statements <- c(.statements,
                     paste0("  ", .stanParName(.n), " ~ ", .p, ";"))
    .notes <- c(.notes,
                paste0("default prior ", .stanParName(.n), " ~ ", .p,
                       " on the inter-occasion magnitude (SD) of '", .n, "'"))
  }
  list(statements = .statements, notes = .notes)
}

#' Re-attach the declared IOV priors to a rewritten ui's prior rows
#'
#' @param iniDf the REWRITTEN ui's `iniDf`
#' @param pri the [rxode2::rxUiPriors()] data frame for that ui
#' @return `pri` with every IOV magnitude theta row replaced by the prior the
#'   user declared on the original occasion eta (dropped entirely when none
#'   was declared -- the generator supplies the default)
#' @noRd
.stanIovRepairPriors <- function(iniDf, pri) {
  .mag <- .stanIovMagnitude(iniDf)
  if (length(.mag) == 0L) return(pri)
  .cap <- .stanIovEnv$priors
  if (!setequal(.mag, names(.cap))) {
    stop("cannot recover the priors declared on the inter-occasion ",
         "parameter(s) ", paste0("'", .mag, "'", collapse = ", "),
         ": nlmixr2bayes did not see this model before nlmixr2est's IOV ",
         "preprocessing rewrote it (fit the model with nlmixr2(..., ",
         "est=\"stan\") rather than reusing a rewritten ui)", call. = FALSE)
  }
  # whatever the rewrite left on the magnitude theta belongs to the theta it
  # was copied from, never to this one
  pri <- pri[!(pri$name %in% .mag & is.na(pri$neta1)), , drop = FALSE]
  for (.n in .mag) {
    .p <- .cap[[.n]]
    if (is.na(.p)) next
    .lk <- .stanPriorLookup(.p)
    if (.lk$kind %in% c("matrix", "multivariate")) {
      stop("prior '", .lk$name, "' on the inter-occasion parameter '", .n,
           "' is a ", .lk$kind, " prior; est=\"stan\" samples an IOV ",
           "magnitude as a single standard-deviation theta, so its prior ",
           "must be univariate", call. = FALSE)
    }
    .r <- iniDf[which(iniDf$name == .n & !is.na(iniDf$ntheta)), , drop = FALSE]
    pri <- rbind(pri,
                 data.frame(name = .n, prior = .p, neta1 = NA_real_,
                            neta2 = NA_real_, lower = .r$lower[1],
                            upper = .r$upper[1], stringsAsFactors = FALSE))
  }
  pri
}

#' The model's priors, with nlmixr2est's IOV rewrite repaired
#'
#' Every nlmixr2bayes read of [rxode2::rxUiPriors()] goes through here, so
#' that an IOV model's declared priors are seen identically by the
#' `est="stan"` assertions and by the generator.
#' @noRd
.stanUiPriors <- function(ui) {
  .stanIovRepairPriors(ui$iniDf, rxode2::rxUiPriors(ui))
}
