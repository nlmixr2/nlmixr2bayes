## Handle registry.
##
## A handle bundles everything a solve needs that is NOT a Stan parameter: the
## compiled rxode2 model (built with calcSens so it carries analytic
## sensitivities), the event table, fixed parameter values, solver controls, and
## which states at which rows Stan wants back.  Stan then passes only the
## opaque integer handle plus the differentiated parameter vector.

#' Register an rxode2 model as a Stan-callable solver
#'
#' Compiles `model` with analytic first-order sensitivities with respect to
#' `sens`, and returns an integer handle that a Stan program can pass to
#' `rx_solve()`.
#'
#' The returned outputs are ODE *states*.  rxode2 emits sensitivities only for
#' states (`rx__sens_<state>_BY_<param>__`), never for derived `lhs`
#' quantities, so any derived prediction (a concentration `center / v`, say)
#' must be formed in the Stan program, where Stan differentiates it itself.
#'
#' @param model rxode2 model source: a string, a `rxode2` expression, or an
#'   already-compiled model.  It is recompiled with `calcSens = sens`.
#' @param events event table (see [rxode2::et()]) or a data frame of events.
#' @param sens character vector of parameter names Stan will supply and
#'   differentiate with respect to, in the exact order Stan will pass them.
#' @param output character vector of state names to return, in order.
#' @param params named numeric vector of any remaining fixed parameters.
#' @param atol,rtol solver tolerances.  Defaults are deliberately tight: HMC
#'   needs the solver error well below the gradient signal or the trajectories
#'   stop being reversible.
#' @param method rxode2 integration method.
#' @param eventSens `"jump"` for analytic dosing-event sensitivities (modeled
#'   lag time, bioavailability, duration, rate), `"fd"` for the legacy
#'   finite-difference behaviour, or `NULL` to leave rxode2's default alone.
#' @param perSubject give every subject its own block of `sens` parameters,
#'   as a population model needs.  Stan then supplies `nsub * length(sens)`
#'   values, block-major (subject 1's parameters first).  Because subjects are
#'   independent the returned Jacobian stays `ny x length(sens)`: each output
#'   only depends on its own subject's block.
#' @param quiet suppress rxode2's per-solve diagnostics.  HMC proposes bad
#'   parameters routinely and each rejected proposal otherwise prints a block of
#'   `lsoda`/`intdy` warnings, which buries everything else over a real fit.
#'   Only that repeated chatter is silenced: a failed solve is still turned into
#'   an error that Stan sees and rejects on, and [rxsSolveStats()] keeps
#'   counting failures either way, so nothing is lost.
#' @param fast use Path A, which re-drives the solve rxode2 has already built
#'   instead of calling [rxode2::rxSolve()] again.  That skips roughly 2.3 ms of
#'   per-call R setup and teardown, which dominates small-to-medium problems.
#'   Any other `rxSolve()` in the session invalidates it; the bridge detects
#'   that and falls back to the slow path for one solve, which re-arms it.
#'   It is switched off automatically when rxode2 reports a non-zero
#'   `needSort` for the sensitivity-augmented model *and* the event table
#'   contains dose rows -- a modelled `alag()`, `f()`, `dur()` or `rate()`,
#'   whose event ordering Path A cannot refresh.  Both conditions are required
#'   because all four quantities act on doses, so an observation-only event
#'   table cannot be re-ordered; rxode2 sets the flag for a delay differential
#'   equation, which has no doses at all.  [print()] on the handle reports when
#'   the fast path was refused.
#' @param ... passed through to [rxode2::rxSolve()] on every solve.
#'
#' @return an integer handle, with attributes describing the layout
#' @author Lukas A. Widmer
#' @export
rxsRegister <- function(model, events, sens, output,
                        params = numeric(0),
                        atol = 1e-8, rtol = 1e-8,
                        method = "liblsoda",
                        eventSens = NULL,
                        perSubject = FALSE,
                        quiet = TRUE,
                        fast = TRUE,
                        ...) {
  stopifnot(is.character(sens), length(sens) > 0, !anyDuplicated(sens))
  stopifnot(is.character(output), length(output) > 0)

  args <- list(model, calcSens = sens)
  if (!is.null(eventSens)) args$eventSens <- eventSens
  m <- do.call(rxode2::rxode2, args)

  states <- rxode2::rxState(m)
  missingOut <- setdiff(output, states)
  if (length(missingOut)) {
    stop("rxsRegister(): output state(s) not in the model: ",
         paste(missingOut, collapse = ", "),
         "\n  available states: ", paste(states, collapse = ", "),
         call. = FALSE)
  }

  ## Column names rxode2 will produce, laid out state-major so that the
  ## flattened Jacobian is contiguous per parameter.
  sensCols <- outer(output, sens,
                    function(s, p) paste0("rx__sens_", s, "_BY_", p, "__"))
  dim(sensCols) <- c(length(output), length(sens))

  h <- list(model = m, events = events, sens = sens, output = output,
            params = params, atol = atol, rtol = rtol, method = method,
            sensCols = sensCols, perSubject = isTRUE(perSubject),
            quiet = isTRUE(quiet), fast = isTRUE(fast), dots = list(...))

  h$ids <- .rxsEventIds(events)
  h$nsub <- length(h$ids)
  if (h$perSubject && h$nsub < 1L) {
    stop("rxsRegister(): perSubject = TRUE needs an `id` column in `events`",
         call. = FALSE)
  }
  nBlock <- length(sens)
  nBlocks <- if (h$perSubject) h$nsub else 1L

  ## One probe solve: validates the model/events/params combination and fixes
  ## the number of returned rows before Stan ever calls in.
  probe <- .rxsRawSolve(h, rep(0.1, nBlock * nBlocks))
  missingCols <- setdiff(c(output, as.vector(sensCols)), names(probe))
  if (length(missingCols)) {
    stop("rxsRegister(): rxode2 did not return expected column(s): ",
         paste(missingCols, collapse = ", "), call. = FALSE)
  }

  h$nobs <- nrow(probe)
  ny <- h$nobs * length(output)

  ## Which parameter block each output row belongs to.  rxSolve stacks rows
  ## subject-major and the outputs repeat that block per requested state.
  ## With a single subject rxSolve drops the `id` column altogether, so fall
  ## back to the one block we know every row belongs to.
  rowBlock <- if (h$perSubject) {
    if (is.null(probe$id)) rep(0L, nrow(probe)) else match(probe$id, h$ids) - 1L
  } else {
    rep(0L, h$nobs)
  }
  if (length(rowBlock) != h$nobs || anyNA(rowBlock)) {
    stop("rxsRegister(): could not map solved rows back to subject ids",
         call. = FALSE)
  }
  outBlock <- rep(rowBlock, times = length(output))

  ## Index maps Path A needs, resolved once against the compiled model.
  states <- rxode2::rxState(m)
  h$sensIdx <- match(sens, rxode2::rxModelVars(m)$params) - 1L
  h$outIdx <- match(output, states) - 1L
  h$sensState <- as.integer(t(matrix(match(as.vector(sensCols), states) - 1L,
                                     nrow = length(output))))
  h$blockOf <- if (h$perSubject) seq_len(h$nsub) - 1L else rep(0L, h$nsub)
  if (anyNA(c(h$sensIdx, h$outIdx, h$sensState))) {
    h$fast <- FALSE
  }

  ## A modelled alag/f/dur/rate changes where a dose lands in the event order.
  ## rxSolve re-sorts events during setup; Path A only re-drives par_solve, so
  ## its cached ordering goes stale as soon as the parameter moves and the
  ## results are silently wrong near the shifted dose.  needSort is a bitmask:
  ## 1 = f, 2 = alag, 4 = dur, 8 = rate.
  ##
  ## All four act on DOSES, so with no dose rows there is nothing to re-order
  ## and the flag is a false alarm.  It fires that way more often than one
  ## would think: rxode2 reports needSort = 3 for a delay differential
  ## equation once it is recompiled with sensitivities, which is what the
  ## bridge actually solves, so every DDE was giving up the fast path for
  ## nothing.
  needSort <- as.integer(rxode2::rxModelVars(m)$needSort)[1]
  fastOff <- NULL
  if (isTRUE(h$fast) && !is.na(needSort) && needSort != 0L &&
      .rxsHasDoses(events)) {
    h$fast <- FALSE
    fastOff <- paste0("model has modelled dosing quantities (needSort = ",
                      needSort, "), whose event ordering Path A cannot refresh")
  }

  handle <- .rxsEnv$nextHandle
  .rxsEnv$nextHandle <- handle + 1L
  .rxsEnv$handles[[as.character(handle)]] <- h
  .Call(C_rxstanSetDims, handle, ny, nBlock, nBlocks, as.integer(outBlock))
  .rxsArmFast(handle, h)

  structure(handle,
            class = "rxsHandle",
            ny = ny, np = nBlock * nBlocks, nBlock = nBlock, nBlocks = nBlocks,
            nobs = h$nobs, nsub = h$nsub, perSubject = h$perSubject,
            fastDisabled = fastOff,
            sens = sens, output = output)
}

## Does the event table contain anything other than observations?  Answering
## "yes" when we cannot tell keeps the needSort refusal conservative: a missed
## dose would make Path A silently wrong, a missed speedup only costs time.
.rxsHasDoses <- function(events) {
  d <- try(as.data.frame(events), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d$evid)) return(TRUE)
  any(as.integer(d$evid) != 0L, na.rm = TRUE)
}

## Subject ids in the order rxode2 will solve and return them.
.rxsEventIds <- function(events) {
  d <- as.data.frame(events)
  if (!"id" %in% names(d)) return(1L)
  unique(d$id)
}

## Hands the live rx_solve structure to the C fast path.  Must be called while
## the solve that just ran is still the one rxode2 holds.
.rxsArmFast <- function(handle, h) {
  if (!isTRUE(h$fast)) return(invisible(FALSE))
  ok <- .Call(C_rxstanFastSetup, as.integer(handle),
              as.integer(h$sensIdx), as.integer(h$outIdx),
              as.integer(h$sensState), as.integer(h$nobs),
              as.integer(h$blockOf), isTRUE(h$quiet))
  invisible(isTRUE(ok))
}

#' How often a handle has solved, and how often that failed
#'
#' Counted regardless of `quiet`, so a silenced run still reports whether the
#' solver was struggling.
#'
#' @param handle handle returned by [rxsRegister()]
#' @return named numeric vector with `solves` and `failures`
#' @author Lukas A. Widmer
#' @export
rxsSolveStats <- function(handle) {
  .Call(C_rxstanStats, as.integer(unclass(handle)))
}

#' Reset a handle's solve counters
#' @param handle handle returned by [rxsRegister()]
#' @author Lukas A. Widmer
#' @export
rxsResetStats <- function(handle) {
  .Call(C_rxstanResetStats, as.integer(unclass(handle)))
  invisible(NULL)
}

#' Force the next solve back through the slow path
#'
#' Call this if something else in the session has solved an rxode2 model and you
#' want the bridge to rebuild its cached solve immediately rather than on the
#' next detection.
#' @author Lukas A. Widmer
#' @export
rxsInvalidateFast <- function() {
  .Call(C_rxstanFastInvalidate)
  invisible(NULL)
}

#' Is Path A currently armed for this handle?
#' @param handle handle returned by [rxsRegister()]
#' @author Lukas A. Widmer
#' @export
rxsFastAvailable <- function(handle) {
  .Call(C_rxstanFastAvailable, as.integer(unclass(handle)))
}

#' @author Lukas A. Widmer
#' @export
print.rxsHandle <- function(x, ...) {
  cat("<rxstan handle ", unclass(x), ">\n", sep = "")
  cat("  parameters (", attr(x, "np"), "): ",
      paste(attr(x, "sens"), collapse = ", "),
      if (attr(x, "perSubject")) paste0(" x ", attr(x, "nBlocks"), " subjects") else "",
      "\n", sep = "")
  cat("  states     (", length(attr(x, "output")), "): ",
      paste(attr(x, "output"), collapse = ", "), "\n", sep = "")
  cat("  rows       : ", attr(x, "nobs"), " obs -> ", attr(x, "ny"),
      " outputs\n", sep = "")
  if (!is.null(attr(x, "fastDisabled"))) {
    cat("  fast path  : off (", attr(x, "fastDisabled"), ")\n", sep = "")
  }
  invisible(x)
}

#' Release a registered handle
#' @param handle handle returned by [rxsRegister()]
#' @author Lukas A. Widmer
#' @export
rxsRelease <- function(handle) {
  key <- as.character(unclass(handle))
  .rxsEnv$handles[[key]] <- NULL
  .Call(C_rxstanClearDims, as.integer(unclass(handle)))
  invisible(NULL)
}

#' Release every registered handle
#' @author Lukas A. Widmer
#' @export
rxsReleaseAll <- function() {
  for (key in names(.rxsEnv$handles)) rxsRelease(as.integer(key))
  invisible(NULL)
}

#' List registered handles
#' @author Lukas A. Widmer
#' @export
rxsHandles <- function() as.integer(names(.rxsEnv$handles))

.rxsRawSolve <- function(h, p) {
  pars <- if (h$perSubject) {
    ## rxode2 takes one row per subject; Stan hands us the blocks end to end.
    pm <- matrix(as.numeric(p), nrow = h$nsub, ncol = length(h$sens),
                 byrow = TRUE, dimnames = list(NULL, h$sens))
    if (length(h$params)) {
      fixed <- matrix(rep(h$params, each = h$nsub), nrow = h$nsub,
                      dimnames = list(NULL, names(h$params)))
      pm <- cbind(pm, fixed)
    }
    pm
  } else {
    x <- h$params
    x[h$sens] <- as.numeric(p)
    x
  }
  call <- c(list(object = h$model, params = pars, events = h$events,
                 returnType = "data.frame", atol = h$atol, rtol = h$rtol,
                 method = h$method, cores = 1L),
            h$dots)
  if (!isTRUE(h$quiet)) return(do.call(rxode2::rxSolve, call))

  .Call(C_rxstanSetSilent, TRUE)
  on.exit(.Call(C_rxstanSetSilent, FALSE), add = TRUE)
  suppressWarnings(do.call(rxode2::rxSolve, call))
}

## Called from C (rxs_solve_sens) on every Stan gradient evaluation.  Must
## return a plain ny x (np + 1) double matrix: column 1 is the value, the
## remaining np columns are dy/dp.  Anything else is reported as an error
## rather than silently reinterpreted.
.rxsSolveOne <- function(handle, p) {
  h <- .rxsEnv$handles[[as.character(handle)]]
  if (is.null(h)) stop("unknown rxstan handle ", handle)

  s <- .rxsRawSolve(h, p)

  ny <- h$nobs * length(h$output)
  if (nrow(s) != h$nobs) {
    stop("rxode2 returned ", nrow(s), " rows, expected ", h$nobs)
  }

  out <- matrix(0, nrow = ny, ncol = length(h$sens) + 1L)
  out[, 1L] <- unlist(s[, h$output, drop = FALSE], use.names = FALSE)
  for (j in seq_along(h$sens)) {
    out[, j + 1L] <- unlist(s[, h$sensCols[, j], drop = FALSE],
                            use.names = FALSE)
  }

  ## The solve just built is the one rxode2 now holds, so this is the moment to
  ## re-arm Path A after whatever invalidated it.
  .rxsArmFast(handle, h)
  out
}

#' Solve a registered handle from R
#'
#' Goes through the same C entry point Stan uses, so the test suite covers the
#' real ABI without compiling a Stan model.
#'
#' @param handle handle returned by [rxsRegister()]
#' @param p numeric vector of parameter values, in `sens` order
#' @param slow force the slow [rxode2::rxSolve()] route, bypassing Path A
#' @return an `ny x (np + 1)` matrix: value, then `dy/dp` per parameter
#' @author Lukas A. Widmer
#' @export
rxsSolve <- function(handle, p, slow = FALSE) {
  fn <- if (slow) C_rxstanSolveSlow else C_rxstanSolve
  .Call(fn, as.integer(unclass(handle)), as.double(p))
}
