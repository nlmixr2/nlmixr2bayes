#' Control options for the `est="stan"` estimation method
#'
#' The class name `stanControl` makes `nlmixr2(model, data, stanControl(...))`
#' infer `est="stan"`.
#'
#' Two things are pinned and deliberately not user-settable: `cores=1` for the
#' Stan chains (the linked likelihood lives in process-wide state on the
#' nlmixr2est side, so chains run sequentially in one process; the
#' parallelism is over subjects inside each likelihood evaluation, see
#' `likCores`), and the ODE-retry machinery (`stickyRecalcN=1`,
#' `maxOdeRecalc=0`) -- a sampler tolerates a poor gradient but is wrong when
#' the density is not a pure function of state.
#'
#' @param chains number of MCMC chains (run sequentially)
#' @param iter total iterations per chain (Stan convention, includes warmup)
#' @param warmup warmup iterations per chain
#' @param thin thinning interval
#' @param seed random seed (fixed by default: reproducibility is a gate)
#' @param init `"ini"` (default) starts every chain from the model's `ini()`
#'   estimates with `initJitterSd` jitter for chains 2+; `"random"` uses
#'   Stan's default; or a list/function passed through to
#'   [rstan::sampling()]
#' @param initJitterSd unconstrained-scale jitter for chains 2+ under
#'   `init="ini"` (chain 1 starts exactly at the `ini()` values)
#' @param adapt_delta target acceptance rate (above Stan's 0.8 default:
#'   PK posteriors are stiff)
#' @param max_treedepth NUTS maximum tree depth
#' @param likelihood the individual likelihood the linked problem evaluates;
#'   `"focei"` (interaction) or `"foce"` -- the two self-consistent
#'   value/gradient pairs
#' @param likCores subject-parallel thread count inside each likelihood
#'   evaluation
#' @param diagOmegaSdPrior character template for the default half-Cauchy
#'   scale on an omega SD without its own prior; `%s` is replaced by the
#'   initial SD estimate (see Details in the vignette)
#' @param lkjEta shape for the default `lkj_corr_cholesky()` prior on a
#'   correlation block without its own prior
#' @param point posterior summary used as the point estimate on the fit
#'   (`"mean"`)
#' @param ofv objective row made current on the fit: `"focei"` evaluates the
#'   FOCEi objective at the posterior point estimate (comparable across
#'   nlmixr2 methods; a plug-in criterion, not a Bayesian one), `"none"`
#'   leaves the Bayesian fit without one
#' @param rhatMax,essBulkMin,essTailMin,maxDivergent diagnostic gates
#' @param onDiagnostic what a failed gate does: `"warn"` (default),
#'   `"error"`, `"message"`, `"none"`
#' @param run when `FALSE`, return the generated Stan program + data instead
#'   of compiling or sampling (no Stan toolchain needed)
#' @param cache cache compiled models by content hash
#' @param cacheDir cache directory (`NULL` = `tools::R_user_dir()`)
#' @param stanFile also write the generated `.stan` program here
#' @param verbose show compiler/sampler output
#' @param sigdig ODE solver accuracy (`atol=rtol=10^-sigdig`); tight by
#'   default so the density is smooth at sampler scale
#' @param rxControl solving options ([rxode2::rxControl()])
#' @param addProp,sumProd,optExpression,literalFix,calcTables,compress,ci,sigdigTable
#'   standard nlmixr2 passthroughs for the output machinery
#' @param ... only `genRxControl`
#' @return a `stanControl` list
#' @export
#' @author Matthew L. Fidler
stanControl <- function(chains = 4L, iter = 2000L, warmup = floor(iter / 2),
                        thin = 1L, seed = 42L,
                        init = c("ini", "random"), initJitterSd = 0.1,
                        adapt_delta = 0.9, max_treedepth = 12L, # nolint: object_name_linter.
                        likelihood = c("focei", "foce"),
                        likCores = rxode2::getRxThreads(),
                        diagOmegaSdPrior = "cauchy(0, %s)",
                        lkjEta = 2,
                        point = c("mean", "median"),
                        ofv = c("focei", "none"),
                        rhatMax = 1.01, essBulkMin = 100, essTailMin = 100,
                        maxDivergent = 0L,
                        onDiagnostic = c("warn", "error", "message", "none"),
                        run = TRUE, cache = TRUE, cacheDir = NULL,
                        stanFile = NULL, verbose = FALSE,
                        sigdig = 8,
                        rxControl = NULL,
                        addProp = c("combined2", "combined1"),
                        sumProd = FALSE, optExpression = TRUE,
                        literalFix = TRUE,
                        calcTables = TRUE, compress = TRUE, ci = 0.95,
                        sigdigTable = NULL, ...) {
  checkmate::assertIntegerish(chains, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(iter, lower = 2, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(warmup, lower = 1, upper = iter - 1, len = 1)
  checkmate::assertIntegerish(thin, lower = 1, len = 1)
  checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(initJitterSd, lower = 0, len = 1)
  checkmate::assertNumeric(adapt_delta, lower = 0, upper = 1, len = 1)
  checkmate::assertIntegerish(max_treedepth, lower = 1, upper = 20, len = 1)
  checkmate::assertIntegerish(likCores, lower = 1, len = 1)
  checkmate::assertCharacter(diagOmegaSdPrior, len = 1, pattern = "%s")
  checkmate::assertNumeric(lkjEta, lower = 0, len = 1)
  checkmate::assertNumeric(rhatMax, lower = 1, len = 1)
  checkmate::assertNumeric(essBulkMin, lower = 0, len = 1)
  checkmate::assertNumeric(essTailMin, lower = 0, len = 1)
  checkmate::assertIntegerish(maxDivergent, lower = 0, len = 1)
  checkmate::assertLogical(run, len = 1, any.missing = FALSE)
  checkmate::assertLogical(cache, len = 1, any.missing = FALSE)
  checkmate::assertLogical(verbose, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(sigdig, lower = 3, len = 1)
  checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)
  checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)
  checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)
  checkmate::assertLogical(compress, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(ci, lower = 0, upper = 1, len = 1)
  if (!is.null(stanFile)) checkmate::assertCharacter(stanFile, len = 1)
  if (!is.null(cacheDir)) checkmate::assertCharacter(cacheDir, len = 1)
  if (is.character(init)) init <- match.arg(init)
  likelihood <- match.arg(likelihood)
  point <- match.arg(point)
  ofv <- match.arg(ofv)
  onDiagnostic <- match.arg(onDiagnostic)
  if (checkmate::testIntegerish(addProp, lower = 1, upper = 2, len = 1)) {
    addProp <- c("combined1", "combined2")[addProp]
  } else {
    addProp <- match.arg(addProp)
  }
  .xtra <- list(...)
  .bad <- setdiff(names(.xtra), "genRxControl")
  if (length(.bad) > 0) {
    stop("unused argument: ", paste0("'", .bad, "'", collapse = ", "),
         call. = FALSE)
  }
  if (!is.null(.xtra$genRxControl)) {
    genRxControl <- .xtra$genRxControl
  } else {
    genRxControl <- FALSE
    .tol <- 10^(-sigdig)
    if (is.null(rxControl)) {
      # pinned for correctness: the retry machinery makes the density depend
      # on evaluation history (the batch entries also reset it per call)
      rxControl <- rxode2::rxControl(rtol = .tol, atol = .tol,
                                     ssRtol = .tol, ssAtol = .tol,
                                     maxsteps = 100000L)
      genRxControl <- TRUE
    } else if (is.list(rxControl) && !inherits(rxControl, "rxControl")) {
      rxControl <- do.call(rxode2::rxControl, rxControl)
    }
    if (!inherits(rxControl, "rxControl")) {
      stop("'rxControl' needs to be ode solving options from rxode2::rxControl()",
           call. = FALSE)
    }
  }
  .ret <- list(chains = as.integer(chains), iter = as.integer(iter),
               warmup = as.integer(warmup), thin = as.integer(thin),
               seed = as.integer(seed), init = init,
               initJitterSd = initJitterSd,
               adapt_delta = adapt_delta,
               max_treedepth = as.integer(max_treedepth),
               likelihood = likelihood, likCores = as.integer(likCores),
               diagOmegaSdPrior = diagOmegaSdPrior, lkjEta = lkjEta,
               point = point, ofv = ofv,
               rhatMax = rhatMax, essBulkMin = essBulkMin,
               essTailMin = essTailMin, maxDivergent = as.integer(maxDivergent),
               onDiagnostic = onDiagnostic,
               run = run, cache = cache, cacheDir = cacheDir,
               stanFile = stanFile, verbose = verbose,
               sigdig = sigdig, rxControl = rxControl,
               genRxControl = genRxControl,
               addProp = addProp, sumProd = sumProd,
               optExpression = optExpression, literalFix = literalFix,
               calcTables = calcTables, compress = compress, ci = ci,
               sigdigTable = sigdigTable)
  class(.ret) <- "stanControl"
  .ret
}

#' @export
rxUiDeparse.stanControl <- function(object, var) {
  .default <- stanControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

#' Validate a stanControl for the nlmixr2 dispatch (mandatory: the default
#' errors before `nlmixr2Est` dispatch)
#'
#' @param control a length-1 list holding the control (the
#'   `getValidNlmixrControl` convention)
#' @return a valid `stanControl`
#' @export
getValidNlmixrCtl.stan <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- stanControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    .ctl <- do.call("stanControl", .ctl)
  }
  if (!inherits(.ctl, "stanControl")) {
    cli::cli_inform(paste0("invalid control for est=\"", .cls,
                           "\", using default"))
    .ctl <- stanControl()
  } else {
    .ctl <- do.call(stanControl, .ctl)
  }
  .ctl
}

#' @export
nmObjHandleControlObject.stanControl <- function(control, env) {
  assign("stanControl", control, envir = env)
}

#' @export
nmObjGetControl.stan <- function(x, ...) {
  .env <- x[[1]]
  if (exists("stanControl", .env)) {
    .control <- get("stanControl", .env)
    if (inherits(.control, "stanControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "stanControl")) return(.control)
  }
  stop("cannot find stan related control object", call. = FALSE)
}

#' Zero-iteration foceiControl that drives the output/table machinery
#' @noRd
.stanControlToFoceiControl <- function(env, assign = TRUE) {
  .c <- env$stanControl
  .f <- nlmixr2est::foceiControl(
                                 rxControl = .c$rxControl,
                                 maxOuterIterations = 0L, maxInnerIterations = 0L,
                                 covMethod = 0L, interaction = 1L, scaleTo = 0,
                                 sumProd = .c$sumProd, optExpression = .c$optExpression,
                                 literalFix = .c$literalFix,
                                 calcTables = .c$calcTables, addProp = .c$addProp,
                                 compress = .c$compress, ci = .c$ci, sigdigTable = .c$sigdigTable,
                                 etaMat = env$etaMat)
  if (assign) env$control <- .f
  invisible(.f)
}
