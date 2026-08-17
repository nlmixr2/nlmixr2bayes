# The NIMBLE generator (phase 1, minimal scope): a custom distribution that
# wraps nlmixr2bayes_cond_batch_theta (through the compiled shim, see
# R/nimbleShim.R) as a NIMBLE potential, plus a BUGS-code builder and a
# driver that reuses stanLinkSetup()/setThetaBase()/setMuRef() verbatim --
# those are already engine-agnostic (R/stanLink.R, R/stanEst.R).
#
# Scope, deliberately narrow for phase 1 (refuse rather than mis-state the
# model, matching the est="stan" side's philosophy in R/stanGen.R):
#   * diagonal omega only (no correlated blocks, no fixed-variance etas)
#   * no finite mixtures, no mu-referenced covariates, no TBS Jacobian
#   * NIMBLE's default derivative-free MCMC only -- gradEta/gradTheta from
#     the linked likelihood are computed (the C entry always fills them) but
#     discarded; NIMBLE's nimbleExternalCall has no documented AD hook, so
#     HMC/NUTS parity with the Stan backend is not attempted here.
# The eta/theta PRIORS are NIMBLE-native normal priors (theta centered at
# ui$iniDf$est on the natural scale; eta ~ dnorm(0, sd=sqrt(omega))) rather
# than a translation of the Stan-syntax statements from stanPriors() --
# reusing those would require re-parsing Stan syntax under different
# parameterization conventions (e.g. NIMBLE's dnorm second argument
# defaults to precision, not sd). Full ini({})-driven prior parity is future
# work, not silently approximated here.

.nimbleGenEnv <- new.env(parent = emptyenv())

#' Register (once per session) the external call + custom distribution that
#' every nimbleLinkedSample() model uses.
#' @noRd
.nimbleCondLikSetup <- function() {
  # Cheap (no compilation) -- run on EVERY call, not just the first, so a
  # cache hit can't bypass it: reported as running "before every
  # nimbleLinkedSample() call" and a cache-gated placement would silently
  # contradict that for the entire rest of a session after the first call.
  .nimbleAssertBuildOk()
  # The cache is only trustworthy if EVERY piece of state it gates is still
  # actually true, not just "we did this once": (a) all three cached R
  # objects exist (a session interrupted between the three end-of-function
  # assigns below, or direct tampering with the internal env, could leave a
  # partial cache -- registerDistributions() was confirmed idempotent-safe
  # to call again, so falling through to a full re-setup on ANY gap is
  # cheap and correct, not just defensive); (b) nimble is still attached
  # (confirmed reproducible: detach("package:nimble") between two
  # nimbleLinkedSample() calls left it detached after the second, since the
  # attach check used to live only in the slow path below); (c) the cached
  # shim object file still exists on disk (tempdir() contents can be
  # cleared externally within one session without the directory itself
  # going away). All confirmed independently, not merely plausible.
  .cacheOk <- !is.null(.nimbleGenEnv$dFoceiCondLik) &&
    !is.null(.nimbleGenEnv$rFoceiCondLik) &&
    !is.null(.nimbleGenEnv$condBatchExternal) &&
    !is.null(.nimbleGenEnv$oFile) && file.exists(.nimbleGenEnv$oFile)
  if (.cacheOk) {
    if (!"package:nimble" %in% search()) {
      suppressPackageStartupMessages(attachNamespace("nimble"))
    }
    # the .GlobalEnv bindings below are a SEPARATE piece of state from the
    # cache above (something as ordinary as the user's own rm(list=ls())
    # wipes them without touching .nimbleGenEnv) -- restore them on every
    # cache-hit call, not just the first, or nimbleModel()'s plain-name
    # lookup breaks for the rest of the session with a cryptic "R function
    # 'dFoceiCondLik' ... does not exist" error.
    for (.nm in c("condBatchExternal", "dFoceiCondLik", "rFoceiCondLik")) {
      if (!exists(.nm, envir = .GlobalEnv, inherits = FALSE)) {
        assign(.nm, .nimbleGenEnv[[.nm]], envir = .GlobalEnv)
      }
    }
    return(invisible(TRUE))
  }
  # no rxode2::rxReq("nimble") here: .nimbleAssertBuildOk() at the top of
  # this function already stopped with a clear message if nimble were
  # missing, so by this point it is guaranteed installed.
  # nimble's model-building internals (e.g. init_isDataEnv() ->
  # getNimbleOption()) call unexported helpers by bare name, resolved via the
  # search path -- namespace-qualified nimble:: calls alone are not enough;
  # the package must be attached, exactly like every nimble user script's
  # `library(nimble)`.
  if (!"package:nimble" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("nimble"))
  }
  .o <- .nimbleShimCompile()
  # nimble's build machinery does an unanchored substitution for a literal
  # ".c" somewhere in its oFile path handling (observed: the cache dir
  # ".cache" gets clipped to "ache", corrupting the path and breaking the
  # make rule). Sidestep by handing it a path guaranteed not to contain
  # ".c" anywhere, not just past the cache dir's content-hash cache.
  #
  # tempfile()'s random component is alphanumeric and never contains a
  # period, so retrying tempfile() cannot change whether the match came
  # from the DIRECTORY portion -- every call in one session shares the same
  # tempdir(). If tempdir() itself contains ".c" (an unusual but real
  # possibility depending on TMPDIR/the OS's temp path), a `while
  # (grepl(...)) tempfile()`-style retry loop never terminates: it keeps
  # generating new filenames inside the same still-offending directory.
  # Check the directory once, up front, and fall back rather than loop.
  .oDir <- tempdir()
  if (grepl(".c", .oDir, fixed = TRUE)) {
    .oDir <- normalizePath(".", mustWork = TRUE)
    if (grepl(".c", .oDir, fixed = TRUE)) {
      stop("cannot find a directory without a literal '.c' substring for ",
           "the nimble shim object file (both tempdir() and the working ",
           "directory contain one) -- nimble's build machinery corrupts ",
           "such paths; set a different TMPDIR or working directory and ",
           "retry", call. = FALSE)
    }
  }
  .oSafe <- tempfile("nlmixr2bayes_nimble_shim_", tmpdir = .oDir, fileext = ".o")
  file.copy(.o, .oSafe, overwrite = TRUE)
  condBatchExternal <- nimble::nimbleExternalCall(
    function(theta = double(1), ntheta = integer(0), eta = double(1),
             nid = integer(0), neta = integer(0), value = double(1),
             gradEta = double(1), gradTheta = double(1)) {},
    Cfun = "nlmixr2bayes_nimble_cond_batch_theta",
    headerFile = .nimbleShimHeader(),
    returnType = integer(0),
    oFile = .oSafe
  )
  # like dFoceiCondLik/rFoceiCondLik below, must be globally visible BEFORE
  # dFoceiCondLik is defined against it -- nimble's compiler resolves a
  # referenced nimbleFunction by plain-name lookup at the referencing
  # function's compile time, not by lexical closure capture.
  assign("condBatchExternal", condBatchExternal, envir = .GlobalEnv)

  # eta arrives PRE-FLATTENED (double(1), length nid*neta, row-major) rather
  # than as a double(2) matrix: nimble's BUGS-code parser collapses a
  # whole-array range slice like eta[1:nid, 1:neta] to a lower rank when
  # neta (or nid) is exactly 1 -- a single eta is a common model shape, not
  # an edge case, and that collapse breaks dFoceiCondLik's declared rank-2
  # parameter. The flattening loop instead lives in .nimbleBuildCode()'s
  # BUGS code as scalar-indexed (unambiguous) assignments; only a 1D range
  # ever crosses the distribution-call boundary.
  dFoceiCondLik <- nimble::nimbleFunction(
    run = function(x = double(0), theta = double(1), eta = double(1),
                   ntheta = integer(0), nid = integer(0), neta = integer(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      value <- numeric(nid)
      gradEta <- numeric(nid * neta)
      gradTheta <- numeric(nid * ntheta)
      rc <- condBatchExternal(theta, ntheta, eta, nid, neta, value,
                              gradEta, gradTheta)
      logProb <- -Inf
      if (rc >= 0) logProb <- sum(value)
      if (log) return(logProb)
      return(exp(logProb))
    }
  )

  rFoceiCondLik <- nimble::nimbleFunction(
    run = function(n = integer(0), theta = double(1), eta = double(1),
                   ntheta = integer(0), nid = integer(0), neta = integer(0)) {
      returnType(double(0))
      return(0)
    }
  )

  # nimbleModel()'s parser re-resolves 'dFoceiCondLik'/'rFoceiCondLik' by
  # plain name lookup at BUILD time (not just once at registerDistributions()
  # time, which only captures them via userEnv=parent.frame() for its own
  # internal metadata) -- from wherever nimbleModel() is later called, which
  # is not this function's frame.  NIMBLE's own examples define custom
  # distributions at the top level (globalenv) for exactly this reason;
  # matching that is the documented-idiom fix, not a workaround.
  assign("dFoceiCondLik", dFoceiCondLik, envir = .GlobalEnv)
  assign("rFoceiCondLik", rFoceiCondLik, envir = .GlobalEnv)

  nimble::registerDistributions(list(
    dFoceiCondLik = list(
      BUGSdist = "dFoceiCondLik(theta, eta, ntheta, nid, neta)",
      Rdist = "dFoceiCondLik(theta, eta, ntheta, nid, neta)",
      types = c("value = double(0)", "theta = double(1)", "eta = double(1)",
               "ntheta = integer(0)", "nid = integer(0)", "neta = integer(0)"),
      discrete = FALSE
    )
  ))

  .nimbleGenEnv$condBatchExternal <- condBatchExternal
  .nimbleGenEnv$dFoceiCondLik <- dFoceiCondLik
  .nimbleGenEnv$rFoceiCondLik <- rFoceiCondLik
  .nimbleGenEnv$oFile <- .oSafe
  invisible(TRUE)
}

#' Refuse model shapes phase 1 does not support (rather than mis-state them)
#'
#' Two of these are NIMBLE limitations, not choices: a whole-array range
#' slice (`theta[1:ntheta]`, `etaFlat[1:(nid*neta)]`) passed at a
#' distribution-call-site collapses to a lower rank when its length is
#' exactly 1, and no restructuring found sidesteps it -- confirmed with a
#' `dimensions=` override, a bare (unsliced) reference, and rebuilding the
#' array via a scalar-indexed for loop (the trick that fixes the analogous
#' 2D eta collapse), all independently, all failing identically. A single
#' free theta or a single-subject/single-eta model are not edge cases for a
#' real nlmixr2 model, so refuse them explicitly here rather than let a
#' cryptic NIMBLE compiler error be the first sign something is wrong.
#' `fix()`ed thetas are refused too: `literalFix=TRUE` (the default) drops
#' them from `foceiLikLoad`'s internal theta count, so the count this
#' function's caller builds (`nrow(map$theta)`, which still includes them)
#' would disagree with the linked problem's -- caught safely by
#' `nimbleLinkedSample()`'s own consistency check either way, but refusing
#' up front gives a clear reason instead of a generic mismatch error.
#'
#' A multivariate `ini({})` prior on population parameters, or any prior on
#' an omega block, is refused here too: `.nimbleBuildCode()` only ever
#' consults a NON-multivariate matching row of `priPop` for a theta's prior,
#' so silently proceeding would mean the user's explicit prior is simply
#' never applied -- exactly the kind of silent mis-handling this function
#' exists to prevent. Omega is not sampled yet (phase 1/2 fixes it at its
#' `ini()` value), so a prior on it can never be honored either.
#' @param priPop `stanPriors(ui)$pop`, or `NULL` to skip the prior checks
#' @param priOmega `stanPriors(ui)$omega`, or `NULL` to skip the check
#' @noRd
.nimbleAssertSupported <- function(map, priPop = NULL, priOmega = NULL) {
  if (map$nMix > 1L) {
    stop("nimble linked sampling does not yet support finite mixtures",
         call. = FALSE)
  }
  if (nrow(map$muRefCov) > 0L) {
    stop("nimble linked sampling does not yet support mu-referenced ",
         "covariates", call. = FALSE)
  }
  if (!is.null(priPop) && nrow(priPop) > 0L) {
    .mv <- unique(priPop$name[priPop$kind == "multivariate" &
                               priPop$name %in% map$theta$name])
    if (length(.mv) > 0L) {
      stop("nimble linked sampling does not yet support multivariate ",
           "priors on population parameter(s) '",
           paste(.mv, collapse = ", "), "'", call. = FALSE)
    }
  }
  if (!is.null(priOmega) && nrow(priOmega) > 0L) {
    stop("nimble linked sampling does not yet support priors on omega ",
         "blocks (omega is fixed at its ini() value, not sampled)",
         call. = FALSE)
  }
  if (any(map$theta$fix)) {
    stop("nimble linked sampling does not yet support fix()ed population ",
         "parameter(s) '", paste(map$theta$name[map$theta$fix], collapse = ", "),
         "'", call. = FALSE)
  }
  if (nrow(map$theta) < 2L) {
    stop("nimble linked sampling needs at least 2 population parameters ",
         "(nimble collapses a length-1 array at a distribution call site; ",
         "this model has ", nrow(map$theta), ")", call. = FALSE)
  }
  for (.b in map$blocks) {
    if (.b$k != 1L) {
      stop("nimble linked sampling supports diagonal omega blocks only; ",
           "block over '", paste(.b$members, collapse = ", "),
           "' is correlated (dimension ", .b$k, ")", call. = FALSE)
    }
    if (any(.b$fix)) {
      stop("nimble linked sampling does not yet support fixed-variance ",
           "eta '", .b$members, "'", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Build the NIMBLE BUGS code for a (phase-1-supported) parameter map
#'
#' Prior hyperparameters (theta mean/sd, fixed-theta values, eta sd) are
#' inlined as LITERALS directly in the generated code text -- mirroring
#' `stanGen.R`'s `.stanNum()` treatment of fixed values -- rather than
#' passed through as indexable `constants=` vectors. nimble infers a
#' length-1 numeric constant as a bare scalar with no array dimension, and
#' a model with a single free theta or a single eta (a common shape, not an
#' edge case) then breaks on `thetaSd[1]`/`etaSd[1]`-style indexed access
#' with no reliable override (`dimensions=` on a constant either conflicts
#' with nimble's own inference or is silently insufficient). Literals sidestep
#' the ambiguity entirely: there is no constant array to mis-infer the rank of.
#' `priPop` (`stanPriors(ui)$pop`) supplies the `ini({})`-declared prior for a
#' free theta when one matches (see `.nimbleThetaPriorText()`,
#' R/nimblePriors.R); a theta without a matching row falls back to the
#' original default (a weakly-informative normal at its own declared
#' bounds). Bounds always apply via NIMBLE's `T()` truncation, whether from
#' the matched prior's support-promoted bounds or the theta's own.
#' @param thetaSdVec free-theta DEFAULT prior SDs (natural scale, used only
#'   when no `priPop` row matches), same order as `map$theta`
#' @param etaSdVec per-eta prior SD (natural scale), same order as `map$eta`
#' @param priPop `stanPriors(ui)$pop`, or `NULL` for the default prior on
#'   every free theta
#' @noRd
.nimbleBuildCode <- function(map, thetaSdVec, etaSdVec, priPop = NULL) {
  .lines <- character(0)
  .thetaNames <- map$theta$name
  for (.i in seq_len(nrow(map$theta))) {
    .r <- map$theta[.i, ]
    if (.r$fix) {
      .lines <- c(.lines, sprintf("theta[%d] <- %s", .i, .stanNum(.r$est)))
      next
    }
    .priRow <- NULL
    if (!is.null(priPop) && nrow(priPop) > 0L) {
      .w <- which(priPop$name == .r$name & priPop$kind != "multivariate")
      if (length(.w) == 1L) .priRow <- priPop[.w, ]
    }
    .distText <- .nimbleThetaPriorText(.r$name, .r$est, thetaSdVec[.i],
                                       .priRow, .r$lower, .r$upper,
                                       .thetaNames)
    .lines <- c(.lines, sprintf("theta[%d] ~ %s", .i, .distText))
  }
  for (.b in map$blocks) {
    .lines <- c(.lines, sprintf(
      "for (i in 1:nid) { eta[i, %d] ~ dnorm(0, sd = %s) }",
      .b$start, .stanNum(etaSdVec[.b$start])))
  }
  # flatten via scalar-indexed (unambiguous) assignments in BUGS code itself
  # -- only a 1D range (nid*neta elements, never exactly 1 for any model
  # with at least one subject and one eta) crosses the distribution-call
  # boundary; see the eta=double(1) comment in .nimbleCondLikSetup().
  .lines <- c(.lines,
             "for (i in 1:nid) { for (j in 1:neta) { etaFlat[(i - 1) * neta + j] <- eta[i, j] } }",
             "zero ~ dFoceiCondLik(theta[1:ntheta], etaFlat[1:(nid * neta)], ntheta, nid, neta)")
  .txt <- paste0("nimble::nimbleCode({\n  ", paste(.lines, collapse = "\n  "), "\n})")
  eval(parse(text = .txt))
}

#' Sample a linked FOCEi conditional likelihood with NIMBLE (phase 1)
#'
#' Diagonal-omega, no-covariate, no-mixture models only; see the file header
#' comment for the full scope statement. Reuses [stanLinkSetup()] (and the
#' same theta-base/mu-reference setup [stanEst.R] uses) verbatim -- that
#' machinery is engine-agnostic.
#'
#' @param ui rxode2 ui
#' @param data estimation data
#' @param niter total MCMC iterations per chain
#' @param nburnin burn-in iterations
#' @param thetaSd prior SD (natural scale) for free thetas; default
#'   `pmax(3 * abs(est), 1)`
#' @param likelihood `"focei"` or `"foce"` (see [stanLinkSetup()])
#' @param rxControl solving options
#' @param cores subject-parallel thread count for each likelihood evaluation
#' @param literalFix inline fixed population parameters as literals (kept
#'   matched with the linked problem, as in [stanLinkSetup()])
#' @param seed passed to [nimble::runMCMC()]
#' @param verbose show NIMBLE's model-compiler output
#' @return the `runMCMC()` samples
#' @export
#' @author Matthew L Fidler
nimbleLinkedSample <- function(ui, data, niter = 2000L, nburnin = 1000L,
                               thetaSd = NULL,
                               likelihood = c("focei", "foce"),
                               rxControl = rxode2::rxControl(),
                               cores = rxode2::getRxThreads(),
                               literalFix = TRUE, seed = NULL,
                               verbose = FALSE) {
  likelihood <- match.arg(likelihood)
  if (!is.null(thetaSd)) {
    # an NA/NaN/Inf here would get baked as a literal into the generated
    # NIMBLE code text (.stanNum(NA) is the string "NA", parseable but not
    # a usable prior) -- refuse it here rather than emit a broken model.
    checkmate::assertNumeric(thetaSd, finite = TRUE, any.missing = FALSE,
                             min.len = 1)
  }
  # no rxode2::rxReq("nimble") here: it used to intercept a missing-nimble
  # session with a generic message before .nimbleCondLikSetup()'s own
  # .nimbleAssertBuildOk() call (unconditional, at its very top) ever got a
  # chance to give the clearer, function-specific one -- let that check run
  # first instead of shadowing it.
  .nimbleCondLikSetup()
  .ui <- rxode2::rxode2(ui)
  .map <- .stanMap(.ui)
  .pri <- stanPriors(.ui)
  .nimbleAssertSupported(.map, .pri$pop, .pri$omega)
  .h <- stanLinkSetup(.ui, data, likelihood = likelihood, rxControl = rxControl,
                      thetaSens = FALSE, literalFix = literalFix, cores = cores)
  on.exit(stanLinkFree(), add = TRUE)
  on.exit(.Call(`_nlmixr2bayes_clearThetaBase`), add = TRUE)
  if (!identical(.h$etaNames, .map$eta$name) ||
        !identical(.h$thetaNames, .map$theta$name)) {
    stop("the linked problem's parameters do not match the model map",
         call. = FALSE) # nocov
  }
  .Call(`_nlmixr2bayes_setThetaBase`, as.double(.h$initPar))
  .Call(`_nlmixr2bayes_setMuRef`, as.integer(.map$muRefIdx))

  .ntheta <- nrow(.map$theta)
  .neta <- nrow(.map$eta)
  .nid <- .h$nid
  if (.nid * .neta < 2L) {
    stop("nimble linked sampling needs at least 2 total (subject, eta) ",
         "combinations (the same length-1 distribution-call-site collapse ",
         "as the theta check; this model has ", .nid, " subject(s) and ",
         .neta, " eta(s))", call. = FALSE)
  }
  .thetaMean <- .map$theta$est
  .thetaSdVec <- if (is.null(thetaSd)) {
    pmax(3 * abs(.thetaMean), 1)
  } else {
    rep(thetaSd, length.out = .ntheta)
  }
  .etaSd <- rep(NA_real_, .neta)
  for (.b in .map$blocks) .etaSd[.b$start] <- sqrt(.b$init[1, 1])

  .code <- .nimbleBuildCode(.map, .thetaSdVec, .etaSd, .pri$pop)
  .constants <- list(ntheta = .ntheta, nid = .nid, neta = .neta)
  .inits <- list(theta = .map$theta$est, eta = matrix(0, .nid, .neta))
  .dataList <- list(zero = 0)
  # `eta` is the one MODEL VARIABLE (a stochastic node, not a plain
  # constant) whose shape genuinely needs pinning: nimble's BUGS-code parser
  # can otherwise collapse a whole-array reference to a lower rank when nid
  # or neta is exactly 1 -- a single eta is a common model shape, not an
  # edge case. Plain constants (thetaMean/thetaSd/etaSd/thetaFix) are NOT
  # listed here: nimble infers their dimension directly and unambiguously
  # from the constants list itself, and an explicit override that disagrees
  # with that inference (as happens for a length-1 constant) errors instead
  # of being applied.
  .dimensions <- list(theta = .ntheta, eta = c(.nid, .neta))

  # calculate=FALSE: nimbleExternalCall only runs in COMPILED nimbleFunctions
  # (documented NIMBLE limitation), so the R-interpreted initial-value pass
  # nimbleModel() would otherwise run here fails on the custom distribution.
  # Initial logProbs are calculated below, AFTER compilation, instead.
  .m <- nimble::nimbleModel(.code, constants = .constants, data = .dataList,
                            inits = .inits, dimensions = .dimensions,
                            calculate = FALSE)
  .cm <- nimble::compileNimble(.m, showCompilerOutput = verbose)
  .cm$calculate()
  .conf <- nimble::configureMCMC(.m)
  .mcmc <- nimble::buildMCMC(.conf)
  .cmcmc <- nimble::compileNimble(.mcmc, project = .m)
  nimble::runMCMC(.cmcmc, niter = niter, nburnin = nburnin, setSeed = seed)
}
