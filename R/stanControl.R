#' Control options for the `est="stan"` estimation method
#'
#' The class name `stanControl` makes `nlmixr2(model, data, stanControl(...))`
#' infer `est="stan"`.
#'
#' Parallelism has two axes.  On unix, the CHAINS run in forked
#' processes (`chainCores`, inheriting the [rxode2::getRxThreads()]
#' budget capped at `chains`) with serial evaluation inside each chain.
#' On Windows, `chainCores > 1` runs one chain per PSOCK worker (each
#' worker builds its own linked likelihood); `cores` remains the
#' subject-parallel setting per worker.  With `chainCores = 1` the chains run sequentially and each
#' likelihood evaluation is subject-parallel instead (`cores`,
#' defaulting to `getRxThreads()`).  `getRxThreads()` honors
#' `OMP_THREAD_LIMIT`, which CRAN sets to 2, so checks stay within
#' CRAN's core policy automatically.
#'
#' The ODE-retry machinery (`maxOdeRecalc`, `stickyRecalcN`) stays at
#' nlmixr2est's defaults and is ENABLED: a hard point retries with relaxed
#' tolerances instead of rejecting, which beats truncating the posterior to
#' `-Inf` where the density is actually finite.  This is safe for a sampler
#' because the C batch entries reset all retry/sticky/tolerance state at the
#' top of every evaluation (per batch and per subject), so the same point
#' reproduces the same retry ladder and the same value -- the density stays
#' a pure function of state, which gate G2 proves bitwise, with the retry
#' machinery live.
#'
#' @param algorithm the Stan inference algorithm: `"NUTS"` (default; the
#'   No-U-Turn HMC sampler -- full Bayes), or Stan's two ADVI variational
#'   approximations `"meanfield"` (independent Gaussians on the
#'   unconstrained scale) and `"fullrank"` (one full-covariance Gaussian),
#'   run through [rstan::vb()], or `"pathfinder"` (Stan's multi-path
#'   Pathfinder, run through the StanEstimators package against the
#'   compiled model's exact log density and analytic gradient --
#'   L-BFGS trajectories from jittered starts, ELBO-selected normal
#'   approximations, PSIS-resampled draws; needs `StanEstimators`
#'   installed since rstan does not expose Pathfinder yet).  ADVI is typically 10-100x faster and is a
#'   fair first look, but it is an APPROXIMATION -- it understates tails
#'   and correlations; the Pareto-k diagnostic (`khat`) replaces
#'   Rhat/divergences and `khat > 0.7` means the approximation is not
#'   reliable and NUTS should be used.  All algorithms run the same
#'   generated program against the same linked rxode2/nlmixr2est
#'   likelihood and return the same complete nlmixr2 fit
#' @param pathfinderPaths number of independent Pathfinder paths
#'   (`algorithm = "pathfinder"`); the multi-path PSIS mixture is what
#'   makes Pathfinder robust to a single path landing badly
#' @param vbIter maximum ADVI iterations (`algorithm` != "NUTS")
#' @param vbTolRelObj ADVI relative-ELBO convergence tolerance
#' @param vbOutputSamples posterior draws taken from the fitted
#'   approximation (these flow into every posterior summary on the fit)
#' @param chains number of MCMC chains
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
#' @param adapt_delta target acceptance rate; matches Stan's own default
#'   (0.8) so an `est="stan"` run and a hand-written rstan run of the
#'   same model behave identically out of the box.  Stiff PK posteriors
#'   that show divergences should raise it toward 0.95-0.99
#' @param max_treedepth NUTS maximum tree depth (Stan's default, 10)
#' @param likelihood the individual likelihood the linked problem evaluates;
#'   `"focei"` (interaction) or `"foce"` -- the two self-consistent
#'   value/gradient pairs
#' @param etaParam eta parameterization: `"noncentered"` (default; Stan
#'   samples `z` with `eta = z L'`, the funnel-robust form) or `"centered"`
#'   (Stan samples `eta` directly with `eta[i] ~ multi_normal_cholesky(0, L)`).
#'   The posteriors are identical by construction (gate G13); the default
#'   usually mixes better with few subjects
#' @param print the familiar nlmixr2est iteration-print cadence: every
#'   `print`-th log-density evaluation prints a row of natural-scale thetas
#'   plus the current omega variances/covariances (`om.<eta>`,
#'   `cov.<eta1>.<eta2>`) and `-2*sum(conditional log-lik)`, and records it
#'   in `$parHistData`.  A NUTS chain makes many evaluations per iteration,
#'   so the default (100) is roughly a row every dozen sampler iterations;
#'   `0` turns printing and history off
#' @param maxOdeRecalc failed-solve tolerance relaxations inside each
#'   likelihood evaluation before the FD fallback (0 disables; the
#'   evaluators reset this state per call so the density stays a pure
#'   function of the parameters either way)
#' @param fallbackFD after the relaxations, retry the prediction model
#'   with shi21 central-difference eta gradients before rejecting the
#'   point with -Inf
#' @param chainCores processes for parallel chains (rstan's `cores`).
#'   Default (`NULL`): on unix, `min(chains, rxode2::getRxThreads())` --
#'   chains run in FORKED processes, each fork receiving a copy-on-write
#'   duplicate of the linked likelihood state (verified: draws are
#'   bit-identical to a sequential run, ~1.6x faster on 2 chains).  When
#'   chains fork, `cores` (subject threads) defaults to 1 -- OpenMP
#'   inside a forked child of an OpenMP-using parent can deadlock -- and
#'   the iteration print is disabled (the ticks would happen in the
#'   forks' memory).  Set `chainCores = 1` for sequential chains with
#'   subject-parallel evaluation and the iteration table.  On Windows
#'   `chainCores > 1` uses PSOCK workers (one linked likelihood per
#'   worker process) and runs one chain per worker; `cores` remains the
#'   subject-parallel setting per worker
#' @param cores subject-parallel thread count inside each likelihood
#'   evaluation (default [rxode2::getRxThreads()]; on CRAN this is capped
#'   at 2 through `OMP_THREAD_LIMIT`).  This is the per-run analogue of the
#'   `cores` argument every other nlmixr2est control carries; Stan's own
#'   chains still run sequentially
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
#' @param rxControl solving options ([rxode2::rxControl()]).  When `NULL`
#'   (default) the solver is the dense AutoSwitch composite
#'   (`method="dop853+ros4", dense=TRUE`) at `sigdig`-tight tolerances:
#'   dense-output Dormand-Prince stepping -- the twin of Stan's own
#'   default ODE solver (`ode_rk45`), 25-40% cheaper per gradient than
#'   lsoda on typical PK models -- with an automatic dense Rosenbrock
#'   fallback when the system turns stiff
#' @param addProp,sumProd,optExpression,literalFix,calcTables,compress,ci,sigdigTable
#'   standard nlmixr2 passthroughs for the output machinery
#' @param ... only `genRxControl`
#' @return a `stanControl` list
#' @export
#' @author Matthew L Fidler
stanControl <- function(chains = 4L, iter = 2000L, warmup = floor(iter / 2),
                        algorithm = c("NUTS", "meanfield", "fullrank",
                                      "pathfinder"),
                        pathfinderPaths = 4L,
                        vbIter = 10000L, vbTolRelObj = 0.01,
                        vbOutputSamples = 1000L,
                        thin = 1L, seed = 42L,
                        init = c("ini", "random"), initJitterSd = 0.1,
                        adapt_delta = 0.8, max_treedepth = 10L, # nolint: object_name_linter.
                        likelihood = c("focei", "foce"),
                        etaParam = c("noncentered", "centered"),
                        cores = rxode2::getRxThreads(),
                        print = 100L,
                        maxOdeRecalc = 3L, fallbackFD = TRUE,
                        chainCores = NULL,
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
  algorithm <- match.arg(algorithm)
  checkmate::assertIntegerish(pathfinderPaths, lower = 1, len = 1,
                              any.missing = FALSE)
  checkmate::assertIntegerish(vbIter, lower = 10, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(vbTolRelObj, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(vbOutputSamples, lower = 10, len = 1,
                              any.missing = FALSE)
  checkmate::assertIntegerish(chains, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(iter, lower = 2, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(warmup, lower = 1, upper = iter - 1, len = 1)
  checkmate::assertIntegerish(thin, lower = 1, len = 1)
  checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(initJitterSd, lower = 0, len = 1)
  checkmate::assertNumeric(adapt_delta, lower = 0, upper = 1, len = 1)
  checkmate::assertIntegerish(max_treedepth, lower = 1, upper = 20, len = 1)
  checkmate::assertIntegerish(cores, lower = 1, len = 1)
  checkmate::assertIntegerish(print, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(maxOdeRecalc, lower = 0, len = 1,
                              any.missing = FALSE)
  checkmate::assertLogical(fallbackFD, len = 1, any.missing = FALSE)
  .coresExplicit <- !missing(cores)
  if (is.null(chainCores)) {
    # default: parallel chains inheriting the core budget (capped at the
    # chain count).  On unix these are forked chains; on Windows they use
    # PSOCK workers with one linked likelihood per worker.
    chainCores <- if (identical(.Platform$OS.type, "unix")) {
      min(as.integer(chains), as.integer(rxode2::getRxThreads()))
    } else {
      1L
    }
  }
  checkmate::assertIntegerish(chainCores, lower = 1, len = 1,
                              any.missing = FALSE)
  if (chainCores > 1L) {
    if (identical(.Platform$OS.type, "unix")) {
      if (.coresExplicit && cores > 1L) {
        warning("subject-parallel OpenMP (cores = ", cores, ") inside ",
                "FORKED chains (chainCores = ", chainCores, ") risks the ",
                "OpenMP-after-fork hazard; cores = 1 is the safe setting",
                call. = FALSE)
      } else if (!.coresExplicit) {
        # the forks are the parallelism; OpenMP inside a forked child of an
        # OpenMP-tainted parent can deadlock, so the default stays serial
        # within each chain
        cores <- 1L
      }
    } else {
      # PSOCK workers are independent processes, so the OpenMP-after-fork
      # hazard does not apply here; keep the requested per-worker thread
      # count and guard with the aggregate-core warning below.
    }
  }
  if (as.integer(chainCores) * as.integer(cores) >
        as.integer(rxode2::getRxThreads())) {
    warning("requested chainCores * cores (", as.integer(chainCores), " * ",
            as.integer(cores), ") exceeds rxode2::getRxThreads() (",
            as.integer(rxode2::getRxThreads()), "); this may oversubscribe CPUs",
            call. = FALSE)
  }
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
  etaParam <- match.arg(etaParam)
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
      # on evaluation history (the batch entries also reset it per call).
      # method: the dense AutoSwitch composite dop853+ros4 -- dense-output
      # Dormand-Prince stepping (the twin of Stan's own default ode_rk45:
      # large internal steps, interpolation at observation times, ~25-40%
      # cheaper per gradient than lsoda on typical PK models) with an
      # automatic dense Rosenbrock (ros4) fallback when the system turns
      # stiff, so stiff models stay correct without user intervention.
      rxControl <- rxode2::rxControl(rtol = .tol, atol = .tol,
                                     ssRtol = .tol, ssAtol = .tol,
                                     method = "dop853+ros4", dense = TRUE,
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
               algorithm = algorithm,
               pathfinderPaths = as.integer(pathfinderPaths),
               vbIter = as.integer(vbIter),
               vbTolRelObj = vbTolRelObj,
               vbOutputSamples = as.integer(vbOutputSamples),
               warmup = as.integer(warmup), thin = as.integer(thin),
               seed = as.integer(seed), init = init,
               initJitterSd = initJitterSd,
               adapt_delta = adapt_delta,
               max_treedepth = as.integer(max_treedepth),
               likelihood = likelihood, etaParam = etaParam,
               cores = as.integer(cores),
               print = as.integer(print),
               maxOdeRecalc = as.integer(maxOdeRecalc),
               fallbackFD = fallbackFD,
               chainCores = as.integer(chainCores),
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

#' @author Matthew L Fidler
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
#' @author Matthew L Fidler
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

#' @author Matthew L Fidler
#' @export
nmObjHandleControlObject.stanControl <- function(control, env) {
  assign("stanControl", control, envir = env)
}

#' @author Matthew L Fidler
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
