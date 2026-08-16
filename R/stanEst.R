# est="stan" (Specs 7, 12, 13, 14): dispatch, assertions, sampling driver,
# and the posterior -> nlmixr2 fit mapping.

#' Suggested prior() lines for a model without any
#' @noRd
.stanSuggestPriors <- function(ui) {
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta) & !.iniDf$fix, , drop = FALSE]
  .lines <- vapply(seq_len(nrow(.th)), function(.i) {
    .r <- .th[.i, ]
    if (!is.na(.r$err) && is.finite(.r$lower) && .r$lower >= 0) {
      paste0("prior(", .r$name, ") ~ dcauchy(0, ",
             .stanNum(signif(5 * max(abs(.r$est), 0.1), 3)), ")")
    } else {
      paste0("prior(", .r$name, ") ~ dnorm(", .stanNum(signif(.r$est, 3)),
             ", ", .stanNum(signif(10 * max(1, abs(.r$est)), 3)), ")")
    }
  }, character(1))
  .et <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2 &
                  !.iniDf$fix, , drop = FALSE]
  if (nrow(.et) > 0L) {
    .lines <- c(.lines, paste0("prior(", .et$name[1], ") ~ invWishart(",
                               nrow(.et) + 3L, ")"))
  }
  .lines
}

#' Does the linked nlmixr2est build the theta forward-sensitivity model with
#' rxode2's analytic event ("jump") sensitivities?  With nlmixr2/nlmixr2est#946
#' (PR #947) the model is compiled per the control's `eventSens` and its event
#' shape installed around the gradient batch, so a dose-handling theta's
#' derivative carries the jump condition at the event.  Both sides ship as
#' version 7.0.3, so the capability is probed structurally: the fix introduced
#' the `eventSens` argument on the model builder.
#' @noRd
.stanHasEventThetaSens <- function() {
  .f <- get0(".impmapThetaSensModel", envir = asNamespace("nlmixr2est"))
  is.function(.f) && "eventSens" %in% names(formals(.f))
}

#' Does the linked nlmixr2est carry the estimated transform-both-sides
#' lambda's conditional sensitivity?  nlmixr2/nlmixr2est#949 (PR #950) added
#' the rx__sens_rx_lambda__BY_THETA columns + the DV-side chain rule; both
#' sides ship as 7.0.3, so probe the builder's body for the marker.
#' @noRd
.stanHasTbsLambdaSens <- function() {
  .f <- get0("rxUiGet.impmapThetaSens", envir = asNamespace("nlmixr2est"))
  is.function(.f) && any(grepl("rx__sens_rx_lambda__BY_THETA",
                               deparse(body(.f)), fixed = TRUE))
}

#' Per-endpoint DV-transform Jacobian statistics for estimated
#' transform-both-sides lambdas (gen$tbsJac): Box-Cox sJ = sum(log y);
#' Yeo-Johnson sJ = sum(log1p(y)) over y >= 0 minus sum(log1p(-y)) over
#' y < 0.  Observations are matched to the lambda's endpoint by CMT when
#' the model has several endpoints.
#' @noRd
.stanTbsJacValues <- function(gen, ui, dataSav) {
  if (length(gen$tbsJac) == 0L) return(gen)
  .obs <- dataSav[dataSav$EVID == 0 & !is.na(dataSav$DV), , drop = FALSE]
  .pd <- ui$predDf
  for (.j in gen$tbsJac) {
    .rows <- .obs
    if (nrow(.pd) > 1L && "CMT" %in% names(.obs)) {
      .cmt <- .pd$cmt[match(.j$condition, .pd$cond)]
      .rows <- .obs[.obs$CMT == .cmt, , drop = FALSE]
    }
    .dv <- .rows$DV
    if (identical(.j$transform, "boxCox")) {
      if (any(.dv <= 0)) {
        stop("boxCox(", .j$theta, ") needs strictly positive DV and ",
             "endpoint '", .j$condition, "' has ", sum(.dv <= 0),
             " non-positive observation(s)", call. = FALSE)
      }
      .s <- sum(log(.dv))
    } else {
      .s <- sum(log1p(.dv[.dv >= 0])) - sum(log1p(-.dv[.dv < 0]))
    }
    gen$data[[.j$name]] <- .s
  }
  gen
}

#' Estimated thetas that enter dose handling (alag/f/dur/rate), transitively
#' through intermediate assignments
#' @noRd
.stanEventThetas <- function(ui) {
  .exprs <- ui$lstExpr
  .lhs <- vapply(.exprs, function(e) {
    if (is.call(e) && length(e) >= 3L) deparse1(e[[2]]) else ""
  }, character(1))
  .isEvent <- grepl("^(alag|lag|f|F|rate|dur)\\(", .lhs)
  if (!any(.isEvent)) return(character(0))
  .need <- unique(unlist(lapply(.exprs[.isEvent], function(e) all.vars(e[[3]]))))
  repeat {
    .add <- character(0)
    for (.i in which(!.isEvent)) {
      if (.lhs[.i] %in% .need) .add <- c(.add, all.vars(.exprs[[.i]][[3]]))
    }
    .new <- setdiff(.add, .need)
    if (length(.new) == 0L) break
    .need <- c(.need, .new)
  }
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta) & !.iniDf$fix, "name"]
  intersect(.need, .th)
}

#' Capability assertions + prior triage for est="stan"
#' @noRd
.stanAssert <- function(ui) {
  rxode2::assertRxUiPrediction(ui, " for est=\"stan\"", .var.name = ui$modelName)
  # population-only (no-eta) models are tier 0 via the nlm C API when the
  # loaded nlmixr2est provides it (nlmixr2/nlmixr2est#953)
  .noEta <- !any(!is.na(ui$iniDf$neta1))
  if (.noEta && !.stanHasNlmApi()) {
    stop("est=\"stan\" needs a mixed model with this nlmixr2est ",
         "(population-only models use the nlm C API; update nlmixr2est, ",
         "nlmixr2/nlmixr2est#953)", call. = FALSE)
  }
  rxode2::assertRxUiRandomOnIdOnly(ui, " for est=\"stan\"", .var.name = ui$modelName)
  rxode2::assertRxUiNoMix(ui, " for est=\"stan\"", .var.name = ui$modelName)
  if (!.noEta) {
    rxode2::assertRxUiMixedOnly(ui, " for est=\"stan\"", .var.name = ui$modelName)
  }
  # (rxode2 >= 5.1.7 always has llik support, so no assertRxUiTransformNormal)
  # transform-both-sides with an ESTIMATED lambda: two lambda-dependent
  # pieces are needed.  The DV-transform Jacobian is pure data-times-lambda
  # algebra the generator can emit Stan-side with exact autodiff (Box-Cox:
  # (lambda-1)*sum(log y) per endpoint), so it is NOT the blocker; the
  # blocker is the linked conditional's own d/dlambda through the
  # transformed residual h(y;lambda) - h(f;lambda), whose theta-sensitivity
  # column is silently zero upstream (nlmixr2/nlmixr2est#949, FD-measured).
  # Refused until that column is real; a FIXED lambda works today (every
  # lambda term is then constant in the sampled parameters).
  .iniDf <- ui$iniDf
  .tbs <- which(.iniDf$err %in% c("boxCox", "yeoJohnson") & !.iniDf$fix)
  if (length(.tbs) > 0L && !.stanHasTbsLambdaSens()) {
    stop("est=\"stan\" cannot estimate the transform-both-sides ",
         "parameter(s) ",
         paste0("'", .iniDf$name[.tbs], "'", collapse = ", "),
         " with this nlmixr2est: the linked conditional's d/dlambda ",
         "sensitivity column is silently zero; fix() the lambda or ",
         "update nlmixr2est (>= the nlmixr2/nlmixr2est#949 fix)",
         call. = FALSE)
  }
  # estimated dose-handling parameters: the derivative through the event
  # needs a jump condition.  nlmixr2est with nlmixr2/nlmixr2est#946 compiles
  # the theta forward-sensitivity model with rxode2's analytic event ("jump")
  # sensitivities, so the column is real (FD-verified); an OLDER nlmixr2est
  # advertises the theta in the sensitivity index but leaves its column
  # silently zero, which a gradient-based sampler experiences as a
  # value/gradient mismatch -- refused rather than sampled wrong.
  .evTh <- .stanEventThetas(ui)
  if (length(.evTh) > 0L && !.stanHasEventThetaSens()) {
    stop("est=\"stan\" cannot estimate dose-handling parameter(s) ",
         paste0("'", .evTh, "'", collapse = ", "),
         " (alag/f/dur/rate) with this nlmixr2est: the linked theta ",
         "sensitivities do not propagate through the event jump; fix() ",
         "them or update nlmixr2est (>= the nlmixr2/nlmixr2est#946 fix)",
         call. = FALSE)
  }
  # D10: no priors at all is an error, with the exact lines to add --
  # Bayesian inference with silently-invented priors produces a
  # publishable-looking wrong answer
  .pri <- rxode2::rxUiPriors(ui)
  if (nrow(.pri) == 0L) {
    stop("est=\"stan\" needs prior distributions in the ini({}) block and ",
         "this model declares none; for example:\n  ",
         paste(.stanSuggestPriors(ui), collapse = "\n  "),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Bayesian estimation via Stan linked to the rxode2/nlmixr2est likelihood
#'
#' Dispatched by `nlmixr2(model, data, est="stan", stanControl(...))`.  The
#' first nlmixr2 estimation method to declare `nlmixr2Priors = "all"`: prior
#' distributions come from the model's `ini({})` block and are required
#' (a model without any errors, printing the lines to add).
#'
#' @param env nlmixr2 estimation environment
#' @param ... passed through
#' @return an nlmixr2 fit (or an `nlmixr2stanCode` when
#'   `stanControl(run=FALSE)`)
#' @keywords internal
#' @export
nlmixr2Est.stan <- function(env, ...) {
  .ui <- env$ui
  .stanAssert(.ui)
  .control <- env$control
  if (!inherits(.control, "stanControl")) {
    .control <- stanControl() # nocov
  }
  .stanFamilyFit(env, .ui, .control, ...)
}
attr(nlmixr2Est.stan, "nlmixr2Priors") <- "all"
attr(nlmixr2Est.stan, "type") <- "Markov chain Monte Carlo"
attr(nlmixr2Est.stan, "description") <-
  "Hamiltonian Monte Carlo (Stan/rstan) linked to the rxode2/nlmixr2est likelihood"
attr(nlmixr2Est.stan, "covPresent") <- TRUE
attr(nlmixr2Est.stan, "mu") <- FALSE
# unbounded=FALSE is load-bearing: the natural lower/upper survive in iniDf
# for the generator to turn into Stan constraints
attr(nlmixr2Est.stan, "unbounded") <- FALSE
#' IOV readiness: everything downstream is wired (the .uiApplyIov hook
#' expansion flows through the generator's fixed-eta blocks and the whole
#' assembled target was FD-checked end to end), but the upstream
#' theta-sensitivity column for the IOV magnitude theta is wrong by a
#' theta-dependent factor (x1.42 at 0.1, x1.14 at 0.2 -- measured, filed
#' as nlmixr2/nlmixr2est#952).  A scaled-wrong gradient is exactly the
#' silent kind a sampler turns into a wrong posterior, so IOV stays OFF
#' until that column FD-agrees.
#' @noRd
.stanHasIovSens <- function() {
  FALSE # flip via a probe when nlmixr2/nlmixr2est#952 lands
}
# IOV via nlmixr2est's preprocessing hook (.uiApplyIov): occasion-level
# random effects become per-occasion FIXED unit-variance etas scaled by an
# estimated sd theta in the model code; the generator's fixed-eta blocks
# (constant L, z still sampled) carry them, and the sd theta rides the
# forward sensitivities like any structural theta.  Gated on the upstream
# sensitivity being right (#952).
attr(nlmixr2Est.stan, "iov") <- function(control) .stanHasIovSens()

#' The orchestration: preprocess -> map -> generate -> (compile -> link ->
#' sample -> finalize)
#' @noRd
.stanFamilyFit <- function(env, ui, control, ...) {
  .ret <- new.env(parent = emptyenv())
  .ret$table <- if (is.null(env$table)) nlmixr2est::tableControl() else env$table
  .ret$stanControl <- control
  nlmixr2est::.foceiPreProcessData(env$data, .ret, ui, control$rxControl)
  .map <- .stanMap(ui)
  .pri <- stanPriors(ui)
  .gen <- .stanGenerate(ui, .map, .pri, control)
  .nid <- length(.ret$idLvl)
  .gen$data$N <- .nid
  # DV-transform Jacobian statistics for any estimated TBS lambda
  .gen <- .stanTbsJacValues(.gen, ui, .ret$dataSav)
  # per-subject mu-referenced covariate values for the coefficient gradient
  # scatter (d/dtheta_p = cov_i * d/deta_k); a time-varying covariate is
  # flagged, not refused -- its coefficient rides the forward-sensitivity
  # model like any other structural theta (nlmixr2est treats a time-varying
  # regressor the same way)
  .cov <- .stanMuRefCovValues(.map, .ret$dataSav)
  for (.n in .gen$notes) cli::cli_inform(paste0("est=\"stan\": ", .n))
  if (!is.null(control$stanFile)) {
    writeLines(.gen$code, control$stanFile)
  }
  if (!isTRUE(control$run)) {
    # $env carries the ui the way a fit does, so nlmixr2's post-dispatch
    # restore (which un-does the literal-fix/zero-omega preprocessing on the
    # returned object's ui) has somewhere to write
    .env <- new.env(parent = emptyenv())
    .env$ui <- ui
    .out <- list(code = .gen$code, data = .gen$data, map = .map,
                 priors = .pri, notes = .gen$notes, ui = ui,
                 control = control, env = .env)
    class(.out) <- "nlmixr2stanCode"
    return(.out)
  }
  rxode2::rxReq("rstan")
  # ---- compile (cached), then link --------------------------------------
  .sm <- stanCompile(.gen$code, cache = control$cache,
                     cacheDir = control$cacheDir, verbose = control$verbose)
  if (isTRUE(.gen$pop)) {
    # ---- tier 0: population-only (no etas) via the nlm C API ------------
    .h <- stanPopLinkSetup(ui, env$data, rxControl = control$rxControl)
    on.exit(stanLinkFree(), add = TRUE)
    if (.h$ntheta != sum(!.map$theta$fix)) {
      stop("the tier-0 problem's parameter count does not match the model ",
           "map (fix() thetas with literalFix=TRUE)", call. = FALSE)
    }
    .sf <- rxode2::rxWithSeed(control$seed, {
      .init <- if (identical(control$init, "ini")) {
        .stanInit(.map, .gen$blockSpecs, .nid, control)
      } else {
        control$init
      }
      rstan::sampling(.sm, data = .gen$data, chains = control$chains,
                      iter = control$iter, warmup = control$warmup,
                      thin = control$thin, seed = control$seed,
                      init = .init, cores = 1,
                      refresh = if (control$verbose) max(1L, control$iter %/% 10L) else 0L,
                      control = list(adapt_delta = control$adapt_delta,
                                     max_treedepth = control$max_treedepth))
    })
    .dx <- .stanDiagnostics(.sf, control)
    # the nlm objective (-2 log-likelihood) at the posterior point estimate,
    # evaluated while the link is still up
    .pf <- if (identical(control$point, "median")) stats::median else mean
    .thPt <- apply(rstan::extract(.sf, pars = "theta")$theta, 2, .pf)
    .popObj <- tryCatch(2 * .popEval(.thPt[!.map$theta$fix])$value,
                        error = function(e) NA_real_)
    stanLinkFree()
    return(.stanFinalizeEnvPop(.ret, ui, env, .sf, .map, .gen, .dx, control,
                               popObj = .popObj))
  }
  # covariate coefficients on a subject-CONSTANT covariate can come from the
  # scatter, so only they are excluded from the sensitivity requirement; a
  # time-varying coefficient needs the forward-sensitivity model
  .mrcConst <- .map$muRefCov$thetaIdx[!.cov$timeVarying]
  .needSens <- any(!.map$theta$fix &
                     !(seq_len(nrow(.map$theta)) %in%
                         c(.map$muRefIdx, .mrcConst)))
  .h <- stanLinkSetup(ui, env$data, likelihood = control$likelihood,
                      rxControl = control$rxControl,
                      thetaSens = .needSens, literalFix = control$literalFix,
                      cores = control$likCores)
  on.exit(stanLinkFree(), add = TRUE)
  on.exit(.Call(`_nlmixr2stan_clearThetaBase`), add = TRUE)
  if (!identical(.h$etaNames, .map$eta$name) ||
        .h$ntheta != nrow(.map$theta)) {
    stop("the linked problem's parameters do not match the model map",
         call. = FALSE) # nocov
  }
  .stanAssertThetaGradCover(.map, .h$thetaSensIdx)
  # tier-2 state: base parameter vector (omega tail fixed at link values) +
  # the mu-reference map for the theta-gradient assembly
  .Call(`_nlmixr2stan_setThetaBase`, as.double(.h$initPar))
  .Call(`_nlmixr2stan_setMuRef`, as.integer(.map$muRefIdx))
  # covariate-coefficient scatter: d/dtheta_p = cov_i * d/deta_k.  Only for
  # SUBJECT-CONSTANT coefficients the theta-sensitivity model does NOT
  # already carry -- when upstream classifies the coefficient as a plain
  # structural theta it gets an exact forward sensitivity, and adding the
  # (equal) scatter on top would double the column.  A TIME-VARYING
  # coefficient must be sensitivity-covered (the scatter identity does not
  # factor); with neither source its gradient would be silently zero, so
  # refuse.  (etaIdx is 1-based in the map, 0-based in the shim; covVal
  # columns follow muRefCov rows.)
  # (a fix()ed coefficient needs no gradient at all -- with literalFix=FALSE
  # it survives into the map, and refusing it would reject a valid model)
  .mrcFree <- !.map$theta$fix[.map$muRefCov$thetaIdx]
  .mrcTvBad <- which(.cov$timeVarying & .mrcFree &
                       !(.map$muRefCov$thetaIdx %in% .h$thetaSensIdx))
  if (length(.mrcTvBad) > 0L) {
    stop("time-varying mu-referenced covariate coefficient(s) ",
         paste0("'", .map$muRefCov$name[.mrcTvBad], "'", collapse = ", "),
         " have no forward sensitivity in the linked model, so their ",
         "gradient would be silently zero", call. = FALSE)
  }
  .mrcS <- which(!.cov$timeVarying & .mrcFree &
                   !(.map$muRefCov$thetaIdx %in% .h$thetaSensIdx))
  if (length(.mrcS) > 0L) {
    .Call(`_nlmixr2stan_setMuRefCov`,
          as.integer(.map$muRefCov$thetaIdx[.mrcS]),
          as.integer(.map$muRefCov$etaIdx[.mrcS] - 1L),
          .cov$val[, .mrcS, drop = FALSE])
  }
  # gradient conditioning: keep nlmixr2est's Omega^-1 commensurate with the
  # model's initial Omega (the conditional value is Omega-free)
  .omInv <- tryCatch(solve(ui$omega), error = function(e) NULL)
  if (!is.null(.omInv)) {
    tryCatch(.linkSetOmegaInv(.omInv), error = function(e) NULL)
  }
  # ---- sample -------------------------------------------------------------
  # the init jitter draws come from R's RNG; rxWithSeed scopes the seed and
  # restores the user's RNG state (touching .Random.seed directly is
  # forbidden by CRAN)
  .sf <- rxode2::rxWithSeed(control$seed, {
    .init <- if (identical(control$init, "ini")) {
      .stanInit(.map, .gen$blockSpecs, .nid, control)
    } else {
      control$init
    }
    rstan::sampling(.sm, data = .gen$data, chains = control$chains,
                    iter = control$iter, warmup = control$warmup,
                    thin = control$thin, seed = control$seed,
                    init = .init, cores = 1,
                    refresh = if (control$verbose) max(1L, control$iter %/% 10L) else 0L,
                    control = list(adapt_delta = control$adapt_delta,
                                   max_treedepth = control$max_treedepth))
  })
  .dx <- .stanDiagnostics(.sf, control)
  # Free the link BEFORE finalize: nlmixr2CreateOutputFromUi (and the
  # setOfv("FOCEi") row) run zero-iteration focei fits whose setup tears down
  # and rebuilds the process-global inner problem -- our link would be stale
  # underneath them either way, and nothing after sampling evaluates it.  The
  # on.exit registrations above stay as no-op safety.
  .Call(`_nlmixr2stan_clearThetaBase`)
  stanLinkFree()
  .stanFinalizeEnv(.ret, ui, env, .sf, .map, .gen, .dx, control)
}

#' Sampler diagnostics with actionable, loud messages
#' @noRd
.stanDiagnostics <- function(sf, control) {
  .sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  .nDiv <- sum(vapply(.sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  .nTree <- sum(vapply(.sp, function(x) {
    sum(x[, "treedepth__"] >= control$max_treedepth)
  }, numeric(1)))
  .sum <- rstan::summary(sf)$summary
  .keep <- !grepl("^(z_|etaP_|eta\\[|omegaOut|logLikSubj|lp__)",
                  rownames(.sum))
  .maxRhat <- suppressWarnings(max(.sum[.keep, "Rhat"], na.rm = TRUE))
  .minEss <- suppressWarnings(min(.sum[.keep, "n_eff"], na.rm = TRUE))
  .msg <- character(0)
  if (.nDiv > control$maxDivergent) {
    .msg <- c(.msg, paste0(.nDiv, " divergent transition(s) after warmup: ",
                           "the posterior is biased; raise adapt_delta (now ",
                           control$adapt_delta, ") toward 0.99"))
  }
  if (is.finite(.maxRhat) && .maxRhat > control$rhatMax) {
    .msg <- c(.msg, paste0("max Rhat ", signif(.maxRhat, 4), " > ",
                           control$rhatMax,
                           ": the chains have not mixed; do not use these estimates"))
  }
  if (is.finite(.minEss) && .minEss < control$essBulkMin) {
    .msg <- c(.msg, paste0("min ESS ", round(.minEss), " < ",
                           control$essBulkMin,
                           ": posterior summaries are not resolved; increase iter"))
  }
  if (.nTree > 0) {
    .msg <- c(.msg, paste0(.nTree, " transition(s) saturated max_treedepth (",
                           control$max_treedepth, "); efficiency, not validity"))
  }
  if (length(.msg) > 0L) {
    .txt <- paste(.msg, collapse = "\n")
    switch(control$onDiagnostic,
           error = stop(.txt, call. = FALSE),
           warn = warning(.txt, call. = FALSE),
           message = message(.txt),
           none = invisible())
  }
  list(nDivergent = .nDiv, nMaxTreedepth = .nTree, maxRhat = .maxRhat,
       minEss = .minEss, messages = .msg)
}

#' Posterior -> nlmixr2 fit for tier 0 (population-only, no etas)
#' @noRd
.stanFinalizeEnvPop <- function(ret, ui, env, sf, map, gen, dx, control,
                                popObj = NA_real_) {
  .pointFun <- if (identical(control$point, "median")) stats::median else mean
  .thDraw <- rstan::extract(sf, pars = "theta")$theta
  .fullTheta <- apply(.thDraw, 2, .pointFun)
  names(.fullTheta) <- map$theta$name
  .free <- !map$theta$fix
  .cov <- stats::cov(.thDraw[, .free, drop = FALSE])
  dimnames(.cov) <- list(map$theta$name[.free], map$theta$name[.free])
  .sum <- rstan::summary(sf)$summary
  .keep <- !grepl("^lp__", rownames(.sum))
  .posteriorSummary <- as.data.frame(.sum[.keep, , drop = FALSE])
  .ui <- rxode2::rxUiDecompress(ui)
  env2 <- ret
  env2$ui <- .ui
  env2$fullTheta <- .fullTheta
  env2$cov <- .cov
  env2$covMethod <- "stan.posterior"
  env2$objective <- popObj
  env2$adjObf <- FALSE
  env2$extra <- paste0(" (", control$chains, " chains x ",
                       control$iter - control$warmup, " draws; max Rhat ",
                       signif(dx$maxRhat, 4), "; min ESS ", round(dx$minEss),
                       "; population-only tier 0)")
  env2$method <- "Stan (HMC)"
  env2$est <- "stan"
  env2$ofvType <- "stan"
  env2$model <- .ui$ebe
  env2$message <- if (dx$nDivergent > 0) {
    paste0(dx$nDivergent, " divergent transitions")
  } else {
    ""
  }
  env2$stanfit <- sf
  env2$stanCode <- gen$code
  env2$stanData <- gen$data
  env2$stanDiagnostics <- dx
  env2$posteriorSummary <- .posteriorSummary
  nlmixr2est::.nlmixr2FitUpdateParams(env2)
  nlmixr2est::nmObjHandleControlObject(env2$stanControl, env2)
  .stanControlToFoceiControl(env2)
  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(env2$ui, data = env2$origData,
                                                control = env2$control,
                                                table = env2$table,
                                                env = env2, est = "stan")
  .env <- .fit$env
  .env$method <- "Stan (HMC)"
  .fit
}

#' Posterior -> nlmixr2 fit (the .nonmemFinalizeEnv template)
#' @noRd
.stanFinalizeEnv <- function(ret, ui, env, sf, map, gen, dx, control) {
  .ex <- rstan::extract(sf, pars = c("theta", "eta", "omegaOut", "logLikSubj"))
  .pointFun <- if (identical(control$point, "median")) stats::median else mean
  .thDraw <- .ex$theta                      # draws x ntheta
  .fullTheta <- apply(.thDraw, 2, .pointFun)
  names(.fullTheta) <- map$theta$name
  .omega <- apply(.ex$omegaOut, c(2, 3), .pointFun)
  dimnames(.omega) <- list(map$eta$name, map$eta$name)
  if (identical(control$point, "median")) {
    .chol <- tryCatch(chol(.omega), error = function(e) NULL)
    if (is.null(.chol)) {
      stop("the elementwise posterior-median omega is not positive definite;",
           " use point=\"mean\"", call. = FALSE) # nocov
    }
  }
  .etaMean <- apply(.ex$eta, c(2, 3), .pointFun) # nid x neta
  .nid <- nrow(.etaMean)
  .etaObf <- data.frame(ID = seq_len(.nid))
  for (.k in seq_len(ncol(.etaMean))) .etaObf[[map$eta$name[.k]]] <- .etaMean[, .k]
  .etaObf$OBJI <- -2 * colMeans(.ex$logLikSubj)
  .free <- !map$theta$fix
  .cov <- stats::cov(.thDraw[, .free, drop = FALSE])
  dimnames(.cov) <- list(map$theta$name[.free], map$theta$name[.free])
  # posterior summary on the monitored parameters (exact quantiles -- the
  # printed $parFixed CI is a Gaussian approximation from $cov)
  .sum <- rstan::summary(sf)$summary
  .keep <- !grepl("^(z_|eta\\[|omegaOut|logLikSubj)", rownames(.sum))
  .posteriorSummary <- as.data.frame(.sum[.keep, , drop = FALSE])

  # ---- the nlmixr2CreateOutputFromUi env contract -------------------------
  .ui <- rxode2::rxUiDecompress(ui)
  env2 <- ret
  env2$ui <- .ui
  env2$fullTheta <- .fullTheta
  env2$omega <- .omega
  env2$etaObf <- .etaObf
  env2$etaMat <- as.matrix(.etaObf[, map$eta$name, drop = FALSE])
  env2$cov <- .cov
  env2$covMethod <- "stan.posterior"
  env2$objective <- NA_real_
  env2$adjObf <- TRUE
  env2$extra <- paste0(" (", control$chains, " chains x ",
                       control$iter - control$warmup, " draws; max Rhat ",
                       signif(dx$maxRhat, 4), "; min ESS ", round(dx$minEss),
                       "; Wald CI from posterior cov)")
  env2$method <- "Stan (HMC)"
  env2$est <- "stan"
  env2$ofvType <- "stan"
  env2$theta <- data.frame(lower = map$theta$lower, theta = .fullTheta,
                           fixed = map$theta$fix, upper = map$theta$upper,
                           row.names = map$theta$name)
  env2$model <- .ui$ebe
  env2$message <- if (dx$nDivergent > 0) {
    paste0(dx$nDivergent, " divergent transitions")
  } else {
    ""
  }
  env2$stanfit <- sf
  env2$stanCode <- gen$code
  env2$stanData <- gen$data
  env2$stanDiagnostics <- dx
  env2$posteriorSummary <- .posteriorSummary
  nlmixr2est::.nlmixr2FitUpdateParams(env2)
  nlmixr2est::nmObjHandleControlObject(env2$stanControl, env2)
  .stanControlToFoceiControl(env2)
  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(env2$ui, data = env2$origData,
                                                control = env2$control,
                                                table = env2$table,
                                                env = env2, est = "stan")
  .env <- .fit$env
  .env$method <- "Stan (HMC)"
  # H3: the automatic FOCEi objective row (fired when CWRES exists) runs a
  # zero-iteration focei whose machinery overwrites $etaObf with the FOCEi
  # EBEs; keep those as $etaObfFocei and restore the posterior etas
  if (!identical(.env$etaObf, .etaObf)) {
    assign("etaObfFocei", .env$etaObf, envir = .env)
    assign("etaObf", .etaObf, envir = .env)
  }
  if (identical(control$ofv, "focei")) {
    # the FOCEi comparability row: a zero-iteration focei evaluation at the
    # posterior point estimate (runs with nlmixr2est's prior gate bypassed
    # internally).  Its machinery overwrites $etaObf with FOCEi EBEs; stash
    # those and restore the posterior etas (H3).
    .ok <- tryCatch({
      nlmixr2est::setOfv(.fit, "FOCEi")
      TRUE
    }, error = function(e) {
      cli::cli_warn(paste0("the FOCEi objective row could not be computed: ",
                           conditionMessage(e)))
      FALSE
    })
    if (.ok && !identical(.env$etaObf, .etaObf)) {
      assign("etaObfFocei", .env$etaObf, envir = .env)
      assign("etaObf", .etaObf, envir = .env)
    }
  }
  # ---- WAIC / LOO (D7), subject-level -------------------------------------
  # logLikSubj is the per-SUBJECT conditional log-likelihood, so these are
  # leave-one-SUBJECT-out quantities -- the natural cross-validation unit
  # for a hierarchical model (subjects are the exchangeable unit), and
  # conditional on the eta factorization (the only one the linkage exposes;
  # per-observation granularity is a possible upstream extension).
  if (requireNamespace("loo", quietly = TRUE)) {
    .llArr <- rstan::extract(sf, pars = "logLikSubj", permuted = FALSE)
    .rEff <- tryCatch(loo::relative_eff(exp(.llArr)), error = function(e) NULL)
    .loo <- tryCatch(suppressWarnings(loo::loo(.llArr, r_eff = .rEff)),
                     error = function(e) NULL)
    .waic <- tryCatch(suppressWarnings(loo::waic(.llArr)),
                      error = function(e) NULL)
    if (!is.null(.loo)) assign("loo", .loo, envir = .env)
    if (!is.null(.waic)) assign("waic", .waic, envir = .env)
    .odf <- get0("objDf", envir = .env)
    if (is.data.frame(.odf)) {
      .add <- function(odf, nm, elpd) {
        .row <- odf[1, , drop = FALSE]
        .row[1, ] <- NA
        if ("OBJF" %in% names(.row)) .row[1, "OBJF"] <- -2 * elpd
        if ("Log-likelihood" %in% names(.row)) {
          .row[1, "Log-likelihood"] <- elpd
        }
        rownames(.row) <- nm
        rbind(odf, .row)
      }
      if (!is.null(.waic)) {
        .odf <- .add(.odf, "WAIC (subject-level)",
                     .waic$estimates["elpd_waic", "Estimate"])
      }
      if (!is.null(.loo)) {
        .odf <- .add(.odf, "LOO (leave-one-subject-out)",
                     .loo$estimates["elpd_loo", "Estimate"])
      }
      assign("objDf", .odf, envir = .env)
    }
  }
  .fit
}

#' @export
print.nlmixr2stanCode <- function(x, ...) {
  cat(x$code, sep = "\n")
  invisible(x)
}
