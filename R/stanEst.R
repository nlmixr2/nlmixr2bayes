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
#' transform-both-sides lambdas (gen$tbsJac).  The lambda-dependent part of
#' log|dh(y; lambda)/dy| is always (lambda - 1) * sJ where sJ is the
#' Box-Cox statistic sum(log t) or the Yeo-Johnson statistic
#' sum(log1p(t), t >= 0) - sum(log1p(-t), t < 0) evaluated on the
#' TRANSFORM-BASE value t: the raw DV for a plain-normal base, the
#' logit/probit of the (bounded, FIXED-bounds) DV for the combined
#' "logit + yeoJohnson" / "probit + yeoJohnson" endpoints -- the inner
#' transform's own Jacobian term is lambda-free and constant, so it drops.
#' (The linked d/dlambda sensitivity handles these combined transforms
#' exactly -- FD-verified at 9e-11 on logitNorm + yeoJohnson.)
#' Observations are matched to the lambda's endpoint by CMT when the model
#' has several endpoints.
#' @noRd
.stanTbsJacValues <- function(gen, ui, dataSav) {
  if (length(gen$tbsJac) == 0L) return(gen)
  .obs <- dataSav[dataSav$EVID == 0 & !is.na(dataSav$DV), , drop = FALSE]
  .pd <- ui$predDf
  for (.j in gen$tbsJac) {
    .rows <- .obs
    .pw <- match(.j$condition, .pd$cond)
    if (nrow(.pd) > 1L && "CMT" %in% names(.obs)) {
      .rows <- .obs[.obs$CMT == .pd$cmt[.pw], , drop = FALSE]
    }
    .dv <- .rows$DV
    .tr <- as.character(.pd$transform[.pw])
    # the transform-base value t the lambda acts on
    .t <- if (.tr %in% c("boxCox", "yeoJohnson")) {
      .dv
    } else if (.tr %in% c("logit + yeoJohnson", "probit + yeoJohnson")) {
      .lo <- .pd$trLow[.pw]
      .hi <- .pd$trHi[.pw]
      if (any(.dv <= .lo | .dv >= .hi)) {
        stop("endpoint '", .j$condition, "' (", .tr, ") needs DV strictly ",
             "inside (", .lo, ", ", .hi, ") and has ",
             sum(.dv <= .lo | .dv >= .hi), " observation(s) outside",
             call. = FALSE)
      }
      .u <- (.dv - .lo) / (.hi - .lo)
      if (startsWith(.tr, "logit")) stats::qlogis(.u) else stats::qnorm(.u)
    } else {
      stop("estimated transform-both-sides lambda on endpoint '",
           .j$condition, "' with transform '", .tr, "' is not supported ",
           "by est=\"stan\" yet; fix() the lambda", call. = FALSE)
    }
    if (identical(.j$transform, "boxCox")) {
      if (any(.t <= 0)) {
        stop("boxCox(", .j$theta, ") needs strictly positive DV and ",
             "endpoint '", .j$condition, "' has ", sum(.t <= 0),
             " non-positive observation(s)", call. = FALSE)
      }
      .s <- sum(log(.t))
    } else {
      .s <- sum(log1p(.t[.t >= 0])) - sum(log1p(-.t[.t < 0]))
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

#' Does the loaded nlmixr2est provide the sampler iteration-print API
#' (foceiLikIterPrintStart/End + FOCEi table entry 8)?
#' @noRd
.stanHasIterPrint <- function() {
  is.function(get0("foceiLikIterPrintStart", envir = asNamespace("nlmixr2est")))
}

#' Initial display vector for the iteration print: natural-scale theta
#' estimates + the initial omega variances/covariances of every estimated
#' block + a trailing 0 placeholder for the raw-evaluation-count slot (the
#' C tick shim overwrites it every call; see stanGen.R's "nEval" slot), in
#' gen$dispNames order
#' @noRd
.stanDispInit <- function(map, gen) {
  .v <- map$theta$est
  for (.sp in gen$blockSpecs) {
    if (identical(.sp$type, "fixed")) next
    .m <- .sp$block$init
    for (.i in seq_len(nrow(.m))) {
      for (.j in seq_len(.i)) {
        .v <- c(.v, .m[.i, .j])
      }
    }
  }
  c(.v, 0)
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
  # finite mixtures (nlmixr2/nlmixr2est#955): supported for K = 2 when the
  # loaded nlmixr2est blesses the component-major batch layout (the nMix
  # table entry exists); the component-conditional rows are log-sum-exped
  # Stan-side with the mixing weight as an autodiff parameter
  .mixP <- tryCatch(ui$mixProbs, error = function(e) NULL)
  .nMixUi <- length(.mixP) + 1L
  if (.nMixUi > 1L) {
    if (identical(.Call(`_nlmixr2bayes_nMix`), -2L)) {
      stop("est=\"stan\" mixture support needs an nlmixr2est whose FOCEi C ",
           "API blesses the component-major layout ",
           "(nlmixr2/nlmixr2est#955); update nlmixr2est", call. = FALSE)
    }
    if (.nMixUi > 2L) {
      stop("est=\"stan\" supports 2-component mixtures for now (this ",
           "model has ", .nMixUi, " components)", call. = FALSE)
    }
  }
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
  .pri <- .stanUiPriors(ui)
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
#' @return an nlmixr2 fit (or an `nlmixr2bayesCode` when
#'   `stanControl(run=FALSE)`)
#' @keywords internal
#' @author Matthew L Fidler
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
# a declared non-normal random effect arrives here already expanded by
# nlmixr2est's pre-processing hook: the latent eta is a FIXED unit
# variance (which the generator already supports), the transform lives
# inside the linked nlmixr2_cond_all2() likelihood, and the copula
# correlation is an ordinary unbounded theta -- so nothing in the
# generated Stan program is special cased for it
attr(nlmixr2Est.stan, "etaDist") <- TRUE
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
  # nlmixr2/nlmixr2est#952 (the IOV magnitude theta's sensitivity column
  # scaled by a theta-dependent factor) was fixed together with an
  # iovXform pass-through on foceiLikLoad -- probe that formal, the same
  # structural-probe pattern as the other capability checks
  "iovXform" %in% names(formals(nlmixr2est::foceiLikLoad))
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
                 control = control, env = .env, dispNames = .gen$dispNames)
    class(.out) <- "nlmixr2bayesCode"
    return(.out)
  }
  rxode2::rxReq("rstan")
  # ---- compile (cached), then link --------------------------------------
  .sm <- stanCompile(.gen$code, cache = control$cache,
                     cacheDir = control$cacheDir, verbose = control$verbose)
  if (isTRUE(.gen$pop)) {
    # ---- tier 0: population-only (no etas) via the nlm C API ------------
    .h <- stanPopLinkSetup(ui, env$data, rxControl = control$rxControl,
                           cores = control$cores, print = control$print)
    on.exit(stanLinkFree(), add = TRUE)
    if (.h$ntheta != sum(!.map$theta$fix)) {
      stop("the tier-0 problem's parameter count does not match the model ",
           "map (fix() thetas with literalFix=TRUE)", call. = FALSE)
    }
    .sf <- .stanRunInference(.sm, .gen, .map, .nid, control)
    .dx <- .stanDiagnostics(.sf, control)
    # the nlm objective (-2 log-likelihood) at the posterior point estimate,
    # evaluated while the link is still up
    .pf <- if (identical(control$point, "median")) stats::median else mean
    .thPt <- apply(.stanExtract(.sf, pars = "theta")$theta, 2, .pf)
    .popObj <- tryCatch(2 * .popEval(.thPt[!.map$theta$fix])$value,
                        error = function(e) NA_real_)
    if (control$print > 0L) {
      .ph <- tryCatch(nlmixr2est::nlmGetParHist(TRUE),
                      error = function(e) NULL)
      if (is.data.frame(.ph)) .ret$parHistData <- .ph
    }
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
                      cores = control$cores,
                      maxOdeRecalc = control$maxOdeRecalc,
                      fallbackFD = control$fallbackFD)
  on.exit(stanLinkFree(), add = TRUE)
  on.exit(.Call(`_nlmixr2bayes_clearThetaBase`), add = TRUE)
  if (!identical(.h$etaNames, .map$eta$name) ||
        .h$ntheta != nrow(.map$theta)) {
    stop("the linked problem's parameters do not match the model map",
         call. = FALSE) # nocov
  }
  .stanAssertThetaGradCover(.map, .h$thetaSensIdx)
  # tier-2 state: base parameter vector (omega tail fixed at link values) +
  # the mu-reference map for the theta-gradient assembly
  .Call(`_nlmixr2bayes_setThetaBase`, as.double(.h$initPar))
  .Call(`_nlmixr2bayes_setMuRef`, as.integer(.map$muRefIdx))
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
    .cv <- .cov$val[, .mrcS, drop = FALSE]
    if (.map$nMix > 1L) {
      # component-major expanded rows share the physical subject's covariate
      .cv <- do.call(rbind, rep(list(.cv), .map$nMix))
    }
    .Call(`_nlmixr2bayes_setMuRefCov`,
          as.integer(.map$muRefCov$thetaIdx[.mrcS]),
          as.integer(.map$muRefCov$etaIdx[.mrcS] - 1L),
          .cv)
  }
  if (.map$nMix > 1L) {
    .nm <- .Call(`_nlmixr2bayes_nMix`)
    if (!identical(.nm, .map$nMix)) {
      stop("the linked problem reports ", .nm, " mixture component(s) but ",
           "the model map expects ", .map$nMix, call. = FALSE) # nocov
    }
  }
  # gradient conditioning: keep nlmixr2est's Omega^-1 commensurate with the
  # model's initial Omega (the conditional value is Omega-free)
  .omInv <- tryCatch(solve(ui$omega), error = function(e) NULL)
  if (!is.null(.omInv)) {
    tryCatch(.linkSetOmegaInv(.omInv), error = function(e) NULL)
  }
  # ---- nlmixr2est-style iteration print + parameter history --------------
  # (the scale.h residency in nlmixr2est; the compiled model's tick calls
  # entry 8 of the FOCEi C table every log-density evaluation, gated to
  # control$print)
  .iterOn <- control$print > 0L && .stanHasIterPrint() &&
    control$chainCores <= 1L
  if (control$print > 0L && .stanHasIterPrint() &&
        control$chainCores > 1L) {
    # parallel chains evaluate in child processes; the print ticks happen in
    # worker-local memory and never reach this process's history
    message("iteration printing needs sequential chains; it is disabled ",
            "with parallel chains (set chainCores = 1 for the iteration ",
            "table)")
  }
  if (control$print > 0L && !.stanHasIterPrint()) {
    warning("this nlmixr2est does not provide the iteration-print API; ",
            "update nlmixr2est for the familiar iteration table",
            call. = FALSE)
  }
  .ph <- NULL
  if (.iterOn) {
    .Call(`_nlmixr2bayes_resetEvalCount`)
    nlmixr2est::foceiLikIterPrintStart(control$print,
                                       .stanDispInit(.map, .gen),
                                       .gen$dispNames)
    on.exit(try(nlmixr2est::foceiLikIterPrintEnd(), silent = TRUE),
            add = TRUE)
  }
  # ---- sample -------------------------------------------------------------
  # the init jitter draws come from R's RNG; rxWithSeed scopes the seed and
  # restores the user's RNG state (touching .Random.seed directly is
  # forbidden by CRAN)
  .usePsock <- identical(.Platform$OS.type, "windows") &&
    identical(control$algorithm, "NUTS") &&
    control$chainCores > 1L
  if (.usePsock) {
    # Windows PSOCK workers must build their OWN link state; free the
    # parent's preflight link before the workers start.
    .Call(`_nlmixr2bayes_clearThetaBase`)
    stanLinkFree()
    .sf <- .stanRunInferencePsock(.sm, .gen, .map, .cov, .needSens, .nid,
                                  control, ui, env$data)
  } else {
    .sf <- .stanRunInference(.sm, .gen, .map, .nid, control)
  }
  if (.iterOn) {
    .ph <- tryCatch(nlmixr2est::foceiLikIterPrintEnd(),
                    error = function(e) NULL)
    if (is.data.frame(.ph)) .ret$parHistData <- .ph
  }
  .dx <- .stanDiagnostics(.sf, control)
  if (inherits(.sf, "nlmixr2bayesPathfinder")) {
    # the generated-quantities pass (rstan::gqs -> nlmixr2_cond_all2)
    # needs the linked likelihood, so it must run BEFORE the teardown
    # below; NUTS/ADVI evaluate their GQ inside rstan while sampling
    .pathfinderPrefetch(.sf)
  }
  # Free the link BEFORE finalize: nlmixr2CreateOutputFromUi (and the
  # setOfv("FOCEi") row) run zero-iteration focei fits whose setup tears down
  # and rebuilds the process-global inner problem -- our link would be stale
  # underneath them either way, and nothing after sampling evaluates it.  The
  # on.exit registrations above stay as no-op safety.
  .Call(`_nlmixr2bayes_clearThetaBase`)
  stanLinkFree()
  .stanFinalizeEnv(.ret, ui, env, .sf, .map, .gen, .dx, control)
}


#' Run the selected Stan inference algorithm on the compiled model
#'
#' NUTS goes through [rstan::sampling()] (unix: chain parallelism by fork;
#' Windows: one chain per PSOCK worker with one linked likelihood per
#' worker); the ADVI variants go through
#' [rstan::vb()] with importance resampling, whose draws then flow through
#' the SAME posterior->fit machinery.  Both return a stanfit.
#' @noRd
.stanInitForControl <- function(control, map, blockSpecs, nid) {
  if (identical(control$init, "ini")) {
    .stanInit(map, blockSpecs, nid, control)
  } else {
    control$init
  }
}

#' Per-chain init payload for one-chain rstan::sampling()
#' @noRd
.stanChainInit <- function(init, chainId, chains) {
  if (!is.list(init) || chains <= 1L) return(init)
  if (length(init) >= chainId && is.list(init[[chainId]])) return(init[[chainId]])
  init
}

#' Build the linked-likelihood state for one process
#' @noRd
.stanLinkSetupForRun <- function(ui, data, map, cov, needSens, control) {
  .h <- stanLinkSetup(ui, data, likelihood = control$likelihood,
                      rxControl = control$rxControl,
                      thetaSens = needSens, literalFix = control$literalFix,
                      cores = control$cores,
                      maxOdeRecalc = control$maxOdeRecalc,
                      fallbackFD = control$fallbackFD)
  if (!identical(.h$etaNames, map$eta$name) || .h$ntheta != nrow(map$theta)) {
    stop("the linked problem's parameters do not match the model map",
         call. = FALSE) # nocov
  }
  .stanAssertThetaGradCover(map, .h$thetaSensIdx)
  .Call(`_nlmixr2bayes_setThetaBase`, as.double(.h$initPar))
  .Call(`_nlmixr2bayes_setMuRef`, as.integer(map$muRefIdx))
  .mrcFree <- !map$theta$fix[map$muRefCov$thetaIdx]
  .mrcTvBad <- which(cov$timeVarying & .mrcFree &
                       !(map$muRefCov$thetaIdx %in% .h$thetaSensIdx))
  if (length(.mrcTvBad) > 0L) {
    stop("time-varying mu-referenced covariate coefficient(s) ",
         paste0("'", map$muRefCov$name[.mrcTvBad], "'", collapse = ", "),
         " have no forward sensitivity in the linked model, so their ",
         "gradient would be silently zero", call. = FALSE)
  }
  .mrcS <- which(!cov$timeVarying & .mrcFree &
                   !(map$muRefCov$thetaIdx %in% .h$thetaSensIdx))
  if (length(.mrcS) > 0L) {
    .cv <- cov$val[, .mrcS, drop = FALSE]
    if (map$nMix > 1L) .cv <- do.call(rbind, rep(list(.cv), map$nMix))
    .Call(`_nlmixr2bayes_setMuRefCov`,
          as.integer(map$muRefCov$thetaIdx[.mrcS]),
          as.integer(map$muRefCov$etaIdx[.mrcS] - 1L),
          .cv)
  }
  if (map$nMix > 1L) {
    .nm <- .Call(`_nlmixr2bayes_nMix`)
    if (!identical(.nm, map$nMix)) {
      stop("the linked problem reports ", .nm, " mixture component(s) but ",
           "the model map expects ", map$nMix, call. = FALSE) # nocov
    }
  }
  .omInv <- tryCatch(solve(ui$omega), error = function(e) NULL)
  if (!is.null(.omInv)) {
    tryCatch(.linkSetOmegaInv(.omInv), error = function(e) NULL)
  }
  invisible(.h)
}

#' Combine one-chain stanfit objects into a multi-chain stanfit
#' @noRd
.stanCombineSflist <- function(sflist) {
  .sfFun <- get0("sflist2stanfit", envir = asNamespace("rstan"),
                 inherits = FALSE)
  if (is.null(.sfFun)) {
    stop("this rstan does not provide sflist2stanfit; update rstan",
         call. = FALSE)
  }
  .sfFun(sflist)
}

#' Windows PSOCK chain parallelism for NUTS
#' @noRd
.stanRunInferencePsock <- function(sm, gen, map, cov, needSens, nid, control,
                                   ui, data) {
  .init <- .stanInitForControl(control, map, gen$blockSpecs, nid)
  .chainIds <- as.list(seq_len(control$chains))
  .cl <- parallel::makePSOCKcluster(control$chainCores)
  on.exit(parallel::stopCluster(.cl), add = TRUE)
  parallel::clusterEvalQ(.cl, {
    requireNamespace("nlmixr2bayes", quietly = TRUE)
    requireNamespace("rstan", quietly = TRUE)
    requireNamespace("rxode2", quietly = TRUE)
    requireNamespace("nlmixr2est", quietly = TRUE)
    NULL
  })
  .worker <- local({
    .ui <- ui
    .data <- data
    .map <- map
    .cov <- cov
    .needSens <- needSens
    .control <- control
    .gen <- gen
    .init <- .init
    .sm <- sm
    function(.cid) {
      nlmixr2bayes:::.stanLinkSetupForRun(.ui, .data, .map, .cov,
                                          .needSens, .control)
      on.exit(nlmixr2bayes:::stanLinkFree(), add = TRUE)
      on.exit(.Call(`_nlmixr2bayes_clearThetaBase`), add = TRUE)
      .i1 <- nlmixr2bayes:::.stanChainInit(.init, .cid, .control$chains)
      rstan::sampling(.sm, data = .gen$data, chains = 1L, chain_id = .cid,
                      iter = .control$iter, warmup = .control$warmup,
                      thin = .control$thin, seed = .control$seed,
                      init = list(.i1), cores = 1L,
                      refresh = if (.control$verbose) {
                        max(1L, .control$iter %/% 10L)
                      } else {
                        0L
                      },
                      control = list(adapt_delta = .control$adapt_delta,
                                     max_treedepth = .control$max_treedepth))
    }
  })
  .sfl <- parallel::parLapply(.cl, .chainIds, .worker)
  .stanCombineSflist(.sfl)
}

.stanRunInference <- function(sm, gen, map, nid, control) {
  rxode2::rxWithSeed(control$seed, {
    .init <- .stanInitForControl(control, map, gen$blockSpecs, nid)
    if (identical(control$algorithm, "pathfinder")) {
      .stanRunPathfinder(sm, gen, map, nid, control, .init)
    } else if (identical(control$algorithm, "NUTS")) {
      rstan::sampling(sm, data = gen$data, chains = control$chains,
                      iter = control$iter, warmup = control$warmup,
                      thin = control$thin, seed = control$seed,
                      init = .init, cores = control$chainCores,
                      refresh = if (control$verbose) {
                        max(1L, control$iter %/% 10L)
                      } else {
                        0L
                      },
                      control = list(adapt_delta = control$adapt_delta,
                                     max_treedepth = control$max_treedepth))
    } else {
      # one variational run; init is a single list -- chain 1 of the "ini"
      # scheme (exactly the ini() estimates, no jitter)
      .i1 <- if (is.list(.init) && length(.init) >= 1L &&
                   is.list(.init[[1L]])) {
        .init[[1L]]
      } else {
        .init
      }
      rstan::vb(sm, data = gen$data, algorithm = control$algorithm,
                iter = control$vbIter, tol_rel_obj = control$vbTolRelObj,
                output_samples = control$vbOutputSamples,
                seed = control$seed, init = .i1,
                importance_resampling = TRUE,
                refresh = if (control$verbose) 100L else 0L)
    }
  })
}

#' Method label for the fit ($method) by algorithm
#' @noRd
.stanMethodLabel <- function(control) {
  switch(control$algorithm,
         NUTS = "Stan (HMC)",
         pathfinder = "Stan (Pathfinder)",
         paste0("Stan (ADVI ", control$algorithm, ")"))
}

#' Is real Pathfinder available?  rstan does not expose the Pathfinder
#' service (StanHeaders 2.32 predates it), and CmdStan cannot work here at
#' all (a separate executable cannot reach the in-process linked
#' likelihood), so Pathfinder runs through StanEstimators' in-process
#' callbacks against the compiled model's log density + analytic gradient.
#' @noRd
.stanHasPathfinder <- function() {
  requireNamespace("StanEstimators", quietly = TRUE)
}

#' Sampler diagnostics with actionable, loud messages
#' @noRd
.stanDiagnostics <- function(sf, control) {
  if (identical(control$algorithm, "pathfinder")) {
    return(.stanDiagnosticsPathfinder(sf, control))
  }
  if (!identical(control$algorithm, "NUTS")) {
    return(.stanDiagnosticsVb(sf, control))
  }
  .sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  .nDiv <- sum(vapply(.sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  .nTree <- sum(vapply(.sp, function(x) {
    sum(x[, "treedepth__"] >= control$max_treedepth)
  }, numeric(1)))
  .sum <- .stanSummaryDf(sf)
  .keep <- !grepl("^(z_|etaP_|eta\\[|omegaOut|logLikSubj|mixProbOut|lp__)",
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
       minEss = .minEss, khat = NA_real_, messages = .msg)
}

#' ADVI diagnostics: the Pareto-k of the importance ratios replaces
#' Rhat/ESS/divergences (there are no chains).  khat <= 0.7 means the
#' variational approximation supports importance correction; above that the
#' approximation is unreliable and NUTS should be used.
#' @noRd
.stanDiagnosticsVb <- function(sf, control) {
  # rstan::vb(importance_resampling=TRUE) stores the PSIS result at
  # sim$diagnostics: $psis$pareto_k is the overall khat; the (unnamed)
  # first element holds the raw log_p__/log_g__ draws as a fallback
  .khat <- tryCatch(as.numeric(sf@sim$diagnostics$psis$pareto_k),
                    error = function(e) NA_real_)
  if (length(.khat) != 1L || !is.finite(.khat)) {
    .khat <- tryCatch({
      .d <- sf@sim$diagnostics[[1L]]
      if (!is.null(.d$log_p__) && !is.null(.d$log_g__) &&
            requireNamespace("loo", quietly = TRUE)) {
        .lw <- .d$log_p__ - .d$log_g__
        .lw <- .lw[is.finite(.lw)]
        suppressWarnings(loo::psis(.lw, r_eff = NA)$diagnostics$pareto_k)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
  }
  .msg <- character(0)
  if (is.finite(.khat) && .khat > 0.7) {
    .msg <- c(.msg, paste0("ADVI Pareto khat ", signif(.khat, 3), " > 0.7: ",
                           "the variational approximation is unreliable; ",
                           "use algorithm=\"NUTS\""))
  }
  if (length(.msg) > 0L) {
    .txt <- paste(.msg, collapse = "\n")
    switch(control$onDiagnostic,
           error = stop(.txt, call. = FALSE),
           warn = warning(.txt, call. = FALSE),
           message = message(.txt),
           none = invisible())
  }
  list(nDivergent = 0L, nMaxTreedepth = 0L, maxRhat = NA_real_,
       minEss = NA_real_, khat = .khat, messages = .msg)
}

#' Pathfinder diagnostics: khat of the pooled importance ratios (same
#' gate as ADVI) plus how many of the requested paths yielded a usable
#' local Gaussian
#' @noRd
.stanDiagnosticsPathfinder <- function(sf, control) {
  .khat <- sf$khat
  .msg <- character(0)
  if (is.finite(.khat) && .khat > 0.7) {
    .msg <- c(.msg, paste0("Pathfinder Pareto khat ", signif(.khat, 3),
                           " > 0.7: the approximation is unreliable; ",
                           "use algorithm=\"NUTS\""))
  }
  if (sf$nPathsOk < control$pathfinderPaths) {
    .msg <- c(.msg, paste0(control$pathfinderPaths - sf$nPathsOk, " of ",
                           control$pathfinderPaths,
                           " Pathfinder path(s) failed to produce a ",
                           "usable local Gaussian"))
  }
  if (length(.msg) > 0L) {
    .txt <- paste(.msg, collapse = "\n")
    switch(control$onDiagnostic,
           error = stop(.txt, call. = FALSE),
           warn = warning(.txt, call. = FALSE),
           message = message(.txt),
           none = invisible())
  }
  list(nDivergent = 0L, nMaxTreedepth = 0L, maxRhat = NA_real_,
       minEss = NA_real_, khat = .khat, messages = .msg)
}

#' Relabel a posterior-summary row-name vector from the Stan-mangled
#' identifiers the generated program declares (dots/other punctuation in an
#' nlmixr2 parameter name are not legal in a Stan identifier, so
#' `.stanParName()` replaces them with `_`, e.g. `add.sd` -> `add_sd`) back
#' to the original nlmixr2 model names, so `$posteriorSummary` reads like the
#' model rather than like its Stan translation.
#'
#' Thetas (including error-model and covariate-coefficient parameters, which
#' are ordinary theta rows) have an exact 1:1 Stan-name -> model-name
#' mapping already carried on `map$theta` -- those are always translated.
#' Omega is different: what a user recognizes as "the model's omega" is the
#' natural variance/covariance entries of `omegaOut` (an neta x neta
#' generated quantity assembled from whatever internal SD/correlation/
#' covariance parameterization the block actually samples), so on the eta
#' side this relabels `omegaOut[i,j]` cells -- using the SAME om./cov.
#' convention as the iteration-print table -- rather than trying to name the
#' internal per-block sampling parameters (a Cholesky factor entry has no
#' single-parameter meaning to translate).
#' @noRd
.stanPosteriorRowLabels <- function(rn, map, blocks = NULL) {
  .lbl <- stats::setNames(map$theta$name, map$theta$par)
  if (!is.null(blocks) && length(blocks) > 0L) {
    for (.blk in blocks) {
      .mem <- .blk$members
      .k <- length(.mem)
      for (.di in seq_len(.k)) {
        for (.dj in seq_len(.di)) {
          .nm <- paste0("omegaOut[", .blk$start - 1L + .di, ",",
                        .blk$start - 1L + .dj, "]")
          .lbl[.nm] <- if (.di == .dj) {
            paste0("om.", .mem[.di])
          } else {
            paste0("cov.", .mem[.di], ".", .mem[.dj])
          }
          # omegaOut is symmetric; the mirrored cell carries the same value
          .lbl[paste0("omegaOut[", .blk$start - 1L + .dj, ",",
                     .blk$start - 1L + .di, "]")] <- .lbl[.nm]
        }
      }
    }
  }
  .m <- unname(.lbl[rn])
  ifelse(is.na(.m), rn, .m)
}

#' The internal per-block omega sampling parameters (sd<id>/Lcorr<id>/
#' omega<id>/L<id>/z<id>/etaP<id> -- see stanGen.R's `.stanBlockId()`), as
#' the base (unindexed) Stan names `.stanFinalizeEnv` declares them under.
#' These have no clean single-parameter meaning to translate (a Cholesky
#' factor entry, a per-subject random-effect draw) and are redundant with
#' the derived, relabeled `omegaOut` cells, so posterior summaries drop them
#' by exact name match here rather than a generic prefix regex (a
#' `_b[0-9]+`-style pattern breaks once block ids are member-name-derived,
#' and a plain prefix match risks swallowing a theta that happens to start
#' the same way, e.g. an error parameter named "sd.something").
#' @noRd
.stanBlockInternalNames <- function(blocks) {
  unlist(lapply(blocks, function(.blk) {
    .id <- .stanBlockId(.blk$members)
    paste0(c("sd", "Lcorr", "omega", "L", "z", "etaP"), .id)
  }), use.names = FALSE)
}

#' Put the control(s) on the fit env
#'
#' The engine always runs a `stanControl`, but an `est="nuts"`/`"advi"`/
#' `"pathfinder"` run also carries the sugar control the user wrote
#' ([nutsControl()]/[adviControl()]/[pathfinderControl()]) -- that is what
#' `nmObjGetControl` and `rxUiDeparse` should report back, the way an
#' `est="foce"` fit reports a `foceControl` rather than the `foceiControl`
#' its engine ran.
#' @param env the estimation environment (holds `stanEstName` + the sugar
#'   control)
#' @param env2 the fit environment
#' @return nothing, called for the side effect
#' @noRd
.stanHandleControlObjects <- function(env, env2) {
  .est <- env$stanEstName
  if (!is.null(.est)) {
    .sugar <- env[[paste0(.est, "Control")]]
    if (!is.null(.sugar)) {
      nlmixr2est::nmObjHandleControlObject(.sugar, env2)
    }
  }
  nlmixr2est::nmObjHandleControlObject(env2$stanControl, env2)
  invisible()
}

#' Posterior -> nlmixr2 fit for tier 0 (population-only, no etas)
#' @noRd
.stanFinalizeEnvPop <- function(ret, ui, env, sf, map, gen, dx, control,
                                popObj = NA_real_) {
  # est="nuts"/"advi"/"pathfinder" sugar records its own name on the fit
  .estName <- if (is.null(env$stanEstName)) "stan" else env$stanEstName
  .pointFun <- if (identical(control$point, "median")) stats::median else mean
  .thDraw <- .stanExtract(sf, pars = "theta")$theta
  .fullTheta <- apply(.thDraw, 2, .pointFun)
  names(.fullTheta) <- map$theta$name
  .free <- !map$theta$fix
  .cov <- stats::cov(.thDraw[, .free, drop = FALSE])
  dimnames(.cov) <- list(map$theta$name[.free], map$theta$name[.free])
  .sum <- .stanSummaryDf(sf)
  # theta[] duplicates the individually-named theta parameters above; drop it
  .keep <- !grepl("^(lp__|theta\\[)", rownames(.sum))
  .posteriorSummary <- as.data.frame(.sum[.keep, , drop = FALSE])
  rownames(.posteriorSummary) <- .stanPosteriorRowLabels(
    rownames(.posteriorSummary), map)
  .ui <- rxode2::rxUiDecompress(ui)
  env2 <- ret
  env2$ui <- .ui
  env2$fullTheta <- .fullTheta
  env2$cov <- .cov
  env2$covMethod <- "stan.posterior"
  env2$objective <- popObj
  env2$adjObf <- FALSE
  env2$extra <- if (identical(control$algorithm, "NUTS")) {
    paste0(" (", control$chains, " chains x ",
           control$iter - control$warmup, " draws; max Rhat ",
           signif(dx$maxRhat, 4), "; min ESS ", round(dx$minEss),
           "; population-only tier 0)")
  } else if (identical(control$algorithm, "pathfinder")) {
    paste0(" (Pathfinder, ", control$pathfinderPaths, " paths, ",
           control$vbOutputSamples, " draws; Pareto khat ",
           signif(dx$khat, 3), "; population-only tier 0)")
  } else {
    paste0(" (ADVI ", control$algorithm, ", ", control$vbOutputSamples,
           " draws; Pareto khat ", signif(dx$khat, 3),
           "; population-only tier 0)")
  }
  env2$method <- .stanMethodLabel(control)
  env2$est <- .estName
  env2$ofvType <- "stan"
  env2$model <- .ui$ebe
  env2$message <- if (dx$nDivergent > 0) {
    paste0(dx$nDivergent, " divergent transitions")
  } else {
    ""
  }
  if (map$nMix > 1L) {
    .mixProb <- apply(.stanExtract(sf, pars = "mixProbOut")$mixProbOut,
                      c(2, 3), .pointFun)
    dimnames(.mixProb) <- list(NULL, paste0("mix", seq_len(map$nMix)))
    env2$mixProb <- data.frame(ID = seq_len(nrow(.mixProb)), .mixProb)
  }
  env2$stanfit <- sf
  env2$stanCode <- gen$code
  env2$stanData <- gen$data
  env2$stanDiagnostics <- dx
  env2$posteriorSummary <- .posteriorSummary
  nlmixr2est::.nlmixr2FitUpdateParams(env2)
  .stanHandleControlObjects(env, env2)
  .stanControlToFoceiControl(env2)
  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(env2$ui, data = env2$origData,
                                                control = env2$control,
                                                table = env2$table,
                                                env = env2, est = .estName)
  .env <- .fit$env
  .env$method <- .stanMethodLabel(control)
  .fit
}

#' Posterior -> nlmixr2 fit (the .nonmemFinalizeEnv template)
#' @noRd
.stanFinalizeEnv <- function(ret, ui, env, sf, map, gen, dx, control) {
  # est="nuts"/"advi"/"pathfinder" sugar records its own name on the fit
  .estName <- if (is.null(env$stanEstName)) "stan" else env$stanEstName
  .ex <- .stanExtract(sf, pars = c("theta", "eta", "omegaOut", "logLikSubj"))
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
  if (map$nMix > 1L) {
    # component-specific etas: report the MEMBERSHIP-WEIGHTED eta per
    # subject (draw-wise, then the point summary), and keep the posterior
    # membership probabilities on the fit
    .mp <- .stanExtract(sf, pars = "mixProbOut")$mixProbOut # draws x N x K
    .nid <- dim(.mp)[2]
    .neta <- dim(.ex$eta)[3]
    .etaW <- array(0, c(dim(.ex$eta)[1], .nid, .neta))
    for (.k in seq_len(map$nMix)) {
      .rows <- seq_len(.nid) + (.k - 1L) * .nid
      for (.j in seq_len(.neta)) {
        .etaW[, , .j] <- .etaW[, , .j] + .mp[, , .k] * .ex$eta[, .rows, .j]
      }
    }
    .etaMean <- apply(.etaW, c(2, 3), .pointFun)
  } else {
    .etaMean <- apply(.ex$eta, c(2, 3), .pointFun) # nid x neta
    .nid <- nrow(.etaMean)
  }
  .etaObf <- data.frame(ID = seq_len(.nid))
  for (.k in seq_len(ncol(.etaMean))) .etaObf[[map$eta$name[.k]]] <- .etaMean[, .k]
  .etaObf$OBJI <- -2 * colMeans(.ex$logLikSubj)
  .free <- !map$theta$fix
  .cov <- stats::cov(.thDraw[, .free, drop = FALSE])
  dimnames(.cov) <- list(map$theta$name[.free], map$theta$name[.free])
  # posterior summary on the monitored parameters (exact quantiles -- the
  # printed $parFixed CI is a Gaussian approximation from $cov).  The raw
  # per-block omega sampling parameterization (sd<id>/Lcorr<id>/omega<id>/
  # L<id>/z<id>/etaP<id>) is dropped by exact name (`.stanBlockInternalNames`)
  # in favor of the derived, natural-scale omegaOut cells, which get
  # relabeled to om./cov. names below -- that is what a user recognizes as
  # "the model's omega", not whichever internal parameterization a given
  # block happens to sample from.
  .sum <- .stanSummaryDf(sf)
  .rn <- rownames(.sum)
  .rnBase <- sub("\\[.*\\]$", "", .rn)
  .keep <- !(.rnBase %in% .stanBlockInternalNames(map$blocks)) &
    !grepl("^(eta\\[|logLikSubj|mixProbOut|theta\\[)", .rn)
  .posteriorSummary <- as.data.frame(.sum[.keep, , drop = FALSE])
  rownames(.posteriorSummary) <- .stanPosteriorRowLabels(
    rownames(.posteriorSummary), map, map$blocks)
  # cross-block omegaOut cells are always exactly 0 (no modeled correlation
  # across blocks) and have no om./cov. label -- drop rather than show noise
  .posteriorSummary <- .posteriorSummary[
    !grepl("^omegaOut\\[", rownames(.posteriorSummary)), , drop = FALSE]

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
  env2$extra <- if (identical(control$algorithm, "NUTS")) {
    paste0(" (", control$chains, " chains x ",
           control$iter - control$warmup, " draws; max Rhat ",
           signif(dx$maxRhat, 4), "; min ESS ", round(dx$minEss),
           "; Wald CI from posterior cov)")
  } else if (identical(control$algorithm, "pathfinder")) {
    paste0(" (Pathfinder, ", control$pathfinderPaths, " paths, ",
           control$vbOutputSamples, " draws; Pareto khat ",
           signif(dx$khat, 3), "; Wald CI from posterior cov)")
  } else {
    paste0(" (ADVI ", control$algorithm, ", ", control$vbOutputSamples,
           " draws; Pareto khat ", signif(dx$khat, 3),
           "; Wald CI from posterior cov)")
  }
  env2$method <- .stanMethodLabel(control)
  env2$est <- .estName
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
  if (map$nMix > 1L) {
    .mixProb <- apply(.stanExtract(sf, pars = "mixProbOut")$mixProbOut,
                      c(2, 3), .pointFun)
    dimnames(.mixProb) <- list(NULL, paste0("mix", seq_len(map$nMix)))
    env2$mixProb <- data.frame(ID = seq_len(nrow(.mixProb)), .mixProb)
  }
  env2$stanfit <- sf
  env2$stanCode <- gen$code
  env2$stanData <- gen$data
  env2$stanDiagnostics <- dx
  env2$posteriorSummary <- .posteriorSummary
  nlmixr2est::.nlmixr2FitUpdateParams(env2)
  .stanHandleControlObjects(env, env2)
  .stanControlToFoceiControl(env2)
  .fit <- nlmixr2est::nlmixr2CreateOutputFromUi(env2$ui, data = env2$origData,
                                                control = env2$control,
                                                table = env2$table,
                                                env = env2, est = .estName)
  .env <- .fit$env
  .env$method <- .stanMethodLabel(control)
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
    if (inherits(sf, "nlmixr2bayesPathfinder")) {
      # resampled draws have no chain structure; r_eff stays NULL
      .ll <- .stanExtract(sf, pars = "logLikSubj")$logLikSubj
      .llArr <- array(.ll, dim = c(nrow(.ll), 1L, ncol(.ll)))
      .rEff <- NULL
    } else {
      .llArr <- rstan::extract(sf, pars = "logLikSubj", permuted = FALSE)
      .rEff <- tryCatch(loo::relative_eff(exp(.llArr)),
                        error = function(e) NULL)
    }
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

#' @author Matthew L Fidler
#' @export
print.nlmixr2bayesCode <- function(x, ...) {
  cat(x$code, sep = "\n")
  invisible(x)
}
