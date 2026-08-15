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

#' Capability assertions + prior triage for est="stan"
#' @noRd
.stanAssert <- function(ui) {
  rxode2::assertRxUiPrediction(ui, " for est=\"stan\"", .var.name = ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(ui, " for est=\"stan\"", .var.name = ui$modelName)
  rxode2::assertRxUiNoMix(ui, " for est=\"stan\"", .var.name = ui$modelName)
  rxode2::assertRxUiMixedOnly(ui, " for est=\"stan\"", .var.name = ui$modelName)
  # (rxode2 >= 5.1.7 always has llik support, so no assertRxUiTransformNormal)
  # transform-both-sides with an ESTIMATED lambda: the conditional value the
  # linkage returns omits the DV-transform Jacobian, which is not constant in
  # an estimated lambda -- the target would be silently wrong
  .iniDf <- ui$iniDf
  .tbs <- which(.iniDf$err %in% c("boxCox", "yeoJohnson") & !.iniDf$fix)
  if (length(.tbs) > 0L) {
    stop("est=\"stan\" cannot estimate the transform-both-sides parameter(s) ",
         paste0("'", .iniDf$name[.tbs], "'", collapse = ", "),
         ": the linked conditional likelihood omits the DV-transform ",
         "Jacobian, which is not constant in an estimated lambda; fix the ",
         "lambda or drop the transform", call. = FALSE)
  }
  # mu-referenced covariates: the covariate-coefficient theta gradient needs
  # the per-subject covariate values, which the linkage does not carry yet
  .mrc <- ui$muRefCovariateDataFrame
  if (is.data.frame(.mrc) && nrow(.mrc) > 0L) {
    stop("est=\"stan\" does not yet support mu-referenced covariates (",
         paste(unique(.mrc$covariateParameter), collapse = ", "),
         "): the covariate-coefficient gradient is not wired",
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
attr(nlmixr2Est.stan, "iov") <- FALSE

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
  for (.n in .gen$notes) cli::cli_inform(paste0("est=\"stan\": ", .n))
  if (!is.null(control$stanFile)) {
    writeLines(.gen$code, control$stanFile)
  }
  if (!isTRUE(control$run)) {
    .out <- list(code = .gen$code, data = .gen$data, map = .map,
                 priors = .pri, notes = .gen$notes, ui = ui,
                 control = control)
    class(.out) <- "nlmixr2stanCode"
    return(.out)
  }
  rxode2::rxReq("rstan")
  # ---- compile (cached), then link --------------------------------------
  .sm <- stanCompile(.gen$code, cache = control$cache,
                     cacheDir = control$cacheDir, verbose = control$verbose)
  .needSens <- any(!.map$theta$fix &
                     !(seq_len(nrow(.map$theta)) %in% .map$muRefIdx))
  .h <- stanLinkSetup(ui, env$data, likelihood = control$likelihood,
                      rxControl = control$rxControl,
                      thetaSens = .needSens, cores = control$likCores)
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
  # gradient conditioning: keep nlmixr2est's Omega^-1 commensurate with the
  # model's initial Omega (the conditional value is Omega-free)
  .omInv <- tryCatch(solve(ui$omega), error = function(e) NULL)
  if (!is.null(.omInv)) {
    tryCatch(.linkSetOmegaInv(.omInv), error = function(e) NULL)
  }
  # ---- sample -------------------------------------------------------------
  # jitter draws for the inits come from R's RNG; leave the user's RNG state
  # as we found it
  .oldSeed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit({
    if (!is.null(.oldSeed)) {
      assign(".Random.seed", .oldSeed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv())) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(control$seed)
  .init <- if (identical(control$init, "ini")) {
    .stanInit(.map, .gen$blockSpecs, .nid, control)
  } else control$init
  .sf <- rstan::sampling(.sm, data = .gen$data, chains = control$chains,
                         iter = control$iter, warmup = control$warmup,
                         thin = control$thin, seed = control$seed,
                         init = .init, cores = 1,
                         refresh = if (control$verbose) max(1L, control$iter %/% 10L) else 0L,
                         control = list(adapt_delta = control$adapt_delta,
                                        max_treedepth = control$max_treedepth))
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
  .keep <- !grepl("^(z_|eta\\[|omegaOut|logLikSubj|lp__)", rownames(.sum))
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
  } else ""
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
  .fit
}

#' @export
print.nlmixr2stanCode <- function(x, ...) {
  cat(x$code, sep = "\n")
  invisible(x)
}
