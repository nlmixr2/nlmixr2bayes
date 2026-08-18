## Generate a Stan program from an nlmixr2 model.
##
## The parameter block Stan hands the bridge is, per subject, the population
## thetas followed by that subject's etas.  Repeating the thetas in every block
## looks wasteful but is exactly right: the operands are the same Stan `var`s
## each time, so the reverse sweep accumulates d/dtheta across subjects by
## itself, and no special casing is needed for shared parameters.
##
## The bridge returns ODE states.  Everything downstream of the states -- the
## parameter transformations and the prediction -- is re-emitted as Stan code,
## because rxode2 has no sensitivities for derived `lhs` values and Stan
## differentiates that part itself.

## Stan identifiers cannot contain dots.
.rxsName <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

.rxsBinary <- c("+", "-", "*", "/")
.rxsFuns <- c(exp = "exp", log = "log", sqrt = "sqrt", abs = "abs",
              sin = "sin", cos = "cos", tan = "tan", logit = "logit",
              expit = "inv_logit")

## R expression -> Stan expression.  Deliberately narrow: anything outside this
## set is rejected rather than mistranslated.
.rxsExprToStan <- function(e) {
  if (is.numeric(e)) {
    return(format(e, digits = 17))
  }
  if (is.name(e)) {
    return(.rxsName(as.character(e)))
  }
  if (!is.call(e)) {
    stop("rxstan codegen: cannot translate ", deparse(e), call. = FALSE)
  }

  op <- as.character(e[[1]])

  if (op == "(") return(paste0("(", .rxsExprToStan(e[[2]]), ")"))

  if (op %in% .rxsBinary) {
    if (length(e) == 2L) return(paste0("(", op, .rxsExprToStan(e[[2]]), ")"))
    return(paste0("(", .rxsExprToStan(e[[2]]), " ", op, " ",
                  .rxsExprToStan(e[[3]]), ")"))
  }
  if (op == "^") {
    return(paste0("pow(", .rxsExprToStan(e[[2]]), ", ",
                  .rxsExprToStan(e[[3]]), ")"))
  }
  if (op %in% names(.rxsFuns)) {
    return(paste0(.rxsFuns[[op]], "(",
                  paste(vapply(as.list(e)[-1], .rxsExprToStan, character(1)),
                        collapse = ", "), ")"))
  }
  stop("rxstan codegen: unsupported function `", op, "` in ", deparse(e),
       call. = FALSE)
}

## Splits the model into parameter transformations, ODEs, prediction
## assignments and the residual error line.  `linCmtRepl` replaces a
## `linCmt()` right-hand side, since Stan sees the expanded system's states.
.rxsSplitModel <- function(ui, linCmtRepl = NULL) {
  states <- rxode2::rxState(rxode2::rxode2(rxode2::rxNorm(ui)))

  pre <- list()
  post <- list()
  err <- list()
  seenState <- character(0)

  for (e in ui$lstExpr) {
    head1 <- as.character(e[[1]])
    if (head1 == "~") {
      ## Keyed by endpoint, since a multi-endpoint model has one line each.
      err[[as.character(e[[2]])]] <- e
      next
    }
    lhs <- e[[2]]
    if (is.call(lhs)) next  # d/dt(...) and alag()/f()/dur()/rate() are rxode2's

    if (!is.null(linCmtRepl) && is.call(e[[3]]) &&
        identical(as.character(e[[3]][[1]]), "linCmt")) {
      e[[3]] <- linCmtRepl
    }

    vars <- all.vars(e[[3]])
    if (any(vars %in% c(states, seenState))) {
      post[[length(post) + 1L]] <- e
      seenState <- c(seenState, as.character(lhs))
    } else {
      pre[[length(pre) + 1L]] <- e
    }
  }

  list(states = states, pre = pre, post = post, err = err)
}

## Residual error line -> how to build the Stan density for one observation.
## `dist` picks the family, `loc` the location expression and `sd` the scale;
## nlmixr2 treats lnorm as a TRANSFORM (predDf$transform == "lnorm") rather
## than a separate errType, so it lands here as its own function name.
.rxsErrorSd <- function(ui, errLine) {
  rhs <- errLine[[3]]
  terms <- if (is.call(rhs) && identical(as.character(rhs[[1]]), "+")) {
    list(rhs[[2]], rhs[[3]])
  } else {
    list(rhs)
  }

  add <- NULL
  prop <- NULL
  lnorm <- NULL
  powPar <- NULL
  for (t in terms) {
    if (!is.call(t)) stop("rxstan codegen: cannot read error model ",
                          deparse(rhs), call. = FALSE)
    fn <- as.character(t[[1]])
    par <- as.character(t[[2]])
    if (fn == "add") add <- par
    else if (fn == "prop") prop <- par
    else if (fn == "lnorm") lnorm <- par
    else if (fn == "pow") {
      powPar <- c(par, as.character(t[[3]]))
    } else {
      stop("rxstan codegen: unsupported residual error '", fn, "'. ",
           "add(), prop(), add() + prop(), lnorm() and pow() are handled.",
           call. = FALSE)
    }
  }

  ## pow() lets the exponent be an existing theta, which nlmixr2 then moves
  ## into the error parameters -- that would silently reshape the theta block,
  ## so only the standalone form is accepted.
  if (!is.null(powPar) && (!is.null(add) || !is.null(prop) || !is.null(lnorm))) {
    stop("rxstan codegen: pow() combined with another error term is not ",
         "supported; use pow() on its own", call. = FALSE)
  }
  if (!is.null(lnorm) && (!is.null(add) || !is.null(prop))) {
    stop("rxstan codegen: lnorm() combined with another error term is not ",
         "supported", call. = FALSE)
  }

  if (!is.null(lnorm)) {
    ## Additive on the log scale.  lognormal_lpdf is the density of dv itself,
    ## so log_lik stays on the observation scale and loo() remains valid.
    return(list(pars = lnorm, dist = "lognormal",
                loc = "log(pred[i])", sd = .rxsName(lnorm)))
  }
  if (!is.null(powPar)) {
    return(list(pars = powPar, dist = "normal", loc = "pred[i]",
                sd = sprintf("%s * pow(pred[i], %s)",
                             .rxsName(powPar[1]), .rxsName(powPar[2]))))
  }
  if (!is.null(add) && !is.null(prop)) {
    list(pars = c(add, prop), dist = "normal", loc = "pred[i]",
         sd = sprintf("sqrt(square(%s) + square(%s * pred[i]))",
                      .rxsName(add), .rxsName(prop)))
  } else if (!is.null(add)) {
    list(pars = add, dist = "normal", loc = "pred[i]", sd = .rxsName(add))
  } else if (!is.null(prop)) {
    list(pars = prop, dist = "normal", loc = "pred[i]",
         sd = sprintf("%s * pred[i]", .rxsName(prop)))
  } else {
    stop("rxstan codegen: no residual error found in ", deparse(rhs),
         call. = FALSE)
  }
}

## One observation's log density.  Emitted into BOTH the model block and
## generated quantities so log_lik cannot drift away from what was fitted.
## Uses target += rather than ~ because loo needs the normalizing constant.
.rxsLogLik <- function(errSpec, cens = FALSE) {
  if (cens) {
    sprintf("rxs_obs_ll_%s(dv[i], %s, %s, cens[i], hasLimit[i], limit[i])",
            errSpec$dist, errSpec$loc, errSpec$sd)
  } else {
    sprintf("%s_lpdf(dv[i] | %s, %s)", errSpec$dist, errSpec$loc, errSpec$sd)
  }
}

## The per-observation likelihood, branching on endpoint when there is more
## than one.  `assign` is the statement around the density, so the model block
## and generated quantities share this and cannot diverge.  A single-endpoint
## model emits exactly what it did before multiple endpoints existed.
.rxsLikLines <- function(errSpecs, cens, assign, indent = "    ") {
  if (length(errSpecs) == 1L) {
    return(paste0(indent, sprintf(assign, .rxsLogLik(errSpecs[[1L]], cens))))
  }
  vapply(seq_along(errSpecs), function(k) {
    kw <- if (k == 1L) "if" else "else if"
    paste0(indent, sprintf("%s (dvid[i] == %d) ", kw, k),
           sprintf(assign, .rxsLogLik(errSpecs[[k]], cens)))
  }, character(1))
}

## Which endpoint's prediction each observation row takes.
.rxsPredLines <- function(predVars, indent = "        ") {
  if (length(predVars) == 1L) {
    return(paste0(indent, sprintf("pred[i] = %s;", predVars)))
  }
  vapply(seq_along(predVars), function(k) {
    kw <- if (k == 1L) "if" else "else if"
    paste0(indent, sprintf("%s (dvid[i] == %d) pred[i] = %s;", kw, k,
                           predVars[k]))
  }, character(1))
}

## The censored likelihood, emitted as a Stan function so the model block and
## generated quantities call one definition rather than repeating the
## branching.  M3: a censored record contributes the probability mass over the
## interval it is known to lie in, not a density at the reported value.
## One function per distinct family, since endpoints need not share one.
.rxsCensFun <- function(errSpecs) {
  dists <- unique(vapply(errSpecs, `[[`, character(1), "dist"))
  unlist(lapply(dists, function(d) {
    c(
      sprintf("  real rxs_obs_ll_%s(real y, real mu, real sigma, int cens, int hasLimit, real limit) {", d),
      sprintf("    if (cens == 0) return %s_lpdf(y | mu, sigma);", d),
      "    if (cens == 1) {  // below the limit of quantification",
      sprintf("      if (hasLimit == 1) return log_diff_exp(%s_lcdf(y | mu, sigma), %s_lcdf(limit | mu, sigma));", d, d),
      sprintf("      return %s_lcdf(y | mu, sigma);", d),
      "    }",
      sprintf("    if (hasLimit == 1) return log_diff_exp(%s_lcdf(limit | mu, sigma), %s_lcdf(y | mu, sigma));", d, d),
      sprintf("    return %s_lccdf(y | mu, sigma);", d),
      "  }")
  }), use.names = FALSE)
}

## nlmixr2 puts censoring in the DATA, not the model: CENS is 0 observed,
## 1 left-censored at dv, -1 right-censored at dv, and LIMIT is the other end
## of the interval when there is one.
.rxsCensCols <- function(obs) {
  nm <- names(obs)
  cens <- if ("cens" %in% tolower(nm)) obs[[nm[match("cens", tolower(nm))]]] else NULL
  limit <- if ("limit" %in% tolower(nm)) obs[[nm[match("limit", tolower(nm))]]] else NULL

  if (is.null(cens)) return(NULL)
  ## Values are not validated here: rxsRegister()'s probe solve has already run
  ## and rxode2 rejects anything outside -1/0/1, naming the offending row.
  cens <- as.integer(cens)
  cens[is.na(cens)] <- 0L
  if (all(cens == 0L)) return(NULL)

  lim <- if (is.null(limit)) rep(NA_real_, length(cens)) else as.numeric(limit)
  hasLimit <- as.integer(!is.na(lim) & is.finite(lim) & cens != 0L)
  lim[hasLimit == 0L] <- 0
  list(cens = cens, hasLimit = hasLimit, limit = lim)
}

## `<lower=,upper=>` for one ini() row, omitting infinite ends.
.rxsBound <- function(lower, upper) {
  parts <- character(0)
  if (is.finite(lower)) parts <- c(parts, sprintf("lower=%.17g", lower))
  if (is.finite(upper)) parts <- c(parts, sprintf("upper=%.17g", upper))
  if (!length(parts)) return("")
  paste0("<", paste(parts, collapse = ", "), ">")
}

## The prior for one parameter, overridable by name through `priors`.
.rxsPrior <- function(name, default, priors) {
  if (!is.null(priors[[name]])) priors[[name]] else default
}

## --- linCmt() ---------------------------------------------------------------
##
## rxode2 will not give us sensitivities for linCmt().  calcSens creates the
## rx__sens_* compartments but nothing ever writes to them, so a linCmt model
## silently solves to all-zero derivatives; its real Jacobians live in a
## separate internal path (linCmtB/sensType) that rxSolve does not expose.
##
## So the closed form is replaced with the equivalent ODE system, which the
## bridge already handles and differentiates correctly.  The expansion is then
## checked against linCmt() itself, so getting the parameterization wrong is a
## loud failure rather than quietly wrong numbers.

.rxsLinCmtNames <- list(
  cl = c("cl", "CL"),
  v = c("v", "v1", "vc", "V", "V1", "Vc"),
  q = c("q", "q1", "Q", "Q1"),
  vp = c("vp", "vp1", "v2", "Vp", "Vp1", "V2"),
  q2 = c("q2", "Q2"),
  vp2 = c("vp2", "v3", "Vp2", "V3"),
  ka = c("ka", "Ka", "KA"))

.rxsPick <- function(candidates, lhs) {
  hit <- candidates[candidates %in% lhs]
  if (length(hit)) hit[1] else NA_character_
}

## Structure of a linCmt() model: how many compartments, is there absorption,
## and which names the user gave the parameters.  Detection is the model-vars
## `linCmt` flag (2 with absorption, 1 without, -100 for an ordinary model);
## predDf$linCmt is a different thing and is FALSE even for linCmt models.
.rxsLinCmtInfo <- function(ui) {
  m <- rxode2::rxode2(rxode2::rxNorm(ui))
  mv <- rxode2::rxModelVars(m)
  if (as.integer(mv$flags[["linCmt"]]) <= 0L) return(NULL)

  lhs <- mv$lhs
  info <- list(ncmt = as.integer(mv$flags[["ncmt"]]),
               hasKa = as.integer(mv$flags[["ka"]]) == 1L,
               states = mv$state)
  for (nm in names(.rxsLinCmtNames)) {
    info[[nm]] <- .rxsPick(.rxsLinCmtNames[[nm]], lhs)
  }
  info
}

## The ODE system equivalent to linCmt(), in the clearance parameterization.
.rxsLinCmtOdes <- function(info) {
  need <- c("cl", "v", if (info$ncmt >= 2L) c("q", "vp"),
            if (info$ncmt >= 3L) c("q2", "vp2"), if (info$hasKa) "ka")
  missing <- need[is.na(unlist(info[need]))]
  if (length(missing)) {
    stop("rxstan codegen: linCmt() expansion needs the clearance ",
         "parameterization; could not find ", paste(missing, collapse = ", "),
         " among the model's assignments", call. = FALSE)
  }

  el <- sprintf("(%s / %s)", info$cl, info$v)
  toC <- character(0)
  depot <- character(0)
  periph <- character(0)

  if (info$hasKa) {
    depot <- sprintf("d/dt(depot) = -%s * depot;", info$ka)
    toC <- c(toC, sprintf("+ %s * depot", info$ka))
  }
  if (info$ncmt >= 2L) {
    k12 <- sprintf("(%s / %s)", info$q, info$v)
    k21 <- sprintf("(%s / %s)", info$q, info$vp)
    toC <- c(toC, sprintf("- %s * central + %s * peripheral1", k12, k21))
    periph <- c(periph,
                sprintf("d/dt(peripheral1) = %s * central - %s * peripheral1;",
                        k12, k21))
  }
  if (info$ncmt >= 3L) {
    k13 <- sprintf("(%s / %s)", info$q2, info$v)
    k31 <- sprintf("(%s / %s)", info$q2, info$vp2)
    toC <- c(toC, sprintf("- %s * central + %s * peripheral2", k13, k31))
    periph <- c(periph,
                sprintf("d/dt(peripheral2) = %s * central - %s * peripheral2;",
                        k13, k31))
  }

  central <- sprintf("d/dt(central) = %s - %s * central;",
                     paste(c("0", toC), collapse = " "), el)

  # Compartment order must match linCmt's, or dosing by cmt number would shift.
  c(depot, central, periph)
}

## The expansion is only trustworthy if it reproduces linCmt() itself, so solve
## both at the ini() values and insist they agree.  This is what turns a wrong
## parameterization into an error instead of plausible-looking wrong numbers.
.rxsCheckLinCmt <- function(ui, modelText, info, tol = 1e-6) {
  ini <- ui$iniDf
  pars <- stats::setNames(ini$est[!is.na(ini$ntheta) & is.na(ini$err)],
                          ini$name[!is.na(ini$ntheta) & is.na(ini$err)])
  etas <- ini$name[!is.na(ini$neta1) & ini$neta1 == ini$neta2]
  pars <- c(pars, stats::setNames(rep(0, length(etas)), etas))

  ev <- rxode2::et(amt = 100, cmt = if (info$hasKa) "depot" else "central")
  ev <- rxode2::et(ev, c(0.25, 1, 4, 12, 24))

  solve1 <- function(m) {
    rxode2::rxSolve(m, params = pars, events = ev, returnType = "data.frame",
                    cores = 1L, atol = 1e-10, rtol = 1e-10)
  }
  closed <- try(solve1(rxode2::rxode2(rxode2::rxNorm(ui))), silent = TRUE)
  expanded <- try(solve1(rxode2::rxode2(modelText)), silent = TRUE)
  if (inherits(closed, "try-error") || inherits(expanded, "try-error")) {
    stop("rxstan codegen: could not verify the linCmt() expansion against the ",
         "closed form", call. = FALSE)
  }

  v <- as.character(ui$predDf$var)
  d <- max(abs(closed[[v]] - expanded[[v]])) / max(1, max(abs(closed[[v]])))
  if (!is.finite(d) || d > tol) {
    stop("rxstan codegen: the linCmt() ODE expansion disagrees with linCmt() ",
         "itself (relative difference ", signif(d, 3), "). The model is ",
         "probably using a parameterization the expansion does not cover; ",
         "write it with explicit ODEs.", call. = FALSE)
  }
  invisible(d)
}

## Connected components of the off-diagonal structure of omega: each is a
## correlated block, a lone eta is a block of one.
.rxsOmegaBlocks <- function(omega) {
  n <- nrow(omega)
  if (!n) return(list())
  seen <- logical(n)
  blocks <- list()
  for (i in seq_len(n)) {
    if (seen[i]) next
    grp <- i
    repeat {
      linked <- which(apply(abs(omega[grp, , drop = FALSE]) > 0, 2, any))
      grp2 <- sort(unique(c(grp, linked)))
      if (identical(grp2, grp)) break
      grp <- grp2
    }
    seen[grp] <- TRUE
    blocks[[length(blocks) + 1L]] <- grp
  }
  blocks
}

#' Generate a Stan program and a bridge handle from an nlmixr2 model
#'
#' Supports a diagonal or block-correlated `omega`, `add()` / `prop()` /
#' `add() + prop()` / `lnorm()` / `pow()` residual error, BLQ censoring, and
#' any number of endpoints.  With more than one endpoint each observation row
#' must say which it belongs to, through a `dvid` column or through `cmt`
#' naming the endpoint.
#'
#' `linCmt()` is rewritten as the equivalent ODE system, because rxode2 emits
#' no usable sensitivities for the closed-form solution.  The expansion is
#' checked against `linCmt()` itself on every call, so a parameterization it
#' does not cover is an error rather than quietly wrong numbers.  One to three
#' compartments, with or without first-order absorption, in the clearance
#' parameterization.
#'
#' `ini()` bounds become Stan constraints, and `fix()` removes the parameter
#' from the `parameters` block entirely -- it becomes a constant in
#' `transformed data`, the same treatment brms gives a `constant()` prior.  A
#' fixed theta is also dropped from the sensitivity system, so rxode2 has one
#' fewer equation to solve.
#'
#' Default priors are weak and centered on the `ini()` estimates, which are
#' starting values rather than genuine prior beliefs; override them through
#' `priors` for anything real.
#'
#' @section Log-scale residual error:
#' `lnorm()` puts the location at `log(pred)`, so it is the *relative* accuracy
#' of the prediction that matters and the analytic sensitivity carries a
#' `1 / pred` factor.  Where the curve has decayed to the solver's absolute
#' tolerance that factor amplifies solver noise without bound, and the gradient
#' stops agreeing with the value even though the chain rule is correct --
#' measured at a relative difference of about `1e-2` for predictions around
#' `1e-11`, against `1e-9` for the same data under `add()`.
#'
#' Tightening `atol` does not rescue it (a prediction of `1e-11` would need
#' `atol` near `1e-17`).  Keep the predictions scaled instead: start the
#' sampler somewhere sensible with [rxsInit()] so it does not wander into
#' regions where the curve underflows, and be wary of observations far into the
#' terminal phase.  [rxsCheckGradient()] is what detects this.
#'
#' @section Solver tolerance:
#' The [rxsRegister()] default of `1e-8` suits a plain PK model but not every
#' system.  A PD compartment with feedback is harder to integrate, and on the
#' PK/PD model in the test suite the gradient check sits at `1.4e-5` at the
#' default, `3.1e-8` at `1e-10`, and no better at `1e-12` -- solver accuracy
#' rather than a wrong derivative, and the finite-difference noise floor below
#' that.  If [rxsCheckGradient()] is worse than about `1e-6`, tighten
#' `atol`/`rtol` before suspecting anything else.
#'
#' @section Covariates:
#' Any symbol in the model that is neither an estimated parameter nor a state
#' is looked for as a column of `data`.  rxode2 reads it straight from the
#' event table, so a covariate costs nothing on the solver side and does not
#' disarm the fast path; Stan is only told about it because the parameter
#' transforms are re-emitted there to build the prediction.
#'
#' Whether a covariate is passed once per subject or once per observation is
#' decided by looking at the data, not by naming convention: a column constant
#' within every subject becomes `vector[nSub]`, one that moves becomes
#' `vector[nObs]` and rxode2 interpolates it during the solve according to
#' `covsInterpolation` (last observation carried forward by default).  Note
#' that a covariate which happens to be constant in the data at hand will be
#' treated as baseline, which is harmless for the fit but changes the
#' generated program.
#'
#' A symbol that matches no column is an error naming the symbol, raised
#' before the model reaches rxode2 -- otherwise it surfaces as a dump of the
#' whole expanded sensitivity system.
#'
#' @param ui an nlmixr2 model (a function, or the `rxUi` object from
#'   [rxode2::rxode2()])
#' @param data a data frame with `id`, `time`, `dv` and the usual event columns
#' @param priors named list of Stan distribution strings overriding the
#'   defaults, keyed by the model's own parameter names, e.g.
#'   `list(tka = "normal(0, 2)", eta.ka = "exponential(1)", add.sd = "cauchy(0, 1)")`.
#'   An eta's entry is the prior on its between-subject SD.
#' @param priorSd standard deviation of the default theta priors.
#' @param lkjEta shape of the `lkj_corr_cholesky` prior on each correlated
#'   `omega` block.  The default of 2 mildly favors weaker correlations; 1 is
#'   uniform over correlation matrices.
#' @param ... passed to [rxsRegister()]
#' @return a list with the generated Stan `code`, the bridge `handle`, the
#'   `standata` to pass to [rstan::sampling()], and the parameter names.  Names
#'   are the model's own; `stanNames` gives the sanitized versions actually used
#'   in the generated program, since Stan identifiers cannot contain dots.
#' @author Lukas A. Widmer
#' @export
rxsStanFromUi <- function(ui, data, priors = list(), priorSd = 10, lkjEta = 2,
                          ...) {
  if (is.function(ui)) ui <- rxode2::rxode2(ui)
  ini <- ui$iniDf

  thetaDf <- ini[!is.na(ini$ntheta) & is.na(ini$err), , drop = FALSE]
  errDf <- ini[!is.na(ini$err), , drop = FALSE]
  ## Only the diagonal rows name an eta; off-diagonals are covariance entries
  ## carrying a composite name like "(eta.ka,eta.cl)".
  etaDf <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
  etaDf <- etaDf[order(etaDf$neta1), , drop = FALSE]

  omegaMat <- if (nrow(etaDf)) {
    as.matrix(ui$omega)[etaDf$name, etaDf$name, drop = FALSE]
  } else {
    matrix(numeric(0), 0, 0)
  }
  blocks <- .rxsOmegaBlocks(omegaMat)
  if (any(etaDf$fix & vapply(seq_len(nrow(etaDf)), function(i) {
        any(lengths(Filter(function(b) i %in% b, blocks)) > 1L)
      }, logical(1)))) {
    stop("rxstan codegen: fixing an eta inside a correlated block is not ",
         "supported", call. = FALSE)
  }
  if (isTRUE(ui$predDf$linCmt)) {
    stop("rxstan codegen: linCmt() models are not supported yet. rxode2 does ",
         "emit sensitivities for them, but against linCmt's own parameters ",
         "(p1, v1, ka) and under a different column naming convention, so the ",
         "bridge needs a separate path. Write the model with explicit ODEs ",
         "for now.", call. = FALSE)
  }
  if (any(errDf$fix)) {
    stop("rxstan codegen: a fixed residual error parameter is not supported ",
         "yet", call. = FALSE)
  }
  unknown <- setdiff(names(priors), c(thetaDf$name, etaDf$name, errDf$name))
  if (length(unknown)) {
    stop("rxstan codegen: priors given for unknown parameter(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }

  ## linCmt() is replaced by its ODE equivalent, because rxode2 gives no
  ## sensitivities for the closed form.
  linInfo <- .rxsLinCmtInfo(ui)
  linRepl <- if (!is.null(linInfo)) {
    str2lang(sprintf("central / %s", linInfo$v))
  }

  parts <- .rxsSplitModel(ui, linCmtRepl = linRepl)
  epVar <- as.character(ui$predDf$var)
  missingErr <- setdiff(epVar, names(parts$err))
  if (length(missingErr)) {
    stop("rxstan codegen: no residual error line for endpoint(s): ",
         paste(missingErr, collapse = ", "), call. = FALSE)
  }
  errSpecs <- lapply(epVar, function(v) .rxsErrorSd(ui, parts$err[[v]]))
  names(errSpecs) <- epVar

  ## fix() means the value leaves the parameter block entirely and becomes a
  ## constant -- the same thing brms does for a constant() prior.  It is never
  ## emulated with a tight prior, which would still sample it.  A fixed theta
  ## also drops out of `sens`, so rxode2 solves one fewer sensitivity equation.
  thetaEst <- thetaDf[!thetaDf$fix, , drop = FALSE]
  thetaFix <- thetaDf[thetaDf$fix, , drop = FALSE]

  thetaNames <- thetaEst$name
  etaNames <- etaDf$name
  block <- c(thetaNames, etaNames)
  if (!length(block)) {
    stop("rxstan codegen: every parameter is fixed, nothing to estimate",
         call. = FALSE)
  }

  fixedParams <- stats::setNames(thetaFix$est, thetaFix$name)

  ## Covariates: symbols in the transforms that are neither parameters nor
  ## states but are columns of `data`.  Checked BEFORE the probe solve,
  ## because rxode2 reports a missing one by dumping the whole expanded
  ## sensitivity system, which buries the single word that matters.
  defined <- c(vapply(parts$pre, function(e) as.character(e[[2]]), character(1)),
               vapply(parts$post, function(e) as.character(e[[2]]), character(1)))
  used <- unique(unlist(lapply(c(parts$pre, parts$post),
                               function(e) all.vars(e[[3]])), use.names = FALSE))
  known <- c(thetaNames, etaNames, thetaFix$name, parts$states, defined, "pi")
  cand <- setdiff(used, known)
  covNames <- intersect(cand, names(data))
  orphan <- setdiff(cand, names(data))
  if (length(orphan)) {
    stop("rxstan codegen: the model uses ", paste(orphan, collapse = ", "),
         ", which is neither an estimated parameter nor a column of `data`",
         call. = FALSE)
  }

  ## The bridge solves with the estimated thetas and etas as rxode2 parameters;
  ## fixed thetas are passed as constants instead.
  modelText <- rxode2::rxNorm(ui)
  if (!is.null(linInfo)) {
    keep <- grep("linCmt\\(", strsplit(modelText, "\n")[[1]], invert = TRUE,
                 value = TRUE)
    modelText <- paste(c(keep, .rxsLinCmtOdes(linInfo),
                         sprintf("%s = central / %s;", ui$predDf$var,
                                 linInfo$v)),
                       collapse = "\n")
    .rxsCheckLinCmt(ui, modelText, linInfo)
    parts$states <- rxode2::rxState(rxode2::rxode2(modelText))
  }
  ids <- unique(data$id)
  handle <- rxsRegister(modelText, events = data, sens = block,
                        output = parts$states, params = fixedParams,
                        perSubject = TRUE, ...)

  nSub <- length(ids)
  nTheta <- length(thetaNames)
  nEta <- length(etaNames)
  nState <- length(parts$states)

  obs <- data[is.na(data$evid) | data$evid == 0, , drop = FALSE]
  nObs <- nrow(obs)
  if (nObs != attr(handle, "nobs")) {
    rxsRelease(handle)
    stop("rxstan codegen: ", nObs, " observation rows in `data` but the solve ",
         "returned ", attr(handle, "nobs"), call. = FALSE)
  }

  ## --- Stan source ---------------------------------------------------------
  ## Parameters are declared under their own names rather than as theta[k], so
  ## the generated program and its output read like the nlmixr2 model.
  sName <- .rxsName(thetaNames)
  eName <- .rxsName(etaNames)
  oName <- paste0("omega_", eName)
  ePar <- .rxsName(unlist(lapply(errSpecs, `[[`, "pars"), use.names = FALSE))

  thetaDecl <- sprintf("  real%s %s;",
                       mapply(.rxsBound, thetaEst$lower, thetaEst$upper),
                       sName)
  etaFixed <- etaDf$fix

  ## One correlated block becomes an SD vector plus a Cholesky correlation
  ## factor; a lone eta stays a single SD.  nlmixr2 states etas as VARIANCES,
  ## Stan uses SDs here, hence the sqrt on the fixed ones.
  blkName <- sprintf("blk%d", seq_along(blocks))
  blkSize <- lengths(blocks)
  corr <- which(blkSize > 1L)
  solo <- which(blkSize == 1L)
  soloIdx <- unlist(blocks[solo], use.names = FALSE)

  omegaDecl <- c(
    if (length(corr)) {
      c(sprintf("  vector<lower=0>[%d] omega_%s;", blkSize[corr], blkName[corr]),
        sprintf("  cholesky_factor_corr[%d] L_%s;", blkSize[corr], blkName[corr]))
    },
    if (length(soloIdx) && any(!etaFixed[soloIdx])) {
      sprintf("  real<lower=0> %s;", oName[soloIdx][!etaFixed[soloIdx]])
    })

  etaLines <- unlist(lapply(seq_along(blocks), function(b) {
    idx <- blocks[[b]]
    if (length(idx) == 1L) {
      return(sprintf("    eta[%d, s] = %s * z[%d, s];", idx, oName[idx], idx))
    }
    zvec <- paste(sprintf("z[%d, s]", idx), collapse = ", ")
    c(sprintf("    {"),
      sprintf("      vector[%d] e = diag_pre_multiply(omega_%s, L_%s) * [%s]';",
              length(idx), blkName[b], blkName[b], zvec),
      sprintf("      eta[%d, s] = e[%d];", idx, seq_along(idx)),
      "    }")
  }), use.names = FALSE)
  errDecl <- sprintf("  real%s %s;",
                     mapply(.rxsBound, errDf$lower, errDf$upper), ePar)

  ## Fixed values live in transformed data: visible in the program, but not
  ## sampled.  nlmixr2 states an eta as a VARIANCE, Stan uses the SD here.
  fixedDecl <- c(
    if (nrow(thetaFix)) sprintf("  real %s = %.17g;  // fixed by ini()",
                                .rxsName(thetaFix$name), thetaFix$est),
    if (length(soloIdx) && any(etaFixed[soloIdx])) {
      k <- soloIdx[etaFixed[soloIdx]]
      sprintf("  real %s = %.17g;  // fixed by ini() (sd of variance %.17g)",
              oName[k], sqrt(etaDf$est[k]), etaDf$est[k])
    })

  ## Thetas are already in scope under their own names (estimated ones as
  ## parameters, fixed ones as transformed data), so only the per-subject etas
  ## need a local -- redeclaring a theta here would shadow it and Stan refuses.
  decl <- sprintf("        real %s = eta[%d, s];", eName, seq_along(eName))

  preLines <- vapply(parts$pre, function(e) {
    sprintf("        real %s = %s;", .rxsName(as.character(e[[2]])),
            .rxsExprToStan(e[[3]]))
  }, character(1))

  stateLines <- sprintf("        real %s = ys[%d * nObs + i];",
                        .rxsName(parts$states),
                        seq_along(parts$states) - 1L)

  postLines <- vapply(parts$post, function(e) {
    sprintf("        real %s = %s;", .rxsName(as.character(e[[2]])),
            .rxsExprToStan(e[[3]]))
  }, character(1))

  predVars <- .rxsName(epVar)
  censCols <- .rxsCensCols(obs)
  hasCens <- !is.null(censCols)

  ## rxode2 reads covariates straight from the event table -- baseline or
  ## time-varying alike, interpolating the latter during the solve -- so
  ## nothing changes on that side.  Stan needs them declared, because the
  ## transforms are re-emitted there to build the prediction.
  ##
  ## A covariate constant within a subject only needs one value per subject;
  ## one that moves needs a value per observation.  Getting this wrong is not
  ## cosmetic -- collapsing a time-varying covariate to one value per subject
  ## silently changes the model.
  covs <- lapply(covNames, function(cv) {
    v <- obs[[cv]]
    if (anyNA(v)) {
      stop("rxstan codegen: covariate '", cv, "' has missing values at ",
           sum(is.na(v)), " observation(s)", call. = FALSE)
    }
    perId <- tapply(v, obs$id, function(x) length(unique(x)))
    baseline <- all(perId == 1L)
    list(name = cv, stan = .rxsName(cv), baseline = baseline,
         values = if (baseline) {
           as.numeric(vapply(split(v, obs$id)[as.character(ids)],
                             function(x) x[1], numeric(1)))
         } else {
           as.numeric(v)
         })
  })
  names(covs) <- covNames

  covDecl <- vapply(covs, function(cv) {
    sprintf("  vector[%s] cov_%s;", if (cv$baseline) "nSub" else "nObs", cv$stan)
  }, character(1))
  ## Aliased to the model's own name so the re-emitted transform needs no
  ## rewriting.
  covLocal <- vapply(covs, function(cv) {
    sprintf("        real %s = cov_%s[%s];", cv$stan, cv$stan,
            if (cv$baseline) "s" else "i")
  }, character(1))

  ## Which endpoint each observation row belongs to.  nlmixr2 identifies it
  ## with dvid, or with cmt naming the endpoint.
  dvidVec <- NULL
  if (length(epVar) > 1L) {
    nm <- tolower(names(obs))
    key <- if ("dvid" %in% nm) {
      as.character(obs[[names(obs)[match("dvid", nm)]]])
    } else if ("cmt" %in% nm) {
      as.character(obs[[names(obs)[match("cmt", nm)]]])
    } else {
      stop("rxstan codegen: ", length(epVar), " endpoints, so `data` needs a ",
           "`dvid` or `cmt` column saying which endpoint each observation is",
           call. = FALSE)
    }
    dvidVec <- match(key, epVar)
    if (anyNA(dvidVec)) {
      stop("rxstan codegen: observation(s) name an endpoint the model does ",
           "not have: ", paste(unique(key[is.na(dvidVec)]), collapse = ", "),
           ". The model has ", paste(epVar, collapse = ", "), call. = FALSE)
    }
  }

  ## Caught here rather than left to Stan: a lognormal residual needs positive
  ## observations, and without this the first bad row surfaces as
  ## "lognormal_lpdf: Random variable is -0.105, but must be nonnegative"
  ## partway through warmup, with nothing pointing at the data.
  for (k in seq_along(errSpecs)) {
    if (!identical(errSpecs[[k]]$dist, "lognormal")) next
    rows <- if (is.null(dvidVec)) seq_len(nrow(obs)) else which(dvidVec == k)
    bad <- rows[!is.na(obs$dv[rows]) & obs$dv[rows] <= 0]
    if (length(bad)) {
      stop("rxstan codegen: lnorm() needs positive observations, but dv has ",
           length(bad), " value(s) <= 0 for endpoint ", epVar[k],
           " (first at row ", bad[1], ", dv = ",
           format(obs$dv[bad[1]]), ")", call. = FALSE)
    }
  }

  ## Defaults are weak and centered on the ini() estimates; `priors` overrides
  ## them by the model's own parameter names.
  thetaPrior <- sprintf("  %s ~ %s;", sName,
                        mapply(.rxsPrior, thetaNames,
                               sprintf("normal(%.17g, %.17g)", thetaEst$est,
                                       priorSd),
                               MoreArgs = list(priors = priors)))
  omegaPrior <- c(
    if (length(soloIdx) && any(!etaFixed[soloIdx])) {
      k <- soloIdx[!etaFixed[soloIdx]]
      sprintf("  %s ~ %s;", oName[k],
              mapply(.rxsPrior, etaNames[k], "normal(0, 1)",
                     MoreArgs = list(priors = priors)))
    },
    if (length(corr)) {
      c(sprintf("  omega_%s ~ normal(0, 1);", blkName[corr]),
        sprintf("  L_%s ~ lkj_corr_cholesky(%.17g);", blkName[corr], lkjEta))
    })
  errPrior <- sprintf("  %s ~ %s;", ePar,
                      mapply(.rxsPrior, errDf$name, "normal(0, 10)",
                             MoreArgs = list(priors = priors)))

  code <- c(
    "// Generated by rxstan::rxsStanFromUi() -- do not edit by hand.",
    "//",
    "// rx_solve() is left undefined here: the rxstan bridge supplies it, and",
    "// returns the ODE states together with rxode2's analytic sensitivities.",
    "// The parameter transformations and the prediction are re-emitted below",
    "// so that Stan differentiates that part itself.",
    "functions {",
    "  vector rx_solve(int handle, vector p);",
    if (!is.null(censCols)) .rxsCensFun(errSpecs),
    "}",
    "data {",
    "  int<lower=1> nSub;",
    "  int<lower=1> nObs;",
    "  int<lower=1> handle;",
    "  vector[nObs] dv;",
    "  array[nObs] int<lower=1, upper=nSub> subj;",
    if (length(epVar) > 1L) sprintf("  array[nObs] int<lower=1, upper=%d> dvid;",
                                    length(epVar)),
    if (length(covDecl)) covDecl,
    if (!is.null(censCols)) "  array[nObs] int<lower=-1, upper=1> cens;",
    if (!is.null(censCols)) "  array[nObs] int<lower=0, upper=1> hasLimit;",
    if (!is.null(censCols)) "  vector[nObs] limit;",
    "}",
    if (length(fixedDecl)) "transformed data {",
    fixedDecl,
    if (length(fixedDecl)) "}",
    "parameters {",
    thetaDecl,
    omegaDecl,
    if (nEta) sprintf("  matrix[%d, nSub] z;", nEta),
    errDecl,
    "}",
    "transformed parameters {",
    sprintf("  vector[%d] theta = [%s]';", nTheta, paste(sName, collapse = ", ")),
    if (nEta) sprintf("  matrix[%d, nSub] eta;", nEta),
    "  vector[nObs] pred;",
    if (nEta) "  for (s in 1 : nSub) {",
    if (nEta) etaLines,
    if (nEta) "  }",
    ## p and ys live in a local block so Stan does not write them out with
    ## every draw; pred is deliberately kept, since generated quantities needs
    ## it and it doubles as the posterior predictive mean.
    "  {",
    sprintf("    vector[%d * nSub] p;", nTheta + nEta),
    "    for (s in 1 : nSub) {",
    sprintf("      p[(s - 1) * %d + 1 : (s - 1) * %d + %d] = theta;",
            nTheta + nEta, nTheta + nEta, nTheta),
    if (nEta) sprintf("      p[(s - 1) * %d + %d : s * %d] = eta[ : , s];",
                      nTheta + nEta, nTheta + 1L, nTheta + nEta),
    "    }",
    sprintf("    vector[%d * nObs] ys = rx_solve(handle, p);", nState),
    "    for (i in 1 : nObs) {",
    "      {",
    "        int s = subj[i];",
    if (length(covLocal)) covLocal,
    decl,
    preLines,
    stateLines,
    postLines,
    .rxsPredLines(predVars),
    "      }",
    "    }",
    "  }",
    "}",
    "model {",
    thetaPrior,
    omegaPrior,
    if (nEta) "  to_vector(z) ~ std_normal();",
    errPrior,
    "  for (i in 1 : nObs) {",
    .rxsLikLines(errSpecs, hasCens, "target += %s;"),
    "  }",
    "}",
    "generated quantities {",
    "  vector[nObs] log_lik;",
    if (length(corr)) {
      sprintf("  matrix[%d, %d] corr_%s = multiply_lower_tri_self_transpose(L_%s);",
              blkSize[corr], blkSize[corr], blkName[corr], blkName[corr])
    },
    "  for (i in 1 : nObs) {",
    .rxsLikLines(errSpecs, hasCens, "log_lik[i] = %s;"),
    "  }",
    "}")

  standata <- list(nSub = nSub, nObs = nObs,
                   handle = as.integer(unclass(handle)),
                   dv = as.numeric(obs$dv),
                   subj = match(obs$id, ids))
  if (!is.null(dvidVec)) standata$dvid <- dvidVec
  for (cv in covs) standata[[paste0("cov_", cv$stan)]] <- cv$values
  if (!is.null(censCols)) {
    standata$cens <- censCols$cens
    standata$hasLimit <- censCols$hasLimit
    standata$limit <- censCols$limit
  }

  ## Point inits at the ini() values.  A 1-cmt oral model is bimodal
  ## (flip-flop: swapping absorption and elimination fits nearly as well), so
  ## diffuse inits routinely send chains to different modes; starting from
  ## ini() is what prevents it.  See [rxsInit()] for jittered multi-chain use.
  inits <- list()
  if (nrow(thetaEst)) inits[sName] <- as.list(thetaEst$est)
  if (length(soloIdx)) {
    k <- soloIdx[!etaFixed[soloIdx]]
    if (length(k)) inits[oName[k]] <- as.list(sqrt(etaDf$est[k]))
  }
  for (b in corr) {
    idx <- blocks[[b]]
    inits[[paste0("omega_", blkName[b])]] <- sqrt(etaDf$est[idx])
    inits[[paste0("L_", blkName[b])]] <- diag(length(idx))
  }
  if (nEta) inits$z <- matrix(0, nrow = nEta, ncol = nSub)
  if (length(ePar)) inits[ePar] <- as.list(errDf$est)

  list(code = paste(code[!vapply(code, is.null, logical(1))], collapse = "\n"),
       handle = handle, standata = standata, inits = inits,
       thetaNames = thetaNames, etaNames = etaNames,
       errNames = unlist(lapply(errSpecs, `[[`, "pars"), use.names = FALSE),
       states = parts$states,
       stanNames = list(theta = sName, eta = eName, err = ePar,
                        theta0 = thetaEst$est))
}
