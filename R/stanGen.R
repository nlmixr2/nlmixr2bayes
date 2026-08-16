# The .stan generator (Spec 9).  Pure functions -- no rstan.
#
# Program shape: `parameters` declares the free thetas with their (support-
# promoted) bounds and the per-block Omega parameterization; `transformed
# parameters` assembles the full theta vector (fixed values inlined as
# literals) and the non-centred etas (eta = z L'); `model` carries the prior
# statements from stanPriors(), z ~ std_normal(), and the external
# conditional likelihood.  The observations never enter Stan's data block --
# they live in the linked nlmixr2est problem.
#
# Omega rule (the one that makes hand Jacobians unnecessary): the declared
# Stan parameter is whatever the prior is written on.  invWishart/wishart ->
# cov_matrix with that prior; lkjCorr(Cholesky) -> cholesky_factor_corr + SD
# vector; no block prior -> LKJ(lkjEta) + a default half-Cauchy on each SD,
# supplied AND announced.  Stan's own constraining transform then supplies
# the only Jacobian needed.

#' Classify one omega block against the omega-prior rows
#' @noRd
.stanOmegaBlockSpec <- function(block, omPri, ctl) {
  .first <- block$members[1L]
  # fixed etas: the variance is a model CONSTANT, not a sampled parameter --
  # the eta is still random (z ~ std_normal), but L is fixed.  This is what
  # the IOV preprocessing hook produces (per-occasion unit-variance etas
  # scaled by an estimated sd theta in the model code), and it serves plain
  # fixed-variance etas identically.  Fixed etas with modeled correlations
  # (a fixed block, k > 1) are out of scope and must error, not be
  # approximated.
  if (any(block$fix)) {
    if (block$k > 1L) {
      stop("the omega block over ", paste(block$members, collapse = ", "),
           " contains fixed eta(s) with modeled correlations; ",
           "est=\"stan\" supports fixed DIAGONAL etas only", call. = FALSE)
    }
    if (any(omPri$name %in% block$members)) {
      stop("eta '", .first, "' is fixed; a prior on a fixed eta is a ",
           "contradiction", call. = FALSE)
    }
    return(list(type = "fixed", default = FALSE, nu = NULL, scale = NULL,
                lkjEta = NULL))
  }
  .row <- omPri[omPri$name == .first, , drop = FALSE]
  if (nrow(.row) == 0L) {
    # no block prior: defaults, announced by the caller
    if (any(omPri$name %in% block$members)) {
      stop("an omega prior inside block '",
           paste(block$members, collapse = ","),
           "' is not on the block's first diagonal element; ",
           "normal (TNPRI-style) priors on omega values are not yet supported",
           call. = FALSE)
    }
    return(list(type = if (block$k == 1L) "sd" else "lkj",
                default = TRUE, nu = NULL, scale = NULL,
                lkjEta = ctl$lkjEta))
  }
  .lk <- .stanPriorLookup(.row$prior)
  if (.lk$kind != "matrix") {
    stop("prior '", .lk$name, "' on omega element '", .first,
         "' is not a matrix prior; normal (TNPRI-style) priors on omega ",
         "values are not yet supported by est=\"stan\"", call. = FALSE)
  }
  if (.lk$stanName %in% c("inv_wishart", "wishart",
                          "inv_wishart_cholesky", "wishart_cholesky")) {
    .nu <- eval(.lk$args[[1]], envir = baseenv())
    .scale <- if (length(.lk$args) >= 2L) .stanEvalMvArg(.lk$args[[2]]) else block$init
    if (!is.matrix(.scale)) .scale <- matrix(.scale, 1, 1)
    return(list(type = "cov", default = FALSE,
                stanName = sub("_cholesky$", "", .lk$stanName),
                nu = .nu, scale = .scale, lkjEta = NULL))
  }
  if (.lk$stanName %in% c("lkj_corr", "lkj_corr_cholesky")) {
    return(list(type = "lkj", default = FALSE, nu = NULL, scale = NULL,
                lkjEta = eval(.lk$args[[1]], envir = baseenv())))
  }
  stop("omega prior '", .lk$name, "' is not supported", call. = FALSE) # nocov
}

#' Generate the Stan program + data + notes for a ui
#'
#' @param ui rxode2 ui
#' @param map the parameter map (`.stanMap`)
#' @param pri `stanPriors(ui)`
#' @param ctl `stanControl()`
#' @return list(code=, data=, notes=, blockSpecs=)
#' @noRd
.stanGenerate <- function(ui, map, pri, ctl) {
  .notes <- character(0)
  .free <- map$theta[!map$theta$fix, , drop = FALSE]
  .prByName <- stats::setNames(as.list(seq_len(nrow(pri$pop))), pri$pop$name)
  # ---- parameters: free thetas -------------------------------------------
  .decl <- character(0)
  for (.i in seq_len(nrow(.free))) {
    .r <- .free[.i, ]
    .con <- .stanConstraint(.r$lower, .r$upper)
    .w <- which(pri$pop$name == .r$name & pri$pop$kind != "multivariate")
    if (length(.w) == 1L) {
      # support-promoted bounds from the prior win (they are only ever
      # tighter than the parameter's own)
      .con <- pri$pop$constraint[.w]
    }
    .decl <- c(.decl, paste0("  real", .con, " ", .r$par, ";"))
  }
  # flat thetas: allowed, but say so loudly
  .flat <- .free$name[!(.free$name %in% pri$pop$name) &
                        !(.free$name %in% unlist(strsplit(stats::na.omit(pri$pop$members), ",")))]
  if (length(.flat) > 0L) {
    .notes <- c(.notes, paste0("flat (improper) priors on: ",
                               paste(.flat, collapse = ", ")))
  }
  # ---- omega blocks -------------------------------------------------------
  .blockSpecs <- list()
  .tpL <- character(0)   # transformed-parameter lines building L_b
  .priorL <- character(0)
  .gqOmega <- character(0)
  for (.b in seq_along(map$blocks)) {
    .blk <- map$blocks[[.b]]
    .sp <- .stanOmegaBlockSpec(.blk, pri$omega, ctl)
    .sp$block <- .blk
    .blockSpecs[[.b]] <- .sp
    .k <- .blk$k
    .id <- paste0("_b", .b)
    if (.sp$type == "fixed") {
      # fixed-variance eta: L is a constant, nothing is declared or given a
      # prior; z is still sampled, so the eta remains a random effect with
      # a KNOWN variance (the IOV hook's unit-variance occasion effects)
      .tpL <- c(.tpL, paste0("  matrix[1,1] L", .id, " = [[",
                             .stanNum(sqrt(.blk$init[1, 1])), "]];"))
    } else if (.sp$type == "cov") {
      .decl <- c(.decl, paste0("  cov_matrix[", .k, "] omega", .id, ";"))
      .priorL <- c(.priorL, paste0("  omega", .id, " ~ ", .sp$stanName, "(",
                                   .stanNum(.sp$nu), ", ",
                                   .stanMatLit(.sp$scale), ");"))
      .tpL <- c(.tpL, paste0("  matrix[", .k, ",", .k, "] L", .id,
                             " = cholesky_decompose(omega", .id, ");"))
    } else if (.sp$type == "sd") {
      .decl <- c(.decl, paste0("  real<lower=0> sd", .id, ";"))
      .sdPrior <- sprintf(ctl$diagOmegaSdPrior,
                          .stanNum(2.5 * sqrt(.blk$init[1, 1])))
      .priorL <- c(.priorL, paste0("  sd", .id, " ~ ", .sdPrior, ";"))
      if (.sp$default) {
        .notes <- c(.notes, paste0("default prior sd", .id, " ~ ", .sdPrior,
                                   " on omega diagonal '", .blk$members, "'"))
      }
      .tpL <- c(.tpL, paste0("  matrix[1,1] L", .id, " = [[sd", .id, "]];"))
    } else { # lkj
      .decl <- c(.decl,
                 paste0("  cholesky_factor_corr[", .k, "] Lcorr", .id, ";"),
                 paste0("  vector<lower=0>[", .k, "] sd", .id, ";"))
      .priorL <- c(.priorL, paste0("  Lcorr", .id, " ~ lkj_corr_cholesky(",
                                   .stanNum(.sp$lkjEta), ");"))
      for (.j in seq_len(.k)) {
        .sdPrior <- sprintf(ctl$diagOmegaSdPrior,
                            .stanNum(2.5 * sqrt(.blk$init[.j, .j])))
        .priorL <- c(.priorL, paste0("  sd", .id, "[", .j, "] ~ ", .sdPrior, ";"))
      }
      if (.sp$default) {
        .notes <- c(.notes, paste0("default priors lkj_corr_cholesky(",
                                   .stanNum(.sp$lkjEta), ") + half-Cauchy SDs",
                                   " on omega block '",
                                   paste(.blk$members, collapse = ","), "'"))
      }
      .tpL <- c(.tpL, paste0("  matrix[", .k, ",", .k, "] L", .id,
                             " = diag_pre_multiply(sd", .id, ", Lcorr", .id, ");"))
    }
    if (identical(ctl$etaParam, "centered")) {
      # centred: eta sampled directly; multi_normal_cholesky supplies the
      # (Omega-dependent) density, so no Jacobian bookkeeping is needed --
      # the declared parameter is the one the prior is written on
      .decl <- c(.decl, paste0("  matrix[N, ", .k, "] etaP", .id, ";"))
      .priorL <- c(.priorL,
                   paste0("  for (i in 1:N) etaP", .id,
                          "[i] ~ multi_normal_cholesky(rep_vector(0, ", .k,
                          "), L", .id, ");"))
      .tpL <- c(.tpL, paste0("  eta[, ", .blk$start, ":", .blk$end,
                             "] = etaP", .id, ";"))
    } else {
      .decl <- c(.decl, paste0("  matrix[N, ", .k, "] z", .id, ";"))
      .priorL <- c(.priorL, paste0("  to_vector(z", .id, ") ~ std_normal();"))
      .tpL <- c(.tpL, paste0("  eta[, ", .blk$start, ":", .blk$end,
                             "] = z", .id, " * L", .id, "';"))
    }
    .gqOmega <- c(.gqOmega, paste0("  omegaOut[", .blk$start, ":", .blk$end,
                                   ", ", .blk$start, ":", .blk$end,
                                   "] = L", .id, " * L", .id, "';"))
  }
  # ---- theta vector -------------------------------------------------------
  .thL <- paste0("  vector[", nrow(map$theta), "] theta;")
  for (.i in seq_len(nrow(map$theta))) {
    .r <- map$theta[.i, ]
    .thL <- c(.thL, paste0("  theta[", .i, "] = ",
                           if (.r$fix) .stanNum(.r$est) else .r$par, ";"))
  }
  # ---- theta prior statements (declared parameters only) ------------------
  .thPrior <- character(0)
  for (.i in seq_len(nrow(pri$pop))) {
    .r <- pri$pop[.i, ]
    if (.r$kind == "multivariate") {
      .members <- strsplit(.r$members, ",")[[1]]
      if (any(map$theta$fix[match(.members, map$theta$name)])) {
        # dropping the statement would silently leave the FREE members with a
        # flat prior; representing it would need the conditional normal given
        # the fixed member -- refuse rather than mis-state the model
        stop("the multivariate prior on (", .r$members, ") includes a fixed ",
             "parameter; est=\"stan\" cannot honor it -- unfix the parameter ",
             "or use univariate priors", call. = FALSE)
      }
    } else if (isTRUE(map$theta$fix[match(.r$name, map$theta$name)])) {
      .notes <- c(.notes, paste0("prior on fixed parameter '", .r$name,
                                 "' skipped"))
      next
    }
    .thPrior <- c(.thPrior, paste0("  ", .r$statement))
  }
  # ---- estimated transform-both-sides lambda: DV-transform Jacobian ------
  # The linked conditional evaluates the TRANSFORMED-scale density and omits
  # the Jacobian sum(log |dh(y; lambda)/dy|).  For Box-Cox/Yeo-Johnson that
  # is (lambda - 1) * sJ with sJ a pure DATA statistic per endpoint
  # (Box-Cox: sum(log y); Yeo-Johnson: sum(log1p(y)) over y >= 0 minus
  # sum(log1p(-y)) over y < 0), so it is emitted STAN-SIDE: lambda's full
  # gradient is then exact -- d(cond)/dlambda from the linked sensitivities
  # (nlmixr2/nlmixr2est#949) plus autodiff of this term.  The est layer
  # fills the statistic into the data list (.stanTbsJacValues).
  .iniDf <- ui$iniDf
  .tbsRows <- which(.iniDf$err %in% c("boxCox", "yeoJohnson") &
                      !.iniDf$fix & !is.na(.iniDf$ntheta))
  .tbsDecl <- character(0)
  .tbsTarget <- character(0)
  .tbsJac <- list()
  for (.i in .tbsRows) {
    .p <- match(.iniDf$name[.i], map$theta$name)
    .dn <- paste0("sumLogJac_", .p)
    .tbsDecl <- c(.tbsDecl,
                  paste0("  real ", .dn,
                         ";  // DV-transform Jacobian statistic"))
    .tbsTarget <- c(.tbsTarget, paste0("  target += (", map$theta$par[.p],
                                       " - 1) * ", .dn, ";"))
    .tbsJac[[length(.tbsJac) + 1L]] <-
      list(name = .dn, theta = .iniDf$name[.i],
           transform = .iniDf$err[.i], condition = .iniDf$condition[.i])
  }
  # ---- assemble -----------------------------------------------------------
  if (nrow(map$eta) == 0L) {
    # tier 0: population-only (no etas) via the nlm C API
    # (nlmixr2/nlmixr2est#953) -- one external scalar carrying the whole
    # data log-likelihood and its analytic theta gradient.  The nlm problem
    # estimates only FREE thetas, while the program's theta vector carries
    # fixed rows as literals -- a length mismatch the C entry would reject
    # at the first evaluation (after the compile), so refuse it here: with
    # literalFix=TRUE (the default) fixed thetas dissolve before this point
    if (any(map$theta$fix)) {
      stop("population-only est=\"stan\" needs fixed thetas dissolved: ",
           "keep literalFix=TRUE (the nlm problem estimates only free ",
           "thetas, and the program's theta vector would not match)",
           call. = FALSE)
    }
    .code <- paste(c(
                     "functions {",
                     "  // external (allow_undefined): log p(y | theta) + analytic",
                     "  // d/d(theta), via nlmixr2est's nlm population likelihood",
                     "  real nlmixr2_pop_ll(vector theta);",
                     "}",
                     "data {",
                     "  int<lower=1> N;",
                     .tbsDecl,
                     "}",
                     "parameters {",
                     .decl,
                     "}",
                     "transformed parameters {",
                     .thL,
                     "}",
                     "model {",
                     "  // priors from ini({}); the external function returns the",
                     "  // COMPLETE population data log-likelihood",
                     .thPrior,
                     "  target += nlmixr2_pop_ll(theta);",
                     .tbsTarget,
                     "}"), collapse = "\n")
    return(list(code = .code, data = list(N = NA_integer_), notes = .notes,
                blockSpecs = .blockSpecs, tbsJac = .tbsJac, pop = TRUE))
  }
  .code <- paste(c(
                   "functions {",
                   "  // external (allow_undefined): value + analytic d/d(eta) and",
                   "  // d/d(theta) supplied by the linked rxode2/nlmixr2est likelihood",
                   "  vector nlmixr2_cond_all2(matrix etaMat, vector theta);",
                   "}",
                   "data {",
                   "  int<lower=1> N;",
                   .tbsDecl,
                   "}",
                   "parameters {",
                   .decl,
                   "}",
                   "transformed parameters {",
                   paste0("  matrix[N, ", nrow(map$eta), "] eta;"),
                   .thL,
                   .tpL,
                   "}",
                   "model {",
                   "  // priors from ini({}); Stan owns the FULL, normalized eta prior --",
                   "  // the external function returns the CONDITIONAL log p(y_i|eta_i)",
                   .thPrior,
                   .priorL,
                   "  target += sum(nlmixr2_cond_all2(eta, theta));",
                   .tbsTarget,
                   "}",
                   "generated quantities {",
                   paste0("  matrix[", nrow(map$eta), ",", nrow(map$eta),
                          "] omegaOut = rep_matrix(0, ", nrow(map$eta), ",",
                          nrow(map$eta), ");"),
                   "  vector[N] logLikSubj;",
                   .gqOmega,
                   "  logLikSubj = nlmixr2_cond_all2(eta, theta);",
                   "}"), collapse = "\n")
  list(code = .code, data = list(N = NA_integer_), notes = .notes,
       blockSpecs = .blockSpecs, tbsJac = .tbsJac, pop = FALSE)
}

#' Per-chain initial values on the constrained scale
#' @noRd
.stanInit <- function(map, blockSpecs, nid, ctl) {
  .free <- map$theta[!map$theta$fix, , drop = FALSE]
  lapply(seq_len(ctl$chains), function(.chain) {
    .jit <- if (.chain == 1L) 0 else ctl$initJitterSd
    .ini <- list()
    for (.i in seq_len(nrow(.free))) {
      .r <- .free[.i, ]
      .v <- .r$est + .jit * stats::rnorm(1) * max(abs(.r$est), 0.1)
      # keep strictly inside the declared bounds
      if (is.finite(.r$lower)) .v <- max(.v, .r$lower + 1e-4 * max(1, abs(.r$lower)))
      if (is.finite(.r$upper)) .v <- min(.v, .r$upper - 1e-4 * max(1, abs(.r$upper)))
      .ini[[.r$par]] <- .v
    }
    for (.b in seq_along(blockSpecs)) {
      .sp <- blockSpecs[[.b]]
      .id <- paste0("_b", .b)
      .k <- .sp$block$k
      if (.sp$type == "fixed") {
        # constant L: nothing to initialize beyond z below
      } else if (.sp$type == "cov") {
        .ini[[paste0("omega", .id)]] <- .sp$block$init
      } else if (.sp$type == "sd") {
        .ini[[paste0("sd", .id)]] <- sqrt(.sp$block$init[1, 1])
      } else {
        .ini[[paste0("sd", .id)]] <- sqrt(diag(.sp$block$init))
        .ini[[paste0("Lcorr", .id)]] <- diag(.k)
      }
      if (identical(ctl$etaParam, "centered")) {
        .ini[[paste0("etaP", .id)]] <- matrix(.jit * stats::rnorm(nid * .k),
                                              nid, .k)
      } else {
        .ini[[paste0("z", .id)]] <- matrix(.jit * stats::rnorm(nid * .k),
                                           nid, .k)
      }
    }
    .ini
  })
}
