# The parameter map (Spec 8): iniDf row order is the contract.  Thetas are
# `!is.na(ntheta)` ordered by ntheta (fixed rows stay in the vector so the
# THETA[] indices match what foceiSetup_ installed, but are not declared as
# Stan parameters); etas are the neta1==neta2 diagonals in neta1 order, with
# omega blocks recovered from the lotri structure (contiguous by
# construction).  muRefIdx[k] is the 1-based theta index mu-referencing eta k
# (0 = none) -- the C shim uses it to fill the mu-referenced theta-gradient
# columns from the eta gradient.

#' Build the Stan parameter map for a ui
#'
#' @param ui rxode2 ui
#' @return list(theta=, eta=, blocks=, muRefIdx=)
#' @noRd
.stanMap <- function(ui) {
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), , drop = FALSE]
  .th <- .th[order(.th$ntheta), , drop = FALSE]
  .theta <- data.frame(name = .th$name, par = .stanParName(.th$name),
                       est = .th$est, lower = .th$lower, upper = .th$upper,
                       fix = .th$fix, ntheta = .th$ntheta,
                       stringsAsFactors = FALSE)
  .et <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, ,
                drop = FALSE]
  .et <- .et[order(.et$neta1), , drop = FALSE]
  .eta <- data.frame(name = .et$name, neta = .et$neta1, est = .et$est,
                     fix = .et$fix, stringsAsFactors = FALSE)
  # omega blocks: split the (block-diagonal) lotri matrix; blocks are
  # contiguous in the neta1 diagonal order
  .om <- ui$omega
  .etaNames <- .eta$name
  .blocks <- list()
  .done <- character(0)
  .lst <- lotri::lotriMatInv(.om)
  for (.b in .lst) {
    .m <- as.matrix(.b)
    .nm <- dimnames(.m)[[1]]
    .idx <- match(.nm, .etaNames)
    .blocks[[length(.blocks) + 1L]] <-
      list(members = .nm, idx = sort(.idx), start = min(.idx),
           end = max(.idx), k = length(.nm),
           init = .m[order(.idx), order(.idx), drop = FALSE])
    .done <- c(.done, .nm)
  }
  # sort blocks by starting eta index; sanity: contiguity
  .blocks <- .blocks[order(vapply(.blocks, function(b) b$start, numeric(1)))]
  for (.b in .blocks) {
    if (!identical(.b$idx, seq(.b$start, .b$end))) {
      stop("omega block over ", paste(.b$members, collapse = ", "),
           " is not contiguous in eta order", call. = FALSE) # nocov
    }
  }
  # mu-reference: theta index per eta (0 = no mu-referenced theta)
  .muRef <- ui$muRefDataFrame
  .muRefIdx <- integer(nrow(.eta))
  if (is.data.frame(.muRef) && nrow(.muRef) > 0L) {
    for (.k in seq_len(nrow(.eta))) {
      .w <- which(.muRef$eta == .eta$name[.k])
      if (length(.w) == 1L) {
        .p <- match(.muRef$theta[.w], .theta$name)
        if (!is.na(.p)) .muRefIdx[.k] <- .p
      }
    }
  }
  # mu-referenced covariate coefficients: for mu_k = theta_p0 + theta_p*cov +
  # eta_k the chain rule gives d/dtheta_p = cov_i * d/deta_k, so each
  # coefficient needs its theta index, the eta it rides on, and (later, from
  # the data) the per-subject covariate value
  .mrc <- ui$muRefCovariateDataFrame
  .muRefCov <- data.frame(thetaIdx = integer(0), etaIdx = integer(0),
                          covariate = character(0), name = character(0),
                          stringsAsFactors = FALSE)
  if (is.data.frame(.mrc) && nrow(.mrc) > 0L) {
    for (.i in seq_len(nrow(.mrc))) {
      .p <- match(.mrc$covariateParameter[.i], .theta$name)
      .w <- which(.muRef$theta == .mrc$theta[.i])
      .k <- if (length(.w) == 1L) match(.muRef$eta[.w], .eta$name) else NA_integer_
      if (is.na(.p) || is.na(.k)) {
        stop("cannot resolve the mu-referenced covariate coefficient '",
             .mrc$covariateParameter[.i], "' (covariate '",
             .mrc$covariate[.i], "') to a theta/eta pair", call. = FALSE) # nocov
      }
      .muRefCov <- rbind(.muRefCov,
                         data.frame(thetaIdx = .p, etaIdx = .k,
                                    covariate = .mrc$covariate[.i],
                                    name = .mrc$covariateParameter[.i],
                                    stringsAsFactors = FALSE))
    }
  }
  list(theta = .theta, eta = .eta, blocks = .blocks, muRefIdx = .muRefIdx,
       muRefCov = .muRefCov)
}

#' Every estimated theta must have a gradient source: mu-reference or the
#' forward-sensitivity index.  A theta with neither would get a silent zero
#' gradient -- refuse instead.
#' @noRd
.stanAssertThetaGradCover <- function(map, thetaSensIdx) {
  .free <- which(!map$theta$fix)
  .muCovered <- map$muRefIdx[map$muRefIdx > 0L]
  .covCovered <- map$muRefCov$thetaIdx
  .sensCovered <- thetaSensIdx
  .uncovered <- setdiff(.free, Reduce(union, list(.muCovered, .covCovered,
                                                  .sensCovered)))
  if (length(.uncovered) > 0L) {
    stop("no gradient source (mu-reference or theta sensitivity) for ",
         "estimated parameter(s) ",
         paste0("'", map$theta$name[.uncovered], "'", collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Per-subject values of the mu-referenced covariates, nid x nrow(muRefCov)
#' in idLvl order.  The scatter identity d/dtheta_p = cov_i * d/deta_k only
#' factors when the covariate is constant within subject, so a time-varying
#' covariate is refused (the rows nlmixr2est's preprocessing adds, e.g. the
#' EVID=9 initialization row, carry NA and are ignored).
#' @noRd
.stanMuRefCovValues <- function(map, dataSav) {
  .mrc <- map$muRefCov
  if (nrow(.mrc) == 0L) return(matrix(0, 0L, 0L))
  .id <- dataSav$ID
  .nid <- length(unique(.id))
  .val <- matrix(0, .nid, nrow(.mrc))
  for (.j in seq_len(nrow(.mrc))) {
    .cn <- .mrc$covariate[.j]
    if (!.cn %in% names(dataSav)) {
      stop("mu-referenced covariate '", .cn, "' is not a data column",
           call. = FALSE) # nocov
    }
    .cv <- dataSav[[.cn]]
    for (.i in seq_len(.nid)) {
      .u <- unique(.cv[.id == .i & !is.na(.cv)])
      if (length(.u) == 0L) {
        stop("mu-referenced covariate '", .cn, "' has no value for ",
             "subject ", .i, call. = FALSE)
      }
      if (length(.u) > 1L) {
        stop("est=\"stan\" needs mu-referenced covariates constant within ",
             "subject and '", .cn, "' varies within subject ", .i,
             ": the coefficient gradient d/d(", .mrc$name[.j],
             ") = ", .cn, " * d/d(eta) only factors for a ",
             "subject-constant covariate", call. = FALSE)
      }
      .val[.i, .j] <- .u
    }
  }
  .val
}
