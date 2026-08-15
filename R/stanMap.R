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
  list(theta = .theta, eta = .eta, blocks = .blocks, muRefIdx = .muRefIdx)
}

#' Every estimated theta must have a gradient source: mu-reference or the
#' forward-sensitivity index.  A theta with neither would get a silent zero
#' gradient -- refuse instead.
#' @noRd
.stanAssertThetaGradCover <- function(map, thetaSensIdx) {
  .free <- which(!map$theta$fix)
  .muCovered <- map$muRefIdx[map$muRefIdx > 0L]
  .sensCovered <- thetaSensIdx
  .uncovered <- setdiff(.free, union(.muCovered, .sensCovered))
  if (length(.uncovered) > 0L) {
    stop("no gradient source (mu-reference or theta sensitivity) for ",
         "estimated parameter(s) ",
         paste0("'", map$theta$name[.uncovered], "'", collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}
