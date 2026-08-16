# Pathfinder (Zhang, Carpenter, Gelman & Vehtari 2022, JMLR 23(306)) in
# pure R against the compiled Stan model's exact log density and analytic
# gradient (rstan::log_prob / grad_log_prob).  rstan does not expose the
# Pathfinder service (StanHeaders 2.32 predates it) and CmdStan cannot
# reach the in-process linked likelihood at all (separate executable), so
# the algorithm itself runs here:
#
#   1. L-BFGS ascent on log p from a jittered start, keeping the whole
#      trajectory and the (s, y) curvature pairs;
#   2. at each iterate, a local Gaussian N(theta_k, H_k) with H_k the
#      dense L-BFGS inverse-Hessian estimate (diagonal init gamma_k I,
#      last m updates applied explicitly -- our unconstrained dimensions
#      are small enough that dense algebra is trivial);
#   3. Monte-Carlo ELBO for each iterate's Gaussian; the best iterate is
#      the path's approximation (stopping BEFORE the mode is what lets
#      Pathfinder avoid collapsing into funnel necks);
#   4. multi-path: draws pooled over paths against the mixture density,
#      Pareto-smoothed importance resampling (loo::psis) with the khat
#      diagnostic.
#
# Everything below .stanRunPathfinder is pure R over fn/gr closures, so
# the algorithm is unit-testable against analytic targets without Stan.

# ---- L-BFGS with trajectory ------------------------------------------------

#' Dense inverse-Hessian estimate from the last m (s, y) pairs
#' (recursive BFGS inverse update from a scaled diagonal), returned as a
#' Cholesky factor; NULL when the curvature is unusable
#' @noRd
.pfInvHessChol <- function(sList, yList, n) {
  .m <- length(sList)
  if (.m == 0L) return(NULL)
  .s <- sList[[.m]]
  .y <- yList[[.m]]
  .sy <- sum(.s * .y)
  if (!is.finite(.sy) || .sy <= 0) return(NULL)
  .gamma <- .sy / sum(.y * .y)
  .h <- diag(.gamma, n)
  for (.k in seq_len(.m)) {
    .s <- sList[[.k]]
    .y <- yList[[.k]]
    .rho <- 1 / sum(.s * .y)
    if (!is.finite(.rho) || .rho <= 0) next
    # H <- (I - rho s y') H (I - rho y s') + rho s s'
    .hy <- .h %*% .y
    .h <- .h - .rho * (.s %*% t(.hy) + .hy %*% t(.s)) +
      (.rho^2 * sum(.y * .hy) + .rho) * (.s %*% t(.s))
  }
  .h <- (.h + t(.h)) / 2
  .ch <- tryCatch(chol(.h), error = function(e) NULL)
  .ch
}

#' log density of N(mu, LL') at the columns of z (L upper-tri chol of the
#' covariance, as chol() returns)
#' @noRd
.pfMvnLogd <- function(z, mu, cholCov) {
  .n <- length(mu)
  .ld <- sum(log(diag(cholCov)))
  .q <- backsolve(cholCov, t(z) - mu, transpose = TRUE)
  -0.5 * colSums(.q^2) - .ld - 0.5 * .n * log(2 * pi)
}

#' draws x n matrix from N(mu, LL')
#' @noRd
.pfMvnDraw <- function(nDraw, mu, cholCov) {
  .n <- length(mu)
  .z <- matrix(stats::rnorm(nDraw * .n), nDraw, .n)
  sweep(.z %*% cholCov, 2, mu, "+")
}

#' One Pathfinder path: L-BFGS trajectory + per-iterate Gaussian + ELBO
#' selection.  fn(x) is the log density (-Inf allowed), gr(x) its
#' gradient.  Returns NULL when no usable approximation was found.
#' @noRd
.pathfinderOne <- function(fn, gr, x0, maxIter = 100L, history = 6L,
                           elboDraws = 25L, tolGrad = 1e-8) {
  .n <- length(x0)
  .x <- x0
  .f <- fn(.x)
  if (!is.finite(.f)) return(NULL)
  .g <- gr(.x)
  if (any(!is.finite(.g))) return(NULL)
  .sList <- list()
  .yList <- list()
  # candidate approximations along the path
  .cand <- list()
  for (.it in seq_len(maxIter)) {
    # ascent direction from the current inverse-Hessian estimate
    .ch <- .pfInvHessChol(.sList, .yList, .n)
    .dir <- if (is.null(.ch)) {
      .g / max(1, sqrt(sum(.g^2)))
    } else {
      # H g via the factor
      as.numeric(t(.ch) %*% (.ch %*% .g))
    }
    if (sum(.dir * .g) <= 0) .dir <- .g # ensure ascent
    # backtracking Armijo line search (maximization)
    .step <- 1
    .ok <- FALSE
    for (.ls in 1:30) {
      .xn <- .x + .step * .dir
      .fn2 <- fn(.xn)
      if (is.finite(.fn2) && .fn2 >= .f + 1e-4 * .step * sum(.g * .dir)) {
        .ok <- TRUE
        break
      }
      .step <- .step / 2
    }
    if (!.ok) break
    .gn <- gr(.xn)
    if (any(!is.finite(.gn))) break
    .s <- .xn - .x
    # y on the NEGATIVE log density scale so s'y > 0 near a maximum
    .y <- .g - .gn
    if (sum(.s * .y) > 1e-12 * sqrt(sum(.s^2)) * sqrt(sum(.y^2))) {
      .sList <- c(.sList, list(.s))
      .yList <- c(.yList, list(.y))
      if (length(.sList) > history) {
        .sList <- .sList[-1L]
        .yList <- .yList[-1L]
      }
    }
    .x <- .xn
    .f <- .fn2
    .g <- .gn
    # local Gaussian candidate at this iterate
    .ch <- .pfInvHessChol(.sList, .yList, .n)
    if (!is.null(.ch)) {
      .z <- .pfMvnDraw(elboDraws, .x, .ch)
      .lq <- .pfMvnLogd(.z, .x, .ch)
      .lp <- apply(.z, 1, fn)
      .keep <- is.finite(.lp)
      if (sum(.keep) >= max(2L, elboDraws %/% 2L)) {
        .elbo <- mean(.lp[.keep] - .lq[.keep])
        .cand[[length(.cand) + 1L]] <- list(mu = .x, chol = .ch,
                                            elbo = .elbo)
      }
    }
    if (sqrt(sum(.g^2)) < tolGrad) break
  }
  if (length(.cand) == 0L) return(NULL)
  .best <- which.max(vapply(.cand, function(c) c$elbo, numeric(1)))
  .cand[[.best]]
}

#' Multi-path Pathfinder over fn/gr: pooled draws against the path
#' mixture, PSIS-smoothed importance resampling.  Returns
#' list(draws, khat, lp, nPathsOk) or NULL when every path failed.
#' @noRd
.pathfinderMulti <- function(fn, gr, x0, paths = 4L, jitterSd = 2,
                             drawsPerPath = 1000L, nDraws = 1000L,
                             maxIter = 100L, history = 6L,
                             elboDraws = 25L) {
  .n <- length(x0)
  .apx <- list()
  for (.j in seq_len(paths)) {
    .xj <- if (.j == 1L) x0 else x0 + stats::rnorm(.n, 0, jitterSd)
    .a <- .pathfinderOne(fn, gr, .xj, maxIter = maxIter,
                         history = history, elboDraws = elboDraws)
    if (!is.null(.a)) .apx[[length(.apx) + 1L]] <- .a
  }
  if (length(.apx) == 0L) return(NULL)
  .nap <- length(.apx)
  # pooled draws + mixture density of the selected Gaussians
  .z <- do.call(rbind, lapply(.apx, function(a) {
    .pfMvnDraw(drawsPerPath, a$mu, a$chol)
  }))
  .lqm <- matrix(0, nrow(.z), .nap)
  for (.k in seq_len(.nap)) {
    .lqm[, .k] <- .pfMvnLogd(.z, .apx[[.k]]$mu, .apx[[.k]]$chol)
  }
  .lq <- apply(.lqm, 1, function(r) {
    .m <- max(r)
    .m + log(mean(exp(r - .m)))
  })
  .lp <- apply(.z, 1, fn)
  .keep <- is.finite(.lp)
  .z <- .z[.keep, , drop = FALSE]
  .lw <- .lp[.keep] - .lq[.keep]
  if (nrow(.z) < 10L) return(NULL)
  # PSIS smoothing + khat (loo is a hard requirement of the pathfinder
  # algorithm; it is in Suggests for the LOO fit rows already)
  .khat <- NA_real_
  .w <- .lw - max(.lw)
  if (requireNamespace("loo", quietly = TRUE)) {
    .ps <- tryCatch(suppressWarnings(loo::psis(.lw, r_eff = NA)),
                    error = function(e) NULL)
    if (!is.null(.ps)) {
      .khat <- .ps$diagnostics$pareto_k
      .w <- as.numeric(stats::weights(.ps, log = TRUE, normalize = TRUE))
    }
  }
  .idx <- sample.int(nrow(.z), size = nDraws, replace = TRUE,
                     prob = exp(.w - max(.w)))
  list(draws = .z[.idx, , drop = FALSE], khat = .khat,
       lp = .lp[.keep][.idx], nPathsOk = .nap)
}

# ---- the est="stan" runner -------------------------------------------------

#' Run Pathfinder against the compiled model and wrap the draws so the
#' shared posterior->fit machinery can consume them
#' @noRd
.stanRunPathfinder <- function(sm, gen, map, nid, control, init) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("algorithm=\"pathfinder\" needs the loo package for ",
         "Pareto-smoothed importance resampling", call. = FALSE)
  }
  # a minimal stanfit to expose log_prob/grad_log_prob (2 iterations from
  # the ini() values; milliseconds)
  .i1 <- if (is.list(init) && length(init) >= 1L && is.list(init[[1L]])) {
    init[[1L]]
  } else {
    init
  }
  .f0 <- rstan::sampling(sm, data = gen$data, chains = 1, iter = 2,
                         warmup = 1, refresh = 0, seed = control$seed,
                         init = list(.i1))
  .up0 <- rstan::unconstrain_pars(.f0, .i1)
  .fn <- function(x) {
    tryCatch(rstan::log_prob(.f0, x, adjust_transform = TRUE,
                             gradient = FALSE),
             error = function(e) -Inf)
  }
  .gr <- function(x) {
    tryCatch(as.numeric(rstan::grad_log_prob(.f0, x,
                                             adjust_transform = TRUE)),
             error = function(e) rep(NaN, length(x)))
  }
  .res <- .pathfinderMulti(.fn, .gr, .up0,
                           paths = control$pathfinderPaths,
                           jitterSd = max(control$initJitterSd, 0.5),
                           drawsPerPath = control$vbOutputSamples,
                           nDraws = control$vbOutputSamples,
                           maxIter = 250L, history = 6L,
                           elboDraws = 25L)
  if (is.null(.res)) {
    stop("Pathfinder failed on every path (no usable local Gaussian); ",
         "use algorithm=\"NUTS\"", call. = FALSE)
  }
  .pf <- list(fit0 = .f0, draws = .res$draws, khat = .res$khat,
              lp = .res$lp, nPathsOk = .res$nPathsOk, gen = gen,
              cache = new.env(parent = emptyenv()))
  class(.pf) <- "nlmixr2stanPathfinder"
  .pf
}

# ---- draw access for the shared finalize ----------------------------------

#' Force the draw caches (constrained draws + generated quantities) while
#' the linked likelihood is still up -- the GQ block calls the external
#' function, which is torn down before finalize
#' @noRd
.pathfinderPrefetch <- function(pf) {
  .pathfinderConstrain(pf)
  if (is.null(pf$cache$gq)) {
    pf$cache$gq <- tryCatch(.pathfinderGqs(pf), error = function(e) NULL)
  }
  invisible(pf)
}

#' rstan::extract-shaped access for a Pathfinder result: constrain each
#' draw for parameters/transformed parameters; run the model's generated
#' quantities (omegaOut, logLikSubj, mixProbOut) through rstan::gqs so
#' the exact generated code produces them
#' @noRd
.pathfinderConstrain <- function(pf) {
  if (!is.null(pf$cache$constrained)) return(pf$cache$constrained)
  .f0 <- pf$fit0
  .nd <- nrow(pf$draws)
  .cp1 <- rstan::constrain_pars(.f0, pf$draws[1L, ])
  .out <- lapply(.cp1, function(v) {
    array(NA_real_, dim = c(.nd, dim(as.array(v))))
  })
  for (.d in seq_len(.nd)) {
    .cp <- rstan::constrain_pars(.f0, pf$draws[.d, ])
    for (.nm in names(.cp)) {
      .v <- as.array(.cp[[.nm]])
      if (length(dim(.v)) == 1L) {
        .out[[.nm]][.d, ] <- .v
      } else {
        .out[[.nm]][.d, , ] <- .v
      }
    }
  }
  pf$cache$constrained <- .out
  .out
}

#' extract() equivalent over a Pathfinder result
#' @noRd
.stanExtract <- function(sf, pars) {
  if (inherits(sf, "nlmixr2stanPathfinder")) {
    .con <- .pathfinderConstrain(sf)
    .miss <- setdiff(pars, names(.con))
    if (length(.miss) > 0L) {
      if (is.null(sf$cache$gq)) sf$cache$gq <- .pathfinderGqs(sf)
      .con <- c(.con, sf$cache$gq)
      .still <- setdiff(pars, names(.con))
      if (length(.still) > 0L) {
        stop("Pathfinder draws lack ", paste(.still, collapse = ", "),
             call. = FALSE) # nocov
      }
    }
    return(.con[pars])
  }
  rstan::extract(sf, pars = pars)
}

#' generated quantities for Pathfinder draws via rstan::gqs (the flat
#' free-parameter draws matrix feeds the model's own GQ block)
#' @noRd
.pathfinderGqs <- function(pf) {
  .f0 <- pf$fit0
  # flat draws of the PARAMETERS block only, in flatname order
  .parNames <- .f0@model_pars
  .con <- .pathfinderConstrain(pf)
  .nd <- nrow(pf$draws)
  # rebuild the flat matrix from the unconstrained draws directly: the
  # constrained parameter values with stan's flat naming
  .flat <- NULL
  .fn <- character(0)
  for (.nm in .parNames) {
    .v <- .con[[.nm]]
    if (is.null(.v)) next # transformed/gq names resolve later
    .dims <- dim(.v)
    if (length(.dims) == 2L) {
      .m <- .v
      .cn <- if (.dims[2] == 1L) .nm else paste0(.nm, "[", seq_len(.dims[2]), "]")
    } else {
      .m <- matrix(.v, .nd, prod(.dims[-1]))
      .cn <- paste0(.nm, "[",
                    apply(expand.grid(lapply(.dims[-1], seq_len)), 1,
                          paste, collapse = ","), "]")
    }
    colnames(.m) <- .cn
    .flat <- if (is.null(.flat)) .m else cbind(.flat, .m)
    .fn <- c(.fn, .cn)
  }
  .gq <- rstan::gqs(rstan::get_stanmodel(.f0), data = pf$gen$data,
                    draws = .flat)
  rstan::extract(.gq)
}

#' summary()$summary equivalent (mean/sd/quantiles; no Rhat/ESS -- there
#' are no chains)
#' @noRd
.stanSummaryDf <- function(sf) {
  if (!inherits(sf, "nlmixr2stanPathfinder")) {
    return(rstan::summary(sf)$summary)
  }
  .con <- .pathfinderConstrain(sf)
  .rows <- list()
  for (.nm in names(.con)) {
    .v <- .con[[.nm]]
    .dims <- dim(.v)
    .m <- if (length(.dims) == 2L) .v else matrix(.v, .dims[1], prod(.dims[-1]))
    .cn <- if (ncol(.m) == 1L) .nm else paste0(.nm, "[", seq_len(ncol(.m)), "]")
    for (.k in seq_len(ncol(.m))) {
      .x <- .m[, .k]
      .rows[[.cn[.k]]] <- c(mean = mean(.x), se_mean = NA_real_,
                            sd = stats::sd(.x),
                            stats::quantile(.x, c(0.025, 0.25, 0.5,
                                                  0.75, 0.975)),
                            n_eff = NA_real_, Rhat = NA_real_)
    }
  }
  do.call(rbind, .rows)
}
