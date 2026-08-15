# G3c: tie the linked conditional to nlmixr2's own exact marginal.  The
# marginal likelihood is computed OUTSIDE Stan by quadrature over eta using
# only the linked conditional value + the normalized eta prior, and compared
# against est="agq" (nAGQ=101) evaluated at the same theta.
#
# Independence boundary, stated honestly: the two MARGINALIZATIONS are
# independent (a naive 401-point trapezoid + R's dnorm prior here vs agq's
# inner optimization + adaptive Gauss-Hermite there), so this catches any
# bug in the link path (shim, batch assembly, theta handling, eta-prior
# bookkeeping) and in either integration.  The INTEGRAND is shared -- both
# sides evaluate likInner0 -- so a defect inside likInner0's density itself
# would cancel here; the absolute anchor for the integrand is
# test-censor-ll.R's hand-computed textbook densities (and upstream's
# NONMEM validation).

test_that("quadrature over the linked conditional reproduces est=\"agq\" (G3c)", {
  skip_on_cran()
  .d <- .linkData()
  h <- stanLinkSetup(.linkMod, .d, cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  .linkSetTheta(h$initPar)
  .om <- 0.1
  # 401-point trapezoid over +/- 6 sd: quadrature error ~1e-9 for this
  # smooth 1-D integrand
  .gr <- seq(-6 * sqrt(.om), 6 * sqrt(.om), length.out = 401)
  .cond <- vapply(.gr, function(e) .condBatch(matrix(e, 4, 1))$value,
                  numeric(4))
  .w <- diff(.gr)[1]
  .marg <- apply(.cond, 1, function(v) {
    .m <- max(v)
    log(sum(exp(v - .m) * stats::dnorm(.gr, 0, sqrt(.om))) * .w) + .m
  })
  stanLinkFree()
  .fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    .linkMod, .d, est = "agq",
    control = nlmixr2est::foceiControl(nAGQ = 101L, maxOuterIterations = 0L,
                                       maxInnerIterations = 100L,
                                       covMethod = "", calcTables = FALSE,
                                       print = 0))))
  # both are proper -2 log-likelihoods (no NONMEM 2*pi offset); measured
  # agreement 3e-7, gated at 1e-5
  expect_lt(abs(-2 * sum(.marg) - .fit$objective), 1e-5)
})
