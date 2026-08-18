## Population models: every subject gets its own block of parameters.  Because
## subjects are independent the Jacobian is block diagonal, so the bridge
## returns ny x length(sens) rather than the mostly-zero dense matrix.
library(testthat)

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
"
sens <- c("lka", "lcl", "lv")

popEvents <- function(nsub, times = c(0.5, 2, 6, 12)) {
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    e <- rxode2::et(amt = 100, cmt = "depot")
    e <- rxode2::et(e, times)
    d <- as.data.frame(e)
    d$id <- i
    d
  }))
}

popHandle <- function(nsub = 3L, perSubject = TRUE, fast = TRUE) {
  rxsRegister(pkModel, events = popEvents(nsub), sens = sens,
              output = "center", perSubject = perSubject, fast = fast,
              atol = 1e-10, rtol = 1e-10)
}

## Subject-major: subject 1's lka, lcl, lv, then subject 2's, ...
popPars <- function(nsub) {
  as.numeric(t(cbind(log(seq(1.0, 1.4, length.out = nsub)),
                     log(seq(3.5, 4.5, length.out = nsub)),
                     log(seq(28, 33, length.out = nsub)))))
}

test_that("a population handle reports a block layout", {
  h <- popHandle(3L)
  on.exit(rxsRelease(h))

  expect_equal(attr(h, "nsub"), 3L)
  expect_equal(attr(h, "nBlock"), 3L)
  expect_equal(attr(h, "nBlocks"), 3L)
  expect_equal(attr(h, "np"), 9L)
  expect_equal(attr(h, "ny"), 12L)
  expect_true(attr(h, "perSubject"))
})

test_that("the Jacobian stays block sized, not dense", {
  h <- popHandle(5L)
  on.exit(rxsRelease(h))

  got <- rxsSolve(h, popPars(5L))
  ## value column + one column per per-subject parameter, NOT per Stan parameter
  expect_equal(ncol(got), length(sens) + 1L)
  expect_equal(nrow(got), attr(h, "ny"))
})

test_that("each subject matches the same subject solved alone", {
  nsub <- 4L
  h <- popHandle(nsub)
  on.exit(rxsRelease(h))

  p <- popPars(nsub)
  pop <- rxsSolve(h, p)
  nt <- attr(h, "nobs") / nsub

  for (s in seq_len(nsub)) {
    hs <- rxsRegister(pkModel, events = popEvents(1L), sens = sens,
                      output = "center", atol = 1e-10, rtol = 1e-10)
    block <- p[((s - 1L) * 3L + 1L):(s * 3L)]
    alone <- rxsSolve(hs, block)
    rows <- ((s - 1L) * nt + 1L):(s * nt)
    expect_equal(pop[rows, ], alone, tolerance = 1e-12, ignore_attr = TRUE)
    rxsRelease(hs)
  }
})

test_that("fast and slow agree for population models", {
  h <- popHandle(4L)
  on.exit(rxsRelease(h))

  set.seed(3)
  for (i in 1:20) {
    p <- popPars(4L) + stats::runif(12, -0.5, 0.5)
    expect_identical(rxsSolve(h, p), rxsSolve(h, p, slow = TRUE))
  }
})

test_that("sensitivities are truly per subject", {
  nsub <- 3L
  h <- popHandle(nsub)
  on.exit(rxsRelease(h))

  p <- popPars(nsub)
  base <- rxsSolve(h, p)[, 1L]

  ## Perturbing subject 2's block must move subject 2's rows and nothing else.
  nt <- attr(h, "nobs") / nsub
  bumped <- p
  bumped[4] <- bumped[4] + 0.1
  moved <- rxsSolve(h, bumped)[, 1L]

  rows2 <- (nt + 1L):(2L * nt)
  expect_gt(max(abs(moved[rows2] - base[rows2])), 1e-6)
  expect_equal(moved[-rows2], base[-rows2], tolerance = 1e-12)
})

test_that("block sensitivities match per-subject finite differences", {
  nsub <- 3L
  h <- popHandle(nsub)
  on.exit(rxsRelease(h))

  p <- popPars(nsub)
  got <- rxsSolve(h, p)
  nt <- attr(h, "nobs") / nsub

  for (j in seq_along(sens)) {
    fd <- numeric(attr(h, "ny"))
    for (s in seq_len(nsub)) {
      k <- (s - 1L) * 3L + j
      step <- 1e-6 * max(1, abs(p[k]))
      pp <- p; pp[k] <- pp[k] + step
      pm <- p; pm[k] <- pm[k] - step
      rows <- ((s - 1L) * nt + 1L):(s * nt)
      fd[rows] <- (rxsSolve(h, pp)[rows, 1L] - rxsSolve(h, pm)[rows, 1L]) /
        (2 * step)
    }
    expect_equal(got[, j + 1L], fd, tolerance = 1e-6)
  }
})

test_that("shared parameters across subjects still work", {
  h <- popHandle(3L, perSubject = FALSE)
  on.exit(rxsRelease(h))

  expect_equal(attr(h, "np"), 3L)
  expect_equal(attr(h, "nBlocks"), 1L)
  expect_equal(attr(h, "ny"), 12L)

  p <- c(log(1.1), log(4), log(30))
  expect_identical(rxsSolve(h, p), rxsSolve(h, p, slow = TRUE))

  ## Same parameters for everyone means identical profiles.
  y <- rxsSolve(h, p)[, 1L]
  expect_equal(y[1:4], y[5:8], tolerance = 1e-12)
})
