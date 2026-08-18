## Path A: re-driving the solve rxode2 already built must give bit-identical
## results to calling rxSolve() again, and must degrade safely when something
## else in the session takes the solve structure away.
library(testthat)

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
"

mkHandle <- function(fast = TRUE, output = "center") {
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by = 1.5))
  rxsRegister(pkModel, events = ev, sens = c("lka", "lcl", "lv"),
              output = output, atol = 1e-10, rtol = 1e-10, fast = fast)
}

p0 <- c(log(1.1), log(4), log(30))

test_that("Path A is armed by registration", {
  h <- mkHandle()
  on.exit(rxsRelease(h))
  expect_true(rxsFastAvailable(h))
})

test_that("fast and slow agree bit for bit", {
  h <- mkHandle()
  on.exit(rxsRelease(h))

  for (p in list(p0, p0 + 0.3, p0 - 0.25, c(0.5, 1.0, 3.0))) {
    fastRes <- rxsSolve(h, p)
    slowRes <- rxsSolve(h, p, slow = TRUE)
    expect_identical(fastRes, slowRes)
  }
})

test_that("multiple output states keep their layout", {
  h <- mkHandle(output = c("depot", "center"))
  on.exit(rxsRelease(h))

  expect_equal(attr(h, "ny"), 2L * attr(h, "nobs"))
  expect_identical(rxsSolve(h, p0), rxsSolve(h, p0, slow = TRUE))
})

test_that("a foreign rxSolve is detected and does not corrupt results", {
  h <- mkHandle()
  on.exit(rxsRelease(h))

  reference <- rxsSolve(h, p0, slow = TRUE)

  ## Something else solves, freeing and rebuilding rxode2's structure.
  other <- rxode2::rxode2("d/dt(x) <- -0.3 * x")
  rxode2::rxSolve(other, params = c(), events = rxode2::et(seq(0, 5, by = 1)),
                  returnType = "data.frame", cores = 1L)
  expect_false(rxsFastAvailable(h))

  ## The next solve must still be right, and must re-arm Path A.
  expect_identical(rxsSolve(h, p0), reference)
  expect_true(rxsFastAvailable(h))
  expect_identical(rxsSolve(h, p0), reference)
})

test_that("interleaving two handles stays correct", {
  h1 <- mkHandle()
  h2 <- mkHandle(output = "depot")
  on.exit({ rxsRelease(h1); rxsRelease(h2) })

  r1 <- rxsSolve(h1, p0, slow = TRUE)
  r2 <- rxsSolve(h2, p0, slow = TRUE)

  for (i in 1:3) {
    expect_identical(rxsSolve(h1, p0), r1)
    expect_identical(rxsSolve(h2, p0), r2)
  }
})

test_that("fast = FALSE never arms Path A", {
  h <- mkHandle(fast = FALSE)
  on.exit(rxsRelease(h))
  expect_false(rxsFastAvailable(h))
  expect_identical(rxsSolve(h, p0), rxsSolve(h, p0, slow = TRUE))
})

test_that("a failing solve on the fast path is still reported, not fatal", {
  h <- mkHandle()
  on.exit(rxsRelease(h))
  rxsSolve(h, p0)
  expect_true(rxsFastAvailable(h))
  expect_error(rxsSolve(h, c(NaN, 1, 3)), "rxs_solve_sens failed")
  expect_equal(sum(1:10), 55L)
})

test_that("a failed solve does not poison later ones", {
  ## A failed solve leaves rxode2's per-subject state dirty, and par_solve()
  ## alone keeps returning that NaN state -- which silently wedged HMC, since
  ## every subsequent proposal was rejected.
  h <- mkHandle()
  on.exit(rxsRelease(h))

  ref <- rxsSolve(h, p0, slow = TRUE)

  for (bad in list(c(3, 30, -20), c(NaN, 1, 3), c(50, -40, -30))) {
    try(rxsSolve(h, bad), silent = TRUE)
    expect_identical(rxsSolve(h, p0), ref)
  }
})

test_that("the fast path survives a long random walk", {
  ## HMC roams far from where the handle was registered, so agreement near p0
  ## proves very little on its own.
  h <- mkHandle()
  on.exit(rxsRelease(h))

  set.seed(7)
  for (i in 1:100) {
    p <- stats::runif(3, -2, 2)
    fastRes <- try(rxsSolve(h, p), silent = TRUE)
    slowRes <- try(rxsSolve(h, p, slow = TRUE), silent = TRUE)
    expect_equal(inherits(fastRes, "try-error"), inherits(slowRes, "try-error"))
    if (!inherits(fastRes, "try-error")) expect_identical(fastRes, slowRes)
  }
})
