## quiet = TRUE silences rxode2's repeated per-solve diagnostics but must not
## hide a failure: the error still propagates and the counters still count.
library(testthat)

pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
"

mk <- function(quiet) {
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by = 1.5))
  rxsRegister(pkModel, events = ev, sens = c("lka", "lcl", "lv"),
              output = "center", atol = 1e-10, rtol = 1e-10, quiet = quiet)
}

good <- c(log(1.1), log(4), log(30))
bad <- c(3, 30, -20)

test_that("quiet suppresses solver chatter", {
  h <- mk(TRUE)
  on.exit(rxsRelease(h))
  out <- utils::capture.output(
    suppressWarnings(try(rxsSolve(h, bad), silent = TRUE)),
    type = "output")
  expect_false(any(grepl("lsoda|intdy", out)))
})

test_that("a failure is still an error when quiet", {
  h <- mk(TRUE)
  on.exit(rxsRelease(h))
  expect_error(suppressWarnings(rxsSolve(h, bad)), "rxs_solve_sens failed")
})

test_that("quiet does not change the numbers", {
  hq <- mk(TRUE)
  hl <- mk(FALSE)
  on.exit({ rxsRelease(hq); rxsRelease(hl) })
  expect_identical(rxsSolve(hq, good), rxsSolve(hl, good))
})

test_that("failures are counted even when silenced", {
  h <- mk(TRUE)
  on.exit(rxsRelease(h))

  rxsResetStats(h)
  expect_equal(unname(rxsSolveStats(h)), c(0, 0))

  rxsSolve(h, good)
  expect_equal(unname(rxsSolveStats(h)["solves"]), 1)
  expect_equal(unname(rxsSolveStats(h)["failures"]), 0)

  suppressWarnings(try(rxsSolve(h, bad), silent = TRUE))
  st <- rxsSolveStats(h)
  expect_equal(unname(st["solves"]), 2)
  expect_gte(unname(st["failures"]), 1)
})
