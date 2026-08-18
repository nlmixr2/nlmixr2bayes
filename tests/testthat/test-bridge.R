library(testthat)

test_that("rxode2's C function-pointer table is installed", {
  probe <- rxsProbeRxode2()
  expect_true(all(probe), info = paste(names(probe)[!probe], collapse = ", "))
})

## A 1-compartment oral PK model on the log scale, which is how Stan will
## parameterize it.  cp = center / v is deliberately NOT the bridge output:
## rxode2 emits sensitivities for states only, so cp is formed in Stan.
pkModel <- "
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
cp <- center / v
"

pkHandle <- function(...) {
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by = 1.5))
  rxsRegister(pkModel, events = ev,
              sens = c("lka", "lcl", "lv"),
              output = "center",
              atol = 1e-10, rtol = 1e-10, ...)
}

test_that("a handle reports the layout it was registered with", {
  h <- pkHandle()
  on.exit(rxsRelease(h))

  expect_equal(attr(h, "np"), 3L)
  expect_equal(attr(h, "ny"), attr(h, "nobs") * 1L)
  expect_equal(attr(h, "sens"), c("lka", "lcl", "lv"))
})

test_that("the C ABI returns rxode2's analytic sensitivities", {
  h <- pkHandle()
  on.exit(rxsRelease(h))

  p <- c(log(1.1), log(4), log(30))
  got <- rxsSolve(h, p)

  expect_equal(dim(got), c(attr(h, "ny"), 4L))

  ## Independent oracle: central differences over whole solves.
  fd <- vapply(seq_along(p), function(j) {
    hh <- 1e-6 * max(1, abs(p[j]))
    pp <- p; pp[j] <- pp[j] + hh
    pm <- p; pm[j] <- pm[j] - hh
    (rxsSolve(h, pp)[, 1L] - rxsSolve(h, pm)[, 1L]) / (2 * hh)
  }, numeric(attr(h, "ny")))

  expect_equal(got[, -1L], fd, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("an unusable handle is an error, not a crash", {
  expect_error(rxsSolve(99999L, c(1, 2, 3)), "unknown rxstan handle")
})

test_that("an R-level error during the solve is caught, not longjmp'd", {
  ## This is the path that, unguarded, would longjmp straight out of Stan's
  ## log_prob() and skip the destructors releasing the autodiff arena.
  ## R_tryCatchError in rxs_solve_sens has to turn it into a return code.
  h <- pkHandle()
  on.exit(rxsRelease(h))

  ## (a) rxode2 itself refuses the parameters.
  expect_error(rxsSolve(h, c(NaN, 1, 3)), "rxs_solve_sens failed")

  ## (b) the R-side registry entry disappears while the C-side dims survive.
  ## Path A never consults the R registry, so it has to be disarmed first for
  ## this to exercise the R error route at all.
  env <- asNamespace("nlmixr2bayes")$.rxsEnv
  key <- as.character(unclass(h))
  saved <- env$handles[[key]]
  env$handles[[key]] <- NULL
  rxsInvalidateFast()
  on.exit(env$handles[[key]] <- saved, add = TRUE, after = FALSE)

  expect_error(rxsSolve(h, c(0, 1, 3)), "unknown rxstan handle")

  ## R must still be healthy afterwards.
  expect_equal(sum(1:10), 55L)
})

test_that("wrong parameter count is rejected", {
  h <- pkHandle()
  on.exit(rxsRelease(h))
  expect_error(rxsSolve(h, c(1, 2)), "expected 3 parameters")
})
