## eventSens = "jump": analytic sensitivities for MODELED dosing quantities --
## lag time, bioavailability, duration, rate.  Nothing in Stan's own ODE
## interface can express these at all, so there is no external oracle; the
## checks are finite differences plus the exact pre-lag behavior.
library(testthat)

lagModel <- "
ka   <- exp(lka)
cl   <- exp(lcl)
tlag <- exp(ltlag)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - cl * center
alag(depot)  <- tlag
"

sens <- c("lka", "lcl", "ltlag")
p0 <- c(log(1.1), log(0.15), log(0.5))

mkHandle <- function(...) {
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, c(1, 2, 4, 8, 12))
  rxsRegister(lagModel, events = ev, sens = sens, output = "center",
              eventSens = "jump", atol = 1e-10, rtol = 1e-10, ...)
}

test_that("a modeled lag time carries analytic sensitivities", {
  h <- mkHandle()
  on.exit(rxsRelease(h))

  got <- rxsSolve(h, p0)
  expect_equal(dim(got), c(attr(h, "ny"), 4L))

  ## The lag column must not be trivially zero, or the jump is not happening.
  expect_gt(max(abs(got[, 4L])), 1e-3)

  fd <- vapply(seq_along(p0), function(j) {
    step <- 1e-5
    pp <- p0; pp[j] <- pp[j] + step
    pm <- p0; pm[j] <- pm[j] - step
    (rxsSolve(h, pp)[, 1L] - rxsSolve(h, pm)[, 1L]) / (2 * step)
  }, numeric(attr(h, "ny")))

  expect_equal(got[, -1L], fd, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("nothing has happened before the dose is released", {
  ## An exact oracle that owes nothing to finite differencing: with a lag of 3
  ## the central compartment is still empty at t = 1 and 2, and insensitive to
  ## every parameter there.
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, c(1, 2, 4, 8))
  h <- rxsRegister(lagModel, events = ev, sens = sens, output = "center",
                   eventSens = "jump", atol = 1e-10, rtol = 1e-10)
  on.exit(rxsRelease(h))

  got <- rxsSolve(h, c(log(1.1), log(0.15), log(3)))
  expect_equal(got[1:2, 1L], c(0, 0), tolerance = 1e-12)
  expect_equal(as.numeric(got[1:2, -1L]), rep(0, 6), tolerance = 1e-12)

  ## And it is definitely non-zero once the dose has been released.
  expect_gt(abs(got[3, 1L]), 1e-6)
})

test_that("Path A is refused for models whose events need re-sorting", {
  ## A modeled lag moves the dose within the event order.  rxSolve re-sorts
  ## during setup; Path A only re-drives par_solve, so its cached ordering goes
  ## stale and results near the shifted dose are silently wrong.  Correctness
  ## wins over the speedup here.
  h <- mkHandle()
  on.exit(rxsRelease(h))

  expect_false(rxsFastAvailable(h))
  expect_match(attr(h, "fastDisabled"), "needSort")
})

test_that("results are right regardless of the direction the lag moves", {
  ## This is what caught the stale ordering: a lag crossing an observation time
  ## agreed going up and disagreed coming back down.
  h <- mkHandle()
  on.exit(rxsRelease(h))

  ref <- lapply(c(0.3, 0.5, 1.5, 3.0), function(lag) {
    rxsSolve(h, c(p0[1:2], log(lag)), slow = TRUE)
  })
  for (i in seq_along(ref)) {
    lag <- c(0.3, 0.5, 1.5, 3.0)[i]
    expect_identical(rxsSolve(h, c(p0[1:2], log(lag))), ref[[i]])
  }
  ## ... and again in the reverse order.
  for (i in rev(seq_along(ref))) {
    lag <- c(0.3, 0.5, 1.5, 3.0)[i]
    expect_identical(rxsSolve(h, c(p0[1:2], log(lag))), ref[[i]])
  }
})

test_that("interleaving an event-sensitivity handle with a plain one is safe", {
  ## rxode2 keeps the event-sensitivity shape in PROCESS GLOBALS, so two
  ## handles with different shapes sharing one cached solve is a real hazard.
  plainModel <- "
ka <- exp(lka)
cl <- exp(lcl)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - cl * center
"
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, c(1, 2, 4, 8, 12))

  hLag <- mkHandle()
  hPlain <- rxsRegister(plainModel, events = ev, sens = c("lka", "lcl"),
                        output = "center", atol = 1e-10, rtol = 1e-10)
  on.exit({ rxsRelease(hLag); rxsRelease(hPlain) })

  refLag <- rxsSolve(hLag, p0, slow = TRUE)
  refPlain <- rxsSolve(hPlain, p0[1:2], slow = TRUE)

  for (i in 1:3) {
    expect_identical(rxsSolve(hLag, p0), refLag)
    expect_identical(rxsSolve(hPlain, p0[1:2]), refPlain)
  }
})
