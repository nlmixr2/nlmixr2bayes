## Path A is refused when rxode2 reports a non-zero needSort, because a
## modelled lag/bioavailability/duration/rate moves a dose within the event
## order and Path A cannot re-sort.  All four act on DOSES, so a dose-free
## event table makes the flag a false alarm -- which is not hypothetical:
## rxode2 reports needSort = 3 for a DDE once it carries sensitivities.
##
## These tests pin both directions.  Getting the refusal wrong in one
## direction costs speed; in the other it returns silently wrong results near
## a shifted dose, so the "still refuses" cases matter more than the rest.
library(testthat)

ddeModel <- "
kin  <- exp(lkin)
kout <- exp(lkout)
tau  <- exp(ltau)
R(0) <- 1
d/dt(R) <- kin - kout * delay(R, tau)
"

lagModel <- "
ka   <- exp(lka)
cl   <- exp(lcl)
tlag <- exp(ltlag)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - cl * center
alag(depot)  <- tlag
"

test_that("a DDE reports needSort but still gets the fast path", {
  ## The flag really is set, otherwise this test proves nothing -- and it is
  ## the SENSITIVITY build that sets it, which is what the bridge solves.
  m <- rxode2::rxode2(ddeModel, calcSens = c("lkin", "lkout", "ltau"))
  expect_gt(as.integer(rxode2::rxModelVars(m)$needSort)[1], 0L)
  expect_equal(
    as.integer(rxode2::rxModelVars(rxode2::rxode2(ddeModel))$needSort)[1], 0L)

  h <- rxsRegister(ddeModel, events = rxode2::et(seq(0.25, 25, by = 0.25)),
                   sens = c("lkin", "lkout", "ltau"), output = "R",
                   atol = 1e-10, rtol = 1e-10, method = "dop853")
  on.exit(rxsRelease(h))

  expect_null(attr(h, "fastDisabled"))
  invisible(rxsSolve(h, log(c(0.7, 0.4, 1.5))))
  expect_true(rxsFastAvailable(h))
})

test_that("the DDE fast path agrees with the slow one exactly", {
  h <- rxsRegister(ddeModel, events = rxode2::et(seq(0.25, 25, by = 0.25)),
                   sens = c("lkin", "lkout", "ltau"), output = "R",
                   atol = 1e-10, rtol = 1e-10, method = "dop853")
  on.exit(rxsRelease(h))

  ## Sweep the way HMC would rather than checking one point.
  set.seed(11)
  worst <- 0
  for (i in 1:25) {
    p <- log(c(0.7, 0.4, 1.5)) + stats::runif(3, -0.5, 0.5)
    a <- rxsSolve(h, p)
    b <- rxsSolve(h, p, slow = TRUE)
    worst <- max(worst, max(abs(a - b)))
  }
  expect_identical(worst, 0)
  expect_true(rxsFastAvailable(h))
})

test_that("a modelled alag with real doses still refuses the fast path", {
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, c(1, 2, 4, 8, 12))
  h <- rxsRegister(lagModel, events = ev,
                   sens = c("lka", "lcl", "ltlag"), output = "center",
                   eventSens = "jump", atol = 1e-10, rtol = 1e-10)
  on.exit(rxsRelease(h))

  expect_match(attr(h, "fastDisabled"), "needSort")
  invisible(rxsSolve(h, c(log(1.1), log(0.15), log(0.5))))
  expect_false(rxsFastAvailable(h))
})

test_that("an ordinary dosing model is unaffected", {
  pk <- "
ka <- exp(lka)
cl <- exp(lcl)
v  <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
"
  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by = 1.5))
  h <- rxsRegister(pk, events = ev, sens = c("lka", "lcl", "lv"),
                   output = "center", atol = 1e-10, rtol = 1e-10)
  on.exit(rxsRelease(h))

  expect_null(attr(h, "fastDisabled"))
  invisible(rxsSolve(h, c(log(1.1), log(4), log(30))))
  expect_true(rxsFastAvailable(h))
})

test_that(".rxsHasDoses fails safe when it cannot read the events", {
  expect_true(nlmixr2bayes:::.rxsHasDoses("not an event table"))
  expect_true(nlmixr2bayes:::.rxsHasDoses(data.frame(time = 1:3)))
  expect_false(nlmixr2bayes:::.rxsHasDoses(data.frame(time = 1:3, evid = c(0, 0, 0))))
  expect_true(nlmixr2bayes:::.rxsHasDoses(data.frame(time = 1:3, evid = c(1, 0, 0))))

  ## And on the event tables users actually pass.
  ev <- rxode2::et(amt = 100, cmt = "depot")
  expect_true(nlmixr2bayes:::.rxsHasDoses(rxode2::et(ev, seq(0.5, 24, by = 1.5))))
  expect_false(nlmixr2bayes:::.rxsHasDoses(rxode2::et(seq(0.5, 24, by = 1.5))))
})
