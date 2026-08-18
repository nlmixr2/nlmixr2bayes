## rstan::stan_model() leaves PKG_CPPFLAGS/PKG_LIBS pointing at the Stan headers
## for the rest of the session, which makes every later rxode2 model build fail:
## the generated C file picks up Stan's C++ Eigen header and dies on <cmath>.
## rxsStanModel() has to put them back.
library(testthat)

test_that("compiling a Stan model leaves rxode2 able to compile", {
  skip_if_not_installed("rstan")
  if (!nzchar(Sys.getenv("RXSTAN_STAN_TESTS"))) {
    skip("set RXSTAN_STAN_TESTS=1 to run the Stan compilation tests")
  }

  before <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS"), unset = NA)

  ev <- rxode2::et(amt = 100, cmt = "depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by = 4))
  h <- rxsRegister("
ka  <- exp(lka)
cl  <- exp(lcl)
v   <- exp(lv)
d/dt(depot)  <- -ka * depot
d/dt(center) <-  ka * depot - (cl / v) * center
", events = ev, sens = c("lka", "lcl", "lv"), output = "center")
  on.exit(rxsRelease(h))

  rxsStanModel(system.file("stan", "pk_1cmt_oral.stan", package = "nlmixr2bayes"),
               modelName = "envcheck")

  expect_equal(Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS"), unset = NA), before)

  ## The real symptom: a fresh rxode2 model must still build.
  m <- rxode2::rxode2("d/dt(a) <- -kk * a")
  expect_equal(rxode2::rxState(m), "a")
})
