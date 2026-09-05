# Declared non-Gaussian random effect distributions under est="nuts"/"stan".
#
# Nothing in the generated Stan program is special cased for these.  By the
# time nlmixr2bayes sees such a model, nlmixr2est's pre-processing hook has
# expanded it into a latent standard normal with a FIXED unit variance --
# which the generator already supports -- the transform lives inside the
# linked nlmixr2_cond_all2() likelihood, and the copula correlation is an
# ordinary unbounded theta.  What is tested here is that the methods say so,
# and that a real run of one works end to end.

.etaDistMod <- function() {
  .f <- function() {
    ini({
      tv <- 3
      lclm  <- log(1)
      lclrv <- log(0.09)
      add.sd <- c(0, 0.5)
      eta.cl ~ 1
      prior(tv) ~ dnorm(3, 2)
      prior(lclm) ~ dnorm(0, 1)
      prior(lclrv) ~ dnorm(-2.4, 1)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- eta.cl
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .b <- body(.f)
  .i <- .b[[2]]
  .i[[2]][[length(.i[[2]]) + 1L]] <-
    str2lang("dist(eta.cl) ~ dgamma(shape=1/exp(lclrv), rate=1/(exp(lclrv)*exp(lclm)))")
  .b[[2]] <- .i
  body(.f) <- .b
  eval(parse(text=paste(deparse(.f), collapse="\n")))
}

test_that("every Stan method declares support", {
  for (.m in c("stan", "nuts", "advi", "pathfinder")) {
    expect_true(isTRUE(attr(utils::getS3method("nlmixr2Est", .m), "etaDist")),
                info=.m)
  }
})

test_that("a declared gamma random effect samples and recovers its parameters", {
  skip_on_cran()
  skip_if_not(requireNamespace("rstan", quietly=TRUE))
  .m <- .etaDistMod()
  set.seed(7)
  .d <- do.call(rbind, lapply(1:12, function(id) {
    data.frame(ID=id, TIME=c(0.5, 1, 2, 4, 8), DV=NA_real_)
  }))
  .sim <- rxode2::rxSolve(.m, .d, nSub=12, returnType="data.frame",
                          addDosing=FALSE)
  ## the simulated clearances ARE gamma draws, not a normal on a log scale
  expect_true(all(.sim$cl > 0))
  .d$DV <- .sim$cp + stats::rnorm(nrow(.d), 0, 0.5)
  .f <- suppressWarnings(
    nlmixr2est::nlmixr2(.m, .d, est="nuts",
                        control=stanControl(iter=400, warmup=200, chains=1,
                                            seed=42, print=0)))
  ## the gamma's own arguments are ordinary thetas, so they are estimated
  ## like any other; truth is log(1) = 0 and log(0.09) = -2.408
  expect_true(abs(.f$fixef[["lclm"]]) < 1)
  expect_true(abs(.f$fixef[["lclrv"]] + 2.408) < 1.5)
  ## the latent random effect and the declared one both come back
  expect_true("z.eta.cl" %in% names(.f))
  expect_true("eta.cl" %in% names(.f))
  expect_true(all(.f$eta.cl > 0))
})
