test_that("stanControl builds, validates, and pins its invariants", {
  .c <- stanControl()
  expect_s3_class(.c, "stanControl")
  expect_equal(.c$chains, 4L)
  expect_equal(.c$likelihood, "focei")
  expect_true(.c$run)
  # atol/rtol pinned tight (sigdig=8) so the density is smooth at sampler scale
  expect_equal(.c$rxControl$atol, 1e-8)
  expect_error(stanControl(chains = 0), "chains")
  # per-run cores like every other nlmixr2est control (default getRxThreads)
  expect_equal(stanControl(cores = 4)$cores, 4L)
  expect_error(stanControl(bogusArgument = 4), "unused argument")
  expect_error(stanControl(likelihood = "focep"), "'arg'")
  # sampler tuning matches rstan's own defaults so est="stan" and a
  # hand-written rstan run behave identically out of the box
  expect_equal(stanControl()$adapt_delta, 0.8)
  expect_equal(stanControl()$max_treedepth, 10L)
  # default solver: dense AutoSwitch dop853+ros4 (the dense Dormand-Prince
  # twin of Stan's ode_rk45, with automatic stiff fallback)
  expect_equal(stanControl()$rxControl$method,
               rxode2::rxControl(method = "dop853+ros4")$method)
  expect_true(stanControl()$rxControl$dense)
  if (identical(.Platform$OS.type, "unix")) {
    # forked parallel chains by default, inheriting the core budget
    # capped at the chain count; subject threads then default to 1
    # (OpenMP-after-fork hazard)
    .c2 <- stanControl(chains = 2L)
    expect_equal(.c2$chainCores,
                 min(2L, as.integer(rxode2::getRxThreads())))
    if (.c2$chainCores > 1L) expect_equal(.c2$cores, 1L)
    # sequential opt-out restores subject-parallel evaluation
    .c1 <- stanControl(chains = 2L, chainCores = 1L)
    expect_equal(.c1$chainCores, 1L)
    expect_equal(.c1$cores, as.integer(rxode2::getRxThreads()))
    # explicit both: honored with the hazard warning
    expect_warning(stanControl(chainCores = 2L, cores = 2L), "fork")
  } else {
    # Windows chain parallelism uses PSOCK workers and keeps the
    # user-requested per-worker subject-thread count.
    .c2 <- suppressWarnings(stanControl(chains = 2L, chainCores = 2L,
                                        cores = 2L))
    expect_equal(.c2$chainCores, 2L)
    expect_equal(.c2$cores, 2L)
  }
})

test_that("getValidNlmixrCtl.stan coerces and re-validates", {
  .v <- getValidNlmixrCtl.stan(list(NULL))
  expect_s3_class(.v, "stanControl")
  .v <- getValidNlmixrCtl.stan(list(list(chains = 2L)))
  expect_equal(.v$chains, 2L)
  .v <- getValidNlmixrCtl.stan(list(stanControl(iter = 500L)))
  expect_equal(.v$iter, 500L)
})

test_that("the control class name makes est= inferrable", {
  # .nlmixr2inferEst regexes ^(.*?)Control$ off class(control)[1]
  expect_equal(sub("Control$", "", class(stanControl())[1]), "stan")
})
