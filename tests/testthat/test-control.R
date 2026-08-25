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

test_that("the sugar controls build, force their algorithm, and are inferrable", {
  # nutsControl()/adviControl()/pathfinderControl() are stanControl() with
  # the distinguishing algorithm forced and re-classed -- the
  # foceControl()-from-foceiControl() pattern
  .n <- nutsControl()
  expect_s3_class(.n, "nutsControl")
  expect_equal(.n$algorithm, "NUTS")
  .a <- adviControl()
  expect_s3_class(.a, "adviControl")
  expect_equal(.a$algorithm, "meanfield")
  expect_equal(adviControl(algorithm = "fullrank")$algorithm, "fullrank")
  .p <- pathfinderControl()
  expect_s3_class(.p, "pathfinderControl")
  expect_equal(.p$algorithm, "pathfinder")
  # class(control)[1] is what .nlmixr2inferEst regexes, so
  # nlmixr2(model, data, nutsControl()) needs no est=
  expect_equal(sub("Control$", "", class(.n)[1]), "nuts")
  expect_equal(sub("Control$", "", class(.a)[1]), "advi")
  expect_equal(sub("Control$", "", class(.p)[1]), "pathfinder")
  # everything else is stanControl(): same contents, same validation
  expect_equal(names(.n), names(stanControl()))
  expect_equal(nutsControl(iter = 500L, warmup = 100L)$iter, 500L)
  expect_equal(pathfinderControl(pathfinderPaths = 8L)$pathfinderPaths, 8L)
  expect_error(nutsControl(bogusArgument = 4), "unused argument")
  expect_error(adviControl(chains = 0), "chains")
})

test_that("sugar-control algorithm contradictions error, not silently override", {
  expect_error(nutsControl(algorithm = "meanfield"), "nuts")
  expect_error(nutsControl(algorithm = "pathfinder"), "nuts")
  expect_error(adviControl(algorithm = "pathfinder"), "advi")
  expect_error(pathfinderControl(algorithm = "fullrank"), "pathfinder")
  # "NUTS" is only ever stanControl()'s default leaking through: est="advi"
  # picks the ADVI default and est="pathfinder" picks Pathfinder
  expect_equal(do.call(adviControl, stanControl())$algorithm, "meanfield")
  expect_equal(do.call(pathfinderControl, stanControl())$algorithm,
               "pathfinder")
  expect_error(do.call(nutsControl, stanControl(algorithm = "meanfield")),
               "nuts")
})

test_that("getValidNlmixrCtl.<sugar> coerces stanControl and re-validates", {
  expect_s3_class(getValidNlmixrCtl.nuts(list(NULL)), "nutsControl")
  expect_s3_class(getValidNlmixrCtl.advi(list(NULL)), "adviControl")
  expect_s3_class(getValidNlmixrCtl.pathfinder(list(NULL)),
                  "pathfinderControl")
  expect_equal(getValidNlmixrCtl.nuts(list(list(chains = 2L)))$chains, 2L)
  # a plain stanControl() still works for est="nuts"/"advi"/"pathfinder"
  .v <- getValidNlmixrCtl.nuts(list(stanControl(iter = 500L)))
  expect_s3_class(.v, "nutsControl")
  expect_equal(.v$iter, 500L)
  expect_equal(getValidNlmixrCtl.advi(list(stanControl()))$algorithm,
               "meanfield")
  # a contradicting sugar control is an error, not a silent conversion
  expect_error(suppressMessages(getValidNlmixrCtl.nuts(list(adviControl()))),
               "nuts")
})

test_that("the sugar controls deparse as themselves", {
  expect_equal(rxode2::rxUiDeparse(nutsControl(), "ctl"),
               str2lang("ctl <- nutsControl()"))
  expect_equal(rxode2::rxUiDeparse(nutsControl(iter = 500L, warmup = 100L),
                                   "ctl"),
               str2lang("ctl <- nutsControl(iter = 500L, warmup = 100L)"))
  expect_equal(rxode2::rxUiDeparse(adviControl(), "ctl"),
               str2lang("ctl <- adviControl()"))
  expect_equal(rxode2::rxUiDeparse(adviControl(algorithm = "fullrank"), "ctl"),
               str2lang("ctl <- adviControl(algorithm = \"fullrank\")"))
  expect_equal(rxode2::rxUiDeparse(pathfinderControl(), "ctl"),
               str2lang("ctl <- pathfinderControl()"))
})
