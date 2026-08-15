test_that("stanControl builds, validates, and pins its invariants", {
  .c <- stanControl()
  expect_s3_class(.c, "stanControl")
  expect_equal(.c$chains, 4L)
  expect_equal(.c$likelihood, "focei")
  expect_true(.c$run)
  # atol/rtol pinned tight (sigdig=8) so the density is smooth at sampler scale
  expect_equal(.c$rxControl$atol, 1e-8)
  expect_error(stanControl(chains = 0), "chains")
  expect_error(stanControl(cores = 4), "unused argument")
  expect_error(stanControl(likelihood = "focep"), "'arg'")
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
