# Smoke stages 1-2 (issue #1, Spec 5): the model compiles and links against
# the injected external function (stage 1), and rstan's grad_log_prob agrees
# with numDeriv on the assembled target (stage 2) -- one comparison that
# catches a sign flip, a row/column-major transposition, a missing or doubled
# prior term, and a precomputed_gradients operand mismatch.

test_that("stan build probes pass", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  .i <- nlmixr2bayesBuildInfo(force = TRUE)
  expect_true(isTRUE(.i$ok))
  expect_true(isTRUE(.i$callShape))
  expect_true(isTRUE(.i$injectionPoint))
})

test_that("stage 1+2: the linked model compiles, runs, and its gradient is right", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("numDeriv")
  skip_on_cran()
  h <- stanLinkSetup(.linkMod, .linkData(), cores = 1L)
  on.exit(stanLinkFree(), add = TRUE)
  nlmixr2bayes:::.linkSetTheta(h$initPar)
  # keep nlmixr2est's internal Omega^-1 aligned with the Omega Stan uses so
  # the gradient assembly is well conditioned (value is Omega-free)
  .om <- matrix(0.1, 1, 1)
  nlmixr2bayes:::.linkSetOmegaInv(solve(.om))

  # ---- stage 1: compiles and links ----------------------------------------
  sm <- stanCompile(cache = TRUE)
  expect_s4_class(sm, "stanmodel")

  .data <- list(N = h$nid, P = h$neta, L = t(chol(.om)))
  sf <- rstan::sampling(sm, data = .data, chains = 1, iter = 2, warmup = 1,
                        refresh = 0, cores = 1, show_messages = FALSE)
  expect_s4_class(sf, "stanfit")

  # ---- stage 2: gradient + decomposition ----------------------------------
  set.seed(21)
  eta <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.25), h$nid, h$neta)
  up <- rstan::unconstrain_pars(sf, list(eta = eta))
  gA <- rstan::grad_log_prob(sf, up)
  attributes(gA) <- NULL
  gN <- numDeriv::grad(function(u) rstan::log_prob(sf, u), up,
                       method = "Richardson")
  expect_equal(gA, gN, tolerance = 1e-4)

  # the target decomposes exactly: differences between two eta points equal
  # the conditional + eta-prior differences (constants cancel, so this is
  # independent of rstan's proportionality conventions)
  set.seed(22)
  eta2 <- matrix(stats::rnorm(h$nid * h$neta, 0, 0.25), h$nid, h$neta)
  up2 <- rstan::unconstrain_pars(sf, list(eta = eta2))
  .lpDiff <- rstan::log_prob(sf, up, adjust_transform = FALSE) -
    rstan::log_prob(sf, up2, adjust_transform = FALSE)
  .cond <- function(e) sum(nlmixr2bayes:::.condBatch(e)$value)
  .prior <- function(e) {
    sum(stats::dnorm(as.numeric(e), 0, sqrt(0.1), log = TRUE))
  }
  .refDiff <- (.cond(eta) + .prior(eta)) - (.cond(eta2) + .prior(eta2))
  expect_equal(.lpDiff, .refDiff, tolerance = 1e-6)
})

test_that("stanCompile does not leak rstan's compiler environment (rxode2 C compiles survive)", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  # rstan::stan_model permanently sets PKG_CPPFLAGS (with a forced C++
  # -include), PKG_LIBS and USE_CXX17; unrestored, every later rxode2 C
  # model compile in the session dies with "compilation terminated".
  # stanCompile snapshots and restores them.
  .before <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS", "USE_CXX17"),
                        unset = NA_character_)
  .sm <- stanCompile(cache = TRUE) # cached: hits the fast path
  # force a REAL compile path once per shape: a trivially different program
  .code <- paste0(.stanPhase0Code(), "\n// leak-check variant")
  .sm2 <- stanCompile(.code, cache = FALSE)
  .after <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS", "USE_CXX17"),
                       unset = NA_character_)
  expect_identical(.before, .after)
})
