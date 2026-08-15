# G7 (truncation-constant invariance) and G9 (reproducibility).
#
# G7: prior-emission branch (b) declares the constraint and drops T[L,U]
# because the normalizer is constant when every distribution argument is a
# literal.  The gate compiles the shipped program AND a variant with the
# explicit T[0,] restored, and shows their log densities differ by exactly
# a constant across 200 random points -- so the posterior is identical and
# the dropped term really was constant.
#
# G9: bitwise reproducibility for a fixed seed within a session, and the
# compile cache: a second stanCompile of the same code is a fast cache hit,
# and a control knob that changes the generated code (lkjEta) changes the
# code text, so the content-hash key re-keys by construction.

.g7Mod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- c(0, 0.5)
    eta.cl ~ 0.1
    prior(tcl) ~ dnorm(1, 2)
    prior(add.sd) ~ dcauchy(0, 2.5)
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

test_that("half-Cauchy branch (b): explicit T[0,] shifts the target by exactly a constant (G7)", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.g7Mod, .linkData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  # golden: the constraint is declared and no truncation is emitted
  expect_true(any(grepl("real<lower=0> add_sd;", .lines, fixed = TRUE)))
  expect_false(any(grepl("T[", .lines, fixed = TRUE)))
  .w <- grep("add_sd ~ cauchy(0, 2.5);", .lines, fixed = TRUE)
  expect_length(.w, 1L)
  .linesT <- .lines
  .linesT[.w] <- "  add_sd ~ cauchy(0, 2.5) T[0,];"
  .smA <- stanCompile(paste(.lines, collapse = "\n"))
  .smB <- stanCompile(paste(.linesT, collapse = "\n"))
  h <- stanLinkSetup(.g7Mod, .linkData(), thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .sfA <- rstan::sampling(.smA, data = .code$data, chains = 1, iter = 2,
                          warmup = 1, refresh = 0, cores = 1,
                          show_messages = FALSE)
  .sfB <- rstan::sampling(.smB, data = .code$data, chains = 1, iter = 2,
                          warmup = 1, refresh = 0, cores = 1,
                          show_messages = FALSE)
  .n <- rstan::get_num_upars(.sfA)
  expect_equal(rstan::get_num_upars(.sfB), .n)
  set.seed(31)
  .diff <- vapply(1:200, function(i) {
    .u <- stats::rnorm(.n, 0, 0.5)
    rstan::log_prob(.sfA, .u, adjust_transform = FALSE) -
      rstan::log_prob(.sfB, .u, adjust_transform = FALSE)
  }, numeric(1))
  expect_lt(stats::sd(.diff), 1e-10)
  # and the constant is the half-Cauchy normalizer: T[0,] subtracts
  # log P(X > 0) = log(1/2), so shipped - truncated = -log(2)
  expect_equal(mean(.diff), -log(2), tolerance = 1e-8)
})

test_that("fixed seed: bitwise-identical draws within a session (G9)", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.g7Mod, .linkData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .sm <- stanCompile(.code$code)
  h <- stanLinkSetup(.g7Mod, .linkData(), thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`, 1L)
  .run <- function() {
    rstan::extract(
      rstan::sampling(.sm, data = .code$data, chains = 1L, iter = 300L,
                      warmup = 150L, seed = 7L, refresh = 0, cores = 1,
                      show_messages = FALSE),
      permuted = FALSE)
  }
  expect_identical(.run(), .run())
})

test_that("compile cache: content-keyed hit; code-changing knobs re-key (G9d)", {
  skip_if_not_installed("rstan")
  skip_on_cran()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.g7Mod, .linkData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .dir <- file.path(tempdir(), "nlmixr2stan-cache-test")
  unlink(.dir, recursive = TRUE)
  .t1 <- system.time(stanCompile(.code$code, cacheDir = .dir))[["elapsed"]]
  expect_length(list.files(.dir, pattern = "^stanmodel-.*rds$"), 1L)
  .t2 <- system.time(stanCompile(.code$code, cacheDir = .dir))[["elapsed"]]
  expect_lt(.t2, max(2, .t1 / 5)) # warm hit, no recompile
  expect_length(list.files(.dir, pattern = "^stanmodel-.*rds$"), 1L)
  unlink(.dir, recursive = TRUE)
  # a control knob that changes the generated program changes the code text,
  # and the cache key is a content hash of that text -- re-keyed by
  # construction (the risk G9d guards: a knob the key ignores)
  .mod2 <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      eta.cl + eta.v ~ c(0.1,
                         0.02, 0.1)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .cA <- suppressMessages(
    nlmixr2est::nlmixr2(.mod2, .linkData(), est = "stan",
                        control = stanControl(run = FALSE)))
  .cB <- suppressMessages(
    nlmixr2est::nlmixr2(.mod2, .linkData(), est = "stan",
                        control = stanControl(run = FALSE, lkjEta = 4)))
  expect_false(identical(.cA$code, .cB$code))
  expect_true(any(grepl("lkj_corr_cholesky(4)",
                        strsplit(.cB$code, "\n")[[1]], fixed = TRUE)))
})
