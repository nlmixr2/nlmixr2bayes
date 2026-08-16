# First-class IOV: nlmixr2est's preprocessing hook (.uiApplyIov, enabled by
# attr(nlmixr2Est.stan, "iov")) expands `iov.x ~ v | OCC` into per-occasion
# FIXED unit-variance etas scaled by an estimated sd theta in the model
# code.  The generator renders a fixed eta as a CONSTANT Cholesky factor
# (z still sampled -- a random effect with known variance), which also
# gives plain fixed-variance etas for free.  Fixed etas with modeled
# correlations are refused.

.iovMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- c(0, 0.5)
    prior(tcl) ~ dnorm(1, 2)
    prior(add.sd) ~ dcauchy(0, 2.5)
    eta.cl ~ 0.1
    iov.cl ~ 0.02 | OCC
  })
  model({
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.iovData <- function() {
  do.call(rbind, lapply(1:4, function(id) {
    rbind(data.frame(ID = id, TIME = c(1, 4, 8), DV = c(2, 1.4, 0.8),
                     AMT = 0, EVID = 0, OCC = 1),
          data.frame(ID = id, TIME = c(25, 28, 32), DV = c(2.1, 1.5, 0.9),
                     AMT = 0, EVID = 0, OCC = 2))
  }))
}

test_that("IOV expands to fixed unit-variance etas + an sd theta", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iovMod, .iovData(), est = "stan",
                        control = stanControl(run = FALSE)))
  expect_s3_class(.code, "nlmixr2stanCode")
  .lines <- strsplit(.code$code, "\n")[[1]]
  # the IOV magnitude is an ordinary theta ...
  expect_true(any(grepl("theta[4] = iov_cl;", .lines, fixed = TRUE)))
  # ... and the occasion etas are fixed unit-variance blocks: constant L,
  # no omega parameter, z still sampled
  expect_true(any(grepl("matrix[1,1] L_b2 = [[1]];", .lines, fixed = TRUE)))
  expect_true(any(grepl("matrix[1,1] L_b3 = [[1]];", .lines, fixed = TRUE)))
  expect_true(any(grepl("to_vector(z_b2) ~ std_normal();", .lines,
                        fixed = TRUE)))
  expect_false(any(grepl("omega_b2|sd_b2", .lines)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .code$code,
                               allow_undefined = TRUE))
  }
})

test_that("the assembled IOV target's gradient FD-agrees end to end", {
  skip_on_cran()
  skip_if_not(nlmixr2stan:::.stanHasIovSens(),
              "IOV gated off: upstream sens column wrong (nlmixr2est#952)")
  skip_if_not_installed("rstan")
  skip_if_not_installed("numDeriv")
  .d <- .iovData()
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.iovMod, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .sm <- stanCompile(.code$code)
  # link the REWRITTEN ui (the hook's expansion) exactly as est="stan" does
  h <- stanLinkSetup(.code$ui, .d, thetaSens = TRUE, cores = 1L)
  on.exit({
    .Call(nlmixr2stan:::`_nlmixr2stan_clearThetaBase`)
    stanLinkFree()
  }, add = TRUE)
  .Call(nlmixr2stan:::`_nlmixr2stan_setThetaBase`, as.double(h$initPar))
  .Call(nlmixr2stan:::`_nlmixr2stan_setMuRef`,
        as.integer(.code$map$muRefIdx))
  .sf <- rstan::sampling(.sm, data = .code$data, chains = 1, iter = 2,
                         warmup = 1, refresh = 0, cores = 1,
                         show_messages = FALSE,
                         init = list(list(tcl = 1, tv = 3, add_sd = 0.5,
                                          iov_cl = sqrt(0.02),
                                          sd_b1 = sqrt(0.1),
                                          z_b1 = matrix(0, 4, 1),
                                          z_b2 = matrix(0, 4, 1),
                                          z_b3 = matrix(0, 4, 1))))
  set.seed(9)
  .pt <- list(tcl = 1.05, tv = 2.95, add_sd = 0.55, iov_cl = 0.2,
              sd_b1 = 0.3,
              z_b1 = matrix(stats::rnorm(4, 0, 0.4), 4, 1),
              z_b2 = matrix(stats::rnorm(4, 0, 0.4), 4, 1),
              z_b3 = matrix(stats::rnorm(4, 0, 0.4), 4, 1))
  .up <- rstan::unconstrain_pars(.sf, .pt)
  .gA <- rstan::grad_log_prob(.sf, .up)
  attributes(.gA) <- NULL
  .gN <- numDeriv::grad(function(u) rstan::log_prob(.sf, u), .up,
                        method = "Richardson")
  # every block: thetas (incl. the IOV sd through the forward
  # sensitivities), the id-level eta, and BOTH occasion-eta blocks
  expect_equal(.gA, .gN, tolerance = 1e-4)
})

test_that("plain fixed-variance etas work; fixed correlated blocks refuse", {
  skip_on_cran()
  .fixEta <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      eta.cl ~ fix(0.04)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .d <- .iovData()[, c("ID", "TIME", "DV", "AMT", "EVID")]
  .code <- suppressMessages(
    nlmixr2est::nlmixr2(.fixEta, .d, est = "stan",
                        control = stanControl(run = FALSE)))
  .lines <- strsplit(.code$code, "\n")[[1]]
  expect_true(any(grepl("matrix[1,1] L_b1 = [[0.2]];", .lines,
                        fixed = TRUE)))
  expect_false(any(grepl("sd_b1", .lines)))
  # the centred parameterization composes with a fixed block: eta sampled
  # directly against the CONSTANT L
  .codeC <- suppressMessages(
    nlmixr2est::nlmixr2(.fixEta, .d, est = "stan",
                        control = stanControl(run = FALSE,
                                              etaParam = "centered")))
  .linesC <- strsplit(.codeC$code, "\n")[[1]]
  expect_true(any(grepl("matrix[1,1] L_b1 = [[0.2]];", .linesC,
                        fixed = TRUE)))
  expect_true(any(grepl("etaP_b1[i] ~ multi_normal_cholesky", .linesC,
                        fixed = TRUE)))
  if (requireNamespace("rstan", quietly = TRUE)) {
    expect_silent(rstan::stanc(model_code = .codeC$code,
                               allow_undefined = TRUE))
  }
  # fixed etas with modeled correlations: out of scope, must error
  .fixBlk <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      eta.cl + eta.v ~ fix(0.04,
                           0.01, 0.04)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.fixBlk, .d, est = "stan",
                          control = stanControl(run = FALSE))),
    "correlations")
})

test_that("IOV models are refused while the upstream sens column is wrong (#952)", {
  skip_on_cran()
  skip_if(nlmixr2stan:::.stanHasIovSens())
  # with the iov attribute gated off, the | occ etas survive to the
  # id-only assertion and are refused rather than sampled with a
  # ~14% wrong gradient on the IOV magnitude theta
  expect_error(
    suppressMessages(
      nlmixr2est::nlmixr2(.iovMod, .iovData(), est = "stan",
                          control = stanControl(run = FALSE))))
})
