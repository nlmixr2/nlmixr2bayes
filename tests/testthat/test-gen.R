# The generator (Spec 9), exercised through the run=FALSE escape -- no Stan
# toolchain needed.

test_that("run=FALSE returns the generated program without touching rstan", {
  skip_on_cran()
  .out <- suppressMessages(
                           nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                                               control = stanControl(run = FALSE)))
  expect_s3_class(.out, "nlmixr2bayesCode")
  .code <- .out$code
  # external declaration + call
  expect_match(.code, "vector nlmixr2_cond_all2\\(matrix etaMat, vector theta\\);")
  # the model block computes llCond once, adds sum(llCond), and ticks the
  # nlmixr2est-style iteration print (a no-op unless armed)
  expect_match(.code, "vector\\[N\\] llCond = nlmixr2_cond_all2\\(eta, theta\\);")
  expect_match(.code, "target \\+= sumLL;")
  expect_match(.code, "nlmixr2_iter_tick\\(parDisp, -2 \\* sumLL\\);")
  # free thetas declared with their (promoted) bounds
  expect_match(.code, "real tcl;")
  expect_match(.code, "real<lower=0> add_sd;")
  # priors emitted from ini({})
  expect_match(.code, "tcl ~ normal\\(1, 2\\);")
  expect_match(.code, "add_sd ~ cauchy\\(0, 2.5\\);")
  # non-centred eta with the default half-Cauchy SD (announced); the block's
  # Stan parameters are named from its member eta(s), not an opaque "_b1"
  expect_match(.code, "real<lower=0> sd_eta_cl;")
  expect_match(.code, "to_vector\\(z_eta_cl\\) ~ std_normal\\(\\);")
  expect_match(.code, "eta\\[, 1:1\\] = z_eta_cl \\* L_eta_cl';")
  # theta vector assembled in iniDf order
  expect_match(.code, "theta\\[1\\] = tcl;")
  expect_match(.code, "theta\\[3\\] = add_sd;")
  # generated quantities carry Omega and the per-subject conditional loglik
  expect_match(.code, "omegaOut")
  expect_match(.code, "logLikSubj = nlmixr2_cond_all2\\(eta, theta\\);")
  expect_equal(.out$data$N, 4L)
})

test_that("an invWishart omega prior declares cov_matrix with that prior", {
  skip_on_cran()
  .f <- function() {
    ini({
      tcl <- 1
      tv <- 3
      add.sd <- c(0, 0.5)
      eta.cl ~ 0.1
      prior(tcl) ~ dnorm(1, 2)
      prior(add.sd) ~ dcauchy(0, 2.5)
      prior(eta.cl) ~ invWishart(4)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  .out <- suppressMessages(
                           nlmixr2est::nlmixr2(.f, .linkData(), est = "stan",
                                               control = stanControl(run = FALSE)))
  expect_match(.out$code, "cov_matrix\\[1\\] omega_eta_cl;")
  expect_match(.out$code, "omega_eta_cl ~ inv_wishart\\(4, \\[\\[0.1\\]\\]\\);")
  expect_match(.out$code, "L_eta_cl = cholesky_decompose\\(omega_eta_cl\\);")
})

test_that("a model without priors errors with the exact lines to add (D10)", {
  skip_on_cran()
  expect_error(
               suppressMessages(
                                nlmixr2est::nlmixr2(.linkMod, .linkData(), est = "stan",
                                                    control = stanControl(run = FALSE))),
               "prior\\(tcl\\) ~ dnorm")
})

test_that("the generated program parses with stanc", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  .out <- suppressMessages(
                           nlmixr2est::nlmixr2(.estMod, .linkData(), est = "stan",
                                               control = stanControl(run = FALSE)))
  .sc <- rstan::stanc(model_code = .out$code, allow_undefined = TRUE,
                      model_name = "gen_parse")
  expect_true(is.list(.sc))
})

test_that("a multivariate prior with a fixed member is refused, not dropped", {
  skip_on_cran()
  .f <- function() {
    ini({
      tcl <- 1
      tv <- fix(3)
      add.sd <- c(0, 0.5)
      eta.cl ~ 0.1
      prior(tcl, tv) ~ multiNormal(c(1, 3), lotri(tcl + tv ~ c(1, 0.01, 1)))
      prior(add.sd) ~ dcauchy(0, 2.5)
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- 100 / v * exp(-cl / v * time)
      cp ~ add(add.sd)
    })
  }
  expect_error(
               suppressMessages(
                                nlmixr2est::nlmixr2(.f, .linkData(), est = "stan",
                                                    control = stanControl(run = FALSE))),
               "fixed")
})
