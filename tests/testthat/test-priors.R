# Prior emission (issue #1, Spec 4).  Pure R -- no Stan toolchain needed.

.priorMod <- function(priors) {
  # a small ui with configurable prior lines
  .body <- paste0("function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- c(0, 0.7)
    eta.ka ~ 0.6
    ", paste(priors, collapse = "\n    "), "
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}")
  suppressMessages(rxode2::rxode2(eval(str2lang(.body))))
}

test_that("half-Cauchy: bound from the parameter, no T[,] for literal args", {
  .ui <- .priorMod(c("prior(tka) ~ dnorm(0, 10)",
                     "prior(add.sd) ~ dcauchy(0, 5)"))
  .p <- stanPriors(.ui)
  expect_equal(nrow(.p$pop), 2L)
  .tka <- .p$pop[.p$pop$name == "tka", ]
  expect_equal(.tka$statement, "tka ~ normal(0, 10);")
  expect_equal(.tka$constraint, "")
  expect_false(.tka$truncated)
  .sd <- .p$pop[.p$pop$name == "add.sd", ]
  # THE half-Cauchy idiom: <lower=0> declaration, plain statement, no T[,]
  # (the dropped normalizer is a constant because the arguments are literals)
  expect_equal(.sd$par, "add_sd")
  expect_equal(.sd$constraint, "<lower=0>")
  expect_equal(.sd$statement, "add_sd ~ cauchy(0, 5);")
  expect_false(.sd$truncated)
})

test_that("support promotion tightens the declaration (lotri's gap)", {
  # dlnorm has support 'positive' but tka is unbounded; lotri only validates
  # against FINITE bounds, so without promotion the target is -Inf at init
  .ui <- .priorMod("prior(tka) ~ dlnorm(0, 1)")
  .p <- stanPriors(.ui)
  expect_equal(.p$pop$lower, 0)
  expect_equal(.p$pop$constraint, "<lower=0>")
  expect_equal(.p$pop$statement, "tka ~ lognormal(0, 1);")
  # dbeta: unit support promotes both bounds
  .ui <- .priorMod("prior(tka) ~ dbeta(2, 2)")
  .p <- stanPriors(.ui)
  expect_equal(.p$pop$lower, 0)
  expect_equal(.p$pop$upper, 1)
  expect_equal(.p$pop$constraint, "<lower=0,upper=1>")
})

test_that("parameter-dependent arguments force T[,] (rule 2c)", {
  # the scale of add.sd's prior references another estimated parameter: the
  # truncation normalizer is parameter-dependent and cannot be dropped
  # (rxode2 requires every ini() parameter in the model, so the
  # hyperparameter here is an ordinary model theta)
  .ui <- .priorMod("prior(add.sd) ~ dcauchy(0, tcl)")
  .p <- stanPriors(.ui)
  .sd <- .p$pop[.p$pop$name == "add.sd", ]
  expect_true(.sd$truncated)
  expect_equal(.sd$statement, "add_sd ~ cauchy(0, tcl) T[0, ];")
})

test_that("discrete priors are refused with the parameter named (rule 4)", {
  .ui <- .priorMod("prior(tka) ~ dpois(3)")
  expect_error(stanPriors(.ui), "dpois.*tka.*real-valued")
})

test_that("multivariate priors dedupe to one statement (rule 3)", {
  .ui <- .priorMod("prior(tcl, tv) ~ multiNormal(c(1, 3.45), lotri(tcl + tv ~ c(1, 0.01, 1)))")
  .p <- stanPriors(.ui)
  # duplicated byte-identically on tcl and tv rows; must emit exactly once
  expect_equal(sum(.p$pop$kind == "multivariate"), 1L)
  .mv <- .p$pop[.p$pop$kind == "multivariate", ]
  expect_equal(.mv$members, "tcl,tv")
  expect_match(.mv$statement,
               "target \\+= multi_normal_lpdf\\(to_vector\\(\\{tcl, tv\\}\\) \\| \\[1, 3.45\\]', \\[\\[1, 0.01\\], \\[0.01, 1\\]\\]\\);")
})

test_that("omega-block priors are classified for the generator, not emitted", {
  .ui <- .priorMod("prior(eta.ka) ~ invWishart(4)")
  .p <- stanPriors(.ui)
  expect_equal(nrow(.p$pop), 0L)
  expect_equal(nrow(.p$omega), 1L)
  expect_equal(.p$omega$stanName, "inv_wishart")
  expect_equal(.p$omega$kind, "matrix")
  expect_equal(.p$omega$neta1, 1)
})

test_that("every non-discrete univariate distribution in the catalogue parses", {
  .tbl <- lotri::lotriPriorDists()
  .uni <- .tbl[.tbl$kind == "univariate", ]
  for (.i in seq_len(nrow(.uni))) {
    .r <- .uni[.i, ]
    .nArgs <- .r$nReq
    .call <- if (.nArgs == 0L) paste0(.r$name, "()") else {
      paste0(.r$name, "(", paste(rep("1", .nArgs), collapse = ", "), ")")
    }
    .lk <- .stanPriorLookup(.call)
    expect_equal(.lk$stanName, .r$stanName, info = .r$name)
    .u <- .stanPriorUnivariate("tka", .lk, -Inf, Inf)
    expect_match(.u$statement, paste0("~ ", .r$stanName, "\\("), info = .r$name)
  }
})

test_that("name mangling maps nlmixr2 names onto Stan identifiers", {
  expect_equal(.stanParName(c("add.sd", "eta.ka", "tka")),
               c("add_sd", "eta_ka", "tka"))
  expect_equal(.stanParName("2fast"), "p_2fast")
})
