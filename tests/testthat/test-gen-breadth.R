# G10 breadth: every non-discrete univariate distribution in the catalogue,
# carried end-to-end -- ini({}) prior -> generator -> a COMPLETE program that
# stanc accepts.  test-priors.R already checks the lookup/statement layer;
# what that cannot catch is an argument-count or type mismatch against
# Stan's real signature, which only stanc sees.  Distributions are grouped
# by support class into four programs (one ui build each).

.breadthSpec <- list(
  real = c(dn = "dnorm(0, 1)", sn = "stdNormal()",
           emn = "expModNormal(0, 1, 1)", skn = "skewNormal(0, 1, 1)",
           st = "studentT(3, 0, 1)", ca = "dcauchy(0, 1)",
           de = "doubleExponential(0, 1)", lo = "dlogis(0, 1)",
           gu = "gumbel(0, 1)", sde = "skewDoubleExponential(0, 1, 0.5)"),
  positive = c(ln = "dlnorm(0, 1)", ch = "dchisq(3)",
               ic = "invChiSquare(3)", sic = "scaledInvChiSquare(3, 1)",
               ex = "dexp(1)", ga = "dgamma(2, 1)", ig = "invGamma(2, 1)",
               we = "dweibull(2, 1)", fr = "frechet(2, 1)",
               ra = "rayleigh(1)", wi = "wiener(1, 0.1, 0.5, 0.5)",
               pa = "pareto(0.1, 2)", p2 = "paretoType2(0, 1, 2)"),
  unit = c(be = "dbeta(2, 2)", bp = "betaProportion(0.5, 2)",
           un = "dunif(0, 1)"),
  circular = c(vm = "vonMises(0, 1)"))

.breadthBounds <- c(real = "0.5", positive = "c(0, 0.5)",
                    unit = "c(0, 0.5, 1)", circular = "0.5")

.breadthModel <- function(supportClass) {
  .pri <- .breadthSpec[[supportClass]]
  .nm <- paste0("px", seq_along(.pri))
  .ini <- c("tcl <- 1",
            paste0(.nm, " <- ", .breadthBounds[[supportClass]]),
            "add.sd <- c(0, 0.5)",
            "eta.cl ~ 0.1",
            "prior(tcl) ~ dnorm(1, 2)",
            paste0("prior(", .nm, ") ~ ", .pri),
            "prior(add.sd) ~ dcauchy(0, 2.5)")
  .mdl <- c("cl <- exp(tcl + eta.cl)",
            paste0("adj <- 0.001 * (", paste(.nm, collapse = " + "), ")"),
            "cp <- exp(adj) * 100 * exp(-cl * time)",
            "cp ~ add(add.sd)")
  eval(parse(text = paste0(
    "function() { ini({", paste(.ini, collapse = "\n"), "})\n",
    "model({", paste(.mdl, collapse = "\n"), "}) }")))
}

test_that("every univariate distribution survives ini -> generator -> stanc (G10)", {
  skip_on_cran()
  .tbl <- lotri::lotriPriorDists()
  for (.cls in names(.breadthSpec)) {
    .code <- suppressMessages(
      nlmixr2est::nlmixr2(.breadthModel(.cls), .linkData(), est = "stan",
                          control = stanControl(run = FALSE)))
    # every distribution's Stan name appears in its program
    for (.p in .breadthSpec[[.cls]]) {
      .rName <- sub("\\(.*", "", .p)
      .stanName <- .tbl$stanName[.tbl$name == .rName |
                                   .tbl$camelName == .rName][1]
      expect_true(any(grepl(paste0("~ ", .stanName, "("),
                            strsplit(.code$code, "\n")[[1]], fixed = TRUE)),
                  label = paste0(.cls, "/", .rName, " emits ", .stanName))
    }
    if (requireNamespace("rstan", quietly = TRUE)) {
      expect_silent(rstan::stanc(model_code = .code$code,
                                 allow_undefined = TRUE))
    }
  }
})
