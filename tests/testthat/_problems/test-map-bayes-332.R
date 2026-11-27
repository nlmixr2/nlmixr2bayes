# Extracted from test-map-bayes.R:332

# prequel ----------------------------------------------------------------------
.mapBayesMod <- function() {
  ini({
    tcl <- fix(log(0.5))
    tv <- fix(log(5))
    eta.cl + eta.v ~ fix(0.09,
                         0.03, 0.04)
    add.sd <- fix(0.2)
    prop.sd <- fix(0.1)
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(centr) <- -cl / v * centr
    cp <- centr / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}
.mapBayesData <- function() {
  data.frame(
    id = 1,
    time = c(0, 2, 11.9, 12, 24, 26),
    amt = c(100, NA, NA, 100, 100, NA),
    dv = c(NA, 17.5, 6.1, NA, NA, 13.9),
    evid = c(1, 0, 0, 1, 1, 0)
  )
}

# test -------------------------------------------------------------------------
skip_on_cran()
skip_if_not_installed("rstan")
.d1 <- .mapBayesData()
.d <- rbind(
  transform(.d1, id = 7),
  transform(.d1, id = 3),
  transform(.d1[.d1$evid == 1, ], id = 5)
)
.fit <- suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
  .mapBayesMod,
  .d,
  est = "nuts",
  control = stanControl(
    chains = 2L,
    iter = 3000L,
    warmup = 1000L,
    seed = 23L,
    chainCores = 1L,
    cores = 1L,
    print = 0L,
    calcTables = FALSE,
    ofv = "none",
    onDiagnostic = "none"
  )
)))
.eta <- rstan::extract(.fit$stanfit, pars = "eta")$eta
.lvl <- .fit$env$idLvl
expect_setequal(.lvl, c("3", "5", "7"))
